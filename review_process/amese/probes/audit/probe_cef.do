* probe_cef — CODE-PATH AUDIT: regimes C (d1c basic ci), E (d4 within ci), F (d3).
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local AUD /Users/hegland/.claude/jobs/640b9294/tmp/amese/audit
do "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh"
do "`AUD'/mkvars.doh"
set emptycells drop

*========== PART C: d1c realization sweep (manual basic+ci), find a blackout ==========
local foundblack 0
local foundclean 0
forval s = 0/15 {
  clear
  set seed `=60000001+`s''
  amese_mkdata d1c
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  local dfr = e(df_r)
  local rk = e(rank)
  mata: d = diagonal(st_matrix("e(V)")); st_local("evmiss", strofreal(colmissing(d))); st_local("evneg", strofreal(sum((d:<0):*(1:-(d:>=.))))); st_local("evzero", strofreal(sum(d:==0))); st_local("evfullmiss", strofreal(sum(rowmissing(st_matrix("e(V)")))))
  cap qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local rc1 = _rc
  local m1 .
  local minv1 .
  if `rc1' == 0 {
    mata: rt = st_matrix("r(table)"); st_local("m1", strofreal(colmissing(rt[2,.]'))); d = diagonal(st_matrix("r(V)")); st_local("minv1", strofreal(min(d), "%21.0g"))
  }
  cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local rc2 = _rc
  local m2 .
  local minv2 .
  if `rc2' == 0 {
    mata: rt = st_matrix("r(table)"); st_local("m2", strofreal(colmissing(rt[2,.]'))); d = diagonal(st_matrix("r(V)")); st_local("minv2", strofreal(min(d), "%21.0g"))
  }
  dis "C_SWEEP|s=`s'|rank=`rk'|dfr=`dfr'|evdiagmiss=`evmiss'|evdiagzero=`evzero'|evdiagneg=`evneg'|evfullmiss=`evfullmiss'|attrc=`rc1'|attSEmiss=`m1'|attminV=`minv1'|amerc=`rc2'|ameSEmiss=`m2'|ameminV=`minv2'"
  if ("`m2'" != "0" | "`m1'" != "0") & `foundblack' == 0 {
    local foundblack 1
    save "`AUD'/d1c_black.dta", replace
    dis "C_SWEEP|saved d1c_black.dta from s=`s'"
  }
  if "`m2'" == "0" & "`m1'" == "0" & `foundclean' == 0 {
    local foundclean 1
    save "`AUD'/d1c_clean.dta", replace
    dis "C_SWEEP|saved d1c_clean.dta from s=`s'"
  }
}

* d1 parallel sweep (same seeds, iid errors): does d1 basic+ci ever black out?
forval s = 0/15 {
  clear
  set seed `=60000001+`s''
  amese_mkdata d1
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  local dfr = e(df_r)
  cap qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m1 .
  if _rc == 0 mata: st_local("m1", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m2 .
  local minv2 .
  if _rc == 0 {
    mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]'))); st_local("minv2", strofreal(min(diagonal(st_matrix("r(V)"))), "%21.0g"))
  }
  dis "C_D1SWEEP|s=`s'|attSEmiss=`m1'|ameSEmiss=`m2'|ameminV=`minv2'"
}

* dissect the blackout realization if found
cap confirm file "`AUD'/d1c_black.dta"
if _rc == 0 {
  use "`AUD'/d1c_black.dta", clear
  noi reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  dis "CBLK_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|rank=" e(rank) "|Nclust=" e(N_clust)
  local dfr = e(df_r)
  mata:
  V = st_matrix("e(V)")
  nm = st_matrixcolstripe("e(V)")
  d = diagonal(V)
  isslope = strpos(nm[.,2], "dose") :> 0
  isom    = strpos(nm[.,2], "o.") :> 0
  live = 1 :- isom
  printf("CBLK_EVDIAG|K=%g|miss=%g|zero=%g|neg=%g|fullmiss=%g\n", rows(d), colmissing(d), sum(d:==0), sum((d:<0):*(1:-(d:>=.))), sum(rowmissing(V)))
  sl = select(d, isslope :& live); du = select(d, (1:-isslope) :& live)
  printf("CBLK_EVDIAG|live slope n=%g min=%21.15g max=%21.15g\n", rows(sl), min(sl), max(sl))
  printf("CBLK_EVDIAG|live dummy n=%g min=%21.15g max=%21.15g\n", rows(du), min(du), max(du))
  cnt = 0
  for (j=1; j<=rows(d); j++) {
    if ((missing(d[j]) | d[j]<0) & cnt<40) { printf("CBLK_BADDIAG|%s|%21.15g\n", nm[j,2], d[j]); cnt++ }
  }
  end
  cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  dis "CBLK_MARG_ATT|rc=" _rc
  cap matrix RTb = r(table)
  cap matrix RVb = r(V)
  cap matrix RJb = r(Jacobian)
  dis "CBLK_MARG_ATT r(table):"
  cap noi matrix list RTb, format(%14.0g)
  mata:
  d = diagonal(st_matrix("RVb"))
  printf("CBLK_MARG_ATT rV diag raw:\n"); d'
  JV = st_matrix("RJb")*st_matrix("e(V)")*st_matrix("RJb")'
  printf("CBLK_MARG_ATT J eV Jt diag:\n"); diagonal(JV)'
  printf("CBLK_MARG_ATT|max|rV-JVJt|=%12.0g\n", max(abs(st_matrix("RVb")-JV)))
  end
  cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  dis "CBLK_MARG_AME|rc=" _rc
  cap matrix RTa = r(table)
  cap matrix RVa = r(V)
  cap matrix RJa = r(Jacobian)
  dis "CBLK_MARG_AME r(table):"
  cap noi matrix list RTa, format(%14.0g)
  mata:
  d = diagonal(st_matrix("RVa"))
  printf("CBLK_MARG_AME rV diag raw:\n"); d'
  JV = st_matrix("RJa")*st_matrix("e(V)")*st_matrix("RJa")'
  printf("CBLK_MARG_AME J eV Jt diag:\n"); diagonal(JV)'
  end
  * fingerprint: Var_cl(phihat) via visible-ringer refit
  qui reghdfe y c.ringer i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t) cluster(i) tol(1e-9)
  mata:
  nm = st_matrixcolstripe("e(b)")
  V = st_matrix("e(V)")
  for (j=1; j<=rows(nm); j++) {
    if (nm[j,2] == "ringer") printf("CBLK_PHIHAT|Var_cl(phihat)=%21.15g\n", V[j,j])
  }
  end
}
else dis "CBLK|no blackout realization found in 16 seeds"

*========== PART E: d4 within ci ==========
clear
set seed 60000001
amese_mkdata d4
save "`AUD'/d4_r1.dta", replace
do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"
use "`AUD'/d4_r1.dta", clear
cap noi wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) contreatwithin safety(off) cluster(i)
dis "E_WOOLDID|rc=" _rc
foreach m in b V b_es_ame V_es_ame b_es_att V_es_att {
  cap confirm matrix e(`m')
  if _rc == 0 dis "E_E|`m'|EXISTS|" rowsof(e(`m')) "x" colsof(e(`m'))
  else dis "E_E|`m'|ABSENT"
}
matrix EWD = e(wdidmainresults)

use "`AUD'/d4_r1.dta", clear
mkvars
set emptycells drop
noi reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "E_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|rank=" e(rank) "|Nclust=" e(N_clust)
local dfr = e(df_r)
mata:
V = st_matrix("e(V)")
nm = st_matrixcolstripe("e(V)")
d = diagonal(V)
isslope = strpos(nm[.,2], "dose") :> 0
isom    = strpos(nm[.,2], "o.") :> 0
live = 1 :- isom
printf("E_EVDIAG|K=%g|miss=%g|zero=%g|neg=%g|fullmiss=%g|nomit=%g\n", rows(d), colmissing(d), sum(d:==0), sum((d:<0):*(1:-(d:>=.))), sum(rowmissing(V)), sum(isom))
sl = select(d, isslope :& live); du = select(d, (1:-isslope) :& live)
printf("E_EVDIAG|live slope n=%g min=%21.15g max=%21.15g\n", rows(sl), min(sl), max(sl))
printf("E_EVDIAG|live dummy n=%g min=%21.15g max=%21.15g\n", rows(du), min(du), max(du))
cnt = 0
for (j=1; j<=rows(d); j++) {
  if ((missing(d[j]) | d[j]<0) & cnt<40) { printf("E_BADDIAG|%s|%21.15g\n", nm[j,2], d[j]); cnt++ }
}
end
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "E_MARG_ATT|rc=" _rc
cap matrix RTb = r(table)
cap matrix RVb = r(V)
dis "E_MARG_ATT r(table):"
cap noi matrix list RTb, format(%14.0g)
mata:
d = diagonal(st_matrix("RVb"))
printf("E_MARG_ATT rV diag raw:\n"); d'
end
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "E_MARG_AME|rc=" _rc
cap matrix RTa = r(table)
cap matrix RVa = r(V)
dis "E_MARG_AME r(table):"
cap noi matrix list RTa, format(%14.0g)
mata:
d = diagonal(st_matrix("RVa"))
printf("E_MARG_AME rV diag raw:\n"); d'
end

* d5 contrast on same seed (regime A twin: expect positive rank-1-style diag)
clear
set seed 60000001
amese_mkdata d5
save "`AUD'/d5_r1.dta", replace
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "E5_MODEL|df_r=" e(df_r) "|rank=" e(rank)
local dfr = e(df_r)
mata:
d = diagonal(st_matrix("e(V)"))
printf("E5_EVDIAG|K=%g|miss=%g|zero=%g|neg=%g\n", rows(d), colmissing(d), sum(d:==0), sum((d:<0):*(1:-(d:>=.))))
end
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
mata:
d = diagonal(st_matrix("r(V)"))
printf("E5_MARG_AME rV diag raw:\n"); d'
end
qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
mata:
d = diagonal(st_matrix("r(V)"))
printf("E5_MARG_ATT rV diag raw:\n"); d'
end

* d4 determinism: 2 more seeds, flags only
forval s = 1/2 {
  clear
  set seed `=60000010+`s''
  amese_mkdata d4
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  local dfr = e(df_r)
  cap qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m1 .
  if _rc == 0 mata: st_local("m1", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m2 .
  if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  dis "E_SEEDX|s=`s'|attSEmiss=`m1'|ameSEmiss=`m2'"
}

*========== PART F: d3 (cohort-constant dose) ==========
clear
set seed 60000001
amese_mkdata d3
save "`AUD'/d3_r1.dta", replace
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "F_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|rank=" e(rank)
local dfr = e(df_r)
mata:
nm = st_matrixcolstripe("e(b)")
isslope = strpos(nm[.,2], "dose") :> 0
isom    = strpos(nm[.,2], "o.") :> 0
printf("F_MODEL|n slope coefs=%g|n omitted slope coefs=%g|n omitted total=%g\n", sum(isslope), sum(isslope :& isom), sum(isom))
end
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "F_MARG_AME|rc=" _rc
cap matrix RTa = r(table)
cap matrix RVa = r(V)
cap matrix RBa = r(b)
dis "F_MARG_AME r(b) raw:"
cap noi matrix list RBa, format(%21.0g)
dis "F_MARG_AME r(V) raw:"
cap noi matrix list RVa, format(%21.0g)
dis "F_MARG_AME r(table):"
cap noi matrix list RTa, format(%14.0g)
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "F_MARG_ATT|rc=" _rc
mata:
d = diagonal(st_matrix("r(V)"))
printf("F_MARG_ATT rV diag raw (expect positive):\n"); d'
end
* robust variant: same omission structure, same ame blackout?
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) vce(robust) tol(1e-9)
local dfr = e(df_r)
cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
local m2 .
if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
dis "F_ROBUST|amerc=" _rc "|ameSEmiss=`m2'"
dis "PROBE_CEF_DONE"
