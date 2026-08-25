* probe_fix2 — completes the audit: verbose reghdfe VCE trace for the zeroing branch,
* d1c clean-realization contrast, regime E (d4/d5), regime F (d3).
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local AUD /Users/hegland/.claude/jobs/640b9294/tmp/amese/audit
do "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh"
do "`AUD'/mkvars.doh"
do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"
set emptycells drop

mata:
void diagreport(string scalar tag)
{
  V = st_matrix("e(V)")
  nm = st_matrixcolstripe("e(V)")
  d = diagonal(V)
  isslope = strpos(nm[.,2], "dose") :> 0
  isom = strpos(nm[.,2], "o.") :> 0
  K = rows(d)
  npos = 0; nzero = 0; nneg = 0; nmiss = 0
  for (j=1; j<=K; j++) {
    if (missing(d[j])) nmiss = nmiss + 1
    else if (d[j] == 0) nzero = nzero + 1
    else if (d[j] < 0) nneg = nneg + 1
    else npos = npos + 1
  }
  printf("%s|K=%g|pos=%g|zero=%g|neg=%g|miss=%g|fullmiss=%g|nomit=%g\n", tag, K, npos, nzero, nneg, nmiss, sum(rowmissing(V)), sum(isom))
  live = 1 :- isom
  sl = select(d, isslope :& live)
  du = select(d, (1:-isslope) :& live)
  if (rows(sl) > 0) printf("%s|live slope n=%g min=%21.15g max=%21.15g\n", tag, rows(sl), min(sl), max(sl))
  if (rows(du) > 0) printf("%s|live dummy n=%g min=%21.15g max=%21.15g\n", tag, rows(du), min(du), max(du))
}
void omitreport(string scalar tag)
{
  nm = st_matrixcolstripe("e(b)")
  for (j=1; j<=rows(nm); j++) {
    if (strpos(nm[j,2], "o.") > 0) printf("%s|OMIT|%s\n", tag, nm[j,2])
  }
}
void margreport(string scalar tag)
{
  rt = st_matrix("r(table)")
  rv = st_matrix("r(V)")
  d = diagonal(rv)
  printf("%s|b row:  ", tag)
  for (j=1; j<=cols(rt); j++) printf(" %12.0g", rt[1,j])
  printf("\n%s|se row: ", tag)
  for (j=1; j<=cols(rt); j++) printf(" %12.0g", rt[2,j])
  printf("\n%s|rV diag:", tag)
  for (j=1; j<=rows(d); j++) printf(" %21.15g", d[j])
  printf("\n%s|semiss=%g|vdiagmiss=%g|vdiagzero=%g|vdiagneg=%g\n", tag, colmissing(rt[2,.]'), colmissing(d), sum(d:==0), sum((d:<0):*(1:-(d:>=.))))
}
end

*===== 2. VERBOSE VCE trace: d1c blackout realization (basic ci) =====
use "`AUD'/d1c_black.dta", clear
foreach v in rel _t targetflag _i null ringer dsd doser dup {
  cap drop `v'
}
mkvars
noi reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9) verbose(1) nosample
dis "VTRACE_CBLK|df_a=" e(df_a) "|df_m=" e(df_m) "|df_r=" e(df_r)

*===== 3. d1c clean-realization contrast =====
local cleanfound 0
foreach s in 60000003 60000004 60000005 60000006 {
  if `cleanfound' == 0 {
    clear
    set seed `s'
    amese_mkdata d1c
    mkvars
    qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
    local dfr = e(df_r)
    mata: d = diagonal(st_matrix("e(V)")); st_local("z", strofreal(sum(d:==0)))
    cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
    local m2 .
    if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
    dis "CCLN_TRY|seed=`s'|evzero=`z'|ameSEmiss=`m2'"
    if "`m2'" == "0" {
      local cleanfound 1
      save "`AUD'/d1c_clean.dta", replace
      dis "CCLN|saved d1c_clean.dta from seed `s'"
    }
  }
}
cap confirm file "`AUD'/d1c_clean.dta"
if _rc == 0 {
  use "`AUD'/d1c_clean.dta", clear
  foreach v in rel _t targetflag _i null ringer {
    cap drop `v'
  }
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  dis "CCLN_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|df_a=" e(df_a) "|rank=" e(rank)
  mata: diagreport("CCLN_EV")
  mata: omitreport("CCLN")
  qui reghdfe y c.ringer i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t) cluster(i) tol(1e-9)
  mata:
  nm = st_matrixcolstripe("e(b)")
  V = st_matrix("e(V)")
  for (j=1; j<=rows(nm); j++) {
    if (nm[j,2] == "ringer") printf("CCLN_PHIHAT|Var_cl(phihat)=%21.15g\n", V[j,j])
  }
  end
}
* phihat point estimate + variance on the BLACKOUT data, and dose forensics
use "`AUD'/d1c_black.dta", clear
foreach v in rel _t targetflag _i null ringer dsd doser dup {
  cap drop `v'
}
mkvars
qui reghdfe y c.ringer i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t) cluster(i) tol(1e-9)
mata:
nm = st_matrixcolstripe("e(b)")
b = st_matrix("e(b)")
V = st_matrix("e(V)")
for (j=1; j<=rows(nm); j++) {
  if (nm[j,2] == "ringer") printf("CBLK_PHIHAT2|phihat=%21.15g|Var_cl=%21.15g\n", b[j], V[j,j])
}
end
mata: diagreport("CBLK_RINGERVIS_EV")
bysort i t: egen double dsd = sd(dose)
qui sum dsd if targetflag == 1
dis "CBLK_DOSESD|min cell-sd of dose (treated)=" r(min) "|max=" r(max)
qui levelsof i if ttre != ., local(tcoh)
local ncohconst 0
foreach g of local tcoh {
  qui sum dose if i == `g'
  if r(sd) < 1e-12 local ++ncohconst
}
dis "CBLK_DOSECOH|n treated cohorts with cohort-constant dose=`ncohconst'"

*===== 4. regime E (d4 within ci) =====
clear
set seed 60000001
amese_mkdata d4
save "`AUD'/d4_r1.dta", replace
use "`AUD'/d4_r1.dta", clear
cap noi wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) contreatwithin safety(off) cluster(i)
dis "E_WOOLDID|rc=" _rc
foreach m in b V b_es_ame V_es_ame b_es_att V_es_att {
  cap confirm matrix e(`m')
  if _rc == 0 dis "E_E|`m'|EXISTS|" rowsof(e(`m')) "x" colsof(e(`m'))
  else dis "E_E|`m'|ABSENT"
}
use "`AUD'/d4_r1.dta", clear
foreach v in rel _t targetflag _i null ringer dsd doser dup {
  cap drop `v'
}
mkvars
noi reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9) verbose(1) nosample
dis "E_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|df_a=" e(df_a) "|rank=" e(rank) "|Nclust=" e(N_clust)
local dfr = e(df_r)
mata: diagreport("E_EV")
mata: omitreport("E")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "E_MARG_ATT|rc=" _rc
mata: margreport("E_MARG_ATT")
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "E_MARG_AME|rc=" _rc
mata: margreport("E_MARG_AME")
* plain-regress discriminator on d4
cap qui regress y i.i i.t c.ringer i._i#i._t#i.targetflag##c.dose, vce(cluster i)
dis "EREG|rc=" _rc "|df_r=" e(df_r) "|rank=" e(rank)
local dfrr = e(df_r)
mata: diagreport("EREG_EV")
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "EREG_MARG_AME|rc=" _rc
mata: margreport("EREG_MARG_AME")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "EREG_MARG_ATT|rc=" _rc
mata: margreport("EREG_MARG_ATT")
* d5 contrast (same seed)
clear
set seed 60000001
amese_mkdata d5
save "`AUD'/d5_r1.dta", replace
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "E5_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|df_a=" e(df_a) "|rank=" e(rank)
local dfr = e(df_r)
mata: diagreport("E5_EV")
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
mata: margreport("E5_MARG_AME")
* d4 determinism: 2 extra seeds
forval s = 1/2 {
  clear
  set seed `=60000010+`s''
  amese_mkdata d4
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  local dfr = e(df_r)
  mata: d = diagonal(st_matrix("e(V)")); st_local("z", strofreal(sum(d:==0)))
  cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m2 .
  if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  dis "E_SEEDX|s=`s'|evzero=`z'|ameSEmiss=`m2'"
}

*===== 5. regime F (d3) =====
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
isom = strpos(nm[.,2], "o.") :> 0
printf("F_MODEL|n slope coefs=%g|n omitted slope=%g|n omitted total=%g\n", sum(isslope), sum(isslope :& isom), sum(isom))
end
mata: diagreport("F_EV")
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "F_MARG_AME|rc=" _rc
mata: margreport("F_MARG_AME")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "F_MARG_ATT|rc=" _rc
mata: margreport("F_MARG_ATT")
* robust variant
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) vce(robust) tol(1e-9)
local dfr = e(df_r)
cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
local m2 .
if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
dis "F_ROBUST|amerc=" _rc "|ameSEmiss=`m2'"


*===== 6. ereturn-post screen mini-experiment: which V property triggers the zeroing? =====
cap program drop postexp
program define postexp, eclass
  args case
  tempname b V
  matrix `b' = (1, 2, 3)
  if "`case'" == "S1" matrix `V' = (1, 0, 0 \ 0, 1e-30, 0 \ 0, 0, 0)
  if "`case'" == "S2" matrix `V' = (1, 0, 0 \ 0, -1e-22, 0 \ 0, 0, 1)
  if "`case'" == "S3" matrix `V' = (1, .5, 0 \ .5 + 1e-12, 1, 0 \ 0, 0, 1)
  if "`case'" == "S4" matrix `V' = (1, .5, 0 \ .5 + 1e-6, 1, 0 \ 0, 0, 1)
  if "`case'" == "S5" matrix `V' = (1, 0, 0 \ 0, 1e-260, 0 \ 0, 0, 1)
  matrix colnames `b' = x1 x2 x3
  matrix colnames `V' = x1 x2 x3
  matrix rownames `V' = x1 x2 x3
  ereturn post `b' `V'
end
foreach c in S1 S2 S3 S4 S5 {
  cap noi postexp `c'
  dis "POSTEXP|`c'|rc=" _rc
  cap noi mata: d = diagonal(st_matrix("e(V)")); printf("POSTEXP|`c'|eV diag: %12.0g %12.0g %12.0g|allzero=%g\n", d[1], d[2], d[3], all(st_matrix("e(V)") :== 0))
}
dis "PROBE_FIX3_DONE"
