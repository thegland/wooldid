* probe_fix4 — final probe: regime E margins/regress/d5/seeds, regime F, POSTEXP.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local AUD /Users/hegland/.claude/jobs/640b9294/tmp/amese/audit
do "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh"
do "`AUD'/mkvars.doh"
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
void margreport(string scalar tag)
{
  rt = st_matrix("r(table)")
  rv = st_matrix("r(V)")
  if (rows(rt) < 2 | rows(rv) < 1) {
    printf("%s|EMPTY r(table)/r(V)\n", tag)
    return
  }
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

*===== regime E margins detail (d4_r1, reghdfe, WITH e(sample)) =====
use "`AUD'/d4_r1.dta", clear
foreach v in rel _t targetflag _i null ringer {
  cap drop `v'
}
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "E_MODEL2|df_r=" e(df_r) "|df_m=" e(df_m) "|df_a=" e(df_a) "|rank=" e(rank)
local dfr = e(df_r)
mata: diagreport("E_EV2")
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

*===== d5 contrast (same seed) =====
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

*===== d4 determinism: 2 extra seeds =====
forval s = 1/2 {
  clear
  set seed `=60000010+`s''
  amese_mkdata d4
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  local dfr = e(df_r)
  mata: d = diagonal(st_matrix("e(V)")); st_local("z", strofreal(sum(d:==0))); st_local("kk", strofreal(rows(d)))
  cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m2 .
  if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  dis "E_SEEDX|s=`s'|evzero=`z' of `kk'|ameSEmiss=`m2'"
}

*===== regime F (d3) =====
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
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) vce(robust) tol(1e-9)
local dfr = e(df_r)
cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
local m2 .
if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
dis "F_ROBUST|amerc=" _rc "|ameSEmiss=`m2'"

*===== POSTEXP: which V property triggers Stata's ereturn-post zeroing? =====
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
dis "PROBE_FIX4_DONE"
