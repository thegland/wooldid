* probe_fix — completes the audit probes after the mata parse aborts in probe_bd/probe_cef.
* Adds the plain-regress discriminators (raw CRVE vs Stata/reghdfe post-processing).
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
  cnt = 0
  for (j=1; j<=K; j++) {
    bad = 0
    if (missing(d[j])) bad = 1
    if (bad == 0) {
      if (d[j] <= 0) {
        if (isom[j] == 0) bad = 1
      }
    }
    if (bad == 1) {
      if (cnt < 25) {
        printf("%s|BADLIVE|%s|%21.15g\n", tag, nm[j,2], d[j])
        cnt = cnt + 1
      }
    }
  }
  if (cnt >= 25) printf("%s|BADLIVE|...truncated at 25\n", tag)
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

*========== PART 1: regime B margins detail (byCohort ci on d1_r1) ==========
use "`AUD'/d1_r1.dta", clear
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.i#c.(ringer)) cluster(i) tol(1e-9)
local dfr = e(df_r)
mata: diagreport("B_EV")
mata: omitreport("B")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "B_MARG_ATT|rc=" _rc
matrix RJ = r(Jacobian)
mata: margreport("B_MARG_ATT")
mata: JV = st_matrix("RJ")*st_matrix("e(V)")*st_matrix("RJ")'; printf("B_MARG_ATT|max|rV-JVJt|=%12.0g\n", max(abs(st_matrix("r(V)")-JV)))
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "B_MARG_AME|rc=" _rc
matrix RJ = r(Jacobian)
matrix RBame = r(b)
mata: margreport("B_MARG_AME")
mata: JV = st_matrix("RJ")*st_matrix("e(V)")*st_matrix("RJ")'; printf("B_MARG_AME|max|rV-JVJt|=%12.0g\n", max(abs(st_matrix("r(V)")-JV)))

*========== PART 2: regime B via plain regress (raw CRVE discriminator) ==========
cap qui regress y i.i i.t i.i#c.ringer i._i#i._t#i.targetflag##c.dose, vce(cluster i)
dis "BREG|rc=" _rc "|df_r=" e(df_r) "|rank=" e(rank)
local dfrr = e(df_r)
mata: diagreport("BREG_EV")
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "BREG_MARG_AME|rc=" _rc
mata: margreport("BREG_MARG_AME")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "BREG_MARG_ATT|rc=" _rc
mata: margreport("BREG_MARG_ATT")

*========== PART 3: regime B second realization (determinism) ==========
clear
set seed 60000002
amese_mkdata d1
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.i#c.(ringer)) cluster(i) tol(1e-9)
local dfr = e(df_r)
mata: d = diagonal(st_matrix("e(V)")); printf("B_SEED2|evzero=%g of %g\n", sum(d:==0), rows(d))
cap qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
local m1 .
if _rc == 0 mata: st_local("m1", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
local m2 .
if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
dis "B_SEED2|attSEmiss=`m1'|ameSEmiss=`m2'"

*========== PART 4: regime A manual on d1_r1 (same-file contrast) ==========
use "`AUD'/d1_r1.dta", clear
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
local dfra = e(df_r)
mata: diagreport("A_EV")
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfra')
mata: margreport("A_MARG_AME")
qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfra')
mata: margreport("A_MARG_ATT")

*========== PART 5: regime D faithfulness (byTime ci on d1_r1) ==========
use "`AUD'/d1_r1.dta", clear
cap noi wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(byTime) cluster(i)
dis "D_WOOLDID|rc=" _rc
matrix DVES = e(V_es_ame)
matrix DBES = e(b_es_ame)
use "`AUD'/d1_r1.dta", clear
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.t#c.(ringer)) cluster(i) tol(1e-9)
dis "D_MODEL|df_r=" e(df_r) "|rank=" e(rank)
local dfrd = e(df_r)
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrd')
matrix DRV = r(V)
matrix DRB = r(b)
matrix DRJ = r(Jacobian)
mata:
RV = st_matrix("DRV")
JV = st_matrix("DRJ")*st_matrix("e(V)")*st_matrix("DRJ")'
printf("D_FAITH|max|rV - J eV Jt| = %12.0g\n", max(abs(RV-JV)))
DV = st_matrix("DVES")
printf("D_FAITH|dim e(V_es_ame)=%g|dim manual rV=%g\n", rows(DV), rows(RV))
if (rows(DV)==rows(RV)) printf("D_FAITH|max|e(V_es_ame) - manual rV| = %12.0g\n", max(abs(DV-RV)))
DB = st_matrix("DBES")
RB = st_matrix("DRB")
if (cols(DB)==cols(RB)) printf("D_FAITH|max|e(b_es_ame) - manual rb| = %12.0g\n", max(abs(DB-RB)))
ev = symeigenvalues(RV)
printf("D_RANK|eigmax=%12.0g|eigmin=%12.0g\n", max(ev), min(ev))
end

*========== PART 6: regime C blackout dissection (d1c_black vs d1c_clean) ==========
use "`AUD'/d1c_black.dta", clear
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "CBLK_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|rank=" e(rank)
local dfr = e(df_r)
mata: diagreport("CBLK_EV")
mata: omitreport("CBLK")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "CBLK_MARG_ATT|rc=" _rc
mata: margreport("CBLK_MARG_ATT")
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "CBLK_MARG_AME|rc=" _rc
mata: margreport("CBLK_MARG_AME")
* raw-CRVE discriminator via plain regress on the SAME blackout data
cap qui regress y i.i i.t c.ringer i._i#i._t#i.targetflag##c.dose, vce(cluster i)
dis "CBLKREG|rc=" _rc "|df_r=" e(df_r) "|rank=" e(rank)
local dfrr = e(df_r)
mata: diagreport("CBLKREG_EV")
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "CBLKREG_MARG_AME|rc=" _rc
mata: margreport("CBLKREG_MARG_AME")
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "CBLKREG_MARG_ATT|rc=" _rc
mata: margreport("CBLKREG_MARG_ATT")
* phihat fingerprint on the blackout data (reghdfe, ringer visible)
qui reghdfe y c.ringer i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t) cluster(i) tol(1e-9)
mata:
nm = st_matrixcolstripe("e(b)")
V = st_matrix("e(V)")
for (j=1; j<=rows(nm); j++) {
  if (nm[j,2] == "ringer") printf("CBLK_PHIHAT|Var_cl(phihat)=%21.15g\n", V[j,j])
}
end
* clean-realization contrast
use "`AUD'/d1c_clean.dta", clear
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "CCLN_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|rank=" e(rank)
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
* dose-variation forensics: is dose within-cell variation degenerate anywhere in the blackout data?
use "`AUD'/d1c_black.dta", clear
mkvars
bysort i t: egen double dsd = sd(dose)
qui sum dsd if targetflag == 1
dis "CBLK_DOSESD|min cell-sd of dose (treated)=" r(min)
qui gen double doser = round(dose, 1e-12)
qui bysort i t doser: gen byte dup = _N
qui sum dup if targetflag == 1
dis "CBLK_DOSEDUP|max within-cell duplicate count=" r(max)

*========== PART 7: regime E (d4 within ci) ==========
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
mkvars
noi reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "E_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|rank=" e(rank) "|Nclust=" e(N_clust)
local dfr = e(df_r)
mata: diagreport("E_EV")
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
* d5 contrast (same seed; regime-A-style expected)
clear
set seed 60000001
amese_mkdata d5
save "`AUD'/d5_r1.dta", replace
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "E5_MODEL|df_r=" e(df_r) "|rank=" e(rank)
local dfr = e(df_r)
mata: diagreport("E5_EV")
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
mata: margreport("E5_MARG_AME")
* d4 determinism: 2 extra seeds, flags only
forval s = 1/2 {
  clear
  set seed `=60000010+`s''
  amese_mkdata d4
  mkvars
  qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
  local dfr = e(df_r)
  mata: d = diagonal(st_matrix("e(V)")); st_local("z", strofreal(sum(d:==0)))
  cap qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m1 .
  if _rc == 0 mata: st_local("m1", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  cap qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
  local m2 .
  if _rc == 0 mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
  dis "E_SEEDX|s=`s'|evzero=`z'|attSEmiss=`m1'|ameSEmiss=`m2'"
}

*========== PART 8: regime F (d3) ==========
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
dis "PROBE_FIX_DONE"
