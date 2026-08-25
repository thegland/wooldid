* BLIND-probe 1: byTime + cluster(i) fingerprint on d1, seed 70710001 (= d1_bytime_ci rep 1)
* Prediction: reported horizon-h ame SE == SE_cl(phi_hat_t(h)) (per-calendar-period absorbed
* baseline dose slope), up to a common dof factor; V_es_ame full rank but entirely the
* phi_t covariance block (within-cell slope noise annihilated).
clear all
set more off
set maxvar 20000
set linesize 240
do /Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado
do /Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh
set seed 70710001
amese_mkdata d1
tempfile dat
save `dat'

wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(byTime) cluster(i) verbose

matrix W = e(wdidmainresults)
matrix list W, format(%12.8f)
matrix V = e(V_es_ame)
mata:
V = st_matrix("V")
printf("reported ES-ame SEs (sqrt diag V_es_ame):\n")
sqrt(diagonal(V))'
printf("eigenvalues of V_es_ame:\n")
symeigenvalues(V)
end

* ---- visible per-period-ringer refit ----
use `dat', clear
gen rel = t - ttre
gen _t = rel + 6
replace _t = 0 if _t == . | _t == 5
gen byte targetflag = _t > 0 & _t < .
gen double _i = i if targetflag
replace _i = 0 if _i == .
gen double ringer = dose
reghdfe y ibn.t#c.ringer i._i#i._t#i.targetflag##c.dose , absorb(i.i i.t) cluster(i) tol(1e-9)
dis as result "visible per-period baseline dose slopes phi_t, cluster(i) SEs:"
forval tt = 1/12 {
  cap dis "t=`tt'  phi_t = " %12.8f _b[`tt'.t#c.ringer] "   SE_cl = " %12.8f _se[`tt'.t#c.ringer]
}
dis as result "horizon map: pre5..pre2 = calendar t 1..4 ; post0..post3 = calendar t 6..9 ; base t=5"
dis as result "compare: reported h SEs (CSV rep1): .1217 .0788 .1156 .1109 | .1719 .0938 .1209 .1265"
exit, clear
