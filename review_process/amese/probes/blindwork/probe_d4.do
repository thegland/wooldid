* BLIND-probe 2: d4 (all-treated staggered, contreatwithin) + cluster(i).
* Rep A = seed 71910003 (SEs post, near-identical ~.2934876): prediction = SE_cl(pooled phi_hat).
* Rep B = seed 71910001 (total SE blackout): capture the caution text + cause.
clear all
set more off
set maxvar 20000
set linesize 240
do /Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado
do /Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh

* ---------- Rep A ----------
set seed 71910003
amese_mkdata d4
tempfile datA
save `datA'
wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) contreatwithin safety(off) cluster(i) verbose
matrix W = e(wdidmainresults)
matrix WA = W[1..., 1..2]
dis as result "REP A e(wdidmainresults) est/se:"
matrix list WA, format(%12.8f)
matrix V = e(V_es_ame)
mata:
V = st_matrix("V")
printf("REP A V_es_ame sqrt-diag / eigenvalues:\n")
sqrt(diagonal(V))'
symeigenvalues(V)
end

* visible pooled-ringer refit on rep A data
use `datA', clear
gen rel = t - ttre
gen _t = rel + 9
replace _t = 0 if _t == . | _t == 8
gen byte targetflag = _t > 0 & _t < .
gen double _i = i if targetflag
replace _i = 0 if _i == .
gen double ringer = dose
reghdfe y c.ringer i._i#i._t#i.targetflag##c.dose , absorb(i.i i.t) cluster(i) tol(1e-9)
dis as result "REP A visible pooled ringer: phi = " %12.8f _b[ringer] "  SE_cl = " %12.8f _se[ringer]
dis as result "compare reported (CSV rep3): all 10 SEs ~ .2934876"

* ---------- Rep B (blackout) ----------
clear
set seed 71910001
amese_mkdata d4
tempfile datB
save `datB'
wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) contreatwithin safety(off) cluster(i)
matrix W = e(wdidmainresults)
matrix WB = W[1..., 1..2]
dis as result "REP B e(wdidmainresults) est/se:"
matrix list WB, format(%12.8f)
dis as result "REP B e(es_ndegen_ame) = " e(es_ndegen_ame)

* visible pooled-ringer refit on rep B data: is phi dropped / SE missing?
use `datB', clear
gen rel = t - ttre
gen _t = rel + 9
replace _t = 0 if _t == . | _t == 8
gen byte targetflag = _t > 0 & _t < .
gen double _i = i if targetflag
replace _i = 0 if _i == .
gen double ringer = dose
reghdfe y c.ringer i._i#i._t#i.targetflag##c.dose , absorb(i.i i.t) cluster(i) tol(1e-9)
dis as result "REP B visible pooled ringer: phi = " %12.8f _b[ringer] "  SE_cl = " %12.8f _se[ringer]
* how many visible-cell coefficients got missing SEs?
mata:
b = st_matrix("e(b)"); Vf = st_matrix("e(V)")
d = diagonal(Vf)
printf("REP B refit: %g coefficients, %g with zero/neg V diagonal\n", rows(d), sum(d :<= 0))
end
exit, clear
