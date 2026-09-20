*=====================================================================
* r4_controls.do -- pin the mechanism of the r3 aborts and dump the analytic
* subgroup values the permute() run contradicts.
*
* ARM V  analytic run: dump e(wdidsgresults_2) VALUES (what permute() must match)
* ARM B1 bootstrap + contreat + ES, NO subgroups          (control)
* ARM B2 bootstrap + contreat + subgroups, NO ES window   (control)
* ARM B3 bootstrap + BINARY + subgroups + ES              (control: BS-26 fix)
* ARM B4 bootstrap + contreat + subgroups + ES            (the failing cell, repeated)
*
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r4_controls.do
*=====================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 31337
set sortseed 31337

clear
set obs 105
gen long unit = _n
gen byte grp  = 1 + mod(_n-1, 5)
gen double ttre = .
replace ttre = 2013 if grp == 1
replace ttre = 2016 if grp == 2
replace ttre = 2019 if grp == 3
gen long cohort = grp
gen double a_i  = rnormal()*0.8
gen double d_i  = 0.4 + 1.1*runiform()
expand 10
bys unit: gen int t = 2009 + _n
gen byte post = (ttre < . & t >= ttre)
gen double dose = 0
replace dose = d_i if post == 1
gen double y = a_i + 0.35*(t-2009) + 0.9*dose*post + rnormal()*0.4
gen byte region = 1
replace region = 2 if ttre == 2019
tempfile P
qui save `P'

*---- ARM V ----
use `P', clear
ereturn clear
dis _n "############ ARM V analytic ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) subgroups(region) esprelength(3) espostlength(2)
dis "V_RC=" _rc
foreach s in 1 2 {
  tempname B
  cap matrix `B' = e(wdidsgresults_`s')
  if _rc == 0 {
    dis "V_sgres_`s'_ROWS: `: rowfullnames `B''"
    matrix list `B', format(%12.0g)
  }
  cap matrix drop `B'
}

*---- ARM B1: no subgroups ----
use `P', clear
ereturn clear
set seed 555
dis _n "############ ARM B1 bootstrap, NO subgroups ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) esprelength(3) espostlength(2) bootstrap(standard) reps(20)
dis "B1_RC=" _rc

*---- ARM B2: subgroups, no ES ----
use `P', clear
ereturn clear
set seed 555
dis _n "############ ARM B2 bootstrap, subgroups, NO ES ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) subgroups(region) bootstrap(standard) reps(20)
dis "B2_RC=" _rc

*---- ARM B3: BINARY + subgroups + ES ----
use `P', clear
ereturn clear
set seed 555
dis _n "############ ARM B3 bootstrap, BINARY, subgroups + ES ############"
cap noi wooldid y cohort t ttre , ///
    fe(unit) cluster(cohort) subgroups(region) esprelength(3) espostlength(2) ///
    bootstrap(standard) reps(20)
dis "B3_RC=" _rc

*---- ARM B4: the failing cell again ----
use `P', clear
ereturn clear
set seed 999
dis _n "############ ARM B4 bootstrap, contreat + subgroups + ES ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) subgroups(region) esprelength(3) espostlength(2) ///
    bootstrap(standard) reps(20)
dis "B4_RC=" _rc

*---- ARM B5: same but robust (not cluster) ----
use `P', clear
ereturn clear
set seed 999
dis _n "############ ARM B5 bootstrap, contreat + subgroups + ES, robust ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) subgroups(region) esprelength(3) espostlength(2) ///
    bootstrap(standard) reps(20)
dis "B5_RC=" _rc

dis "DONE_R4"
