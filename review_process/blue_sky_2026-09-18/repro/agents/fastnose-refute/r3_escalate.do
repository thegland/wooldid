*=====================================================================
* r3_escalate.do -- does the FN-03 divergence escalate?
*
* The BS-26 comment at ado:1974-1981 says the BINARY subgroup-by-horizon call was
* given back its VCE under fastnose because otherwise "a structural-zero horizon
* was reported (and ranked by permute()) and the observed and replication column
* sets diverged, aborting bootstrap()/permute() permanalyticse with r(459) at rep 0."
* The CONTINUOUS twin at ado:3290 was never given the same treatment.  So on a
* contreat() + subgroups() + event-study design:
*   ARM P  permute() fast-mode default -> do the structural-zero horizons get
*          permutation p-values and shift intervals?
*   ARM PA permute() + permanalyticse  -> observed run analytic (drops rows), reps
*          fastnose (keep them): does it abort?
*   ARM BT bootstrap()                 -> reps run fastnose internally (ado:322).
*
* Same panel as r1_fn03.do.
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r3_escalate.do
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

*---------------- ARM P: permute fast-mode default -------------------
use `P', clear
ereturn clear
set seed 555
dis _n "############ ARM P: permute(shuffle) fast default ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) subgroups(region) esprelength(3) espostlength(2) ///
    permute(shuffle) permreps(20)
dis "P_RC=" _rc
foreach s in 1 2 {
  tempname B
  cap matrix `B' = e(wdidsgresults_`s')
  if _rc == 0 {
    dis "P_sgres_`s'_ROWS: `: rowfullnames `B''"
    dis "P_sgres_`s'_COLS: `: colnames `B''"
    matrix list `B', format(%12.0g)
  }
  else dis "P_sgres_`s'_ABSENT"
  cap matrix drop `B'
}

*---------------- ARM PA: permute + permanalyticse -------------------
use `P', clear
ereturn clear
set seed 555
dis _n "############ ARM PA: permute + permanalyticse ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) subgroups(region) esprelength(3) espostlength(2) ///
    permute(shuffle) permreps(20) permanalyticse
dis "PA_RC=" _rc
dis "PA_perm_r_done=" e(perm_reps_done)
foreach s in 1 2 {
  tempname B
  cap matrix `B' = e(wdidsgresults_`s')
  if _rc == 0 dis "PA_sgres_`s'_ROWS: `: rowfullnames `B''"
  else dis "PA_sgres_`s'_ABSENT"
  cap matrix drop `B'
}

*---------------- ARM BT: bootstrap ----------------------------------
use `P', clear
ereturn clear
set seed 555
dis _n "############ ARM BT: bootstrap(standard) ############"
cap noi wooldid y cohort t ttre , contreat(dose) contreatcontrols(basic) ///
    fe(unit) cluster(cohort) subgroups(region) esprelength(3) espostlength(2) ///
    bootstrap(standard) reps(20)
dis "BT_RC=" _rc
foreach s in 1 2 {
  tempname B
  cap matrix `B' = e(wdidsgresults_`s')
  if _rc == 0 {
    dis "BT_sgres_`s'_ROWS: `: rowfullnames `B''"
    matrix list `B', format(%12.0g)
  }
  else dis "BT_sgres_`s'_ABSENT"
  cap matrix drop `B'
}

dis "DONE_R3"
