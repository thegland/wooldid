*=====================================================================
* repro_FN01_FN02.do   SELF-CONTAINED.  No project data, no project generator.
*
* FN-01: permute() + permanalyticse + fastnose is ACCEPTED, computes no analytic
*        standard error, and then PRINTS that the analytic standard errors are in the
*        table and in e(V).  e(V) does not exist and `test' refuses.
* FN-02: under fastnose the header says "SEs and related materials are not estimated in
*        fast mode" and every Main Results inference cell is ".", yet the gteffects table
*        prints full se / p / 95% CI, byte-identical to the analytic run's.
*
* RUN:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FN01_FN02.do
*=====================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 4242
set sortseed 4242

clear
set obs 120
gen long unit = _n
gen long cohort = 1 + mod(_n-1, 4)
gen double ttre = .
replace ttre = 2003 if cohort == 1
replace ttre = 2005 if cohort == 2
replace ttre = 2007 if cohort == 3
gen double a_i = rnormal()
expand 10
bys unit: gen int t = 1999 + _n
gen byte post = (ttre < . & t >= ttre)
gen double y = a_i + 0.3*(t-1999) + 1.0*post + rnormal()*0.4
tempfile P
qui save `P'

dis _n "================ FN-02 : gteffects, ANALYTIC ================"
use `P', clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) gteffects esprelength(2) espostlength(3)
dis "FN02_AN_RC=" _rc
tempname GA
cap matrix `GA' = e(gteffects)
if _rc == 0 matrix list `GA', format(%20.14g)

dis _n "================ FN-02 : gteffects, FASTNOSE ================"
use `P', clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) gteffects esprelength(2) espostlength(3) fastnose
dis "FN02_FN_RC=" _rc
tempname GF
cap matrix `GF' = e(gteffects)
if _rc == 0 matrix list `GF', format(%20.14g)
* the two matrices should be identical even though the run says SEs were not estimated
cap {
  tempname DD
  matrix `DD' = `GA' - `GF'
  dis "FN02_maxabsdiff=" mreldif(`GA', `GF')
}
dis "FN02_diff_rc=" _rc

dis _n "================ FN-01 : permute + permanalyticse, NO fastnose (control) ======"
use `P', clear
ereturn clear
set seed 999
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(40) ///
    permanalyticse esprelength(2) espostlength(2)
dis "FN01_CTRL_RC=" _rc
dis "FN01_CTRL_vcvstate=" e(wdid_vcvstate)
cap confirm matrix e(V)
dis "FN01_CTRL_hasV_rc=" _rc
cap noi test main_att = 0
dis "FN01_CTRL_test_rc=" _rc

dis _n "================ FN-01 : permute + permanalyticse + FASTNOSE (the defect) ====="
use `P', clear
ereturn clear
set seed 999
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(40) ///
    permanalyticse fastnose esprelength(2) espostlength(2)
dis "FN01_RC=" _rc
dis "FN01_vcvstate=" e(wdid_vcvstate)
dis "FN01_perm_analyticse=" e(perm_analyticse)
dis "FN01_header=" e(wdid_vcespec)
cap confirm matrix e(V)
dis "FN01_hasV_rc=" _rc
cap noi test main_att = 0
dis "FN01_test_rc=" _rc
dis "  ^ the header above claims se_analytic is the analytic cluster(cohort) standard"
dis "    error; the printed se_analytic column is entirely missing and e(V) is absent."
dis "DONE_REPRO_FN01_FN02"
