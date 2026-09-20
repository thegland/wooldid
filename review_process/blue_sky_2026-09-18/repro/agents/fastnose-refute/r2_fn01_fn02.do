*=====================================================================
* r2_fn01_fn02.do -- SKEPTIC's INDEPENDENT test of FASTNOSE-01 and FASTNOSE-02.
*
* FN-01 CLAIM: permute() + permanalyticse + fastnose is accepted, computes NO
*   analytic SE at all (byte-identical to the plain permute() fast default), yet
*   prints the permanalyticse header + footer asserting e(V)/e(V_es_*) hold the
*   analytic covariance, and posts e(perm_analyticse)=1.
* FN-02 CLAIM: fastnose + gteffects prints a full se/p/ci gteffects table under a
*   header saying no SEs were estimated, and "the gteffects paragraph never
*   mentions fastnose".
*
* My own binary panel, different from the prober's.
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r2_fn01_fn02.do
*=====================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 90210
set sortseed 90210

clear
set obs 96
gen long unit = _n
gen byte grp = 1 + mod(_n-1, 6)
gen double ttre = .
replace ttre = 2004 if grp == 1
replace ttre = 2006 if grp == 2
replace ttre = 2008 if grp == 3
replace ttre = 2010 if grp == 4
gen long cohort = grp
gen double a_i = rnormal()
expand 12
bys unit: gen int t = 2000 + _n
gen byte post = (ttre < . & t >= ttre)
gen double y = a_i + 0.2*(t-2000) + 1.4*post + rnormal()*0.5
tempfile P
qui save `P'

*=====================================================================
* PART 1 (FN-01): three arms, RNG pinned identically before each call
*=====================================================================
* ARM A: permute + permanalyticse (no fastnose)
use `P', clear
ereturn clear
set seed 4242
dis _n "############ ARM A: permute + permanalyticse ############"
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) ///
    permute(shuffle) permreps(40) permanalyticse esprelength(2) espostlength(2)
dis "A_RC=" _rc
dis "A_vcvstate=" "`e(wdid_vcvstate)'"
dis "A_vcespec=" "`e(wdid_vcespec)'"
dis "A_perm_analyticse=" e(perm_analyticse)
cap confirm matrix e(V)
dis "A_eV_rc=" _rc
cap noi test main_att = 0
dis "A_test_rc=" _rc
tempname AP
cap matrix `AP' = e(perm_p)
if _rc == 0 matrix list `AP', format(%21x)

* ARM B: permute + permanalyticse + fastnose
use `P', clear
ereturn clear
set seed 4242
dis _n "############ ARM B: permute + permanalyticse + fastnose ############"
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) ///
    permute(shuffle) permreps(40) permanalyticse fastnose esprelength(2) espostlength(2)
dis "B_RC=" _rc
dis "B_vcvstate=" "`e(wdid_vcvstate)'"
dis "B_vcespec=" "`e(wdid_vcespec)'"
dis "B_perm_analyticse=" e(perm_analyticse)
cap confirm matrix e(V)
dis "B_eV_rc=" _rc
cap noi test main_att = 0
dis "B_test_rc=" _rc
cap noi lincom main_att
dis "B_lincom_rc=" _rc
tempname BP
cap matrix `BP' = e(perm_p)
if _rc == 0 matrix list `BP', format(%21x)
tempname BM
cap matrix `BM' = e(wdidmainresults)
if _rc == 0 {
  dis "B_mainresults_COLS: `: colnames `BM''"
  matrix list `BM', format(%12.0g)
}

* ARM C: plain permute (fast-mode default)
use `P', clear
ereturn clear
set seed 4242
dis _n "############ ARM C: permute plain (fast default) ############"
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) ///
    permute(shuffle) permreps(40) esprelength(2) espostlength(2)
dis "C_RC=" _rc
dis "C_vcvstate=" "`e(wdid_vcvstate)'"
dis "C_vcespec=" "`e(wdid_vcespec)'"
dis "C_perm_analyticse=" e(perm_analyticse)
cap confirm matrix e(V)
dis "C_eV_rc=" _rc
tempname CP
cap matrix `CP' = e(perm_p)
if _rc == 0 matrix list `CP', format(%21x)
tempname CM
cap matrix `CM' = e(wdidmainresults)
if _rc == 0 {
  dis "C_mainresults_COLS: `: colnames `CM''"
  matrix list `CM', format(%12.0g)
}

*=====================================================================
* PART 2 (FN-02): gteffects analytic vs fastnose
*=====================================================================
use `P', clear
ereturn clear
dis _n "############ GT ANALYTIC ############"
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) gteffects esprelength(2) espostlength(2)
dis "GTA_RC=" _rc
tempname GA
cap matrix `GA' = e(gteffects)
if _rc == 0 {
  dis "GTA_dim=" rowsof(`GA') "x" colsof(`GA')
  dis "GTA_cols: `: colnames `GA''"
  matrix list `GA', format(%21x)
}

use `P', clear
ereturn clear
dis _n "############ GT FASTNOSE ############"
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) gteffects esprelength(2) espostlength(2) fastnose
dis "GTF_RC=" _rc
tempname GF
cap matrix `GF' = e(gteffects)
if _rc == 0 {
  dis "GTF_dim=" rowsof(`GF') "x" colsof(`GF')
  dis "GTF_cols: `: colnames `GF''"
  matrix list `GF', format(%21x)
}

dis "DONE_R2"
