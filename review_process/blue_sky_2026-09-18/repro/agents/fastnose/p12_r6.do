* p12 -- reach the SI-A4/R6 "d3" class: every treated cohort-period DOSE-SLOPE term
* dropped for multicollinearity while the cohort-period DUMMIES survive.  On that class
* the analytic run names the situation and reports the ame results as MISSING; the screen
* is skipped under fastnose (ado:11958), so fast mode should report the regression's
* structural ZEROS as if they were average marginal effects.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fastnose"
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gN 90
local gT 10
local gonsets "2004 2006 2008"
local gysd 0.30
include "`GEN'/gen_core.do"
* the STATE-POLICY-INTENSITY design: the dose is a property of the cohort-period cell,
* not of the unit.  30 units per cohort, so the cohort-period DUMMY is well identified.
gen double cdose = 0
qui replace cdose = 0.5 + 0.4*cohort + 0.2*max(0,h) if post == 1
qui replace y = y + 0.8*cdose*(post==1)
tabstat cdose , by(cohort) stat(min max n)
qui save "`WD'/p12panel.dta", replace

local a1 "y cohort t ttre , contreat(cdose) contreatcontrols(basic) fe(unit) esprelength(2) espostlength(3)"
local a2 "y cohort t ttre , contreat(cdose) contreatcontrols(byTime) fe(unit) esprelength(2) espostlength(3)"
local a3 "y cohort t ttre , contreat(cdose) contreatcontrols(byCohort) fe(unit) esprelength(2) espostlength(3)"
local a4 "y cohort t ttre , contreat(cdose) contreatcontrols(untreatedZero) fe(unit) esprelength(2) espostlength(3)"
local a5 "y cohort t ttre , contreat(cdose) contreatcontrols(basic) esprelength(2) espostlength(3)"

forval k = 1/5 {
  foreach mode in AN FN {
    local ex
    if "`mode'" == "FN" local ex fastnose
    dis "########## A`k' `mode' : `a`k'' `ex'"
    use "`WD'/p12panel.dta", clear
    ereturn clear
    cap noi wooldid `a`k'' `ex'
    dis "A`k'_`mode'_RC=" _rc
    if _rc == 0 {
      tempname MR
      cap matrix `MR' = e(wdidmainresults)
      if _rc == 0 {
        dis "A`k'_`mode'_MAIN:"
        matrix list `MR', format(%12.0g)
        cap matrix drop `MR'
      }
      cap confirm matrix e(b_es_ame)
      dis "A`k'_`mode'_b_es_ame_rc=" _rc
      tempname BA
      cap matrix `BA' = e(b_es_ame)
      if _rc == 0 matrix list `BA', format(%12.0g)
      cap matrix drop `BA'
    }
  }
}
dis "DONE_P12"
