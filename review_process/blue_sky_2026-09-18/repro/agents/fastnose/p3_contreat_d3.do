* p3 -- the "state policy intensity" design: the dose is constant inside every
* cohort-period cell.  wooldid's R6 screen (ado:11958) voids the ame block there,
* because the regression's ame entries are exact structural ZEROS, not estimates.
* That screen is skipped under fastnose.  Does fast mode report the zeros as effects?
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

local gN 60
local gT 10
local gonsets "2004 2006 2008"
local gysd 0.30
local gdose cellconst
include "`GEN'/gen_core.do"
qui save "`WD'/p3panel.dta", replace
dis "NOBS=" _N
tabstat dose , by(cohort) stat(min max n)

foreach arm in A B {
  if "`arm'" == "A" local spec "y unit t ttre , contreat(dose) contreatcontrols(untreatedZero) esprelength(2) espostlength(3)"
  if "`arm'" == "B" local spec "y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) esprelength(2) espostlength(3)"

  dis "########## ARM `arm' ANALYTIC : `spec'"
  use "`WD'/p3panel.dta", clear
  ereturn clear
  cap noi wooldid `spec'
  dis "ARM`arm'_AN_RC=" _rc
  if _rc == 0 {
    dis "ARM`arm'_AN_MAIN:"
    mat list e(wdidmainresults), format(%12.0g)
    cap confirm matrix e(b_es_ame)
    dis "ARM`arm'_AN_b_es_ame_rc=" _rc
  }

  dis "########## ARM `arm' FASTNOSE : `spec' fastnose"
  use "`WD'/p3panel.dta", clear
  ereturn clear
  cap noi wooldid `spec' fastnose
  dis "ARM`arm'_FN_RC=" _rc
  if _rc == 0 {
    dis "ARM`arm'_FN_MAIN:"
    mat list e(wdidmainresults), format(%12.0g)
    cap confirm matrix e(b_es_ame)
    dis "ARM`arm'_FN_b_es_ame_rc=" _rc
    cap mat list e(b_es_ame), format(%12.0g)
  }
}
dis "DONE_P3"
