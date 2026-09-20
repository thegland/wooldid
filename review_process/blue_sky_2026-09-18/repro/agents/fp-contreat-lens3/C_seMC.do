* LENS3 probe C: which reported SE, if any, is honest?  8 draws, all SE arms.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
run "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat-lens3/mkfix.do"
tempname P
postfile `P' int rep double aHC1 double sHC1 double aOLS double sOLS double aUNC double sUNC using Cres , replace
forval r = 1/8 {
  local s = 70000 + `r'
  set seed `s'
  set sortseed `s'
  mkfix . 0
  tempfile G
  qui save `G'
  local a1 = .
  local s1 = .
  local a2 = .
  local s2 = .
  local a3 = .
  local s3 = .
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) customvce(robust)
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a1 = M[2,1]
    local s1 = M[2,2]
  }
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) customvce(ols)
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a2 = M[2,1]
    local s2 = M[2,2]
  }
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) unconditionalse robust
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a3 = M[2,1]
    local s3 = M[2,2]
  }
  dis "C rep`r' HC1 a=" %12.5f `a1' " se=" %12.5f `s1' " | OLS a=" %12.5f `a2' " se=" %12.5f `s2' " | UNC a=" %12.5f `a3' " se=" %12.5f `s3'
  post `P' (`r') (`a1') (`s1') (`a2') (`s2') (`a3') (`s3')
}
postclose `P'
use Cres, clear
foreach v in HC1 OLS UNC {
  qui sum s`v', detail
  dis "CSUM `v' reported se: n=" r(N) " mean=" %12.5f r(mean) " p50=" %12.5f r(p50) " min=" %12.5f r(min) " max=" %12.5f r(max)
  qui gen byte cov`v' = (a`v' - 1.96*s`v' <= 0.4) & (a`v' + 1.96*s`v' >= 0.4) if s`v'<.
  qui count if cov`v'==1
  local nc = r(N)
  qui count if cov`v'<.
  dis "CSUM `v' 95% CI covers 0.400 in `nc' of " r(N) " draws"
}
qui sum aHC1
dis "CSUM ame across C draws: mean=" %12.5f r(mean) " sd=" %12.5f r(sd)
dis "DONE_C"
