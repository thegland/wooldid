* LENS3 probe B: realized sampling distribution of the lattice ame, fast arms.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
run "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat-lens3/mkfix.do"

* -- B0: does the prober's extra gquantiles line change the draw? -----------
set seed 121212
set sortseed 121212
mkfix . 1
cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) customvce(robust)
mat M = e(wdidmainresults)
dis "B0_WITHLATU rc=" _rc " att=" %14.8f M[1,1] " ame=" %14.8f M[2,1] " se=" %14.8g M[2,2]
set seed 121212
set sortseed 121212
mkfix . 0
cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) customvce(robust)
mat M = e(wdidmainresults)
dis "B0_NOLATU   rc=" _rc " att=" %14.8f M[1,1] " ame=" %14.8f M[2,1] " se=" %14.8g M[2,2]

* -- B1: does the all-zero VCV also happen under the iid default? ----------
set seed 4242
set sortseed 4242
mkfix . 0
tempfile Z
qui save `Z'
foreach sp in "lattice(all 3)" "lattice(all 3) robust" "lattice(all 3) cluster(st)" {
  use `Z', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) `sp'
  local rc = _rc
  mat M = e(wdidmainresults)
  dis "B1 SPEC=`sp' rc=`rc' ame=" %14.8f M[2,1] " se=" %14.8g M[2,2]
}

* -- B2: 25 draws, fast reghdfe arms only -----------------------------------
tempname P
postfile `P' int rep double amelat double selat double amenol double senol using Bres , replace
forval r = 1/25 {
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
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) robust
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a1 = M[2,1]
    local s1 = M[2,2]
  }
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) robust
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a2 = M[2,1]
    local s2 = M[2,2]
  }
  dis "B2 rep`r' latame=" %14.6f `a1' " latse=" %12.6f `s1' "  nolatame=" %14.6f `a2' " nolatse=" %12.6f `s2'
  post `P' (`r') (`a1') (`s1') (`a2') (`s2')
}
postclose `P'
use Bres, clear
qui sum amelat, detail
dis "B2SUM LAT ame: n=" r(N) " mean=" %12.6f r(mean) " sd=" %12.6f r(sd) " p50=" %12.6f r(p50) " min=" %12.6f r(min) " max=" %12.6f r(max)
qui gen double dev = abs(amelat-0.4)
qui sum dev, detail
dis "B2SUM LAT |ame-0.4|: p50=" %12.6f r(p50) " p90=" %12.6f r(p90)
qui sum amenol, detail
dis "B2SUM NOLAT ame: mean=" %12.6f r(mean) " sd=" %12.6f r(sd) " p50=" %12.6f r(p50)
qui sum senol
dis "B2SUM NOLAT reported se: mean=" %12.6f r(mean)
qui count if selat < .
dis "B2SUM LAT draws with a NON-MISSING default-path se: " r(N) " of 25"
dis "DONE_B"
