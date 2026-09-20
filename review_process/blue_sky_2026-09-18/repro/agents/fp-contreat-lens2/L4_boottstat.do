* LENS 2 part 4: does the zero-V spec leak a NUMBER anywhere? Two routes a real
* user would take after reading wooldid's own caution ("switch to bootstrap SEs,
* IM p-values, or randomization inference"). boottstat divides by the analytic
* SE, which on this spec is zero.  NOTE: this exact command terminated the L2
* Stata run mid-execution; this file isolates it.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242
which wooldid
clear
set obs 12
gen int st = _n
gen int ncty = 5 + mod(_n*7, 8)
expand ncty
gen long cty = _n
expand 10
bysort cty: gen int year = _n
gen double ttre = .
qui replace ttre = 4 if st <= 3
qui replace ttre = 6 if st > 3 & st <= 6
qui replace ttre = 8 if st > 6 & st <= 8
gen byte post = (ttre < . & year >= ttre)
bysort cty (year): gen double dbase = 2*runiform() if _n == 1
bysort cty (year): replace dbase = dbase[1]
gen double dose = dbase + rnormal(0, 0.5)
gen double aunit = 0
bysort cty (year): replace aunit = 3*runiform() if _n == 1
bysort cty (year): replace aunit = aunit[1]
gen double y = 5 + aunit + 0.25*year + 0.7*dose + 0.4*dose*post + 0.5*post + rnormal(0,1)
tempfile F
qui save `F'

dis as result "=== L4a: plain bootstrap(standard), the remedy the caution names ==="
use `F', clear
cap noi qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust bootstrap(standard) reps(25) rseed(99)
dis "L4a_RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "L4a attb=" M[1,1] " attse=" M[1,2] " amese=" M[2,2]
}
dis as result "=== L4b: boottstat (studentizes by the analytic SE, which is 0 here) ==="
use `F', clear
cap noi qui wooldid y st year ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) robust bootstrap(standard) reps(25) boottstat rseed(99)
dis "L4b_RC=" _rc
cap {
  mat M = e(wdidmainresults)
  dis "L4b attb=" M[1,1] " attse=" M[1,2] " amese=" M[2,2]
  mat li M
}
dis "DONE_L4"
