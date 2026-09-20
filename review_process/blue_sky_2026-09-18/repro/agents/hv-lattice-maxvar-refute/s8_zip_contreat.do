*=====================================================================
* s8_zip_contreat.do  -- SELF-CONTAINED.  ~4 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s8_zip_contreat.do
*
* Realism anchor #2 for sk-maxvar-fp, on the CONTINUOUS-treatment path and at a
* geography scale ordinary in administrative work (5,000 ZIP codes x 12 years,
* 60 of them exposed, continuous exposure intensity).  Stata's UNTOUCHED DEFAULT
* c(maxvar) = 5000.  No lattice, no polynomial, no interactive controls.
*
* si8demand = cells*2 + ilev + tlev + 200.  ilev = 5000 dominates it, although
* 5,000 absorbed unit levels cost zero design-matrix columns on the reghdfe path.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 500012
set sortseed 500012
which wooldid
dis as text "UNTOUCHED default c(maxvar) = " c(maxvar)

clear
set obs 5000
gen long zip = _n
gen int ttre = cond(_n <= 60, 2012 + mod(_n-1,4), .)
gen double zfe = rnormal()
gen double exposure = exp(rnormal()*0.5)
expand 12
sort zip
by zip: gen int t = 2008 + _n - 1
gen double yfe = 0.2*cos(t/3)
gen byte post = (t >= ttre) & ttre < .
gen double y = zfe + yfe + 0.15*exposure + 0.60*exposure*post + rnormal()*0.5
dis as result "PANEL: 5000 ZIPs x 12 years = " _N " rows; 60 exposed"
tempfile D
qui save `D'

*--- ARM 1: DEFAULT safety, DEFAULT maxvar ---
dis as text "=== ARM 1: att contreat(exposure) contreatcontrols(basic), DEFAULT safety ==="
cap noi wooldid y zip t ttre , att contreat(exposure) contreatcontrols(basic)
dis as result "ARM1_RC=" _rc

*--- ARM 2: identical + safety(off), identical DEFAULT maxvar ---
use "`D'", clear
dis as text "=== ARM 2: identical command + safety(off), same DEFAULT maxvar ==="
timer clear 1
timer on 1
cap noi wooldid y zip t ttre , att contreat(exposure) contreatcontrols(basic) safety(off)
local rc2 = _rc
timer off 1
dis as result "ARM2_RC=" `rc2'
timer list 1
if `rc2' == 0 {
  mat A = e(wdidmainresults)
  dis as result "ARM2 row1 = " A[1,1] "  se = " A[1,2] "  N = " e(N)
}

dis "S8_DONE"
