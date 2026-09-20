*=====================================================================
* s7_county_scale.do  -- SELF-CONTAINED.  ~4 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s7_county_scale.do
*
* Realism anchor for the false-positive finding (sk-maxvar-fp): the refusal is
* reachable at US-COUNTY scale (3,143 counties) with an ordinary continuous-dose
* county policy, at Stata's UNTOUCHED DEFAULT c(maxvar) = 5000, with no lattice,
* no polynomial, no interactive controls -- nothing exotic at all.
*
* 3,143 counties x 18 years, 100 of them treated with staggered onsets and a
* continuous exposure intensity.  i = county (unit-level cohorts: the canonical
* staggered design).
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31431
set sortseed 31431
which wooldid
dis as text "UNTOUCHED default c(maxvar) = " c(maxvar)

clear
set obs 3143
gen long county = _n
gen int ttre = cond(_n <= 100, 2008 + mod(_n-1,5), .)
gen double cfe = rnormal()
gen double exposure = exp(rnormal()*0.5)
expand 18
sort county
by county: gen int t = 2004 + _n - 1
gen double yfe = 0.2*cos(t/3)
gen byte post = (t >= ttre) & ttre < .
gen double y = cfe + yfe + 0.15*exposure + 0.60*exposure*post + rnormal()*0.5
dis as result "PANEL: 3143 counties x 18 years = " _N " rows; 100 treated"
tempfile D
qui save `D'

*--- ARM 1: continuous dose, DEFAULT safety, DEFAULT maxvar ---
dis as text "=== ARM 1: att contreat(exposure) contreatcontrols(basic), DEFAULT safety ==="
cap noi wooldid y county t ttre , att contreat(exposure) contreatcontrols(basic)
dis as result "ARM1_RC=" _rc

*--- ARM 2: identical, safety(off), identical DEFAULT maxvar ---
use "`D'", clear
dis as text "=== ARM 2: identical command + safety(off), same DEFAULT maxvar ==="
timer clear 1
timer on 1
cap noi wooldid y county t ttre , att contreat(exposure) contreatcontrols(basic) safety(off)
local rc2 = _rc
timer off 1
dis as result "ARM2_RC=" `rc2'
timer list 1
if `rc2' == 0 {
  mat A = e(wdidmainresults)
  dis as result "ARM2 ame/att row 1 = " A[1,1] "  se = " A[1,2] "  N = " e(N)
}

*--- ARM 3: binary version, DEFAULT safety (si8demand halves) ---
use "`D'", clear
dis as text "=== ARM 3: plain binary att, DEFAULT safety, DEFAULT maxvar ==="
cap noi wooldid y county t ttre , att
dis as result "ARM3_RC=" _rc

dis "S7_DONE"
