*=====================================================================
* s5_family.do  -- SELF-CONTAINED.  ~4 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s5_family.do
*
* MY OWN check of seeded claim hv-lat-02: is the si8demand blind spot a FAMILY
* (interactivecontrols(), coarsecohortcontrols(), ccc_absorb()) rather than a
* lattice special case?  Own data, own shape, maxvar 2048 (Stata/BE's ceiling).
*
* si8demand must be IDENTICAL on every arm (it reads only cells, ilev, tlev,
* and whether contreat() is present), while the realized column count moves.
*=====================================================================
clear all
set more off
set processors 2
set maxvar 2048
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909111
set sortseed 909111
which wooldid
dis as text "c(maxvar) = " c(maxvar)

clear
set obs 30
gen long state = _n
gen int ttre = cond(_n <= 26, 2006 + mod(_n-1,6), .)
gen double sfe = rnormal()
gen byte reg = 1 + mod(_n-1,12)
expand 8
sort state
by state: gen int cnum = _n
gen long county = state*100 + cnum
gen double exposure = exp(rnormal()*0.6)
gen double cfe = rnormal()*0.5
expand 14
sort county
by county: gen int t = 2002 + _n - 1
gen byte post = (t >= ttre) & ttre < .
forval j = 1/8 {
  gen double x`j' = rnormal()
}
gen double y = sfe + cfe + 0.2*exposure + 0.7*exposure*post + 0.1*x1 + 0.1*x2 + rnormal()*0.45
dis as result "PANEL rows=" _N
tempfile D
qui save `D'

*--- predicted si8demand (identical on every arm below) ---
gen int rel = t - ttre
tempvar tg
qui egen byte `tg' = tag(state rel) if post
qui count if `tg'==1
local cells = r(N)
qui levelsof state, local(LS)
local ilev : word count `LS'
qui levelsof t, local(LT)
local tlev : word count `LT'
local pred = `cells'*2 + `ilev' + `tlev' + 200
dis as result "si8demand on EVERY arm below = " `pred' "  (cells=" `cells' ", ilev=" `ilev' ", tlev=" `tlev' ")"
dis as result "c(maxvar)=" c(maxvar) " -> guard silent on every arm? " (`pred' <= c(maxvar))

*--- ARM 1: baseline ---
use "`D'", clear
dis as text "=== ARM 1: att contreat(exposure) contreatcontrols(basic) ==="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
dis as result "ARM1_RC=" _rc

*--- ARM 2: + interactivecontrols(x1..x8), NO lattice ---
use "`D'", clear
dis as text "=== ARM 2: + interactivecontrols(x1 x2 x3 x4 x5 x6 x7 x8) ==="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) interactivecontrols(x1 x2 x3 x4 x5 x6 x7 x8)
dis as result "ARM2_RC=" _rc

*--- ARM 3: + coarsecohortcontrols(i.reg), NO lattice ---
use "`D'", clear
dis as text "=== ARM 3: + coarsecohortcontrols(i.reg) (12-level region) ==="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) coarsecohortcontrols(i.reg)
dis as result "ARM3_RC=" _rc

*--- ARM 4: BINARY mode + interactivecontrols -- no contreat at all ---
* si8demand halves (cells*1), so the guard is even further from the truth.
use "`D'", clear
dis as text "=== ARM 4: BINARY att + interactivecontrols(x1..x8), no contreat ==="
dis as result "si8demand here = " `cells'*1 + `ilev' + `tlev' + 200
cap noi wooldid y state t ttre , att interactivecontrols(x1 x2 x3 x4 x5 x6 x7 x8)
dis as result "ARM4_RC=" _rc

dis "S5_DONE"
