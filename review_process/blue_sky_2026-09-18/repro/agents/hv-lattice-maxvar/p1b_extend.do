*=====================================================================
* p1b_extend.do -- (a) validate my reconstruction of si8demand against
*   wooldid's OWN printed number, (b) narrow the bracket, (c) locate the
*   death point, (d) show contreatpoly() is omitted by the same formula,
*   (e) check P6 (preflight still runs under contreatwithin).
* SELF-CONTAINED.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-lattice-maxvar"

capture program drop mkpanel
program define mkpanel
  syntax , nstates(int) ntreat(int) ncounty(int) nyear(int) seed(int)
  clear
  set seed `seed'
  set sortseed `seed'
  set obs `nstates'
  gen long state = _n
  gen int ttre = cond(_n <= `ntreat', 2005 + mod(_n-1,8), .)
  expand `ncounty'
  sort state
  by state: gen int cnum = _n
  gen long county = state*1000 + cnum
  gen double exposure = exp(rnormal()*0.7)
  gen double cfe      = rnormal()
  expand `nyear'
  sort county
  by county: gen int t = 1999 + _n
  gen double yfe = 0.30*sin(t/3)
  gen byte post = (t >= ttre) & ttre < .
  gen double y = cfe + yfe + 0.25*exposure + 0.80*exposure*post + rnormal()*0.40
end

* ---- build every dataset FIRST (clear all drops user programs) ----
mkpanel , nstates(25) ntreat(18) ncounty(20) nyear(20) seed(20260916)
save "`WD'/p1A.dta", replace
mkpanel , nstates(120) ntreat(110) ncounty(4) nyear(20) seed(777)
save "`WD'/p1B.dta", replace

* recon for B
gen int h = t - ttre
tempvar tg
egen `tg' = tag(state h) if post
qui count if `tg' == 1
local cellsB = r(N)
qui levelsof state, local(LB)
local ilevB : word count `LB'
qui levelsof t, local(LB2)
local tlevB : word count `LB2'
dis "RECON_B cells=" `cellsB' " ilev=" `ilevB' " tlev=" `tlevB' " si8demand_pred=" `cellsB'*2 + `ilevB' + `tlevB' + 200

* fvexpand counts on A for the contreatpoly arms
use "`WD'/p1A.dta", clear
set emptycells drop
gen int h = t - ttre
gen double _ii = cond(post, state, 0)
egen long _hh = group(h) if post
qui replace _hh = 0 if missing(_hh)
gen byte _tf = post
qui gquantiles _lat3 = exposure , xtile nquantiles(3)
qui fvexpand i._ii#i._hh#i._tf##i._lat3##c.exposure
dis "FVX_A_LAT3_POLY1=" wordcount("`r(varlist)'")
qui fvexpand i._ii#i._hh#i._tf##i._lat3##c.exposure##c.exposure##c.exposure
dis "FVX_A_LAT3_POLY3=" wordcount("`r(varlist)'")

*=====================================================================
* (a) TRIP the preflight: design B at maxvar 2048 -> wooldid's own number
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p1B.dta", clear
dis "---- PREFLIGHT TRIP, design B ----"
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
dis "TRIP_RC=" _rc

*=====================================================================
* (e) P6: does the preflight still run under contreatwithin?
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p1B.dta", clear
dis "---- PREFLIGHT TRIP under contreatwithin ----"
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) contreatwithin
dis "TRIP_CTW_RC=" _rc

*=====================================================================
* (d) contreatpoly() is omitted too: lattice(all 3) passes at 2048,
*     lattice(all 3)+contreatpoly(3) does not.  si8demand is 671 for BOTH.
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p1A.dta", clear
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 3)
dis "A_LAT3_POLY1_MAXVAR2048_RC=" _rc

clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p1A.dta", clear
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 3) contreatpoly(3)
dis "A_LAT3_POLY3_MAXVAR2048_RC=" _rc

*=====================================================================
* (c) WHERE does it die?  verbose at 3000 on design A.
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 3000
use "`WD'/p1A.dta", clear
dis "---- VERBOSE DEATH LOCATION ----"
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5) verbose
dis "VERBOSE3000_RC=" _rc

*=====================================================================
* (b) narrow the bracket
*=====================================================================
foreach M in 3500 4000 {
  clear all
  adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
  set maxvar `M'
  use "`WD'/p1A.dta", clear
  cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
  dis "BRACKET2_MAXVAR`M'_RC=" _rc
}
dis "P1B_DONE"
