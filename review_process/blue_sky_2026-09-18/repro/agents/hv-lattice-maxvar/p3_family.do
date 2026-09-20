*=====================================================================
* p3_family.do -- SELF-FOUND, not seeded.
* Is lattice() the ONLY cell-multiplying term si8demand omits?
* ado:1138  coarsecohortcontrols -> i.`_i'#i.`_t'#(vars)
* ado:1171  interactivecontrols  -> i.`_i'#i.`_t'#1.`targetflag'#c.(vars)
* Neither is in the ado:1371 formula either.  Same design, same
* si8demand (671), no lattice at all.
* SELF-CONTAINED.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-lattice-maxvar"

clear
set seed 20260916
set sortseed 20260916
set obs 25
gen long state = _n
gen int ttre = cond(_n <= 18, 2005 + mod(_n-1,8), .)
expand 20
sort state
by state: gen int cnum = _n
gen long county = state*1000 + cnum
gen double exposure = exp(rnormal()*0.7)
gen double cfe      = rnormal()
gen byte reg = 1 + mod(_n-1, 12)
expand 20
sort county
by county: gen int t = 1999 + _n
gen double yfe = 0.30*sin(t/3)
gen byte post = (t >= ttre) & ttre < .
forval j = 1/8 {
  gen double x`j' = rnormal()
}
gen double y = cfe + yfe + 0.25*exposure + 0.80*exposure*post + 0.1*x1 + rnormal()*0.40
save "`WD'/p3.dta", replace

* si8demand recon (identical to p1A: 213 cells, 25 states, 20 periods)
gen int h = t - ttre
tempvar tg
egen `tg' = tag(state h) if post
qui count if `tg' == 1
local cells = r(N)
dis "P3_RECON cells=" `cells' " si8demand_contreat=" `cells'*2 + 25 + 20 + 200 ///
    " si8demand_binary=" `cells' + 25 + 20 + 200

set emptycells drop
gen double _ii = cond(post, state, 0)
egen long _hh = group(h) if post
qui replace _hh = 0 if missing(_hh)
gen byte _tf = post
qui fvexpand i._ii#i._hh#i._tf##c.exposure
dis "P3_FVX_BASE=" wordcount("`r(varlist)'")
qui fvexpand i._ii#i._hh#i._tf##c.exposure i._ii#i._hh#1._tf#c.(x1 x2 x3 x4 x5 x6 x7 x8) i.t#c.(x1 x2 x3 x4 x5 x6 x7 x8)
dis "P3_FVX_INTERCTRL8=" wordcount("`r(varlist)'")
qui fvexpand i._ii#i._hh#i._tf##c.exposure i._ii#i._hh#(i.reg)
dis "P3_FVX_CCC_FACTOR12=" wordcount("`r(varlist)'")

*--- arm 1: 8 interactive controls, NO lattice, maxvar 2048 -------------------
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p3.dta", clear
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) interactivecontrols(x1 x2 x3 x4 x5 x6 x7 x8)
dis "P3_INTERCTRL8_MAXVAR2048_RC=" _rc

*--- arm 2: coarsecohortcontrols with a 12-level factor, NO lattice ------------
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p3.dta", clear
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) coarsecohortcontrols(i.reg)
dis "P3_CCCFACTOR12_MAXVAR2048_RC=" _rc

*--- arm 3: control -- the same run with neither, to prove the base passes -----
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p3.dta", clear
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
dis "P3_BASE_MAXVAR2048_RC=" _rc
dis "P3_DONE"
