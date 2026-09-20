*=====================================================================
* p2b_default_maxvar.do -- THE HEADLINE, rerun.
* 50-state county panel, continuous exposure, lattice(all 5), at Stata's
* DEFAULT maxvar (5000).  Order matters: the two wooldid runs happen FIRST
* at the untouched default; only afterwards is maxvar raised to count the
* columns and to confirm the spec is estimable.
* SELF-CONTAINED.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-lattice-maxvar"

dis "DEFAULT_MAXVAR=" c(maxvar) "  MAXKTH=" c(max_k_theory)

clear
set seed 424242
set sortseed 424242
set obs 50
gen long state = _n
gen int ttre = cond(_n <= 45, 2005 + mod(_n-1,8), .)
expand 15
sort state
by state: gen int cnum = _n
gen long county = state*1000 + cnum
gen double exposure = exp(rnormal()*0.7)
gen double cfe      = rnormal()
expand 20
sort county
by county: gen int t = 1999 + _n
gen double yfe = 0.30*sin(t/3)
gen byte post = (t >= ttre) & ttre < .
gen double y = cfe + yfe + 0.25*exposure + 0.80*exposure*post + rnormal()*0.40
desc, short
save "`WD'/p2.dta", replace

gen int h = t - ttre
tempvar tg
egen `tg' = tag(state h) if post
qui count if `tg' == 1
local cells = r(N)
qui levelsof state, local(L)
local ilev : word count `L'
qui levelsof t, local(L2)
local tlev : word count `L2'
dis "P2_RECON cells=" `cells' " ilev=" `ilev' " tlev=" `tlev' " si8demand=" `cells'*2 + `ilev' + `tlev' + 200

* ---------- run 1: no lattice, UNTOUCHED default maxvar ----------
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use "`WD'/p2.dta", clear
dis "MAXVAR NOW=" c(maxvar)
timer clear 1
timer on 1
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
timer off 1
dis "P2_NOLAT_DEFAULT_RC=" _rc
timer list 1

* ---------- run 2: THE PROBE ----------
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use "`WD'/p2.dta", clear
dis "MAXVAR NOW=" c(maxvar)
dis "=====BEGIN USER-VISIBLE OUTPUT (DEFAULT MAXVAR)====="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
dis "=====END USER-VISIBLE OUTPUT====="
dis "P2_LAT5_DEFAULT_RC=" _rc

* ---------- column counts, now that maxvar is allowed to be raised ----------
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 20000
dis "RAISED maxvar=" c(maxvar) " maxkth=" c(max_k_theory)
use "`WD'/p2.dta", clear
set emptycells drop
gen int h = t - ttre
gen double _ii = cond(post, state, 0)
egen long _hh = group(h) if post
qui replace _hh = 0 if missing(_hh)
gen byte _tf = post
qui gquantiles _lat = exposure , xtile nquantiles(5)
cap noi qui fvexpand i._ii#i._hh#i._tf##c.exposure
dis "P2_FVX_NOLAT=" wordcount("`r(varlist)'") " rc=" _rc
cap noi qui fvexpand i._ii#i._hh#i._tf##i._lat##c.exposure
dis "P2_FVX_LAT5=" wordcount("`r(varlist)'") " rc=" _rc

* ---------- confirmation that the spec IS estimable ----------
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 20000
use "`WD'/p2.dta", clear
timer clear 2
timer on 2
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
timer off 2
dis "P2_LAT5_MAXVAR20000_RC=" _rc
timer list 2
if _rc == 0 dis "P2_ATT=" el(e(wdidmainresults),1,1)
dis "P2_DONE"
