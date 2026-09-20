*=====================================================================
* p1_bisect.do -- does si8demand (the maxvar preflight) account for lattice()?
* SELF-CONTAINED: builds its own data, no repo data, no repo generator.
*
* Shape: county panel, i = STATE (the cohort variable), staggered state
* policy onsets, continuous county-level EXPOSURE (Bartik/shift-share
* style, defined on every row and time-invariant).  lattice(all 5) bins
* the exposure.  This is an ordinary continuous-treatment DiD.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-lattice-maxvar"

dis "FLAVOR=" c(flavor) " MAXVAR=" c(maxvar) " MAXKTH=" c(max_k_theory)
which wooldid

*---------------------------------------------------------------- DATA A (small)
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
  label var y "outcome"
end

mkpanel , nstates(25) ntreat(18) ncounty(20) nyear(20) seed(20260916)
desc, short
save "`WD'/p1A.dta", replace

*------------------------------------------------- reconstruct si8demand by hand
* si8cells = distinct (_i,_t) among treated post rows ; _i=state on post rows
* si8ilev  = distinct i ; si8tlev = distinct t ; +200 slack ; x2 under contreat
gen int h = t - ttre
tempvar tg
egen `tg' = tag(state h) if post
qui count if `tg' == 1
local cells = r(N)
qui levelsof state, local(L)
local ilev : word count `L'
qui levelsof t, local(L2)
local tlev : word count `L2'
local si8demand = `cells'*2 + `ilev' + `tlev' + 200
dis "RECON_A cells=" `cells' " ilev=" `ilev' " tlev=" `tlev' " si8demand=" `si8demand'
drop `tg'

*------------------------------------------- structural demand via fvexpand
set emptycells drop
gen double _ii = cond(post, state, 0)
egen long _hh = group(h) if post
qui replace _hh = 0 if missing(_hh)
gen byte _tf = post
qui gquantiles _lat = exposure , xtile nquantiles(5)
tab _lat

qui fvexpand i._ii#i._hh#i._tf
local n_nolat : word count `r(varlist)'
dis "FVX_NOLATTICE_cells_only=" `n_nolat'

qui fvexpand i._ii#i._hh#i._tf##c.exposure
local n_cont : word count `r(varlist)'
dis "FVX_CONTREAT_NOLAT=" `n_cont'

qui fvexpand i._ii#i._hh#i._tf##i._lat##c.exposure
local n_lat : word count `r(varlist)'
dis "FVX_CONTREAT_LAT5=" `n_lat'

qui fvexpand i._ii#i._hh#i._tf##i._lat##c.exposure##c.exposure
local n_lat2 : word count `r(varlist)'
dis "FVX_CONTREAT_LAT5_POLY2=" `n_lat2'

*=====================================================================
* 1. BASELINE: does the spec work at all?  (big maxvar)
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 12000
dis "SET maxvar=" c(maxvar) " maxkth=" c(max_k_theory)
use "`WD'/p1A.dta", clear
timer clear 1
timer on 1
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
local rc12000 = _rc
timer off 1
dis "RUN_MAXVAR12000_RC=" `rc12000'
timer list 1
if `rc12000' == 0 {
  dis "ATT_b=" el(e(wdidmainresults),1,1)
  dis "e_rank=" e(rank) " e_df_m=" e(df_m) " e_N=" e(N)
}

*=====================================================================
* 2. SAME SPEC, NO LATTICE, at Stata/BE's fixed ceiling 2048
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p1A.dta", clear
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
dis "RUN_NOLAT_MAXVAR2048_RC=" _rc

*=====================================================================
* 3. THE PROBE: lattice(all 5) at 2048.  Preflight predicts ~671.
*    NO verbose -- this is exactly what a real user sees.
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
use "`WD'/p1A.dta", clear
dis "=====BEGIN USER-VISIBLE OUTPUT AT MAXVAR 2048====="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
local rcL = _rc
dis "=====END USER-VISIBLE OUTPUT====="
dis "RUN_LAT_MAXVAR2048_RC=" `rcL'
dis "RC907_TEXT_FOLLOWS:"
cap noi error `rcL'
dis "DATASET AFTER FAILURE:"
desc, short

*=====================================================================
* 4. bracket the true demand: climb maxvar until it survives
*=====================================================================
foreach M in 2500 3000 {
  clear all
  adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
  set maxvar `M'
  use "`WD'/p1A.dta", clear
  cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
  dis "BRACKET_MAXVAR`M'_RC=" _rc
}

*=====================================================================
* 5. VALIDATE the reconstruction: a design big enough to TRIP the preflight
*    at 2048, so wooldid prints its own si8demand and si8cells.
*=====================================================================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 2048
mkpanel , nstates(120) ntreat(110) ncounty(4) nyear(20) seed(777)
gen int h = t - ttre
gen byte post2 = (t >= ttre) & ttre < .
tempvar tg2
egen `tg2' = tag(state h) if post2
qui count if `tg2' == 1
local cellsB = r(N)
qui levelsof state, local(LB)
local ilevB : word count `LB'
qui levelsof t, local(LB2)
local tlevB : word count `LB2'
dis "RECON_B cells=" `cellsB' " ilev=" `ilevB' " tlev=" `tlevB' ///
    " predicted_si8demand_contreat=" `cellsB'*2 + `ilevB' + `tlevB' + 200 ///
    " predicted_si8demand_binary=" `cellsB'*1 + `ilevB' + `tlevB' + 200
drop `tg2' post2 h
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
dis "PREFLIGHT_TRIP_RC=" _rc

dis "P1_DONE"
