* r2_pre_and_missing.do
*  (A) the SECONDARY claim: is "the pretreatment att_i row weights cohorts
*      differently from the main att_i row" cwm-specific, or generic to any
*      time-varying weight?
*  (B) MY OWN angle: a customweightmultiplier() with missing values.
*      sthlp:218 says the multiplier "enters the margins weights only and is
*      deliberately not folded into the regression weight".  ado:704 marks the
*      row OUT of the estimation screen when the multiplier is missing.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771102
set sortseed 771102
which wooldid

clear
set obs 14
gen int county = _n
gen int onset = .
replace onset = 2011 if inrange(county,1,3)
replace onset = 2013 if inrange(county,4,6)
replace onset = 2015 if inrange(county,7,8)
gen double pop = 1
replace pop = 40 if inlist(county,1,5)
replace pop = 18 if inlist(county,2,7)
replace pop = 6  if inlist(county,4,10,12)
replace pop = 3  if inlist(county,3,9)
gen double cfe = rnormal(0,1)
gen double teff = 1.2 + 0.055*pop + rnormal(0,0.3)
expand 10
bysort county: gen int year = 2007 + _n
expand 4
bysort county year: gen int estab = _n
gen double sw = exp(rnormal(1.1,0.55))
gen double yfe = 0.12*(year-2007)
gen byte post = (onset<. & year>=onset)
gen double y = 3 + cfe + yfe + teff*post + rnormal(0,0.8)
gen double popg = pop*(1.15^(year-2007))
gen double swt  = sw*(1.15^(year-2007))
* a realistic administrative gap: population missing for a scattered 6% of county-years
gen double u = runiform()
gegen double cellu = min(u), by(county year)
gen double popmiss = pop
replace popmiss = . if cellu < 0.06
gen double popzero = pop
replace popzero = 0 if cellu < 0.06
gen byte gapcell = cellu < 0.06
qui count if gapcell
dis "NGAPROWS=" r(N)
compress
save "panel_r2.dta", replace

capture program drop blockw
program define blockw
  syntax , TAG(string) WVAR(string)
  matrix GTx = e(gteffects)
  preserve
  clear
  svmat double GTx, names(col)
  qui sum `wvar' if relyr >= 0
  local Tpost = r(sum)
  qui sum `wvar' if relyr < 0
  local Tpre = r(sum)
  qui levelsof onset if onset < ., local(ons)
  foreach o of local ons {
    qui sum `wvar' if onset == `o' & relyr >= 0
    local ap = r(sum)
    local sp = r(sum)/`Tpost'
    qui sum `wvar' if onset == `o' & relyr < 0
    local aq = r(sum)
    local sq = r(sum)/`Tpre'
    dis "`tag'|`wvar'| onset=`o' POSTtot=" %12.6f `ap' " PREtot=" %12.6f `aq' /*
      */ " POSTshr=" %10.6f `sp' " PREshr=" %10.6f `sq' " reldif=" %10.6f reldif(`sp',`sq')
  }
  restore
end

local ESOPT esprelength(3) espostlength(3)

*=============== E1: no cwm, constant-in-time sample weight ===============
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_i att_it `ESOPT' gteffects
dis "E1_RC=" _rc
blockw , tag(E1_nocwm_sw) wvar(w_att_i)
blockw , tag(E1_nocwm_sw) wvar(w_att)

*=============== E2: cwm(pop) -- TIME-CONSTANT multiplier ===============
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_i att_it `ESOPT' gteffects customweightmultiplier(pop)
dis "E2_RC=" _rc
blockw , tag(E2_cwmpop) wvar(w_att_i)

*=============== E3: cwm(popg) -- TIME-VARYING multiplier ===============
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_i att_it `ESOPT' gteffects customweightmultiplier(popg)
dis "E3_RC=" _rc
blockw , tag(E3_cwmpopg) wvar(w_att_i)

*=============== E4: NO cwm at all, TIME-VARYING SAMPLE WEIGHT ===============
* the control for the secondary claim: does the att pre row weight the cohorts
* differently from the att main row with no multiplier anywhere in sight?
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=swt], att att_i att_it `ESOPT' gteffects
dis "E4_RC=" _rc
blockw , tag(E4_nocwm_swt) wvar(w_att)
blockw , tag(E4_nocwm_swt) wvar(w_att_i)

*=============== M-arms: missing customweightmultiplier() ===============
* M0 baseline: cwm(pop), nothing missing
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it att_i customweightmultiplier(pop) verbose
dis "M0_RC=" _rc
dis "M0_eN=" e(N)
gen byte s0 = e(sample)
qui count if s0
dis "M0_NESAMP=" r(N)
cap est restore ___wooldidfullmodel
dis "M0_REST_RC=" _rc
dis "M0_REGN=" e(N)

* M1: same run but the multiplier is MISSING on the gap cells
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it att_i customweightmultiplier(popmiss) verbose
dis "M1_RC=" _rc
dis "M1_eN=" e(N)
gen byte s1 = e(sample)
qui count if s1
dis "M1_NESAMP=" r(N)
cap est restore ___wooldidfullmodel
dis "M1_REST_RC=" _rc
dis "M1_REGN=" e(N)

* M2: the same cells carry a multiplier of ZERO instead of missing
use "panel_r2.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it att_i customweightmultiplier(popzero) verbose
dis "M2_RC=" _rc
dis "M2_eN=" e(N)
gen byte s2 = e(sample)
qui count if s2
dis "M2_NESAMP=" r(N)
cap est restore ___wooldidfullmodel
dis "M2_REST_RC=" _rc
dis "M2_REGN=" e(N)

* M3: the gap rows physically deleted, multiplier = pop
use "panel_r2.dta", clear
drop if gapcell
cap noi wooldid y county year onset [aw=sw], att att_it att_i customweightmultiplier(pop) verbose
dis "M3_RC=" _rc
dis "M3_eN=" e(N)

* M4: cwm missing ONLY on never-treated (control) rows -- a plausible user error:
*     "the population multiplier only applies to the treated"
use "panel_r2.dta", clear
gen double poptr = pop
replace poptr = . if onset >= .
cap noi wooldid y county year onset [aw=sw], att att_it customweightmultiplier(poptr)
dis "M4_RC=" _rc

dis "R2_DONE"
