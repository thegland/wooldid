* probe12: F1 in its most realistic clothing -- an "adoption wave" fixed effect,
* plus: does anything in the OBSERVED output warn the user the term was dropped?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

use "dataB.dta", clear
gen byte wave = (coh <= 2)     // "early vs late/never adopting wave": a COARSENING of i
gen byte region = cond(inlist(state,1,3,5),1,2)   // genuinely cross-cutting: NOT collinear
save "dataB2.dta", replace

* --- observed run, shown in full: is there any signal that fe(wave) was dropped?
use "dataB2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fastnose
dis "OBS_PLAIN_B=" %21x _b[main_att]
use "dataB2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fastnose fe(wave)
dis "OBS_FEWAVE_B=" %21x _b[main_att]
use "dataB2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fastnose fe(coh)
dis "OBS_FECOH_B=" %21x _b[main_att]
use "dataB2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fastnose fe(region)
dis "OBS_FEREGION_B=" %21x _b[main_att]

* --- permutation arms
foreach sp in PLAIN FEWAVE FECOH FEREGION {
  use "dataB2.dta", clear
  local opt
  if "`sp'" == "FEWAVE"   local opt fe(wave)
  if "`sp'" == "FECOH"    local opt fe(coh)
  if "`sp'" == "FEREGION" local opt fe(region)
  cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts `opt'
  dis "PERM12 `sp' rc=" _rc " b=" %21x _b[main_att] " p=" e(perm_p_att) " sd=" e(perm_sd_att) " lb=" e(perm_lb_att) " ub=" e(perm_ub_att)
}
