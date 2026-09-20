* probe15: the full F1 arm set -- every repopts covariate slot that can carry an
* ORIGINAL-cohort-keyed term into a rebuilt-i replication
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

use "dataB.dta", clear
gen byte wave  = (coh <= 2)            // adoption wave: a coarsening of i
gen double cohw = 1 + 0.4*coh          // cohort-keyed weight multiplier
gen double cohx = 10*coh + 0.5*coh*coh // cohort-keyed control
save "dataB3.dta", replace

* observed estimates first: which arms are inert?
foreach sp in PLAIN CTRL FEWAVE ICTRL CCC CCCABS CWM {
  use "dataB3.dta", clear
  local opt
  if "`sp'" == "CTRL"   local opt controls(cohx)
  if "`sp'" == "FEWAVE" local opt fe(wave)
  if "`sp'" == "ICTRL"  local opt interactivecontrols(cohx)
  if "`sp'" == "CCC"    local opt coarsecohortcontrols(c.cohx)
  if "`sp'" == "CCCABS" local opt ccc_absorb(cohx)
  if "`sp'" == "CWM"    local opt customweightmultiplier(cohw)
  cap noi qui wooldid y coh t onset, att cluster(state) fastnose `opt'
  local bb = .
  cap local bb = _b[main_att]
  dis "OBS15 `sp' rc=" _rc " b=" %21x `bb'
}
foreach sp in PLAIN CTRL FEWAVE ICTRL CCC CCCABS CWM {
  use "dataB3.dta", clear
  local opt
  if "`sp'" == "CTRL"   local opt controls(cohx)
  if "`sp'" == "FEWAVE" local opt fe(wave)
  if "`sp'" == "ICTRL"  local opt interactivecontrols(cohx)
  if "`sp'" == "CCC"    local opt coarsecohortcontrols(c.cohx)
  if "`sp'" == "CCCABS" local opt ccc_absorb(cohx)
  if "`sp'" == "CWM"    local opt customweightmultiplier(cohw)
  cap noi qui wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts `opt'
  local rc0 = _rc
  local bb = .
  cap local bb = _b[main_att]
  dis "PERM15 `sp' rc=`rc0' b=" %21x `bb' " p=" e(perm_p_att) " sd=" e(perm_sd_att) " done=" e(perm_completed)
}
