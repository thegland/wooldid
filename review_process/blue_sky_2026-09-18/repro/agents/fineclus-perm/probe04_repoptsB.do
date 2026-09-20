* probe04: repopts contamination on a RAGGED panel where the options actually bite
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* ---- first: do these options move the OBSERVED estimate at all? (instrument check) ----
foreach sp in "NONE" "FESTATE" "FECOHCOPY" "CTRLCOH" "CTRLSTATE" {
  use "dataB.dta", clear
  local opt
  if "`sp'" == "FESTATE"   local opt fe(state)
  if "`sp'" == "FECOHCOPY" local opt fe(cohcopy)
  if "`sp'" == "CTRLCOH"   local opt controls(cohctrl)
  if "`sp'" == "CTRLSTATE" local opt controls(statectrl)
  cap noi qui wooldid y coh t onset, att cluster(state) fastnose `opt'
  dis "INSTR `sp' rc=" _rc " b=" %21x _b[main_att]
}

* ---- now the permutation runs ----
foreach sp in "NONE" "FESTATE" "FECOHCOPY" "CTRLCOH" "CTRLSTATE" {
  use "dataB.dta", clear
  local opt
  if "`sp'" == "FESTATE"   local opt fe(state)
  if "`sp'" == "FECOHCOPY" local opt fe(cohcopy)
  if "`sp'" == "CTRLCOH"   local opt controls(cohctrl)
  if "`sp'" == "CTRLSTATE" local opt controls(statectrl)
  cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(B_`sp') replace `opt'
  dis "PERM `sp' rc=" _rc " M=" e(perm_npossible) " done=" e(perm_completed) " fail=" e(perm_failed) " b=" %21x _b[main_att] " p=" e(perm_p_att) " sd=" e(perm_sd_att)
}

* ---- side-by-side ----
use "B_NONE.dta", clear
rename main_att b_none
foreach sp in FESTATE FECOHCOPY CTRLCOH CTRLSTATE {
  merge 1:1 _rep using "B_`sp'.dta", nogen
  rename main_att b_`sp'
}
gen double d_FESTATE   = b_FESTATE   - b_none
gen double d_FECOHCOPY = b_FECOHCOPY - b_none
gen double d_CTRLCOH   = b_CTRLCOH   - b_none
gen double d_CTRLSTATE = b_CTRLSTATE - b_none
list _rep b_none b_FECOHCOPY d_FECOHCOPY d_FESTATE d_CTRLCOH d_CTRLSTATE in 1/16, sep(0)
foreach v in FESTATE FECOHCOPY CTRLCOH CTRLSTATE {
  qui sum d_`v' if _rep > 0
  dis "DELTA `v': maxabs over reps = " %12.4e max(abs(r(min)),abs(r(max)))
  qui sum d_`v' if _rep == 0
  dis "DELTA `v': at rep 0 (observed) = " %12.4e r(mean)
}
