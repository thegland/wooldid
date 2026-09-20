* p15: the five state surfaces the earlier probes missed
*  A. verbose + ENGINE: does ___wooldidfullmodel hold the OBSERVED model (the
*     2026-09-05 regression: it used to hold the LAST replication's model)?
*  B. a genuine POST-POSTING failure (saveplots to an existing file, no replace)
*     on the plain path AND on the deferred (bootstrap) path
*  C. level() on the engine path (the shell restores c(level) at ado:70 BEFORE the
*     wiring and display run)
*  D. the WOOLDID_CLEANUP=keep leftover the shell's second wooldid_validate should
*     produce under verbose + an engine (ado:312 / ado:4947)
*  E. safety(off) on the all-controls-missing-outcome data -- does the remedy the
*     error message offers produce a NUMBER?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ A1: verbose + bootstrap, cell-weighted estimand (physical re-id leg,
* whose replications have a DIFFERENT N from the observed run -- so e(N) of the
* retained full model is a clean discriminator)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att_i verbose bootstrap(standard) reps(20) dots(0) rseed(7)
dis "A1_RC=" _rc
wdline A1_post
capture noisily estimates dir
capture est restore ___wooldidfullmodel
dis "A1_RESTORE_RC=" _rc
dis "A1_FULLMODEL_N=" e(N)
dis "A1_FULLMODEL_CMD=" "`e(cmd)'"
capture est describe ___wooldidrepmodel
dis "A1_REPMODEL_STILL_THERE=" (_rc == 0)

* ============ A2: verbose + permute
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose permute(shuffle) permreps(20) dots(0) rseed(7)
dis "A2_RC=" _rc
wdline A2_post
capture noisily estimates dir
capture est restore ___wooldidfullmodel
dis "A2_RESTORE_RC=" _rc
dis "A2_FULLMODEL_N=" e(N)
capture est describe ___wooldidrepmodel
dis "A2_REPMODEL_STILL_THERE=" (_rc == 0)

* ============ A3: reference -- verbose on a plain analytic run
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att_i verbose
dis "A3_RC=" _rc
capture est restore ___wooldidfullmodel
dis "A3_RESTORE_RC=" _rc
dis "A3_FULLMODEL_N=" e(N)
capture noisily estimates drop _all

* ============ B1: a POST-POSTING failure on the plain path.
* saveplots writes <name>_ES_att.gph; a second run to the same name with no
* `replace' must fail at graph save, i.e. AFTER ereturn post and AFTER set level.
capture erase p15dup_ES_att.gph
capture erase p15dup_ES.dta
set level 95
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att level(80) esprelength(2) espostlength(2) makeplots saveplots(p15dup)
dis "B1_RC=" _rc
dis "B1_LEVEL_AFTER=" c(level)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att level(80) esprelength(2) espostlength(2) makeplots saveplots(p15dup)
dis "B2_RC=" _rc
dis "B2_LEVEL_AFTER=" c(level)
wdline B2_post
capture noisily wooldid
dis "B2_REPLAY_RC=" _rc
capture noisily estimates dir
capture noisily frames dir

* ============ B3: the same on the DEFERRED (bootstrap) plot path, ado:259-300
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att level(80) esprelength(2) espostlength(2) makeplots saveplots(p15dup) bootstrap(standard) reps(20) dots(0) rseed(7)
dis "B3_RC=" _rc
dis "B3_LEVEL_AFTER=" c(level)
dis "B3_FRAME_AFTER=" "`c(frame)'"
wdline B3_post
capture est describe ___wooldidfullmodel
dis "B3_FULLMODEL_LEFT=" (_rc == 0)
capture noisily wooldid
dis "B3_REPLAY_RC=" _rc

* clear the marker with a good run
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "B4_RC=" _rc
wdline B4_post

* ============ C: level() on the engine paths
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(7) level(80)
dis "C1_RC=" _rc
dis "C1_e_level=" e(level) "|boot_level=" e(boot_level)
dis "C1_LEVEL_AFTER=" c(level)
capture confirm matrix e(wdidmainresults)
if _rc == 0 dis "C1_COLNAMES=" "`: colnames e(wdidmainresults)'"
if _rc == 0 {
  tempname M
  matrix `M' = e(wdidmainresults)
  dis "C1_LB=" %18.12g `M'[1,5] "|UB=" %18.12g `M'[1,6]
}
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(40) dots(0) rseed(7)
dis "C2_RC=" _rc
dis "C2_e_level=" e(level) "|boot_level=" e(boot_level)
capture confirm matrix e(wdidmainresults)
if _rc == 0 dis "C2_COLNAMES=" "`: colnames e(wdidmainresults)'"
if _rc == 0 {
  tempname M2
  matrix `M2' = e(wdidmainresults)
  dis "C2_LB=" %18.12g `M2'[1,5] "|UB=" %18.12g `M2'[1,6]
}
* and the session-level variant: set level 80, engine run with NO level()
wdgen
set level 80
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(40) dots(0) rseed(7)
dis "C3_RC=" _rc
dis "C3_e_level=" e(level)
dis "C3_LEVEL_AFTER=" c(level)
capture confirm matrix e(wdidmainresults)
if _rc == 0 dis "C3_COLNAMES=" "`: colnames e(wdidmainresults)'"
set level 95

* ============ D: WOOLDID_CLEANUP leftover under verbose + engine
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose bootstrap(standard) reps(20) dots(0) rseed(7)
dis "D1_RC=" _rc
wdline D1_post
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose permute(shuffle) permreps(20) dots(0) rseed(7)
dis "D2_RC=" _rc
wdline D2_post
* a plain verbose run for contrast (no shell-side validate call)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose
dis "D3_RC=" _rc
wdline D3_post

* ============ E: safety(off) on the all-controls-missing-outcome data --
* the remedy the r459 message itself offers
wdgen
qui replace y = . if missing(ttre)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E1_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att safety(off)
dis "E2_SAFETYOFF_RC=" _rc
dis "E2_e_N=" e(N)
capture confirm matrix e(wdidmainresults)
if _rc == 0 {
  tempname ME
  matrix `ME' = e(wdidmainresults)
  dis "E2_ATT=" %18.12g `ME'[1,1] "|SE=" %18.12g `ME'[1,2]
}

dis "P15_DONE"
