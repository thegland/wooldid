* p10: the stored-estimates surface -- Stata's 300-store limit, store/restore and
* save/use round trips, and what replay does across them
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: estimates store / restore round trip, then replay
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "A_RC=" _rc
dis "A_replayv=" "`e(wdid_replayv)'"
wddumpe est_orig
est store wdsaved
qui regress y post
capture noisily est restore wdsaved
dis "A_RESTORE_RC=" _rc
dis "A_replayv_after_restore=" "`e(wdid_replayv)'"
wddumpe est_restored
capture noisily wooldid
dis "A_REPLAY_AFTER_RESTORE_RC=" _rc

* ============ ARM B: estimates save / use round trip (to disk), then replay
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "B_RC=" _rc
capture noisily estimates save p10est, replace
dis "B_SAVE_RC=" _rc
qui regress y post
capture noisily estimates use p10est
dis "B_USE_RC=" _rc
dis "B_replayv_after_use=" "`e(wdid_replayv)'"
dis "B_cmd_after_use=" "`e(cmd)'"
wddumpe est_used
capture noisily wooldid
dis "B_REPLAY_AFTER_USE_RC=" _rc

* ============ ARM C: Stata's 300-stored-estimates ceiling
wdgen
qui regress y post
forvalues k = 1/295 {
  capture est store zz`k'
}
capture noisily estimates dir
dis "C_STORED=" 295
capture noisily wooldid y unit year ttre, cluster(state) att
dis "C_RC=" _rc
wdline C_post
* push right up to the ceiling
forvalues k = 296/299 {
  capture est store zz`k'
}
capture noisily wooldid y unit year ttre, cluster(state) att
dis "C2_RC=" _rc
wdline C2_post
capture noisily estimates drop _all

* ============ ARM D: what verbose leaves, and whether the NEXT plain run removes it
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose
dis "D_RC=" _rc
capture noisily estimates dir
capture est describe ___wooldidfullmodel
dis "D_FULLMODEL_AFTER_VERBOSE_RC=" _rc
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "D2_RC=" _rc
capture est describe ___wooldidfullmodel
dis "D_FULLMODEL_AFTER_PLAIN_RC=" _rc
capture noisily estimates dir

* ============ ARM E: a user estimate called ___wooldidfullmodel is destroyed
wdgen
qui regress y post
est store ___wooldidfullmodel
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E_RC=" _rc
capture est describe ___wooldidfullmodel
dis "E_USER_MODEL_SURVIVED_RC=" _rc

dis "P10_DONE"
