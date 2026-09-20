* p08: frames, plots, files in the cwd, replay-with-no-data, and the `if' surface
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: run from a NON-DEFAULT frame, success path
frame create work
frame change work
wdgen
dis "A_FRAME_BEFORE=" "`c(frame)'"
capture noisily wooldid y unit year ttre, cluster(state) att
dis "A_RC=" _rc
dis "A_FRAME_AFTER=" "`c(frame)'"
capture noisily frames dir

* ============ ARM B: non-default frame, FAILURE path
wdgen
qui drop if missing(ttre)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "B_RC=" _rc
dis "B_FRAME_AFTER=" "`c(frame)'"
capture noisily frames dir
frame change default
capture noisily frame drop work

* ============ ARM C: plots -- graphs left in memory, files left in the cwd
dis "C_FILES_BEFORE:"
capture noisily dir "p08*"
capture noisily graph drop _all
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) makeplots saveplots(p08plots) replace
dis "C_RC=" _rc
dis "C_GRAPHS_IN_MEMORY:"
capture noisily graph dir
dis "C_FILES_AFTER:"
capture noisily dir "p08*"
capture noisily frames dir
wdline C_post

* ============ ARM D: replay with NO dataset in memory (sthlp:1459 claims this works)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att subgroups(sg) esprelength(2) espostlength(2)
dis "D_RC=" _rc
clear
dis "D_NOBS=" c(N)
capture noisily wooldid
dis "D_REPLAY_RC=" _rc

* ============ ARM E: replay after the DATA in memory changed underneath
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E_RC=" _rc
qui replace y = y + 100
capture noisily wooldid
dis "E_REPLAY_RC=" _rc

* ============ ARM F: replay when the last e-class command was NOT wooldid
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
qui regress y post
capture noisily wooldid
dis "F_REPLAY_AFTER_REGRESS_RC=" _rc

* ============ ARM G: the `if' surface -- is the expression evaluated on the
* user's row order, or after wooldid has re-sorted the data?
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if keepflag == 1, cluster(state) att
dis "G1_RC=" _rc
dis "G1_N=" e(N)
wddumpe if_flagvar
wdgen
qui gen byte keepflag = (_n <= 360)
capture noisily wooldid y unit year ttre if _n <= 360, cluster(state) att
dis "G2_RC=" _rc
dis "G2_N=" e(N)
wddumpe if_underscore_n

* ============ ARM H: a NONDETERMINISTIC `if' -- estimation sample vs. margins subpop
wdgen
set seed 31415
capture noisily wooldid y unit year ttre if runiform() < 0.7, cluster(state) att
dis "H_RC=" _rc
dis "H_N=" e(N)
dis "H_N_dropped=" e(N_dropped)
wddumpe if_runiform

* ============ ARM I: does a plain run change the user's variable list / types?
wdgen
wdsnap I_pre
capture noisily wooldid y unit year ttre, cluster(state) att
dis "I_RC=" _rc
wdsnap I_post

dis "P08_DONE"
