* p03: what a FAILED wooldid leaves behind, and whether it changes the NEXT run
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: parse-stage failure
wdgen
wdline A_pre
capture noisily wooldid y unit year ttre, cluster(state) thisisnotanoption
dis "A_RC=" _rc
wdline A_post
* replay must refuse
capture noisily wooldid
dis "A_REPLAY_RC=" _rc

* ============ ARM B: estimation-stage failure (no never-treated), with level(80)
wdgen
qui drop if missing(ttre)
set level 95
wdline B_pre
capture noisily wooldid y unit year ttre, cluster(state) level(80)
dis "B_RC=" _rc
wdline B_post
dis "B_LEVEL_RESTORED=" (c(level) == 95)

* ============ ARM C: after those two failures, does a GOOD run reproduce?
wdgen
capture noisily wooldid y unit year ttre, cluster(state)
dis "C_RC=" _rc
wddumpe fail_S
wdline C_post

* ============ ARM D: poisson failure path -- is _ppmlhdfe_d left in the data?
wdgen
qui replace y = exp(y/3)
capture noisily wooldid y unit year ttre, cluster(state) poisson
dis "D_RC=" _rc
capture confirm variable _ppmlhdfe_d
dis "D_PPMLHDFE_D_LEFT=" (_rc == 0)
wdline D_post

* ============ ARM D2: a poisson run that fails INSIDE ppmlhdfe, before wooldid
* registers the _ppmlhdfe_d cleanup obligation (ado:1490-1493)
wdgen
qui replace y = exp(y/3)
qui replace y = -1 in 1
capture noisily wooldid y unit year ttre, cluster(state) poisson
dis "D2_RC=" _rc
capture confirm variable _ppmlhdfe_d
dis "D2_PPMLHDFE_D_LEFT=" (_rc == 0)
wdline D2_post

* ============ ARM D3: poisson with perfect separation (documented r322)
wdgen
qui replace y = exp(y/3)
qui replace y = 0 if !missing(ttre) & ttre == 2005 & year == 2005
capture noisily wooldid y unit year ttre, cluster(state) poisson
dis "D3_RC=" _rc
capture confirm variable _ppmlhdfe_d
dis "D3_PPMLHDFE_D_LEFT=" (_rc == 0)
wdline D3_post

* ============ ARM E: a post-posting failure -- pretrendpower on a run with no pre window
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(0) espostlength(2) pretrendpower(0.5)
dis "E_RC=" _rc
wdline E_post
capture noisily wooldid
dis "E_REPLAY_RC=" _rc

* ============ ARM F: repeated failures inside capture (the loop a user actually writes)
wdgen
qui drop if missing(ttre)
forvalues rr = 1/3 {
  capture wooldid y unit year ttre, cluster(state)
  dis "F_ITER`rr'_RC=" _rc
  wdline F_iter`rr'
}

* ============ ARM G: after everything, a clean run again
wdgen
capture noisily wooldid y unit year ttre, cluster(state)
dis "G_RC=" _rc
wddumpe fail_S2
wdline G_post

dis "P03_DONE"
