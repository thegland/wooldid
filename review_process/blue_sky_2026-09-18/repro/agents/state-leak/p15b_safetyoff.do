* p15b: the remaining p15 arms -- (D) globals after a plain verbose run, and
* (E) what safety(off) -- the remedy the r459 message itself offers -- does on
* data whose outcome is missing on every never-treated row.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ---- D: plain verbose run, globals afterwards
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose
dis "D3_RC=" _rc
wdline D3_post
capture noisily estimates drop _all

* ---- E: the all-controls-missing-outcome case
wdgen
qui count if missing(ttre)
dis "E0_NEVERTREATED_ROWS=" r(N)
qui replace y = . if missing(ttre)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E1_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att safety(off)
dis "E2_SAFETYOFF_RC=" _rc
dis "E2_e_N=" e(N)
dis "E2_e_N_dropped=" e(N_dropped)
capture confirm matrix e(wdidmainresults)
if _rc == 0 {
  tempname ME
  matrix `ME' = e(wdidmainresults)
  dis "E2_ATT=" %18.12g `ME'[1,1] "|SE=" %18.12g `ME'[1,2] "|P=" %18.12g `ME'[1,4]
}
* the honest comparison: the SAME data with the outcome intact
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E3_RC=" _rc
dis "E3_e_N=" e(N)
capture confirm matrix e(wdidmainresults)
if _rc == 0 {
  tempname MF
  matrix `MF' = e(wdidmainresults)
  dis "E3_ATT=" %18.12g `MF'[1,1] "|SE=" %18.12g `MF'[1,2] "|P=" %18.12g `MF'[1,4]
}
* and: what a user who instead DROPS the useless rows would get
wdgen
qui replace y = . if missing(ttre)
qui drop if missing(y)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E4_RC=" _rc

dis "P15B_DONE"
