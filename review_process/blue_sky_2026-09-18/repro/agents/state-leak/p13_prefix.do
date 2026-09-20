* p13: Stata's own replication PREFIXES around wooldid -- the documented
* bootstrap-prefix route (sthlp:bootstrap entry) and what it leaves behind.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: the documented prefix call
wdgen
wdline A_pre
capture noisily bootstrap _b, reps(30) seed(42) cluster(state) idcluster(nc): wooldid y unit year ttre, att fastnose
dis "A_RC=" _rc
dis "A_e_cmd=" "`e(cmd)'"
dis "A_e_N_reps=" e(N_reps)
wdline A_post
capture noisily estimates dir
capture noisily frames dir

* ============ ARM B: session state after the prefix -- a plain run must still
* give exactly the reference numbers
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "B_RC=" _rc
wddumpe prefix_then_plain

* ============ ARM C: svy: must fail loudly (documented r101), and leave nothing
wdgen
qui gen double sw = 1
capture noisily svyset unit [pw=sw]
dis "C_SVYSET_RC=" _rc
capture noisily svy: wooldid y unit year ttre, att
dis "C_SVY_RC=" _rc
wdline C_post

* ============ ARM D: after all that, the reference spec again
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "D_RC=" _rc
wddumpe prefix_final
wdline D_post

* ============ ARM E: statsby, another prefix a user would reach for
wdgen
capture noisily statsby _b, by(sg) clear: wooldid y unit year ttre, att
dis "E_RC=" _rc
dis "E_NOBS=" c(N)

dis "P13_DONE"
