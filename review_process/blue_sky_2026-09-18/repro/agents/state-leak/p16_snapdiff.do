* p16: the full before/after SNAPSHOT DIFF -- the method's headline deliverable.
* Six arms, each snapshotted immediately before and immediately after one wooldid
* call, so `diff snap_<arm>_pre.txt snap_<arm>_post.txt' is the whole answer.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ---- 1 plain
wdgen
wdsnap s1_pre
capture noisily wooldid y unit year ttre, cluster(state) att
dis "S1_RC=" _rc
wdsnap s1_post

* ---- 2 nofastcollapse (the path with NO preserve/restore around the margins work)
wdgen
wdsnap s2_pre
capture noisily wooldid y unit year ttre, cluster(state) att nofastcollapse
dis "S2_RC=" _rc
wdsnap s2_post

* ---- 3 verbose (documented retentions)
wdgen
wdsnap s3_pre
capture noisily wooldid y unit year ttre, cluster(state) att verbose
dis "S3_RC=" _rc
wdsnap s3_post
capture noisily estimates drop _all

* ---- 4 a failing run
wdgen
qui drop if missing(ttre)
wdsnap s4_pre
capture noisily wooldid y unit year ttre, cluster(state) att
dis "S4_RC=" _rc
wdsnap s4_post

* ---- 5 permute (clear the LASTRC marker first with a good run)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
wdgen
wdsnap s5_pre
capture noisily wooldid y unit year ttre, cluster(state) att permute(shuffle) permreps(20) dots(0) rseed(7)
dis "S5_RC=" _rc
wdsnap s5_post

* ---- 6 bootstrap + plots + saveplots
wdgen
wdsnap s6_pre
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) makeplots saveplots(p16plots) replace bootstrap(standard) reps(20) dots(0) rseed(7)
dis "S6_RC=" _rc
wdsnap s6_post

* ---- 7 poisson (documented _ppmlhdfe_d handling)
wdgen
qui replace y = exp(y/3)
wdsnap s7_pre
capture noisily wooldid y unit year ttre, cluster(state) att poisson
dis "S7_RC=" _rc
wdsnap s7_post

dis "P16_DONE"
