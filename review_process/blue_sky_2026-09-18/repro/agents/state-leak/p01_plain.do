clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"

* -- assert we are on the pinned build, not the v1.3 in PLUS
which wooldid

wdgen
wdsnap p01pre

capture noisily wooldid y unit year ttre, cluster(state)
dis "P01_RC=" _rc

wdsnap p01post
wddumpe p01

* -- what a user's next runiform() would be
set seed 20260916
gen double __probe1 = .
replace __probe1 = runiform() in 1
dis "P01_POSTRUN_UNIFORM_after_reseed=" %21x __probe1[1]

dis "P01_DONE"
