clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local SHAPE canon
include "`GENDIR'/gen_shapes.do"
dis "SHAPEDESC=`SHAPEDESC'"
desc, short
capture noisily wooldid y unit t ttre, att robust
dis "SMOKE_RC=" _rc
mat l e(wdidmainresults)
ereturn list
