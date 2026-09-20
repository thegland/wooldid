clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local SHAPE canon
include "`GENDIR'/gen_shapes.do"
dis "SHAPEDESC=`SHAPEDESC'"
desc
list in 1/5
cap noisily wooldid y i t ttre, att
dis "C_RC=" _rc
