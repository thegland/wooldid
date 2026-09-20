* p10 -- pin the si8 diagnosis as ENGINE-CONDITIONAL.  Same shape, maxvar 2048,
* 2300 units: the reg path (customvce) should hit a genuine r(907); the default
* reghdfe path should run, because absorb(i.i i.t) never enters the design matrix.
clear all
set more off
set processors 2
set maxvar 2048
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local gcounts 10 10
local gnevern 2280
local gT 8
local gt0 2000
local gonsets 2004 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui levelsof unit, local(uu)
dis "NUNITS=" wordcount("`uu'")
dis "MAXVAR=" c(maxvar)
dis "=== default safety (the preflight) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "PRE_RC=" _rc
dis "=== reghdfe path, safety(off) ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "REGHDFE_RC=" _rc
dis "REGHDFE_N=" e(N)
dis "=== reg path (customvce), safety(off) ==="
capture noisily wooldid ync unit t ttre, att customvce(cluster unit) safety(off)
dis "REGPATH_RC=" _rc
dis "P10_DONE"
