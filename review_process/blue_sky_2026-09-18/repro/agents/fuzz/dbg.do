clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/assert.do"
use "`D'/ds_canon.dta", clear
wdfz_snap
local n0 `r(nobs)'
local v0 `r(nvar)'
local s0 `r(sort)'
local g0 `r(sig)'
dis "SNAP nobs=`n0' nvar=`v0' sort=`s0' sig=`g0'"
qui wooldid y unit t ttre, att robust
set trace on
set tracedepth 1
capture noisily wdfz_assert, id(dbg) nobs0(`n0') nvar0(`v0') sort0(`s0') sig0(`g0')
set trace off
dis "ARC=" _rc
dis "VIOL=`r(viol)'"
