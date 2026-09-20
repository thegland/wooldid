clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
use "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/architecture/base.dta", clear
gen long origord = _n
sort tt id
dis "PRE sortedby: `: sortedby'"
dis "PRE nvars=" c(k) " nobs=" c(N)
local rng0 = c(rngstate)
qui est clear
dis "=== RUN 1 plain ==="
cap noi wooldid y id tt ttre , att
dis "STATE r1 rc=" _rc
dis "POST sortedby: `: sortedby'"
dis "POST nvars=" c(k) " nobs=" c(N)
dis "POST rngsame=" ("`c(rngstate)'" == "`rng0'")
dis "POST level=" c(level) " emptycells=" "`c(emptycells)'" " type=" "`c(type)'"
qui est dir
dis "POST estdir: `r(names)'"
qui macro dir
dis "--- macros ---"
macro dir
dis "--- mata ---"
mata: mata describe

dis "=== RUN 2 with permute (RNG) ==="
local rng1 = c(rngstate)
cap noi wooldid y id tt ttre , att cluster(id) permute(shuffle) permreps(25) rseed(12345)
dis "STATE r2 rc=" _rc
dis "POST2 rngsame=" ("`c(rngstate)'" == "`rng1'")
dis "POST2 sortedby: `: sortedby'"
dis "POST2 nvars=" c(k)

dis "=== RUN 3 late failure: saveplots into a bad dir ==="
local rng2 = c(rngstate)
cap noi wooldid y id tt ttre , att makeplots saveplots(/nonexistentdir_zz/plot)
dis "STATE r3 rc=" _rc
dis "POST3 nvars=" c(k) " nobs=" c(N) " sortedby: `: sortedby'"
dis "POST3 ecmd=" "`e(cmd)'"
dis "--- macros after failure ---"
macro dir
dis "=== RUN 4 replay after failure ==="
cap noi wooldid
dis "STATE r4 rc=" _rc

dis "=== RUN 5 parse failure, then replay ==="
cap noi wooldid y id tt ttre , att esprelength(-3)
dis "STATE r5 rc=" _rc
cap noi wooldid
dis "STATE r5replay rc=" _rc
dis "STATEDONE"
