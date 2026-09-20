* p08_si8_confirm.do -- confirm the si8 preflight false positive at Stata's DEFAULT
* maxvar, on an ordinary school-district-sized panel, and pin the diagnosis: the
* + si8ilev + si8tlev terms are demands of the reg/poisson path, not of the reghdfe
* path the default engine uses (absorb(i.i i.t) -- ado:1476).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

dis "DEFAULT_MAXVAR=" c(maxvar)

* 4800 districts, 8 years, 20 ever treated
local gcounts 10 10
local gnevern 4780
local gT 8
local gt0 2000
local gonsets 2004 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui count
dis "NROWS=" r(N)
dis "MAXVAR_NOW=" c(maxvar)
save "p08.dta", replace

dis "=== A default safety at DEFAULT maxvar ==="
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "A_RC=" _rc
dis "=== B safety(off) at the same maxvar ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "B_RC=" _rc
mat MB = e(wdidmainresults)
mat list MB, format(%10.7f)
dis "B_N=" e(N)

dis "=== C diagnosis: the reg path (customvce) DOES materialize i.i -- at maxvar 2048 ==="
use "p07.dta", clear
set maxvar 2048
dis "MAXVAR_C=" c(maxvar)
capture noisily wooldid ync unit t ttre, att customvce(cluster unit) safety(off)
dis "C_RC=" _rc
dis "=== D and the reghdfe path at the same maxvar 2048 ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "D_RC=" _rc
dis "D_N=" e(N)
dis "P08_DONE"
