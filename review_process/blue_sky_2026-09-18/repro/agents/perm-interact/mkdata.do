* ---- shared data builder for the perm-interact probes ----
clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916
local gseed   20260916
local gN      40
local gT      8
local gt0     2000
local gonsets 2003 2005
local gysd    0.30
local gtsd    0.40
local gaisd   1
local gtaumode hetero
local gdose   cellconst
local gwtype  unit
local gpois   1
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
* a covariate that is strongly predictive so controls()/fe() move the point estimates
gen double zbig = 3*x2 + 2*x1 + rnormal()*0.1
* a coarse grouping crossing cohorts, usable in fe() / coarsecohortcontrols()
gen byte reg5 = region
* a plausible "if" restriction: drop one region
gen byte keepreg = (region != 1)
compress
save "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact/base.dta", replace
