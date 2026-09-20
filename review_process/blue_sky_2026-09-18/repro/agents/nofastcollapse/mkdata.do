clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"

*=== dC : i = cohort (coarse id). many units per cohort => fat cells. realistic. ===
local gN 150
local gT 10
local gonsets 2004 2006 2008
local gnever 0.35
local gysd 0.30
local gdose unitconst
local gwtype unit
local gwspread 2
local gpois 1
include "`GD'/gen_core.do"
gen byte sgc   = mod(cohort,2)
gen long clnest = ceil(cohort/2)
label var sgc "cohort-constant subgroup"
label var clnest "cluster NESTING i(=cohort)"
gen double cwm = 1 + 0.5*mod(unit,3)
gen double cwm0 = cond(mod(unit,10)==0, 0, 1 + 0.5*mod(unit,3))
gen double xcohort = a_i
qui bysort cohort t: replace xcohort = xcohort[1]
save "`DD'/dC.dta", replace
tab cohort ttre, missing

*=== dU : i = unit, repeated cross-sections (4 rows per unit-period) ===
clear
local gseed 771103
local gN 30
local gT 9
local gonsets 2003 2005 2007
local gnever 0.35
local gysd 0.30
local gpercell 4
local gsub 2
local gdose unitconst
local gwtype row
local gwspread 2
local gpois 1
include "`GD'/gen_core.do"
gen byte sgc = mod(unit,2)
gen long clnest = region
gen double cwm = 1 + 0.5*mod(_n,3)
gen double cwm0 = cond(mod(_n,10)==0, 0, 1 + 0.5*mod(_n,3))
save "`DD'/dU.dta", replace
qui count
dis "dU N = " r(N)

*=== dCm : dC shape but with messy missings / unbalance / gaps ===
clear
local gseed 5150
local gN 150
local gT 11
local gonsets 2004 2006 2008
local gnever 0.35
local gysd 0.30
local gdose unitconst
local gwtype unit
local gwspread 2
local gentry 0.15
local gexit 0.15
local ggap 0.08
local gymiss 0.05
local gxmiss 0.05
local gpois 1
include "`GD'/gen_core.do"
gen byte sgc = mod(cohort,2)
gen long clnest = ceil(cohort/2)
gen double cwm0 = cond(mod(unit,10)==0, 0, 1 + 0.5*mod(unit,3))
save "`DD'/dCm.dta", replace
qui count
dis "dCm N = " r(N)
dis "MKDATA_DONE"
