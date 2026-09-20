clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local OUT    "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
* no unit FE AND no within-cohort effect heterogeneity -> exact fit AT i=cohort too,
* so any control option must leave the estimands untouched
local gseed 833001
local gysd 0
local gaisd 0
local gtauhet 0
local gtaumode hetero
local gN 60
local gT 9
local gonsets 2003 2005
local gentry 0.35
local gexit 0.30
local ggap 0.08
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_hom.dta", replace
dis "EX3BUILT hom " _N
