clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

dis "STATA_FLAVOR=" c(flavor)
dis "MAXVAR_DEFAULT=" c(maxvar)
dis "MAX_K_THEORY=" c(max_k_theory)
which wooldid

local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
* small realistic continuous-treatment panel: 12 cohort groups, T=12, thick cells
local gseed 20260916
local gN    240
local gT    12
local gt0   2000
local gonsets "2004 2005 2006"
local gdose skew
local gysd  0.30
include "`GENDIR'/gen_core.do"

desc, short
sum dose, detail
tab cohort t if _n<1

* ---- run a small lattice spec verbosely so we can see the regression string
cap noi wooldid y cohort t ttre , att contreat(dose) contreatcontrols(basic) lattice(all 3) verbose
dis "P0_LATTICE_RC=" _rc
