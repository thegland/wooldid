clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid
include "`WD'/gen_rare.do"
dis "NOBS=" _N
save rare20, replace

*==================================================================
* R1: long ES window -- post5 is identified by state 1 ALONE, so it is
*     missing from every replication that does not draw state 1.  The band
*     is LISTWISE over horizons; the pointwise percentile interval is not.
*==================================================================
dis "=== R1 bootstrap + band, rare long horizon ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) rseed(2026) att esprelength(2) espostlength(5) simultaneousbands dots(0)
dis "R1_RC=" _rc
dis "R1_ENGINE=" "`e(esband_engine)'"
dis "R1_FAMS=" "`e(esband_families)'"
dis "R1_CRIT=" e(esband_crit_att)
dis "R1_BANDN=" e(esband_n_att) "  COMPLETED=" e(boot_completed)
mat NAV = e(boot_navail)
mat list NAV
local BFAM att
local BTAG R1
include "`WD'/bandcmp.do"

*==================================================================
* R2: same, boottstat -> studentized band
*==================================================================
dis "=== R2 boottstat + band ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) rseed(2026) boottstat att esprelength(2) espostlength(5) simultaneousbands dots(0)
dis "R2_RC=" _rc
dis "R2_ENGINE=" "`e(esband_engine)'"
dis "R2_FAMS=" "`e(esband_families)'"
dis "R2_CRIT=" e(esband_crit_att)
dis "R2_BANDN=" e(esband_n_att) "  COMPLETED=" e(boot_completed)
mat list e(boot_navail)
local BFAM att
local BTAG R2
include "`WD'/bandcmp.do"
dis "P06_DONE"
