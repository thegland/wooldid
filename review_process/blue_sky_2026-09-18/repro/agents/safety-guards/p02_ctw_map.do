* p02_ctw_map.do -- contreatwithin + contreat() rewrites the WORKING safety value to
* off (ado:4829).  Which guards does that silently disable, and does a number move?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

* exact DGP; i will be the COHORT so that cohort-period cells hold within-cell dose variation
local gN 90
local gT 10
local gt0 2000
local gonsets 2004 2005 2006
local gdose unitconst
local gysd 0
local gaisd 0
local gtrendsd 0
include "`GEN'/gen_core.do"

* a PURE RELABEL of the cohort identifier: 0-indexed ids, the shape a pandas/merge
* export routinely produces.  cohort 1 (a treated cohort) becomes id 0.
gen long cohort0 = cohort - 1
tab cohort cohort0
save "p02.dta", replace

dis "=== A. binary, default safety, 0-indexed i : the guard at ado:4959 ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0)
dis "A_RC=" _rc

dis "=== B. binary, safety(off), 0-indexed i : documented bypass ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0) safety(off)
dis "B_RC=" _rc
dis "B_ATT=" %21x e(wdid_att)
dis "B_N=" e(N)

dis "=== C. binary, default safety, 1-indexed i : the reference answer ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort)
dis "C_RC=" _rc
dis "C_ATT=" %21x e(wdid_att)
dis "C_N=" e(N)

dis "=== D. contreat + contreatwithin, NO safety() typed, 1-indexed i ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatwithin
dis "D_RC=" _rc
dis "D_ATT=" %21x e(wdid_att)
dis "D_AME=" %21x e(wdid_ame)
dis "D_N=" e(N)

dis "=== E. contreat + contreatwithin, NO safety() typed, 0-indexed i ==="
dis "=== the ado:4959 assert is skipped because the WORKING safety value was rewritten ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0) contreat(dose) contreatwithin
dis "E_RC=" _rc
dis "E_ATT=" %21x e(wdid_att)
dis "E_AME=" %21x e(wdid_ame)
dis "E_N=" e(N)

dis "=== F. does an explicit safety(off) reproduce E exactly? (1-indexed) ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatwithin safety(off)
dis "F_RC=" _rc
dis "F_AME=" %21x e(wdid_ame)

dis "=== G. non-integer i under contreatwithin ==="
gen double cohorthalf = cohort/2
capture noisily wooldid ync cohorthalf t ttre, att cluster(cohort) contreat(dose) contreatwithin
dis "G_RC=" _rc

dis "=== H. float 2^24 caution (ado:4973) : with and without contreatwithin ==="
use "p02.dta", clear
* an 8-digit YYYYMMDD calendar code held in float -- the exact shape the caution exists for
gen float tdate  = 20000101 + (t - 2000)*10000
gen float ttredate = 20000101 + (ttre - 2000)*10000
capture noisily wooldid ync cohort tdate ttredate, att cluster(cohort)
dis "H1_RC=" _rc
capture noisily wooldid ync cohort tdate ttredate, att cluster(cohort) contreat(dose) contreatwithin
dis "H2_RC=" _rc
dis "P02_DONE"
