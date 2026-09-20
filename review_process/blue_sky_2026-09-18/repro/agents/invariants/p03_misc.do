clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/invariants"
do "`WD'/wdgen.do"
cap erase "`WD'/p03.txt"
local OPT "att att_i att_it att_itime esprelength(3) espostlength(3) cluster(clus)"

* ---- determinism control: same command twice
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL A1 rc=" _rc
wdsnap using "`WD'/p03.txt", tag(A1) mats(escomp b_es_att)
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL A2 rc=" _rc
wdsnap using "`WD'/p03.txt", tag(A2) mats(escomp b_es_att)

* ---- subgroups with exactly ONE realized level vs no subgroups
use "`WD'/d1.dta", clear
gen byte sg1 = 1
cap noi wooldid y i t ttre , `OPT' subgroups(sg1)
di "SENTINEL SG1 rc=" _rc
wdsnap using "`WD'/p03.txt", tag(SG1) mats(escomp b_es_att)

* ---- subgroup LEVEL RELABEL: {0,1} -> {7,3}
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT' subgroups(sg)
di "SENTINEL SGA rc=" _rc
wdsnap using "`WD'/p03.txt", tag(SGA) mats(escomp)
matrix list e(b), format(%21x)
use "`WD'/d1.dta", clear
gen byte sg2 = cond(sg==0, 7, 3)
cap noi wooldid y i t ttre , `OPT' subgroups(sg2)
di "SENTINEL SGB rc=" _rc
wdsnap using "`WD'/p03.txt", tag(SGB) mats(escomp)
matrix list e(b), format(%21x)

* ---- simultaneousbands: family independence + RNG neutrality
use "`WD'/d1.dta", clear
set seed 4242
cap noi wooldid y i t ttre , att esprelength(3) espostlength(3) cluster(clus) simultaneousbands
di "SENTINEL SB1 rc=" _rc
di "SB1 crit_att = " %21x e(esband_crit_att)
di "SB1 rng_after = " c(rngstate)
use "`WD'/d1.dta", clear
set seed 4242
cap noi wooldid y i t ttre , att att_it esprelength(3) espostlength(3) cluster(clus) simultaneousbands
di "SENTINEL SB2 rc=" _rc
di "SB2 crit_att = " %21x e(esband_crit_att)
di "SB2 crit_att_it = " %21x e(esband_crit_att_it)
di "SB2 rng_after = " c(rngstate)
use "`WD'/d1.dta", clear
set seed 4242
cap noi wooldid y i t ttre , att esprelength(3) espostlength(3) cluster(clus)
di "SENTINEL SB0 rc=" _rc
di "SB0 rng_after = " c(rngstate)

* ---- gteffects present?
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , att esprelength(3) espostlength(3) cluster(clus) gteffects
di "SENTINEL GT rc=" _rc
cap matrix list e(gteffects), format(%12.6g)

di "SENTINEL DONE"
