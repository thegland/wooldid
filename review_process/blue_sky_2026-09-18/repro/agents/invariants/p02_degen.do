clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/invariants"
do "`WD'/wdgen.do"
cap erase "`WD'/p02.txt"
local OPT "att att_i att_it att_itime esprelength(3) espostlength(3) cluster(clus)"

* ---- K/L: default collapse vs nofastcollapse (binary, with a within-cell control)
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT' controls(x)
di "SENTINEL K rc=" _rc
wdsnap using "`WD'/p02.txt", tag(K_fcol) mats(escomp b_es_att)
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT' controls(x) nofastcollapse
di "SENTINEL L rc=" _rc
wdsnap using "`WD'/p02.txt", tag(L_nofcol) mats(escomp b_es_att)

* ---- M/N: subgroups with ONE realized level vs no subgroups
use "`WD'/d1.dta", clear
gen byte sg1 = 1
cap noi wooldid y i t ttre , `OPT' subgroups(sg1)
di "SENTINEL M rc=" _rc
wdsnap using "`WD'/p02.txt", tag(M_sg1) mats(escomp b_es_att)

* ---- O/P: anticipation(2) vs hand-shifted ttre
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT' anticipation(2)
di "SENTINEL O rc=" _rc
wdsnap using "`WD'/p02.txt", tag(O_antic2) mats(escomp b_es_att)
use "`WD'/d1.dta", clear
replace ttre = ttre - 2
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL P rc=" _rc
wdsnap using "`WD'/p02.txt", tag(P_handshift) mats(escomp b_es_att)

* ---- Q/R: if-qualifier vs physically dropping the rows
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre if t <= 2006 , `OPT'
di "SENTINEL Q rc=" _rc
wdsnap using "`WD'/p02.txt", tag(Q_if) mats(escomp b_es_att)
use "`WD'/d1.dta", clear
keep if t <= 2006
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL R rc=" _rc
wdsnap using "`WD'/p02.txt", tag(R_drop) mats(escomp b_es_att)

* ---- S/T: t -> 2t with NO event-study window (should be exactly invariant)
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , att att_i att_it att_itime cluster(clus)
di "SENTINEL S rc=" _rc
wdsnap using "`WD'/p02.txt", tag(S_noes_base)
use "`WD'/d1.dta", clear
replace t = 2*t
replace ttre = 2*ttre
cap noi wooldid y i t ttre , att att_i att_it att_itime cluster(clus)
di "SENTINEL T rc=" _rc
wdsnap using "`WD'/p02.txt", tag(T_noes_x2)

* ---- U/V: expand-one-cohort vs aweight doubling of that cohort
use "`WD'/d1.dta", clear
expand 2 if i <= 4
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL U rc=" _rc
wdsnap using "`WD'/p02.txt", tag(U_expand4) mats(escomp b_es_att)
use "`WD'/d1.dta", clear
gen double ww = cond(i<=4, 2, 1)
cap noi wooldid y i t ttre [aw=ww] , `OPT'
di "SENTINEL V rc=" _rc
wdsnap using "`WD'/p02.txt", tag(V_aw4) mats(escomp b_es_att)

* ---- W/X: aw scaled by an exact power of two (must be bit-identical)
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre [aw=w] , `OPT'
di "SENTINEL W rc=" _rc
wdsnap using "`WD'/p02.txt", tag(W_aw) mats(escomp b_es_att)
use "`WD'/d1.dta", clear
replace w = w*1024
cap noi wooldid y i t ttre [aw=w] , `OPT'
di "SENTINEL X rc=" _rc
wdsnap using "`WD'/p02.txt", tag(X_aw1024) mats(escomp b_es_att)

di "SENTINEL DONE"
