clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/invariants"
do "`WD'/wdgen.do"
cap erase "`WD'/p01.txt"

* ---------- build one dataset and save it ----------
wdgen_d1
save "`WD'/d1.dta", replace
di "SENTINEL DGP_OK N=" _N

local OPT "att att_i att_it att_itime esprelength(3) espostlength(3) cluster(clus)"

* ---- A: baseline
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL A rc=" _rc
wdsnap using "`WD'/p01.txt", tag(A_base) mats(escomp b_es_att b_es_att_it)

* ---- B: random row permutation
use "`WD'/d1.dta", clear
set seed 771
gen double _u = runiform()
sort _u
drop _u
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL B rc=" _rc
wdsnap using "`WD'/p01.txt", tag(B_roworder) mats(escomp b_es_att b_es_att_it)

* ---- C: relabel i (order-reversing, still positive integers)
use "`WD'/d1.dta", clear
replace i = 1000 - 7*i
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL C rc=" _rc
wdsnap using "`WD'/p01.txt", tag(C_relabel_i) mats(escomp b_es_att b_es_att_it)

* ---- D: duplicate EVERY row
use "`WD'/d1.dta", clear
expand 2
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL D rc=" _rc
wdsnap using "`WD'/p01.txt", tag(D_dup_all) mats(escomp b_es_att b_es_att_it)

* ---- E: affine pure shift of time t -> t + 500
use "`WD'/d1.dta", clear
replace t = t + 500
replace ttre = ttre + 500
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL E rc=" _rc
wdsnap using "`WD'/p01.txt", tag(E_tshift) mats(escomp b_es_att b_es_att_it)

* ---- F: affine scale of time t -> 2t (ttre mapped consistently)
use "`WD'/d1.dta", clear
replace t = 2*t
replace ttre = 2*ttre
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL F rc=" _rc
wdsnap using "`WD'/p01.txt", tag(F_tscale2) mats(escomp b_es_att b_es_att_it)

* ---- G: y -> y + 100
use "`WD'/d1.dta", clear
replace y = y + 100
cap noi wooldid y i t ttre , `OPT'
di "SENTINEL G rc=" _rc
wdsnap using "`WD'/p01.txt", tag(G_yshift) mats(escomp b_es_att b_es_att_it)

* ---- H: aw weights scaled by 1000
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre [aw=w], `OPT'
di "SENTINEL H rc=" _rc
wdsnap using "`WD'/p01.txt", tag(H_aw) mats(escomp b_es_att b_es_att_it)
use "`WD'/d1.dta", clear
replace w = w*1000
cap noi wooldid y i t ttre [aw=w], `OPT'
di "SENTINEL I rc=" _rc
wdsnap using "`WD'/p01.txt", tag(I_aw1000) mats(escomp b_es_att b_es_att_it)

* ---- J: pw vs aw on identical weights
use "`WD'/d1.dta", clear
cap noi wooldid y i t ttre [pw=w], `OPT'
di "SENTINEL J rc=" _rc
wdsnap using "`WD'/p01.txt", tag(J_pw) mats(escomp b_es_att b_es_att_it)

di "SENTINEL DONE"
