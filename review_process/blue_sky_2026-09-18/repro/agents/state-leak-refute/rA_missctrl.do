* rA_missctrl.do -- SKEPTIC repro for state-leak-01 (missing outcome on all never-treated rows)
* Self-contained: generates its own panel, no project test suite, no external data.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

*---------------------------------------------------------------- own DGP
* 40 "states" x 12 years. 3 treated cohorts of 10, 10 never-treated states.
clear
set obs 40
gen int st = _n
gen int onset = .
replace onset = 2009 if st <= 10
replace onset = 2011 if st > 10 & st <= 20
replace onset = 2013 if st > 20 & st <= 30
* st 31..40 never treated
gen double sfe = rnormal()
expand 12
bys st: gen int yr = 2004 + _n
gen double tfe = 0.03*(yr-2004)
gen byte post = (yr >= onset) & onset < .
gen double y = sfe + tfe + 0.40*post + rnormal()*0.5
gen double wt = 1 + runiform()
label var y "outcome"
save "rA_base.dta", replace

local S1 "OK"

*---------------------------------------------------------------- LEG A1 intact
use "rA_base.dta", clear
cap noi wooldid y st yr onset, cluster(st) att
dis "A1_RC=" _rc
dis "A1_N=" e(N) "  A1_NDROP=" e(N_dropped)
dis "A1_ATT=" %21x e(wdid_att)
dis "A1_SE=" %21x e(wdid_att_se)

*---------------------------------------------------------------- LEG A2 y missing on ALL never-treated rows
use "rA_base.dta", clear
replace y = . if onset >= .
count if y >= .
dis "A2_NMISS=" r(N)
dis "=====A2 BEGIN OUTPUT====="
cap noi wooldid y st yr onset, cluster(st) att
dis "A2_RC=" _rc
dis "=====A2 END OUTPUT====="

*---------------------------------------------------------------- LEG A3 same + safety(off)
use "rA_base.dta", clear
replace y = . if onset >= .
dis "=====A3 BEGIN OUTPUT====="
cap noi wooldid y st yr onset, cluster(st) att safety(off)
dis "A3_RC=" _rc
dis "=====A3 END OUTPUT====="
dis "A3_N=" e(N) "  A3_NDROP=" e(N_dropped)
dis "A3_ATT=" %21x e(wdid_att)
dis "A3_SE=" %21x e(wdid_att_se)
scalar A3att = e(wdid_att)
scalar A3se  = e(wdid_att_se)
scalar A3n   = e(N)

*---------------------------------------------------------------- LEG A4 physically drop never-treated rows, default safety
use "rA_base.dta", clear
drop if onset >= .
dis "=====A4 BEGIN OUTPUT====="
cap noi wooldid y st yr onset, cluster(st) att
dis "A4_RC=" _rc
dis "=====A4 END OUTPUT====="

*---------------------------------------------------------------- LEG A5 physically drop never-treated rows + safety(off)  [the advisor's D2]
use "rA_base.dta", clear
drop if onset >= .
dis "=====A5 BEGIN OUTPUT====="
cap noi wooldid y st yr onset, cluster(st) att safety(off)
dis "A5_RC=" _rc
dis "=====A5 END OUTPUT====="
dis "A5_N=" e(N) "  A5_NDROP=" e(N_dropped)
dis "A5_ATT=" %21x e(wdid_att)
dis "A5_SE=" %21x e(wdid_att_se)
dis "MATCH_A3_A5_ATT=" (e(wdid_att)==A3att)
dis "MATCH_A3_A5_SE="  (e(wdid_att_se)==A3se)
dis "A3_N_vs_A5_N=" A3n " " e(N)

*---------------------------------------------------------------- LEG A6 partial missing on controls (realistic)
use "rA_base.dta", clear
replace y = . if onset >= . & runiform() < 0.30
cap noi wooldid y st yr onset, cluster(st) att
dis "A6_RC=" _rc
dis "A6_N=" e(N) "  A6_NDROP=" e(N_dropped)

*---------------------------------------------------------------- LEG A7 zero weight on never-treated (realistic: pop weight 0)
use "rA_base.dta", clear
replace wt = 0 if onset >= .
dis "=====A7 BEGIN OUTPUT====="
cap noi wooldid y st yr onset [aw=wt], cluster(st) att
dis "A7_RC=" _rc
dis "=====A7 END OUTPUT====="

*---------------------------------------------------------------- LEG A8 i missing on never-treated rows
use "rA_base.dta", clear
gen int st2 = st
replace st2 = . if onset >= .
dis "=====A8 BEGIN OUTPUT====="
cap noi wooldid y st2 yr onset, cluster(st) att
dis "A8_RC=" _rc
dis "=====A8 END OUTPUT====="

dis "ALLDONE"
