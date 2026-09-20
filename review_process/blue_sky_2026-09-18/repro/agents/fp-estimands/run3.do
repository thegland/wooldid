clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260918
set sortseed 20260918
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
cd "`WD'"
do wddump.do

* =====================================================================
* PROBE 1: semielasticity -- which prediction is the denominator?
* =====================================================================
do genE.do
export delimited using "dataE.csv", replace
save "dataE.dta", replace

* E0: plain level effects, reghdfe path (reference point)
use "dataE.dta", clear
cap noi wooldid y i t ttre, att gteffects esprelength(3) espostlength(4)
dis "E0_RC=" _rc
wdd_run , file("resE0.csv")

* E1: semielasticity, DEFAULT engine (reghdfe, FE absorbed)
use "dataE.dta", clear
cap noi wooldid y i t ttre, att semielasticity gteffects esprelength(3) espostlength(4)
dis "E1_RC=" _rc
wdd_run , file("resE1.csv")

* E2: semielasticity, EXPLICIT-DUMMY engine (customvce -> regress with i.i i.t)
use "dataE.dta", clear
cap noi wooldid y i t ttre, att semielasticity customvce(robust) gteffects esprelength(3) espostlength(4)
dis "E2_RC=" _rc
wdd_run , file("resE2.csv")

* E3: level effects on the explicit-dummy engine (must equal E0 point estimates)
use "dataE.dta", clear
cap noi wooldid y i t ttre, att customvce(robust) gteffects esprelength(3) espostlength(4)
dis "E3_RC=" _rc
wdd_run , file("resE3.csv")

* E4: hand-built check.  Fit the same model with explicit dummies and ask margins
*     for eydx directly, then again after reghdfe, and print both.
use "dataE.dta", clear
qui gen byte trt = (t >= ttre) & ttre < .
qui gen double relyr = t - ttre
qui egen long cell = group(i relyr) if trt == 1
qui replace cell = 0 if cell == .
qui reg y i.i i.t i.cell
qui margins , dydx(cell) subpop(if trt==1) noestimcheck nose
dis "HAND_REG_ok"
* per-row fitted value with and without the absorbed FE
qui predict double yhat_full, xb
qui reghdfe y i.cell, absorb(i.i i.t) tol(1e-12)
qui predict double yhat_hdfe, xb
qui sum yhat_full if trt==1
dis "YHATFULL_mean=" %21x r(mean)
qui sum yhat_hdfe if trt==1
dis "YHATHDFE_mean=" %21x r(mean)
qui corr yhat_full yhat_hdfe
dis "CORR_full_hdfe=" r(rho)

* =====================================================================
* PROBE 2: the no-pre-period unit drop -- is anything said, is it counted?
* =====================================================================
use "dataC.dta", clear
count
dis "C_ROWS=" r(N)
qui levelsof i, local(allc)
dis "C_UNITS=" wordcount("`allc'")
cap noi wooldid y i t ttre, att att_i gteffects esprelength(3) espostlength(5)
dis "NDROP_RC=" _rc
dis "NDROP_N=" e(N)
dis "NDROP_NDROPPED=" e(N_dropped)

* =====================================================================
* PROBE 3: duplication invariance.  Physically duplicating every row of ONE
* cohort must leave all four estimands and every horizon unchanged.
* =====================================================================
use "dataA.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "DUP0_RC=" _rc
wdd_run , file("resDUP0.csv")

use "dataA.dta", clear
expand 2 if i <= 8
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "DUP1_RC=" _rc
wdd_run , file("resDUP1.csv")

* and with sample weights halved on the duplicated rows, which must ALSO be a no-op
use "dataA.dta", clear
gen byte orig = 1
expand 2 if i <= 8
bysort i t orig: gen byte dupno = _n
replace w = w/2 if i <= 8
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "DUP2_RC=" _rc
wdd_run , file("resDUP2.csv")

use "dataA.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "DUP3_RC=" _rc
wdd_run , file("resDUP3.csv")

dis "ALLDONE_RUN3"
