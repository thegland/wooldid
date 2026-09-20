clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260920
set sortseed 20260920
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
do wddump.do

* ============ Design F: gaps, onset-period gap, out-of-window cohort ============
do genF.do
export delimited using "dataF.csv", replace
save "dataF.dta", replace

use "dataF.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(4)
dis "F1_RC=" _rc
wdd_run , file("resF1.csv")

use "dataF.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esprelength(3) espostlength(4)
dis "F2_RC=" _rc
wdd_run , file("resF2.csv")

use "dataF.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects anticipation(1) esprelength(3) espostlength(4)
dis "F3_RC=" _rc
wdd_run , file("resF3.csv")

use "dataF.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects oostreatedcontrols(3) esprelength(3) espostlength(4)
dis "F4_RC=" _rc
wdd_run , file("resF4.csv")

use "dataF.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(4)
dis "F5_RC=" _rc
wdd_run , file("resF5.csv")

* ============ Design G: poisson engine split (is poisexpresults immune?) ============
do genG.do
export delimited using "dataG.csv", replace
save "dataG.dta", replace

use "dataG.dta", clear
cap noi wooldid y i t ttre, att poisson poisexpresults
dis "G1_RC=" _rc
dis "G1_att=" %21x e(b)[1,1]
use "dataG.dta", clear
cap noi wooldid y i t ttre, att poisson poisexpresults customvce(robust)
dis "G2_RC=" _rc
dis "G2_att=" %21x e(b)[1,1]

use "dataG.dta", clear
cap noi wooldid y i t ttre, att poisson poisexpresults semielasticity
dis "G3_RC=" _rc
dis "G3_att=" %21x e(b)[1,1]
use "dataG.dta", clear
cap noi wooldid y i t ttre, att poisson poisexpresults semielasticity customvce(robust)
dis "G4_RC=" _rc
dis "G4_att=" %21x e(b)[1,1]

use "dataG.dta", clear
cap noi wooldid y i t ttre, att poisson
dis "G5_RC=" _rc
dis "G5_att=" %21x e(b)[1,1]
use "dataG.dta", clear
cap noi wooldid y i t ttre, att poisson customvce(robust)
dis "G6_RC=" _rc
dis "G6_att=" %21x e(b)[1,1]

* ============ reach of the semielasticity FE split: subgroups + horizons ============
use "dataE.dta", clear
gen byte sg = mod(i,2)
cap noi wooldid y i t ttre, att semielasticity subgroups(sg) esprelength(2) espostlength(3)
dis "H1_RC=" _rc
wdd_run , file("resH1.csv")
use "dataE.dta", clear
gen byte sg = mod(i,2)
cap noi wooldid y i t ttre, att semielasticity subgroups(sg) customvce(robust) esprelength(2) espostlength(3)
dis "H2_RC=" _rc
wdd_run , file("resH2.csv")

dis "ALLDONE_RUN5"
