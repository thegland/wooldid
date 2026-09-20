clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
cd "`WD'"
do wddump.do

* ================= Design B: multi-row-per-cell admin shape =================
do genB.do
export delimited using "dataB.csv", replace
save "dataB.dta", replace
count
dis "B_NOBS=" r(N)

use "dataB.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esprelength(2) espostlength(4)
dis "B1_RC=" _rc
wdd_run , file("resB1.csv")

use "dataB.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esprelength(2) espostlength(4)
dis "B2_RC=" _rc
wdd_run , file("resB2.csv")

use "dataB.dta", clear
cap noi wooldid y i t ttre [pw=w], att att_it att_i att_itime gteffects esprelength(2) espostlength(4)
dis "B3_RC=" _rc
wdd_run , file("resB3.csv")

use "dataB.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects customweightmultiplier(cwm) esprelength(2) espostlength(4)
dis "B4_RC=" _rc
wdd_run , file("resB4.csv")

use "dataB.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esfixedbaseperiod esrelativeto(-1) esprelength(2) espostlength(4)
dis "B5_RC=" _rc
wdd_run , file("resB5.csv")

* ================= Design C: cohort entering at its own onset =================
do genC.do
export delimited using "dataC.csv", replace
save "dataC.dta", replace

use "dataC.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "C1_RC=" _rc
wdd_run , file("resC1.csv")

use "dataC.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "C2_RC=" _rc
wdd_run , file("resC2.csv")

* C3: DEFAULT options -- no esprelength at all; this is the silent branch (ado:12086)
use "dataC.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects
dis "C3_RC=" _rc
wdd_run , file("resC3.csv")

* C4: esfixedbaseperiod -- the comparison regime (cohort 3 has no reference period,
*     so it should be dropped entirely)
use "dataC.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esfixedbaseperiod esrelativeto(-1) esprelength(3) espostlength(5)
dis "C4_RC=" _rc
wdd_run , file("resC4.csv")

dis "ALLDONE_RUN2"
