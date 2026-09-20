clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
cd "`WD'"
do wddump.do

which wooldid
dis "BUILDCHECK done"

* ---- confirm fw is refused ----
do genA.do
cap noi wooldid y i t ttre [fw=w], att esprelength(3) espostlength(5)
dis "FWRC=" _rc

* ---- Design A ----
do genA.do
export delimited using "dataA.csv", replace
save "dataA.dta", replace

* A1 unweighted
use "dataA.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "A1_RC=" _rc
wdd_run , file("resA1.csv")

* A2 aweight
use "dataA.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "A2_RC=" _rc
wdd_run , file("resA2.csv")

* A3 pweight
use "dataA.dta", clear
cap noi wooldid y i t ttre [pw=w], att att_it att_i att_itime gteffects esprelength(3) espostlength(5)
dis "A3_RC=" _rc
wdd_run , file("resA3.csv")

* A4 aweight + cwm
use "dataA.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects customweightmultiplier(cwm) esprelength(3) espostlength(5)
dis "A4_RC=" _rc
wdd_run , file("resA4.csv")

* A5 esfixedbaseperiod
use "dataA.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects esfixedbaseperiod esrelativeto(-1) esprelength(3) espostlength(5)
dis "A5_RC=" _rc
wdd_run , file("resA5.csv")

* A6 unweighted, cwm only (no sample weights)
use "dataA.dta", clear
cap noi wooldid y i t ttre, att att_it att_i att_itime gteffects customweightmultiplier(cwm) esprelength(3) espostlength(5)
dis "A6_RC=" _rc
wdd_run , file("resA6.csv")

dis "ALLDONE_RUN1"
