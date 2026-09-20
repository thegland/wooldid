clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260921
set sortseed 20260921
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
do wddump.do
do genI.do
export delimited using "dataI.csv", replace
save "dataI.dta", replace
count
dis "I_NOBS=" r(N)

local base "att att_it att_i att_itime gteffects esprelength(2) espostlength(3)"

use "dataI.dta", clear
cap noi wooldid y i t ttre [aw=w], `base' controls(x) cluster(i)
dis "I1_RC=" _rc
wdd_run , file("resI1.csv")

use "dataI.dta", clear
cap noi wooldid y i t ttre [aw=w], `base' controls(x) cluster(i) nofastcollapse
dis "I2_RC=" _rc
wdd_run , file("resI2.csv")

use "dataI.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it gteffects esprelength(2) espostlength(3) controls(x) cluster(cl) impvalues
dis "I3_RC=" _rc
wdd_run , file("resI3.csv")

use "dataI.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it gteffects esprelength(2) espostlength(3) controls(x) cluster(cl) impvalues nofastcollapse
dis "I4_RC=" _rc
wdd_run , file("resI4.csv")

use "dataI.dta", clear
cap noi wooldid y i t ttre [aw=w], `base' controls(x) cluster(i) subgroups(sg) jointtests
dis "I5_RC=" _rc
wdd_run , file("resI5.csv")

use "dataI.dta", clear
cap noi wooldid y i t ttre [aw=w], `base' controls(x) cluster(i) subgroups(sg) jointtests nofastcollapse
dis "I6_RC=" _rc
wdd_run , file("resI6.csv")

* fe() varying inside the cell
use "dataI.dta", clear
gen int fev = 1 + mod(_n,4)
cap noi wooldid y i t ttre [aw=w], `base' fe(fev) cluster(i)
dis "I7_RC=" _rc
wdd_run , file("resI7.csv")
use "dataI.dta", clear
gen int fev = 1 + mod(_n,4)
cap noi wooldid y i t ttre [aw=w], `base' fe(fev) cluster(i) nofastcollapse
dis "I8_RC=" _rc
wdd_run , file("resI8.csv")

dis "ALLDONE_RUN6"
