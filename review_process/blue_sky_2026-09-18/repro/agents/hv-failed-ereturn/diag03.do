clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
clear
set obs 48
gen int unit = _n
gen byte grp  = mod(_n,4)
gen double onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i = rnormal()
expand 10
bys unit: gen int year = 1999 + _n
gen byte post = !missing(onset) & year >= onset
gen double x1 = rnormal()
gen double y = a_i + 0.10*(year-2000) + 0.8*post*(year-onset+1) + 0.4*x1 + 0.30*rnormal()
replace y = a_i + 0.10*(year-2000) + 0.4*x1 + 0.30*rnormal() if missing(onset)
count if missing(onset)
dis "--- V1 double onset, i=unit ---"
cap noisily wooldid y unit year onset, att
dis "V1_RC=" _rc
dis "--- V2 double onset, i=unit, esfixedbaseperiod ---"
cap noisily wooldid y unit year onset, att esfixedbaseperiod
dis "V2_RC=" _rc
dis "--- V3 shifted time to 1..10 ---"
gen int yr2 = year-1999
gen double ons2 = onset-1999
cap noisily wooldid y unit yr2 ons2, att
dis "V3_RC=" _rc
dis "--- V4 coarse i (grp+1) ---"
gen byte gi = grp+1
cap noisily wooldid y gi year onset, att
dis "V4_RC=" _rc
dis "--- V5 safety off ---"
cap noisily wooldid y unit year onset, att safety(off)
dis "V5_RC=" _rc
dis "--- V6 canon generator, i=unit ---"
clear
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local SHAPE canon
include "`GENDIR'/gen_shapes.do"
cap noisily wooldid y unit t ttre, att
dis "V6_RC=" _rc
dis "--- V7 canon generator, i=cohort ---"
cap noisily wooldid y cohort t ttre, att
dis "V7_RC=" _rc
sum ttre
tab t if missing(ttre)
dis "DIAG03 DONE"
