* probe1: assert the pinned build, then map what each weight type does at the parse surface.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid
dis "BUILDCHECK_START"
* the build must accept nofastcollapse (v1.3 does not)
dis "BUILDCHECK_END"

*--------------------------------------------------------------------
* generator: 24 units x 8 periods x 3 rows per (i,t); 6 never treated
*--------------------------------------------------------------------
clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 6
replace ttre = 5 if i > 6  & i <= 12
replace ttre = 6 if i > 12 & i <= 18
gen double ufe = rnormal()
expand 8
bysort i: gen t = _n
gen double tfe = 0.15*t
expand 3
bysort i t: gen r = _n
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ufe + tfe + 1.25*post + rnormal()*0.5
* integer frequency-style weight, 1..4
gen int k = 1 + mod(_n, 4)
compress
save "fw_base.dta", replace

dis "GEN_DONE N=" _N

*--------------------------------------------------------------------
* A. weight-type surface
*--------------------------------------------------------------------
cap noi wooldid y i t ttre [fw=k], att
dis "FW_PLAIN_RC=" _rc

cap noi wooldid y i t ttre [aw=k], att
dis "AW_PLAIN_RC=" _rc

cap noi wooldid y i t ttre [pw=k], att
dis "PW_PLAIN_RC=" _rc

* non-integer fw : does Stata's own syntax refuse it before wooldid does?
gen double knoni = k + 0.5
cap noi wooldid y i t ttre [fw=knoni], att
dis "FW_NONINT_RC=" _rc

* fw under oldsyntax (a second parse path)
cap noi wooldid y i t ttre [fw=k], att oldsyntax
dis "FW_OLDSYN_RC=" _rc

* what happens if WOOLDID_REPMODE is set in the user's session (leak simulation)
global WOOLDID_REPMODE 1
cap noi wooldid y i t ttre [fw=k], att level(90)
dis "FW_REPMODE_RC=" _rc
dis "FW_REPMODE_LEVEL=" e(level)
global WOOLDID_REPMODE
dis "GLOBAL_AFTER=[${WOOLDID_REPMODE}]"

* and a normal run with the leaked global but no weight: does level() get dropped?
global WOOLDID_REPMODE 1
cap noi wooldid y i t ttre, att level(90)
dis "LEAK_NOWT_RC=" _rc
dis "LEAK_NOWT_LEVEL=" e(level)
global WOOLDID_REPMODE

* baseline for comparison
cap noi wooldid y i t ttre, att level(90)
dis "CLEAN_NOWT_RC=" _rc
dis "CLEAN_NOWT_LEVEL=" e(level)

* zero / missing / negative weights
gen double kz = k
replace kz = 0 in 1/20
cap noi wooldid y i t ttre [aw=kz], att
dis "AW_ZERO_RC=" _rc
dis "AW_ZERO_N=" e(N)

gen double km = k
replace km = . in 1/20
cap noi wooldid y i t ttre [aw=km], att
dis "AW_MISS_RC=" _rc
dis "AW_MISS_N=" e(N)

gen double kn = k
replace kn = -1 in 1/20
cap noi wooldid y i t ttre [aw=kn], att
dis "AW_NEG_RC=" _rc

dis "PROBE1_DONE"
