*==============================================================================
* SELF-CONTAINED repro: a bare `wooldid' typed while another command owns e()
* arms the replay-staleness marker, which then refuses a valid stored model.
* Run:  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_SK_REPLAYGUARD.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 60
gen long uid = _n
gen byte grp = mod(_n,3)
gen ttre = .
replace ttre = 5 if grp==1
replace ttre = 8 if grp==2
gen double ufe = rnormal()
expand 12
bysort uid: gen int t = _n
gen double y = ufe + 0.10*t + rnormal()
replace y = y + 0.8*(t - ttre + 1) if ttre < . & t >= ttre
gen double x1 = rnormal()

macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
set level 95

qui wooldid y uid t ttre, att robust
est store GOOD
dis "## A good wooldid: e(cmd)=[`e(cmd)'] e(N)=" e(N) " LASTRC=[$WOOLDID_LASTRC]"

qui regress y x1
dis "## B another command owns e(): e(cmd)=[`e(cmd)'] LASTRC=[$WOOLDID_LASTRC]"

dis "## C now type a bare wooldid (no arguments) to redisplay the wooldid results:"
cap noi wooldid
dis "## C rc=" _rc " LASTRC=[$WOOLDID_LASTRC]   <- the replay ATTEMPT armed the marker"

est restore GOOD
dis "## D estimates restore GOOD: e(cmd)=[`e(cmd)'] e(N)=" e(N) "  (a perfectly valid stored model)"
cap noi wooldid
dis "## D rc=" _rc "   <- refused, blaming an estimation failure that never happened"

*---- recovery: any successful wooldid clears it
qui wooldid y uid t ttre, att robust
est restore GOOD
cap noi wooldid
dis "## E after a refit, replay of the stored model works again: rc=" _rc

macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
dis "## REPRO_SK_REPLAYGUARD_DONE"
