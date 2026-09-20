*==============================================================================
* SELF-CONTAINED repro: replay is not byte-for-byte, and drops the cautions
* that disown the numbers it reprints.
* Run:  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_SK_REPLAYDIFF.do
* Then, in a shell, in this directory:
*   for f in RJ RB RC; do echo "== $f"; diff ${f}_live.txt ${f}_replay.txt; done
* `qui log using' writes only wooldid's own output to these files, so the diff
* is a literal byte comparison of the live display against the replay display.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
set linesize 100
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
gen double dose = 0
replace dose = 1.5 if grp==1
replace dose = 3.0 if grp==2
tempfile D
save `D'

*---- RJ: jointtests with cluster() at the level of i -----------------------
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui log using "RJ_live.txt", text replace
cap noi wooldid y uid t ttre, att cluster(uid) jointtests esprelength(3) espostlength(3)
qui log close
qui log using "RJ_replay.txt", text replace
cap noi wooldid
qui log close
dis "## RJ done"

*---- RB: bootstrap(standard), 25 reps --------------------------------------
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui log using "RB_live.txt", text replace
cap noi wooldid y uid t ttre, att robust bootstrap(standard) reps(25) rseed(777) dots(0)
qui log close
qui log using "RB_replay.txt", text replace
cap noi wooldid
qui log close
dis "## RB done"

*---- RC: continuous treatment with a cohort-constant dose ------------------
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui log using "RC_live.txt", text replace
cap noi wooldid y uid t ttre, att contreat(dose) contreatcontrols(basic) cluster(uid)
qui log close
qui log using "RC_replay.txt", text replace
cap noi wooldid
qui log close
dis "## RC done"

*---- RP: permute(), 25 reps -----------------------------------------------
use `D', clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
qui log using "RP_live.txt", text replace
cap noi wooldid y uid t ttre, att cluster(uid) permute(shuffle) permreps(25) rseed(777) dots(0)
qui log close
qui log using "RP_replay.txt", text replace
cap noi wooldid
qui log close
dis "## RP done"

dis "## REPRO_SK_REPLAYDIFF_DONE"
