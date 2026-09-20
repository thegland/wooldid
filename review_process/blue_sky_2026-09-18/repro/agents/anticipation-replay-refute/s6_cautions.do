clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 100
which wooldid

use "dat.dta", clear
gen double dose = 0
replace dose = 1.5 if grp==1
replace dose = 3.0 if grp==2

*=========== CASE J: jointtests + cluster at the level of i =================
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
log using "J_live.txt", text replace
cap noi wooldid y uid t ttre, att cluster(uid) jointtests esprelength(3) espostlength(3)
log close
dis "@@J rc=" _rc
log using "J_replay.txt", text replace
cap noi wooldid
log close
dis "@@Jr rc=" _rc
qui ereturn list, all
dis "@@Jscalars: `e(scalars)'"
dis "@@Jmacros: `e(macros)'"

*=========== CASE P: permute(), 25 reps -> resolution floor caution =========
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
log using "P_live.txt", text replace
cap noi wooldid y uid t ttre, att robust permute(shuffle) permreps(25) rseed(777) dots(0)
log close
dis "@@P rc=" _rc
log using "P_replay.txt", text replace
cap noi wooldid
log close
dis "@@Pr rc=" _rc

*=========== CASE B: bootstrap(), 25 reps ==================================
clear all
use "dat.dta", clear
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
log using "B_live.txt", text replace
cap noi wooldid y uid t ttre, att robust bootstrap(standard) reps(25) rseed(777) dots(0)
log close
dis "@@B rc=" _rc
log using "B_replay.txt", text replace
cap noi wooldid
log close
dis "@@Br rc=" _rc

*=========== CASE C: continuous treatment, cohort-constant dose ============
clear all
use "dat.dta", clear
gen double dose = 0
replace dose = 1.5 if grp==1
replace dose = 3.0 if grp==2
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC WOOLDID_CLEANUP
log using "C_live.txt", text replace
cap noi wooldid y uid t ttre, att contreat(dose) cluster(uid)
log close
dis "@@C rc=" _rc
log using "C_replay.txt", text replace
cap noi wooldid
log close
dis "@@Cr rc=" _rc
qui ereturn list, all
dis "@@Cmacros: `e(macros)'"

dis "S6_DONE"
