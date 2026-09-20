*==============================================================================
* SELF-CONTAINED REPRO -- FP-NOPRE-SILENT
*
* CLAIM: under the pooled (default) reference regime, every unit whose effective
* onset is at or before its first in-sample period is dropped from the estimation
* sample entirely (wooldid_reltime, ado:12047-12055).  Nothing is printed, no e()
* entry names the screen, and the help never mentions it.  The only trace is a
* smaller e(N) and an e(N_dropped) whose help text (sthlp:1411) attributes drops
* to the `if' rule.
*
* Every other sample change wooldid makes is announced: the `if' rule prints a
* caution naming the count (sthlp:61); an anticipation-induced version of THIS
* SAME drop prints one (ado:12036-12046); ppmlhdfe separation prints one (ado:1513).
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_nopre_silent.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

* A plain state-year panel, t = 1..8.  Six units are "early adopters" whose record
* begins in their own adoption year -- an entirely ordinary administrative shape.
clear
set obs 36
gen int i = _n
gen double ui = rnormal()*2
gen int ttre = .
replace ttre = 4 if i <= 9                 // ordinary cohort
replace ttre = 6 if i > 9  & i <= 18       // ordinary cohort
replace ttre = 3 if i > 18 & i <= 24       // early adopters, see the drop below
expand 8
bysort i: gen int t = _n
drop if i > 18 & i <= 24 & t < 3           // they enter the panel in their onset year
gen double tau = 0
replace tau = 1 if t >= ttre & ttre < . & i <= 18
replace tau = 3 if t >= ttre & ttre < . & i > 18 & i <= 24   // a DIFFERENT effect
gen double y = ui + 0.25*t + tau + rnormal()*0.4

qui count
local rows = r(N)
qui levelsof i, local(all)
local units : word count `all'
qui levelsof i if ttre == 3, local(early)
local nearly : word count `early'
dis _n "rows in data = `rows'   units = `units'   early-adopter units = `nearly'"

dis _n "===== the run -- watch for ANY message about the early adopters ====="
cap noi wooldid y i t ttre, att att_i gteffects esprelength(3) espostlength(5)
dis _n "RC              = " _rc
dis "e(N)            = " e(N)     "   (rows the estimation actually used)"
dis "rows in data    = `rows'"
dis "e(N_dropped)    = " e(N_dropped)
dis "units in data   = `units'   of which early adopters = `nearly'"
dis _n "The six early-adopter cohorts, whose true effect is 3 against the others' 1,"
dis "contributed NOTHING, and no line above says so."

dis _n "===== contrast: the same drop, caused by anticipation(), IS announced ====="
cap noi wooldid y i t ttre if i <= 18 | i > 24, att anticipation(3) esprelength(2) espostlength(3)
dis "ANTIC_RC=" _rc
dis "REPRO_DONE"
