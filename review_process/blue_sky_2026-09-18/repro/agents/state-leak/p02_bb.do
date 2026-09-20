* p02: attribution control + back-to-back identical spec + RNG bleed on a plain run
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: attribution control -- reghdfe/gtools alone, no wooldid
wdgen
wdline A_pre
global WDH_RNGA1 "`c(rngstate)'"
qui reghdfe y post, absorb(unit year) cluster(state)
capture noisily gegen double __ctl = mean(y), by(unit)
qui drop __ctl
wdline A_post
global WDH_RNGA2 "`c(rngstate)'"
dis "A_RNG_UNCHANGED_BY_REGHDFE_GTOOLS=" ("$WDH_RNGA1" == "$WDH_RNGA2")

* ============ ARM B: plain wooldid -- does it consume the RNG stream?
wdgen
wdline B_pre
global WDH_RNGB1 "`c(rngstate)'"
capture noisily wooldid y unit year ttre, cluster(state)
dis "B_RC=" _rc
global WDH_RNGB2 "`c(rngstate)'"
dis "B_RNG_UNCHANGED_BY_WOOLDID=" ("$WDH_RNGB1" == "$WDH_RNGB2")
wdline B_post
wddumpe bb_S1

* ============ ARM C: SAME spec again, on the data wooldid just left behind
capture noisily wooldid y unit year ttre, cluster(state)
dis "C_RC=" _rc
wdline C_post
wddumpe bb_S2

* ============ ARM D: same spec on freshly regenerated (byte-identical) data
wdgen
capture noisily wooldid y unit year ttre, cluster(state)
dis "D_RC=" _rc
wddumpe bb_S3

* ============ ARM E: coeftabresults / emptycells program-scoping claims
* the worker does `set coeftabresults on' (ado:552) and `set emptycells drop'
* (ado:602) and relies on Stata's program-scoping to undo both.  Test the EFFECT,
* not a c() name: with coeftabresults off, regress leaves no r(table).
wdgen
capture noisily set coeftabresults off
dis "E_SETOFF_RC=" _rc
qui regress y post
capture confirm matrix r(table)
dis "E_RTABLE_PRESENT_BEFORE_WOOLDID=" (_rc == 0)
capture noisily wooldid y unit year ttre, cluster(state)
dis "E_RC=" _rc
qui regress y post
capture confirm matrix r(table)
dis "E_RTABLE_PRESENT_AFTER_WOOLDID=" (_rc == 0)
set coeftabresults on
wdline E_post

* ============ ARM F: emptycells effect test.  With emptycells drop, a factor
* interaction with an empty cell yields no coefficient for that cell.
wdgen
set emptycells keep
dis "F_EMPTYCELLS_BEFORE=" "`c(emptycells)'"
capture noisily wooldid y unit year ttre, cluster(state)
dis "F_RC=" _rc
dis "F_EMPTYCELLS_AFTER=" "`c(emptycells)'"

capture noisily creturn list

dis "P02_DONE"
