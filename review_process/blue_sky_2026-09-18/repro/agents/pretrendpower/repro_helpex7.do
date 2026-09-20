* ================================================================
* repro_helpex7.do -- SELF-CONTAINED.  The help's only worked example
* of pretrendpower() (wooldid.sthlp line 1546, "Example 7") is run
* verbatim, with the cohort variable actually named i as the help
* writes it.  The run answers with its own caution saying the joint
* pre-trend test the example demonstrates is not calibrated.
* Run: stata-mp -e do repro_helpex7.do ; read repro_helpex7.log
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 80
gen long i = _n
gen int ttre = .
replace ttre = 5 if mod(i,5)==1
replace ttre = 7 if mod(i,5)==2
replace ttre = 9 if mod(i,5)==3
expand 10
bysort i: gen int t = _n
gen byte D = (ttre<. & t>=ttre)
gen double y = 0.1*t + 0.3*sin(t) + 1.0*D + 0.3*rnormal()

* ---- verbatim from wooldid.sthlp:1546 ----
wooldid y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod anticipation(1) pretrendpower(50 80)
dis "EX7_RC=" _rc
dis "DONE_REPRO_HELPEX7"
