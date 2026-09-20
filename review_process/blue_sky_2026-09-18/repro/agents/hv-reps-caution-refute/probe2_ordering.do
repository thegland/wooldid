*==============================================================================
* SELF-CONTAINED probe 2 (skeptic's own hunting, adjacent to hv-reps-01/02).
*
* Q1 (r() clobber): the reps() caution block was moved INTO the cluster branch,
*     between `gstats sum` and the `if r(max) < 10` few-cluster caution. Does
*     executing it destroy r(max) and silence the cluster-count cautions?
*     J1 = reps(10) (the caution prints), J2 = reps(60) (it does not).
*     Both must still print "only 6 clusters" and "only 3 treated cluster(s)".
*
* Q2 (the RF-2 rule, applied only inside the subprogram): commit 813c414 moved
*     the caution "below every aborting validation". wooldid_bootctl is called
*     at ado:4941 and the identifying-assumption checks (r459: no never-treated)
*     run AFTER it in the caller. So a clustered low-reps run on a panel with no
*     never-treated units should STILL be advised about the p-value resolution
*     of a run that then aborts.
*     K  = all units eventually treated, cluster(unit), bootstrap reps(10)
*     K2 = same panel, no bootstrap (control: the r459 alone)
*
* Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe2_ordering.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

dis "=== BUILDCHECK ==="
which wooldid

*---------------------------------------------------------------
* Panel 1: 6 states x 12 units, 12 years. Treatment assigned at
* state level: 3 treated states, 3 never-treated states.
*---------------------------------------------------------------
clear
set obs 72
gen long unit = _n
gen byte state = mod(_n,6) + 1
gen ttre = .
replace ttre = 2010 if state==1
replace ttre = 2012 if state==2
replace ttre = 2014 if state==3
gen double ufe = rnormal()
expand 12
bysort unit: gen int year = 2005 + _n - 1
gen byte post = (year >= ttre) & !missing(ttre)
gen double y = 2 + ufe + 0.05*(year-2005) + 0.30*post + rnormal()*0.8
qui levelsof state if !missing(ttre), local(ts)
dis "treated states: `ts'"
save "panel1.dta", replace

dis "=== LEG_J1_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(10) rseed(201) dots(0) cluster(state)
dis "LEG_J1_RC=" _rc
dis "=== LEG_J1_END ==="

dis "=== LEG_J2_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(60) rseed(202) dots(0) cluster(state)
dis "LEG_J2_RC=" _rc
dis "=== LEG_J2_END ==="

*---------------------------------------------------------------
* Panel 2: every unit eventually treated -> the documented r(459)
*---------------------------------------------------------------
clear
set obs 60
gen long unit = _n
gen ttre = 2009 + 2*mod(_n,3)
gen double ufe = rnormal()
expand 12
bysort unit: gen int year = 2005 + _n - 1
gen byte post = (year >= ttre)
gen double y = 2 + ufe + 0.05*(year-2005) + 0.30*post + rnormal()*0.8
qui count if missing(ttre)
dis "never_treated_rows=" r(N)
save "panel2.dta", replace

dis "=== LEG_K_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(10) rseed(203) dots(0) cluster(unit)
dis "LEG_K_RC=" _rc
dis "=== LEG_K_END ==="

dis "=== LEG_K2_START ==="
cap noisily wooldid y unit year ttre, cluster(unit)
dis "LEG_K2_RC=" _rc
dis "=== LEG_K2_END ==="

*---------------------------------------------------------------
* Leg L: the same premature-caution shape on the permute engine,
* with a legal cluster() this time -- so the only abort left is
* the caller's r459.
*---------------------------------------------------------------
dis "=== LEG_L_START ==="
cap noisily wooldid y unit year ttre, permute(shuffle) permreps(20) rseed(204) dots(0) cluster(unit)
dis "LEG_L_RC=" _rc
dis "=== LEG_L_END ==="

dis "=== ALLDONE ==="
