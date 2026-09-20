*=======================================================================
* SELF-CONTAINED REPRO -- fuzz agent
* CLAIM: a continuous-treatment run that asks for a semi-elasticity
*        (`semielasticity', or `contreatelasticitytype(eydx|eyex)')
*        dies inside `margins' with a RAW Mata traceback naming a
*        wooldid TEMPVAR THAT NO LONGER EXISTS -- "variable __00001M not
*        found / st_data(): 3598 Stata returned error" -- and no wooldid
*        `Error:' line at all.  wooldid's own Caution, printed one line
*        earlier, states that under `semielasticity' it "evaluates the
*        margins calls on the full estimation sample instead"; the
*        tempvar that the collapse frame owned has been dropped by the
*        time margins' numerical-derivative pass calls st_data() on it.
*        Under permute() the same root cause surfaces as
*        "Error: the observed comparison statistic for att could not be
*        read from the rep-0 run", which does not name the real problem.
* RUN:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_fuzz_semielasticity.do
*=======================================================================
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*--------------- a plain staggered panel with a continuous dose ---------------
set seed 20260917
set sortseed 20260917
clear
set obs 60
gen long unit   = _n
gen long cohort = 1 + mod(_n, 6)                 // 6 cohort codes
gen byte region = 1 + mod(_n, 5)
gen double a_i  = rnormal()
gen ttre = .
replace ttre = 2003 if inlist(cohort,1,2)
replace ttre = 2005 if inlist(cohort,3,4)
gen double dose_i = runiform()*3                 // unit-level intensity
replace  dose_i = 0 if runiform() < .4           // a spike at zero
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t  = 0.30*(t-2000)
gen byte post   = !missing(ttre) & t >= ttre
gen double dose = cond(post, dose_i, 0)
* NOTE: the outcome is CENTRED NEAR ZERO on purpose.  A semi-elasticity divides by the
* prediction, so what matters is whether the FITTED values straddle zero.  Rebuild with
* "50 +" in front and every call below returns rc 0 -- that is the control, and it is what
* shows the feature itself is fine; what is missing is a guard and a message.
gen double ync  = a_i + g_t + 0.8*dose + 0.3*rnormal()
* centre it so the FITTED values straddle zero -- that is the trigger: a semi-elasticity
* divides by the prediction.  Add a large positive constant instead and every call below
* returns rc 0; that control is what shows the feature itself works.
summarize ync, meanonly
replace ync = ync - r(mean)
summarize ync
di "OUTCOME RANGE: " r(min) " .. " r(max)
gen double w    = 0.5 + runiform()

*--------------- (1) the crash: analytic path, no engine ----------------------
di _n "=== (1) i=cohort, cluster(cohort), semielasticity -- EXPECT rc 0 ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) semielasticity
di "RC_1 = " _rc

di _n "=== (2) the same with nofastcollapse (the Caution's own remedy) ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) semielasticity nofastcollapse
di "RC_2 = " _rc

di _n "=== (3) contreatelasticitytype(eydx) instead ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) contreatelasticitytype(eydx)
di "RC_3 = " _rc

di _n "=== (4) control: the same run WITHOUT the elasticity request ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic)
di "RC_4 = " _rc

di _n "=== (5) i=unit + safety(off), the shape the fuzzer first hit ==="
capture noisily wooldid ync unit t ttre [aw=w], att robust contreat(dose) contreatcontrols(basic) semielasticity safety(off)
di "RC_5 = " _rc

di _n "=== (6) under permute(): the misleading downstream message ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) semielasticity permute(shuffle) permreps(20)
di "RC_6 = " _rc

di _n "EXPECTED: (4) rc 0 (control).  (1)(3)(5) either rc 0 or a wooldid 'Error:' line"
di "          that names the option combination and a remedy."
di "OBSERVED: a raw Mata traceback 'variable __00001x not found / st_data(): 3598'"
di "          with no wooldid message, and under permute() a rep-0 message that"
di "          does not name the real cause."
