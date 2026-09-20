*==============================================================================
* repro_reps_caution.do  -- SELF-CONTAINED.
*
* CLAIM: the low-reps() p-value-resolution Caution that wooldid documents at
* sthlp:535 ("wooldid prints a Caution naming that grid whenever it is coarser
* than .05") and sthlp:634-636 ("reps() and permreps() both print a caution
* when the floor they imply sits above .05") is printed ONLY when cluster() is
* given.  Source: the block at ado:6473-6478 sits inside `if "`clus'" != "" {`
* which opens at ado:6459 and closes at ado:6501.
*
* Run:
*   cd <this directory>
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_reps_caution.do
*   grep -n "WDLEG\|Caution: with reps(\|resolution of the rank-based" repro_reps_caution.log
*
* PASS/FAIL of the claim: the caution fragment "Caution: with reps(" must appear
* after every WDLEG_*_CLUS sentinel and after NONE of the WDLEG_*_NOCLUS ones.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

dis "WDLEG_BUILDCHECK"
which wooldid

*------------------------------------------------------------------ data -----
* 100 units x 10 annual periods, three treated cohorts (onsets 2003/2005/2007),
* 40 never-treated units.  Ordinary balanced staggered panel; nothing exotic.
clear
set obs 100
gen long unit = _n
gen double _u = runiform()
sort _u unit
gen long _rk = _n
gen double ttre = .
qui replace ttre = 2003 if _rk <= 20
qui replace ttre = 2005 if _rk > 20 & _rk <= 40
qui replace ttre = 2007 if _rk > 40 & _rk <= 60
drop _u _rk
sort unit
gen long region = mod(unit - 1, 10) + 1
gen double a_i = rnormal()
expand 10
bysort unit: gen int t = 1999 + _n
gen byte post = (ttre < . & t >= ttre)
* y0: NO treatment effect (placebo).  Keeps the bootstrap p-values off 0, so the
*     conditional footer wooldid_pfloor (ado:10782) cannot mask the demonstration.
* y1: a large treatment effect, to show what the pfloor footer does instead.
gen double y0 = a_i + 0.10*(t - 2000) +                rnormal()*0.5
gen double y1 = a_i + 0.10*(t - 2000) + 1.5*post     + rnormal()*0.5
label var y0 "placebo outcome"
label var y1 "outcome with a real effect"
compress
dis "WDLEG_DATA_DONE N=" _N

*------------------------------------------------------------- leg 1 ---------
dis "WDLEG_1_PCTILE_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0)
dis "WDRC_1=" _rc
dis "WDCOMPLETED_1=" e(boot_completed)
cap mat drop M1
cap matrix M1 = e(wdidmainresults)
cap noi mat list M1, format(%9.4f)

*------------------------------------------------------------- leg 2 ---------
dis "WDLEG_2_PCTILE_CLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0) cluster(unit)
dis "WDRC_2=" _rc

*------------------------------------------------------------- leg 3 ---------
dis "WDLEG_3_BOOTTSTAT_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0) boottstat
dis "WDRC_3=" _rc

*------------------------------------------------------------- leg 4 ---------
dis "WDLEG_4_BOOTTSTAT_CLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0) boottstat cluster(unit)
dis "WDRC_4=" _rc

*------------------------------------------------------------- leg 5 ---------
dis "WDLEG_5_WALD_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0) bootdisplay(wald) jointtests
dis "WDRC_5=" _rc

*------------------------------------------------------------- leg 6 ---------
dis "WDLEG_6_WALD_CLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0) bootdisplay(wald) jointtests cluster(unit)
dis "WDRC_6=" _rc

*------------------------------------------------------------- leg 7 ---------
dis "WDLEG_7_BAYES_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(bayes) reps(10) rseed(101) dots(0)
dis "WDRC_7=" _rc

*------------------------------------------------------------- leg 8 ---------
* Strong effect, unclustered, low reps: every percentile p prints 0, so the
* CONDITIONAL footer wooldid_pfloor does fire.  This is the only thing that
* speaks on the unclustered leg, and only in this corner.
dis "WDLEG_8_PCTILE_NOCLUS_STRONG"
cap noi wooldid y1 unit t ttre, bootstrap(standard) reps(10) rseed(101) dots(0)
dis "WDRC_8=" _rc

*------------------------------------------------------------- leg 9 ---------
* Control: is the hazard the caution describes actually present unclustered?
* reps(10), percentile family -> grid 2/10 = .2.  Same run at reps(400).
dis "WDLEG_9_REFERENCE_400"
cap noi wooldid y1 unit t ttre, bootstrap(standard) reps(400) rseed(101) dots(0)
dis "WDRC_9=" _rc

dis "WDLEG_ALLDONE"
