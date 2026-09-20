*==============================================================================
* SELF-CONTAINED refutation probe for hv-reps-01 / hv-reps-02.
*
* CLAIM UNDER TEST (hv-reps-01): the low-reps() p-value-resolution Caution in
* wooldid_bootctl never prints on a bootstrap run that has no cluster(), at any
* reps(), because the block sits inside `if "`clus'" != "" {`.
* CLAIM UNDER TEST (hv-reps-02): wooldid_permctl prints its permreps()
* resolution caution BEFORE the refusal that permute() needs a single-variable
* cluster(), so a run that never executes is advised about its resolution.
*
* Own data, own legs. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do refute_reps_caution.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

* ---- assert the pinned build, not the v1.3 in PLUS ----
dis "=== BUILDCHECK ==="
which wooldid

* ---- my own panel: 70 counties, 12 years, three staggered cohorts + 20 never ----
clear
set obs 70
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 2010 if grp==1
replace ttre = 2012 if grp==2
replace ttre = 2014 if grp==3
gen double ufe = rnormal()
expand 12
bysort unit: gen int year = 2005 + _n - 1
gen double yfe = 0.05*(year-2005)
gen byte post = (year >= ttre) & !missing(ttre)
gen double y = 2 + ufe + yfe + 0.30*post + rnormal()*0.8
label var y "outcome"
qui count
dis "N_rows=" r(N)
qui levelsof ttre, local(cohs)
dis "cohorts: `cohs'"

*------------------------------------------------------------------
* LEG A : percentile bootstrap, reps(10), NO cluster()
*------------------------------------------------------------------
dis "=== LEG_A_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(10) rseed(101) dots(0)
dis "LEG_A_RC=" _rc
dis "=== LEG_A_END ==="

*------------------------------------------------------------------
* LEG B : identical + cluster(unit)
*------------------------------------------------------------------
dis "=== LEG_B_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(10) rseed(101) dots(0) cluster(unit)
dis "LEG_B_RC=" _rc
dis "=== LEG_B_END ==="

*------------------------------------------------------------------
* LEG C : bootdisplay(wald) + jointtests, reps(18), NO cluster()
*         joint p floor = 1/19 = .0526 > .05
*------------------------------------------------------------------
dis "=== LEG_C_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(18) bootdisplay(wald) jointtests esprelength(2) espostlength(2) rseed(102) dots(0)
dis "LEG_C_RC=" _rc
dis "=== LEG_C_END ==="

*------------------------------------------------------------------
* LEG D : identical to C + cluster(unit)
*------------------------------------------------------------------
dis "=== LEG_D_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(18) bootdisplay(wald) jointtests esprelength(2) espostlength(2) rseed(102) dots(0) cluster(unit)
dis "LEG_D_RC=" _rc
dis "=== LEG_D_END ==="

*------------------------------------------------------------------
* LEG E : boottstat, reps(10), NO cluster()   (floor 2/11 = .182)
*------------------------------------------------------------------
dis "=== LEG_E_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(10) boottstat rseed(103) dots(0)
dis "LEG_E_RC=" _rc
dis "=== LEG_E_END ==="

*------------------------------------------------------------------
* LEG F : reps(39), NO cluster()  -- sthlp:535/538 says 40 is the line,
*         so 39 (grid 2/39 = .0513 > .05) must still caution
*------------------------------------------------------------------
dis "=== LEG_F_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(39) rseed(104) dots(0)
dis "LEG_F_RC=" _rc
dis "=== LEG_F_END ==="

*------------------------------------------------------------------
* LEG F2: reps(39) + cluster(unit) -- positive control on the same reps()
*------------------------------------------------------------------
dis "=== LEG_F2_START ==="
cap noisily wooldid y unit year ttre, bootstrap(standard) reps(39) rseed(104) dots(0) cluster(unit)
dis "LEG_F2_RC=" _rc
dis "=== LEG_F2_END ==="

*------------------------------------------------------------------
* LEG G : permute() with NO cluster() -- caution then abort?
*------------------------------------------------------------------
dis "=== LEG_G_START ==="
cap noisily wooldid y unit year ttre, permute(shuffle) permreps(20) rseed(105) dots(0)
dis "LEG_G_RC=" _rc
dis "=== LEG_G_END ==="

*------------------------------------------------------------------
* LEG H : permreps(19) -- sthlp:636 says "neither engine has a hard minimum"
*------------------------------------------------------------------
dis "=== LEG_H_START ==="
cap noisily wooldid y unit year ttre, permute(shuffle) permreps(19) rseed(106) dots(0) cluster(unit)
dis "LEG_H_RC=" _rc
dis "=== LEG_H_END ==="

*------------------------------------------------------------------
* LEG I : permute() WITH cluster(), permreps(20) -- the caution's healthy path
*------------------------------------------------------------------
dis "=== LEG_I_START ==="
cap noisily wooldid y unit year ttre, permute(shuffle) permreps(20) rseed(107) dots(0) cluster(unit)
dis "LEG_I_RC=" _rc
dis "=== LEG_I_END ==="

dis "=== ALLDONE ==="
