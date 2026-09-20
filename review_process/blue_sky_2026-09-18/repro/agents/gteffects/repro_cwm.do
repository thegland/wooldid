*==========================================================================
* repro_cwm.do -- SELF-CONTAINED.
* CLAIM: with customweightmultiplier(), att_it is no longer "the simple,
* unweighted mean of all of the cohort-period-specific estimates"
* (sthlp att_it) and att_i is no longer "the simple, unweighted mean of the
* cohort-specific treatment effects" (sthlp att_i).  Cells/cohorts are instead
* weighted by their own multiplier mass, because the normaliser W_it / W_ipost
* is built from the sample weight ALONE (ado:1619-1620) while the numerator
* carries the multiplier (ado:1631/1637/1650).  No caution is printed.
* The proof needs no oracle: e(gteffects) reports the realised cell weights,
* and they are not equal.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_cwm.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- a plain staggered panel: 90 units, 9 years, two treated cohorts -------
clear
set obs 90
gen long unit = _n
gen byte coh = 1                       // 1 = never treated
replace coh = 2 if _n > 40 & _n <= 65
replace coh = 3 if _n > 65
gen double ttre = .
replace ttre = 2004 if coh == 2
replace ttre = 2006 if coh == 3
gen double ai = rnormal()
expand 9
bysort unit: gen int year = 1999 + _n
tempvar ts
bysort year: gen double `ts' = rnormal() if _n == 1
bysort year: replace `ts' = `ts'[1]
gen double te = 0
replace te = 1 + 0.25*(year - ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`ts' + te + rnormal()*0.3

* a multiplier a real user would supply: a cell-level population count.
* It is CONSTANT WITHIN each cohort-by-year cell, so it cannot be blamed on
* within-cell variation.
bysort coh year (unit): gen double pop = 1 + 4*runiform() if _n == 1
bysort coh year (unit): replace pop = pop[1]

*--------------------------------------------------------------------------
dis as result _n "=== (A) NO multiplier: cells carry weight 1, cohorts carry 1 ==="
quietly wooldid y coh year ttre, att att_it att_i att_itime gteffects ///
        esprelength(2) espostlength(2)
preserve
clear
quietly svmat double e(gteffects), names(col)
quietly keep if relyr >= 0
summarize w_att_it
dis as result "A: per-cell w_att_it   min=" %12.9g r(min) "  max=" %12.9g r(max)
bysort i: egen double sw = total(w_att_i)
bysort i: keep if _n == 1
summarize sw
dis as result "A: per-cohort w_att_i  min=" %12.9g r(min) "  max=" %12.9g r(max)
restore

*--------------------------------------------------------------------------
dis as result _n "=== (B) WITH customweightmultiplier(pop) -- watch for a caution ==="
noisily wooldid y coh year ttre, att att_it att_i att_itime gteffects ///
        esprelength(2) espostlength(2) customweightmultiplier(pop)
preserve
clear
quietly svmat double e(gteffects), names(col)
quietly keep if relyr >= 0
summarize w_att_it
dis as result "B: per-cell w_att_it   min=" %12.9g r(min) "  max=" %12.9g r(max)
dis as result "B: ratio max/min = " %12.9g r(max)/r(min)
bysort i: egen double sw = total(w_att_i)
bysort i: keep if _n == 1
summarize sw
dis as result "B: per-cohort w_att_i  min=" %12.9g r(min) "  max=" %12.9g r(max)
restore

dis as result _n "EXPECTED under the att_it / att_i help text: every number above = 1."
dis as result "OBSERVED in (B): they are not, and nothing in the output says so."
