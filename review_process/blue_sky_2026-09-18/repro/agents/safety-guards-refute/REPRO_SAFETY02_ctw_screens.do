*==============================================================================
* SELF-CONTAINED REPRO -- SAFETY-02 (contreatwithin's working-value rewrite)
* CLAIM UNDER TEST: contreatwithin + contreat() rewrites the WORKING safety value
* to off (ado:4829) and thereby silently disables (a) the float-above-2^24
* data-integrity caution (ado:4973) and (b) the positive-integer `i` assert
* (ado:4959) -- although sthlp:1199 says "every data-integrity assert" gates on
* the value the USER typed and is untouched by that rewrite.
* Also measures whether any NUMBER moves (it does not) and whether F13, the one
* data-integrity assert the help names in the contreatwithin entry, survives.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*------- (a) the float 2^24 caution: 8-digit calendar codes stored as float ----
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int yr = 2000 + _n - 1
gen double ttreyr = .
replace ttreyr = 2003 + mod(unit,3) if tg==1
gen float  tfl    = 20000101 + (yr-2000)*10000
gen float  ttrefl = cond(ttreyr<., 20000101 + (ttreyr-2000)*10000, .)
gen double dose = runiform()*2
gen double ui = rnormal()
bysort unit (yr): replace ui = ui[1]
gen double y = 1 + 0.2*(yr-2000) + ui + rnormal() + 1.2*dose*(yr>=ttreyr & ttreyr<.)

dis "--- A1 binary, default safety: EXPECT the float caution TWICE (tfl, ttrefl)"
cap noi wooldid y unit tfl ttrefl, att cluster(unit)
dis "A1_RC=" _rc
dis "--- A2 contreat WITHOUT contreatwithin: EXPECT the float caution TWICE"
cap noi wooldid y unit tfl ttrefl, att cluster(unit) contreat(dose) contreatcontrols(basic)
dis "A2_RC=" _rc
dis "--- A3 contreat WITH contreatwithin, no safety() typed: OBSERVED ZERO cautions"
cap noi wooldid y unit tfl ttrefl, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "A3_RC=" _rc
dis "(grep the log between A2 and A3 for '16,777,216' to count the cautions)"

*------- (b) the positive-integer i assert, and whether any number moves -------
clear
set obs 120
gen long unit = _n
gen byte tg = unit <= 90
gen long coh0 = cond(tg==1, mod(unit,3), 9)
gen long coh1 = coh0 + 1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2003 + mod(unit,3) if tg==1
gen double dose = runiform()*2
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t>=ttre & ttre<.)
gen double cohfrac = coh1 + 0.5

dis "--- B1 0-indexed cohort, binary default safety: EXPECT the named r(411)"
cap noi wooldid y coh0 t ttre, att cluster(unit)
dis "B1_RC=" _rc
dis "--- B2 0-indexed cohort, contreatwithin, NO safety() typed: assert skipped"
cap noi wooldid y coh0 t ttre, att att_i cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "B2_RC=" _rc
if _rc==0 {
  mat b2=e(b)
  dis "B2_ATT_HEX=" %21x b2[1,1]
  dis "B2_AME_HEX=" %21x b2[1,3]
}
dis "--- B3 1-indexed control: the numbers must be compared against B2"
cap noi wooldid y coh1 t ttre, att att_i cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "B3_RC=" _rc
if _rc==0 {
  mat b3=e(b)
  dis "B3_ATT_HEX=" %21x b3[1,1]
  dis "B3_AME_HEX=" %21x b3[1,3]
  dis "ATT_MOVED=" cond(b3[1,1]==b2[1,1],"NO","YES")
  dis "AME_MOVED=" cond(b3[1,3]==b2[1,3],"NO","YES")
}
dis "--- B4 NON-integer cohort under contreatwithin: the named r(411) is gone"
cap noi wooldid y cohfrac t ttre, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "B4_RC=" _rc

*------- (c) control: F13 onset-constancy DOES survive contreatwithin ----------
gen double ttre_brk = ttre if t >= ttre & ttre < .
dis "--- C1 broken (onset-onward-only) ttre + contreatwithin: EXPECT r(459)"
cap noi wooldid y coh1 t ttre_brk, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "C1_RC=" _rc
dis "REPRO_SAFETY02_DONE"
