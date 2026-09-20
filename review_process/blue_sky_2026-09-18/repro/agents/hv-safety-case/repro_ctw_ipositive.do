* repro_ctw_ipositive.do  -- SELF-CONTAINED.  No other file needed.
* HV-CTW-IPOS: contreatwithin + contreat() silently skips the ado:4959 refusal
* "i must consist of positive, non-zero integers" (error 411), because the CS-98
* rewrite at ado:4829 sets the WORKING safety value to off ABOVE that check.
* sthlp:1208 lists this refusal, and sthlp:1195 says the rewrite leaves every
* data-integrity assert alone.  Pre-CS-98 (git e26b3dd^) it did fire.
* Consequence: the clean named refusal becomes a raw Stata r(452) naming an
* internal tempvar.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case"
which wooldid

clear
set obs 60
gen long unit = _n
gen double ttre = .
replace ttre = 5 if unit <= 15
replace ttre = 8 if unit > 15 & unit <= 30
gen double ufe = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort unit: gen double t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
save "`WD'/repro_ipos.dta", replace

foreach bad in noninteger negative zero {
  foreach ctw in "" "contreatwithin" {
    use "`WD'/repro_ipos.dta", clear
    gen double i = unit
    if "`bad'" == "noninteger" qui replace i = i + 0.5 if unit == 3
    if "`bad'" == "negative"   qui replace i = -9       if unit == 3
    if "`bad'" == "zero"       qui replace i = i - 1
    dis "=== bad=`bad'  ctw=[`ctw'] ==="
    cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) `ctw'
    dis "RC_`bad'_`ctw'=" _rc
    cap noi dis "ATT_`bad'_`ctw'=" %21x el(e(b),1,1)
  }
}

dis "EXPECTED: every leg refuses with r411 'i must consist of positive, non-zero"
dis "          integers, but does not' -- the screen sthlp:1208 lists."
dis "OBSERVED: without contreatwithin, r411 as documented.  WITH contreatwithin:"
dis "          noninteger -> r452 '__000005: factor variables may not contain"
dis "          noninteger values'; negative -> r452 '... may not contain negative"
dis "          values'; zero -> rc 0, the run completes with no message."
dis "ALLDONE"
