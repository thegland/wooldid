* repro_safety_case.do  -- SELF-CONTAINED.  No other file needed.
* HV-SAFETY-CASE-1: safety() is validated at ado:5473, ~500 lines after ado:4617
* first consults it, and the test is case-sensitive.  This shows exactly what each
* casing does, and shows the one observable bite: a safety-keyed ERROR fires first
* and tells the user to "specify safety(off)" when they already typed safety(OFF).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid                 // MUST print the .../bluesky/build path

clear
set obs 60
gen long unit = _n
gen double i = unit
gen double ttre = .
replace ttre = 5 if unit <= 15
replace ttre = 7 if unit > 15 & unit <= 30
replace ttre = 9 if unit > 30 & unit <= 45
gen double ufe = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort unit: gen double t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = ufe + 0.1*t + 2.0*post*dose + rnormal()*0.5
gen double wt = 1 + runiform()

dis "--- accepted spellings ---"
foreach v in off on {
  cap noi wooldid y i t ttre, att safety(`v')
  dis "RC_safety(`v')=" _rc
}
cap noi wooldid y i t ttre, att safety( off )
dis "RC_safety(_off_)=" _rc

dis "--- rejected spellings (r198, 'safety() must be on or off.') ---"
foreach v in OFF Off oFF ON On {
  cap noi wooldid y i t ttre, att safety(`v')
  dis "RC_safety(`v')=" _rc
}

dis "--- THE BITE: the safety-keyed error at ado:4617 fires ~500 lines before validation ---"
cap noi wooldid y i t ttre, att_i unconditionalse safety(OFF)
dis "RC_bite=" _rc
dis "EXPECTED: r198 'safety() must be on or off.'"
dis "OBSERVED: r198 'unconditional standard errors may be unreliable with non-att"
dis "          estimand types.  Set option safety(off) to disregard this error.'"
dis "          -- the casing is never mentioned."

dis "--- a safety-keyed CAUTION also prints inside the window ---"
cap noi wooldid y i t ttre [aw=wt], att unconditionalse robust safety(OFF)
dis "RC_caution=" _rc
dis "ALLDONE"
