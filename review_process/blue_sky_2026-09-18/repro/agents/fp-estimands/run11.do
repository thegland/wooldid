clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260925
set sortseed 20260925
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
* B: which term does the explicit arm drop?
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
expand 7
bysort i: gen int t = _n
expand 14
gen double ui   = 100 + 20*runiform()
gen double dose = 0.2 + 2.5*runiform()
gen byte  post  = (t >= ttre) & ttre < .
gen double d2   = dose*post
gen double y    = ui + 0.5*t + 3.0*d2 + rnormal()*0.3
save "dataN.dta", replace
dis "=== EXPLICIT arm, verbose ==="
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) customvce(robust) verbose
dis "VB_FULL_RC=" _rc
use "dataN.dta", clear
dis "=== DEFAULT arm, verbose ==="
cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(untreatedZero) verbose
dis "VB_HDFE_RC=" _rc

* A: confirm the late-entering (not always-treated) unit is the one dropped
use "dataC.dta", clear
qui levelsof t, local(tt)
dis "SAMPLE_T_RANGE=" "`: word 1 of `tt''" " to " "`: word `: word count `tt'' of `tt'''"
qui sum ttre if i > 18 & i <= 24
dis "COHORT3_TTRE=" r(mean)
qui sum t if i > 18 & i <= 24
dis "COHORT3_FIRST_T=" r(min)
dis "-> cohort 3's onset lies INSIDE the sample window; it is dropped anyway."
dis "ALLDONE_RUN11"
