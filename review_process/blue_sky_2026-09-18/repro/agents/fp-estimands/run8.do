clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260923
set sortseed 20260923
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"

* continuous-treatment panel, dispersed positive outcome
clear
set obs 30
gen int i = _n
gen double ui = 120 - 110 + 220*runiform()
gen int ttre = .
replace ttre = 4 if i <= 10
replace ttre = 6 if i > 10 & i <= 20
expand 8
bysort i: gen int t = _n
gen double dose = 0
replace dose = 0.5 + 2*runiform() if t >= ttre & ttre < .
bysort i (t): replace dose = dose[_N] if t >= ttre & ttre < .
gen double y = ui + 2.0*dose + rnormal()*1.5
save "dataK.dta", replace

foreach et in dydx dyex eydx eyex {
  use "dataK.dta", clear
  cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(`et')
  dis "CT_`et'_HDFE_RC=" _rc
  if _rc == 0 {
    local nm : colnames e(b)
    dis "CT_`et'_HDFE_names=`nm'"
    dis "CT_`et'_HDFE_ame=" %21x e(b)[1, 2]
  }
  use "dataK.dta", clear
  cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(`et') customvce(robust)
  dis "CT_`et'_FULL_RC=" _rc
  if _rc == 0 dis "CT_`et'_FULL_ame=" %21x e(b)[1, 2]
}
dis "ALLDONE_RUN8"
