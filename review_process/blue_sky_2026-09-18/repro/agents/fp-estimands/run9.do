clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260924
set sortseed 20260924
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
* i = COHORT (coarse), several units per cohort-period carrying different doses, so the
* dose varies inside the cell and the ame is identified.
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
expand 7
bysort i: gen int t = _n
expand 8
bysort i t: gen int u = _n
gen double ui = 20 + 200*runiform()
gen double dose = 0
replace dose = 0.2 + 2.5*runiform() if t >= ttre & ttre < .
gen double y = ui + 3.0*dose + rnormal()*1.5
save "dataL.dta", replace
foreach et in dydx eydx {
  use "dataL.dta", clear
  cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(`et')
  dis "L_`et'_HDFE_RC=" _rc
  if _rc == 0 {
    local nm : colnames e(b)
    dis "L_`et'_HDFE_names=`nm'"
    matrix BB = e(b)
    forval c = 1/`=colsof(BB)' {
      local n2 : word `c' of `nm'
      dis "L_`et'_HDFE  `n2' = " %21x BB[1,`c']
    }
  }
  use "dataL.dta", clear
  cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(`et') customvce(robust)
  dis "L_`et'_FULL_RC=" _rc
  if _rc == 0 {
    local nm : colnames e(b)
    matrix BB = e(b)
    forval c = 1/`=colsof(BB)' {
      local n2 : word `c' of `nm'
      dis "L_`et'_FULL  `n2' = " %21x BB[1,`c']
    }
  }
}
dis "ALLDONE_RUN9"
