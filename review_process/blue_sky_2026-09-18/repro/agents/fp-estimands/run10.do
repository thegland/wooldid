clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260925
set sortseed 20260925
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"
* well-identified continuous-treatment design: coarse i, 14 units per cell, low noise
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
expand 7
bysort i: gen int t = _n
expand 14
gen double ui = 100 + 20*runiform()
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double dz = dose*post
gen double y = ui + 0.5*t + 3.0*dz + rnormal()*0.3
save "dataM.dta", replace
foreach cc in untreatedZero basic {
  foreach et in dydx eydx {
    use "dataM.dta", clear
    if "`cc'" == "untreatedZero" gen double d2 = dz
    if "`cc'" == "basic"         gen double d2 = dose
    cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(`cc') contreatelasticitytype(`et')
    dis "M_`cc'_`et'_HDFE_RC=" _rc
    if _rc == 0 {
      dis "M_`cc'_`et'_HDFE_att=" %21x e(b)[1,1]
      dis "M_`cc'_`et'_HDFE_ame=" %21x e(b)[1,2]
    }
    use "dataM.dta", clear
    if "`cc'" == "untreatedZero" gen double d2 = dz
    if "`cc'" == "basic"         gen double d2 = dose
    cap noi wooldid y i t ttre, att contreat(d2) contreatcontrols(`cc') contreatelasticitytype(`et') customvce(robust)
    dis "M_`cc'_`et'_FULL_RC=" _rc
    if _rc == 0 {
      dis "M_`cc'_`et'_FULL_att=" %21x e(b)[1,1]
      dis "M_`cc'_`et'_FULL_ame=" %21x e(b)[1,2]
    }
  }
}
dis "ALLDONE_RUN10"
