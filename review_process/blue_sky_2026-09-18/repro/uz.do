clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242
* exact-fit log-linear: true dose semi-elasticity 0.25 EVERYWHERE (homogeneous)
local G 5
local U 8
local T 9
set obs `=`G'*`U''
gen long unit = _n
gen byte grp = mod(_n-1,`G')+1
gen ttre = .
replace ttre = 4 if grp==1
replace ttre = 6 if grp==2
replace ttre = 8 if grp==3
expand `T'
bys unit: gen t = _n
gen double aunit = 0.10*mod(unit,7)
gen double gt    = 0.05*t
* dose EXACTLY zero outside treated cohort-periods
gen double dose = 0
replace dose = exp(rnormal(0.2,0.5)) if t>=ttre & ttre<.
gen double eta = aunit + gt + 0.30*(t>=ttre & ttre<.) + 0.25*dose
gen double ycount = exp(eta)
dis "TRUE dose semi-elasticity = 0.25 (homogeneous, exact fit)"
foreach s in "poisson" "poisson robust" "poisson cluster(unit)" {
  dis "=== SPEC: att contreat(dose) contreatcontrols(untreatedZero) `s' ==="
  cap noi wooldid ycount grp t ttre, att contreat(dose) contreatcontrols(untreatedZero) `s'
  dis "RC=" _rc
  cap mat R = e(wdidmainresults)
  cap mat list R, format(%14.10f)
}
dis "=== CONTROL: default reghdfe (linear), truth for the linear twin ==="
gen double ylin = aunit + gt + 0.30*(t>=ttre & ttre<.) + 0.25*dose
cap noi wooldid ylin grp t ttre, att contreat(dose) contreatcontrols(untreatedZero) cluster(unit)
cap mat R2 = e(wdidmainresults)
cap mat list R2, format(%14.10f)
dis "UZ_DONE"
