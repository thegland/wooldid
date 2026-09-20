clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
* TRUE EXACT FIT: every term lies in the span of {i.i, i.t, cell dummies, cell#dose}
set seed 616161
set sortseed 616161
set obs 32
gen long unit = _n
gen byte i0 = mod(unit,4)+1
gen ttre = .
replace ttre = 4 if i0==1
replace ttre = 6 if i0==2
replace ttre = 8 if i0==3
gen byte irev = 9 - i0
expand 10
bys unit: gen t = _n
gen byte trt = (t>=ttre & ttre<.)
* dose varies WITHIN each treated cohort-period cell (8 units per cohort), 0 elsewhere
gen double dose = 0
replace dose = 0.5 + mod(unit,7)*0.3 if trt==1
gen double beta = 0.55 + 0.20*i0 + 0.10*(t-ttre) if trt==1
qui sum beta if trt==1
local TRUE = r(mean)
* NO unit-level term: cohort + time + cell effects only, all absorbed by the saturated model
gen double y = 0.7*i0 + 0.15*t + (0.3*i0 + 0.05*t)*trt + cond(trt==1, beta*dose, 0)
dis "TRUE_AME = " %14.10f `TRUE'
foreach cc in untreatedZero basic {
  foreach iv in i0 irev {
    cap noi wooldid y `iv' t ttre, att contreat(dose) contreatcontrols(`cc')
    mat M = e(wdidmainresults)
    dis "RESULT cc=`cc' i=`iv' R2=" %6.4f e(r2) "  att=" %14.10f M[1,1] "  ame=" %14.10f M[2,1]
  }
}
dis "C5B_DONE"
