clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

* tiny clean staggered panel: 40 units, 8 periods, cohorts 4/6/never
clear
set obs 40
gen long id = _n
gen double ttre = .
replace ttre = 4 if inrange(id,1,13)
replace ttre = 6 if inrange(id,14,26)
expand 8
bysort id: gen int tt = _n
gen double y = 0.3*id/40 + 0.1*tt + rnormal() + 0.8*(tt>=ttre & ttre<.)
gen double wt = 1
save "base.dta", replace

dis "=== A: clean baseline ==="
cap noi wooldid y id tt ttre , att
dis "SENT A rc=" _rc " att=" _b[main_att]

dis "=== B: one row with y<0 sentinel and MISSING id, poisson ==="
use "base.dta", clear
gen double ycount = ceil(exp(y))
set obs `=_N+1'
replace ycount = -1 in L
replace id = . in L
replace tt = 3 in L
cap noi wooldid ycount id tt ttre , att poisson
dis "SENT B rc=" _rc

dis "=== C: y<0 sentinel on a ZERO-WEIGHT row, poisson ==="
use "base.dta", clear
gen double ycount = ceil(exp(y))
replace ycount = -1 in 1
replace wt = 0 in 1
cap noi wooldid ycount id tt ttre [aw=wt] , att poisson
dis "SENT C rc=" _rc

dis "=== D: id==0 on a row whose y is MISSING (binary) ==="
use "base.dta", clear
set obs `=_N+1'
replace id = 0 in L
replace tt = 3 in L
replace y = . in L
replace wt = 1 in L
cap noi wooldid y id tt ttre , att
dis "SENT D rc=" _rc

dis "=== E: id==0 on a ZERO-WEIGHT row (binary) ==="
use "base.dta", clear
set obs `=_N+1'
replace id = 0 in L
replace tt = 3 in L
replace y = 1 in L
replace wt = 0 in L
cap noi wooldid y id tt ttre [aw=wt] , att
dis "SENT E rc=" _rc

dis "=== F: same as E but excluded by the user if() ==="
use "base.dta", clear
set obs `=_N+1'
replace id = 0 in L
replace tt = 3 in L
replace y = 1 in L
replace wt = 1 in L
cap noi wooldid y id tt ttre if id > 0 , att
dis "SENT F rc=" _rc

dis "SAMPLECONCEPT DONE"
