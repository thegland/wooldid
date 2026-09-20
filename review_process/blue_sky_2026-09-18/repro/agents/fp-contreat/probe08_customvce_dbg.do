clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat"

clear
set obs 400
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/10)
gen double ttre = .
qui replace ttre = 5 if i == 1
qui replace ttre = 7 if i == 2
bysort unit (t): gen double ufe = 30 + 20*runiform() if _n==1
bysort unit (t): replace ufe = ufe[1]
gen double dose = 1 + 2*runiform()
gen double kappatrue = 0
gen double betatrue  = 0
qui replace kappatrue = 1.0 + 0.2*(t-ttre) if ttre < . & t >= ttre
qui replace betatrue  = 2.0 - 0.3*(t-ttre) if ttre < . & t >= ttre
gen double y = ufe + 0.5*t + 0.9*dose + kappatrue + betatrue*dose
save "`D'/fix2.dta", replace

dis "==================== A: fe(unit) + customvce ===================="
use "`D'/fix2.dta", clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) fe(unit) customvce(cluster unit) regtol(14) verbose
dis "RC_A=" _rc
dis "R2_A=" e(r2)

dis "==================== B: fe(unit) + cluster (reghdfe) ===================="
use "`D'/fix2.dta", clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) fe(unit) cluster(unit) regtol(14)
dis "RC_B=" _rc
dis "R2_B=" e(r2)

dis "==================== C: NO fe(), customvce ===================="
use "`D'/fix2.dta", clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) customvce(cluster unit) regtol(14)
dis "RC_C=" _rc
dis "R2_C=" e(r2)

dis "==================== D: NO fe(), cluster (reghdfe) ===================="
use "`D'/fix2.dta", clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) cluster(unit) regtol(14)
dis "RC_D=" _rc
dis "R2_D=" e(r2)
dis "DONE"
