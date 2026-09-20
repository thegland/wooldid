clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31415
set sortseed 31415
clear
set obs 320
gen long unit = ceil(_n/8)
gen int  t    = _n - (unit-1)*8
gen byte i    = ceil(unit/10)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
gen double dose = 0
qui replace dose = 1 + 2*runiform() if ttre<. & t>=ttre
gen double b1 = 0
gen double k0 = 0
qui replace k0 = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre<. & t>=ttre
qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre<. & t>=ttre
gen double y = 2 + 0.6*i + 0.35*t + k0 + b1*dose
tempfile F
qui save `F'
* the cell whose slope reg drops (last in the expansion): i=2, relyr=2
qui sum b1 if i==2 & t-ttre==2
dis "BETA_LASTCELL = " %16.10f r(mean)
qui sum b1 if ttre<. & t>=ttre
dis "MEAN_BETA     = " %16.10f r(mean)
dis "PREDICTED reg-path ame = mean(beta) - beta_lastcell"

dis "########## A: customvce(robust), NOISILY"
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) customvce(robust) regtol(14)
dis "A_RC=" _rc

dis "########## B: poisson robust, NOISILY"
use `F', clear
qui gen double xbi = 0.5 + 0.10*i + 0.04*t
qui replace xbi = xbi + 0.10 + 0.20*dose if ttre<. & t>=ttre
qui gen double ypo = exp(xbi + rnormal(0,0.15))
cap noisily wooldid ypo i t ttre , att contreat(dose) contreatcontrols(untreatedZero) poisson robust regtol(9)
dis "B_RC=" _rc

dis "########## C: default reghdfe robust, NOISILY (control)"
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) robust regtol(14)
dis "C_RC=" _rc
dis "DONE24"
