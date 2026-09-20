* ============================================================================
* R9 (NEW FINDING, self-contained): poisson + contreat() + contreatcontrols
* that reduce to the POOLED ringer set always trips the ppmlhdfe rc-3200
* fallback, because the term wooldid itself puts inside absorb() is
*     i.<a variable that is identically 1> ## c.(ringer)
* i.e. a one-level categorical fixed effect - exactly the condition the
* fallback's own caution names and blames on the user's fe().
*
* Consequences: (a) a caution that misattributes wooldid's own construction,
* and whose closing line sends the user to inspect "the fe() specified" when
* they specified none; (b) every such run silently takes the slower
* explicit-regressor path; (c) that path is the collinearity-drop-order
* sensitive one, which is the mechanism behind the untreatedZero defect.
*
* CONTROL: contreatcontrols(byTime) makes `ccringertail' non-empty, which
* clears the one-level head (ado:1301-1303), and the fallback does NOT fire.
* CONTROL 2: a plain binary poisson run (no contreat) does NOT fire it.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R9_poisson_fallback.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 24680
set sortseed 24680

clear
set obs 720
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/12)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
qui replace ttre = 8 if i == 3
gen byte post = (ttre < . & t >= ttre)
gen double dose = 1 + 2*runiform()
gen double xb = 0.5 + 0.10*i + 0.04*t + 0.15*dose
qui replace xb = xb + 0.10 + 0.20*dose if post == 1
gen double ycount = rpoisson(exp(xb))
qui replace ycount = 0 if ycount >= .
* untreatedZero twin
gen double dose0 = 0
qui replace dose0 = 0.5 + 2*runiform() if post == 1
gen double xb0 = 0.5 + 0.10*i + 0.04*t
qui replace xb0 = xb0 + 0.10 + 0.20*dose0 if post == 1
gen double ycount0 = rpoisson(exp(xb0))
qui replace ycount0 = 0 if ycount0 >= .
tempfile F
qui save `F'

dis as result "===== A: poisson + contreat + contreatcontrols(basic) -- EXPECT the fallback ====="
use `F', clear
cap noi wooldid ycount i t ttre , att contreat(dose) contreatcontrols(basic) poisson robust regtol(9)
dis "RC_A=" _rc

dis as result "===== B: poisson + contreat + contreatcontrols(untreatedZero) -- EXPECT the fallback ====="
use `F', clear
cap noi wooldid ycount0 i t ttre , att contreat(dose0) contreatcontrols(untreatedZero) poisson robust regtol(9)
dis "RC_B=" _rc

dis as result "===== C (CONTROL): poisson + contreat + contreatcontrols(basic byTime) -- EXPECT NO fallback ====="
use `F', clear
cap noi wooldid ycount i t ttre , att contreat(dose) contreatcontrols(basic byTime) poisson robust regtol(9)
dis "RC_C=" _rc

dis as result "===== D (CONTROL): plain binary poisson, no contreat -- EXPECT NO fallback ====="
use `F', clear
cap noi wooldid ycount i t ttre , att poisson robust regtol(9)
dis "RC_D=" _rc

dis "DONE_R9"
