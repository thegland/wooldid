* ============================================================================
* SKEPTIC FOLLOW-UP 2 for hv-noparent-F3: the POLE test.
* One model, one panel, ONE additive shift c applied to y. For a linear model
* the marginal effect dy/dx is invariant to c and every fitted value moves by
* exactly c, so margins' eydx = mean_i( dydx / (yhat_i + c) ) must have a POLE
* wherever the shifted fitted values cross zero.
* If wooldid reports a finite, unremarked number on both sides of that pole,
* then (a) the fitted values in the evaluated subpopulation really do straddle
* zero there, and (b) the reported "semi-elasticity" is a mean of ratios with
* opposite signs -- a complete inference row for a quantity with no meaning.
* The level ATT is printed alongside at every c as the invariance control.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f3c_pole.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31337
set sortseed 31337

which wooldid

clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = 8 + rnormal()*0.5
expand 10
bysort unit: gen int t = _n
gen double ybase = ufe + 0.05*t + rnormal()*0.3
replace ybase = ybase + 0.8*(t>=ttre) if ttre<.
tempfile P
save `P', replace

use `P', clear
qui su ybase if t>=ttre & ttre<.
dis "BASE_MIN_TP=" r(min)
dis "BASE_MAX_TP=" r(max)
dis "BASE_MEAN_TP=" r(mean)

* one level run to pin the invariant marginal effect
use `P', clear
gen double y = ybase
capture noisily wooldid y unit t ttre, att robust
dis "LEVEL_RC=" _rc
capture noisily dis "LEVEL_ATT=" el(e(wdidmainresults),1,1)

* sweep the shift through the fitted-value support
foreach c in 0 -6 -8 -8.4 -8.7 -9 -10 -20 {
  use `P', clear
  gen double y = ybase + (`c')
  qui su y if t>=ttre & ttre<.
  dis "SHIFT `c' : MIN_TP=" r(min) " MAX_TP=" r(max) " MEAN_TP=" r(mean)
  capture noisily wooldid y unit t ttre, att robust semielasticity
  dis "SHIFT_`c'_RC=" _rc
  capture noisily dis "SHIFT_`c'_SEMI=" el(e(wdidmainresults),1,1)
  capture noisily dis "SHIFT_`c'_SE=" el(e(wdidmainresults),1,2)
  * level control at the same shift: must not move at all
  use `P', clear
  gen double y = ybase + (`c')
  capture noisily wooldid y unit t ttre, att robust
  capture noisily dis "SHIFT_`c'_LEVELATT=" el(e(wdidmainresults),1,1)
}

dis "DONE_F3C"
