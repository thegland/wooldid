clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

set obs 50
gen i = _n
gen region = mod(_n,5)+1
gen ufe = 40 + 8*rnormal()
gen coh = mod(_n,4)
gen ttre = .
replace ttre = 5 if coh==1
replace ttre = 7 if coh==2
replace ttre = 9 if coh==3
expand 10
bysort i: gen t = _n
gen post = t>5
gen zone = region*10 + post
gen zfe = 2.7*mod(zone,7) - 1.8*region
gen tfe = 0.4*t
gen dose = 1 + 2*runiform()
gen trtpost = (!missing(ttre) & t>=ttre)
gen eps = 0.5*rnormal()
gen y = ufe + tfe + zfe + 0.35*dose + 1.3*dose*trtpost + eps
sum y
dis "YMEAN=" r(mean)

* --- name discovery (once) ---
cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(dydx) robust
dis "N0_RC=" _rc
mat b0 = e(b)
local nm0 : colfullnames b0
dis "NAMES0: `nm0'"

*======================================================================
* A. ENGINE PARITY: default (reghdfe) vs customvce(robust) (plain regress)
*======================================================================
foreach et in dydx dyex eydx eyex {
  cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust
  dis "A_DEF_`et'_RC=" _rc
  mat bD = e(b)
  cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') customvce(robust)
  dis "A_CVE_`et'_RC=" _rc
  mat bC = e(b)
  local cD = colnumb(bD,"main:ame")
  local cC = colnumb(bC,"main:ame")
  local aD = colnumb(bD,"main:att")
  local aC = colnumb(bC,"main:att")
  dis "A_`et'_AME  DEF=" %21.12f bD[1,`cD'] "  CVE=" %21.12f bC[1,`cC'] "  RELGAP=" %14.10f (bD[1,`cD']-bC[1,`cC'])/bC[1,`cC']
  dis "A_`et'_ATT  DEF=" %21.12f bD[1,`aD'] "  CVE=" %21.12f bC[1,`aC'] "  RELGAP=" %14.10f (bD[1,`aD']-bC[1,`aC'])/bC[1,`aC']
  dis "A_`et'_HEXAMEDEF=" %21x bD[1,`cD'] "  HEXAMECVE=" %21x bC[1,`cC']
}

*======================================================================
* B. PLACEMENT INVARIANCE, reghdfe BOTH times:
*    fe(i.zone) -> absorb();  controls(i.zone) -> main equation.
*    Identical column space, identical fit.
*======================================================================
foreach et in dydx dyex eydx eyex {
  cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust fe(i.zone)
  dis "B_FE_`et'_RC=" _rc
  mat bF = e(b)
  local r2F = e(r2)
  cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust controls(i.zone)
  dis "B_CT_`et'_RC=" _rc
  mat bT = e(b)
  local r2T = e(r2)
  dis "B_`et'_R2  FE=" %14.10f `r2F' "  CT=" %14.10f `r2T'
  local cF = colnumb(bF,"main:ame")
  local cT = colnumb(bT,"main:ame")
  local aF = colnumb(bF,"main:att")
  local aT = colnumb(bT,"main:att")
  dis "B_`et'_ATT FE=" %21.12f bF[1,`aF'] "  CT=" %21.12f bT[1,`aT']
  dis "B_`et'_AME FE=" %21.12f bF[1,`cF'] "  CT=" %21.12f bT[1,`cT'] "  RELGAP=" %14.10f (bF[1,`cF']-bT[1,`cT'])/bT[1,`cT']
}

*======================================================================
* C. BINARY semielasticity: default vs customvce; and plain att control
*======================================================================
cap noisily wooldid y i t ttre, robust semielasticity
dis "C_DEF_RC=" _rc
mat sD = e(b)
cap noisily wooldid y i t ttre, customvce(robust) semielasticity
dis "C_CVE_RC=" _rc
mat sC = e(b)
local x1 = colnumb(sD,"main:att")
local x2 = colnumb(sC,"main:att")
dis "C_SEMI DEF=" %21.12f sD[1,`x1'] "  CVE=" %21.12f sC[1,`x2'] "  RELGAP=" %14.10f (sD[1,`x1']-sC[1,`x2'])/sC[1,`x2']

cap noisily wooldid y i t ttre, robust
mat aD = e(b)
cap noisily wooldid y i t ttre, customvce(robust)
mat aC = e(b)
local y1 = colnumb(aD,"main:att")
local y2 = colnumb(aC,"main:att")
dis "C_ATT  DEF=" %21.12f aD[1,`y1'] "  CVE=" %21.12f aC[1,`y2'] "  RELGAP=" %14.10f (aD[1,`y1']-aC[1,`y2'])/aC[1,`y2']

* binary semielasticity placement invariance too
cap noisily wooldid y i t ttre, robust semielasticity fe(i.zone)
mat sF = e(b)
cap noisily wooldid y i t ttre, robust semielasticity controls(i.zone)
mat sT = e(b)
local z1 = colnumb(sF,"main:att")
local z2 = colnumb(sT,"main:att")
dis "C_SEMIPLACE FE=" %21.12f sF[1,`z1'] "  CT=" %21.12f sT[1,`z2'] "  RELGAP=" %14.10f (sF[1,`z1']-sT[1,`z2'])/sT[1,`z2']

*======================================================================
* D. ORACLE: the help's own definition, computed by hand on the SAME fit.
*    eydx(att-weighted ame) = weighted mean over treated post rows of
*    (d yhat / d dose) / yhat_full.  Use the plain-regress fit, whose
*    prediction IS the full conditional mean.
*======================================================================
qui reg y i.i i.t i.zone c.dose ib0.trtpost#c.dose, vce(robust)
predict double xbfull, xb
gen double dydx_hand = _b[c.dose] + _b[1.trtpost#c.dose]*trtpost
gen double eydx_hand = dydx_hand/xbfull
gen double eyex_hand = dydx_hand*dose/xbfull
sum dydx_hand if trtpost==1
dis "D_DYDX_HAND=" %21.12f r(mean)
sum eydx_hand if trtpost==1
dis "D_EYDX_HAND=" %21.12f r(mean)
sum eyex_hand if trtpost==1
dis "D_EYEX_HAND=" %21.12f r(mean)

* what reghdfe's xb looks like vs the full prediction
qui reghdfe y c.dose ib0.trtpost#c.dose, absorb(i.i i.t i.zone) vce(robust)
predict double xbpart, xb
sum xbfull
dis "D_XBFULL mean=" %14.8f r(mean) " sd=" %14.8f r(sd)
sum xbpart
dis "D_XBPART mean=" %14.8f r(mean) " sd=" %14.8f r(sd)
sum y
dis "D_Y      mean=" %14.8f r(mean) " sd=" %14.8f r(sd)

dis "DONE_L1B"
