clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

*----------------------------------------------------------------------
* Realistic panel: 50 units, 10 periods, 3 staggered onsets + never-treated.
* Unit FEs carry most of the LEVEL of y, as in real administrative data.
* dose is always defined (a policy intensity measured in every period);
* it only MOVES y inside treated post cells.
* zone = region x post : a time-varying control that a user may write
* either as fe(i.zone) [absorbed] or controls(i.zone) [main equation].
*----------------------------------------------------------------------
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
sum ufe

*----------------------------------------------------------------------
* A. engine parity on the SAME model: default(reghdfe) vs customvce(regress)
*----------------------------------------------------------------------
foreach et in dydx eydx eyex dyex {
  cap noisily wooldid y i t ttre, contreat(dose) contreatelasticitytype(`et') robust verbose
  dis "A_DEF_`et'_RC=" _rc
  if _rc==0 {
    mat bD = e(b)
    dis "A_DEF_`et'_B:"
    mat list bD, format(%21x)
  }
  cap noisily wooldid y i t ttre, contreat(dose) contreatelasticitytype(`et') customvce(robust) verbose
  dis "A_CVE_`et'_RC=" _rc
  if _rc==0 {
    mat bC = e(b)
    dis "A_CVE_`et'_B:"
    mat list bC, format(%21x)
  }
  cap {
    local nD = colnumb(bD,"ame")
    local nC = colnumb(bC,"ame")
  }
  cap noisily dis "A_`et'_AMEDEF=" %21.12f bD[1,colnumb(bD,"ame")] "  AMECVE=" %21.12f bC[1,colnumb(bC,"ame")] "  RELGAP=" %12.8f (bD[1,colnumb(bD,"ame")]-bC[1,colnumb(bC,"ame")])/bC[1,colnumb(bC,"ame")]
  dis "A_`et'_GAPRC=" _rc
}

*----------------------------------------------------------------------
* B. PLACEMENT INVARIANCE, reghdfe BOTH times:
*    fe(i.zone) -> absorb();  controls(i.zone) -> main equation.
*    Same column space, same fit. dydx must match; a well-defined
*    elasticity must match too.
*----------------------------------------------------------------------
foreach et in dydx eydx eyex {
  cap noisily wooldid y i t ttre, contreat(dose) contreatelasticitytype(`et') robust fe(i.zone)
  dis "B_FE_`et'_RC=" _rc
  if _rc==0 mat bF = e(b)
  cap noisily wooldid y i t ttre, contreat(dose) contreatelasticitytype(`et') robust controls(i.zone)
  dis "B_CT_`et'_RC=" _rc
  if _rc==0 mat bT = e(b)
  cap noisily dis "B_`et'_ATT_FE=" %21.12f bF[1,colnumb(bF,"att")] "  ATT_CT=" %21.12f bT[1,colnumb(bT,"att")]
  cap noisily dis "B_`et'_AME_FE=" %21.12f bF[1,colnumb(bF,"ame")] "  AME_CT=" %21.12f bT[1,colnumb(bT,"ame")] "  RELGAP=" %12.8f (bF[1,colnumb(bF,"ame")]-bT[1,colnumb(bT,"ame")])/bT[1,colnumb(bT,"ame")]
  dis "B_`et'_GAPRC=" _rc
}

*----------------------------------------------------------------------
* C. BINARY-side semielasticity: default vs customvce (no contreat at all)
*----------------------------------------------------------------------
cap noisily wooldid y i t ttre, robust semielasticity
dis "C_DEF_RC=" _rc
if _rc==0 mat sD = e(b)
cap noisily wooldid y i t ttre, customvce(robust) semielasticity
dis "C_CVE_RC=" _rc
if _rc==0 mat sC = e(b)
cap noisily dis "C_SEMI_DEF=" %21.12f sD[1,colnumb(sD,"att")] "  CVE=" %21.12f sC[1,colnumb(sC,"att")] "  RELGAP=" %12.8f (sD[1,colnumb(sD,"att")]-sC[1,colnumb(sC,"att")])/sC[1,colnumb(sC,"att")]
dis "C_GAPRC=" _rc

* and the plain att for the same two paths: must be identical
cap noisily wooldid y i t ttre, robust
if _rc==0 mat aD = e(b)
cap noisily wooldid y i t ttre, customvce(robust)
if _rc==0 mat aC = e(b)
cap noisily dis "C_ATT_DEF=" %21.12f aD[1,colnumb(aD,"att")] "  CVE=" %21.12f aC[1,colnumb(aC,"att")]
dis "C_ATTGAPRC=" _rc

dis "DONE_L1"
