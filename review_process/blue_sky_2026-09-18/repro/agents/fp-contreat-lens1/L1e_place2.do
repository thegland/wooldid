clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*--- 20 states x 5 counties x 10 years; region CROSSED with cohort ---
set obs 20
gen i = _n
gen region = mod(floor((_n-1)/5),4)+1
gen sfe = 40 + 8*rnormal()
gen coh = mod(_n,4)
gen ttre = .
replace ttre = 5 if coh==1
replace ttre = 7 if coh==2
replace ttre = 9 if coh==3
expand 5
bysort i: gen cty = _n
gen cfe = 2*rnormal()
expand 10
bysort i cty: gen t = _n
* wave: a randomly assigned 4-level county-year program indicator with a real
* effect on y. Orthogonal to treatment by construction, so it can be written
* either as fe(wave) [absorbed] or controls(i.wave) [main equation] with the
* SAME column space and the SAME fit.
gen wave = 1+int(4*runiform())
gen wfe  = 3.0*wave
gen tfe = 0.4*t
gen dose = 1 + 2*runiform()
gen trtpost = (!missing(ttre) & t>=ttre)
gen eps = 0.5*rnormal()
gen y = sfe + cfe + tfe + wfe + 0.35*dose + 1.3*dose*trtpost + eps
qui sum y
dis "YMEAN=" %12.6f r(mean) "  YSD=" %12.6f r(sd)

*======================================================================
* B. PLACEMENT INVARIANCE, reghdfe BOTH times, same column space.
*======================================================================
foreach et in dydx dyex eydx eyex {
  cap noisily qui wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust fe(wave)
  dis "B_FE_`et'_RC=" _rc
  mat bF = e(b)
  local r2F = e(r2)
  cap noisily qui wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust controls(i.wave)
  dis "B_CT_`et'_RC=" _rc
  mat bT = e(b)
  local r2T = e(r2)
  dis "B_`et'_R2  FE=" %18.14f `r2F' "  CT=" %18.14f `r2T'
  dis "B_`et'_ATT FE=" %21.12f bF[1,colnumb(bF,"main_att")] "  CT=" %21.12f bT[1,colnumb(bT,"main_att")] "  RELGAP=" %14.10f (bF[1,colnumb(bF,"main_att")]-bT[1,colnumb(bT,"main_att")])/bT[1,colnumb(bT,"main_att")]
  dis "B_`et'_AME FE=" %21.12f bF[1,colnumb(bF,"main_ame")] "  CT=" %21.12f bT[1,colnumb(bT,"main_ame")] "  RELGAP=" %14.10f (bF[1,colnumb(bF,"main_ame")]-bT[1,colnumb(bT,"main_ame")])/bT[1,colnumb(bT,"main_ame")]
}

*--- binary side: plain att and semielasticity, same placement test ---
cap noisily qui wooldid y i t ttre, robust fe(wave)
dis "C3_ATTFE_RC=" _rc
mat aF = e(b)
cap noisily qui wooldid y i t ttre, robust controls(i.wave)
dis "C3_ATTCT_RC=" _rc
mat aT = e(b)
dis "C3_ATT  FE=" %21.12f aF[1,colnumb(aF,"main_att")] "  CT=" %21.12f aT[1,colnumb(aT,"main_att")] "  RELGAP=" %14.10f (aF[1,colnumb(aF,"main_att")]-aT[1,colnumb(aT,"main_att")])/aT[1,colnumb(aT,"main_att")]

cap noisily qui wooldid y i t ttre, robust semielasticity fe(wave)
dis "C3_SEMIFE_RC=" _rc
mat sF = e(b)
cap noisily qui wooldid y i t ttre, robust semielasticity controls(i.wave)
dis "C3_SEMICT_RC=" _rc
mat sT = e(b)
dis "C3_SEMI FE=" %21.12f sF[1,colnumb(sF,"main_att")] "  CT=" %21.12f sT[1,colnumb(sT,"main_att")] "  RELGAP=" %14.10f (sF[1,colnumb(sF,"main_att")]-sT[1,colnumb(sT,"main_att")])/sT[1,colnumb(sT,"main_att")]

*--- and the engine-parity legs on THIS fixture, for the record ---
foreach et in dydx eydx eyex {
  qui cap wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust
  mat bD = e(b)
  qui cap wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') customvce(robust)
  mat bC = e(b)
  dis "F_`et'_AME DEF=" %21.12f bD[1,colnumb(bD,"main_ame")] "  CVE=" %21.12f bC[1,colnumb(bC,"main_ame")] "  RELGAP=" %14.10f (bD[1,colnumb(bD,"main_ame")]-bC[1,colnumb(bC,"main_ame")])/bC[1,colnumb(bC,"main_ame")]
}
qui cap wooldid y i t ttre, robust semielasticity
mat s1 = e(b)
qui cap wooldid y i t ttre, customvce(robust) semielasticity
mat s2 = e(b)
dis "F_semi DEF=" %21.12f s1[1,colnumb(s1,"main_att")] "  CVE=" %21.12f s2[1,colnumb(s2,"main_att")] "  RELGAP=" %14.10f (s1[1,colnumb(s1,"main_att")]-s2[1,colnumb(s2,"main_att")])/s2[1,colnumb(s2,"main_att")]
qui cap wooldid y i t ttre, robust
mat a1 = e(b)
qui cap wooldid y i t ttre, customvce(robust)
mat a2 = e(b)
dis "F_att  DEF=" %21.12f a1[1,colnumb(a1,"main_att")] "  CVE=" %21.12f a2[1,colnumb(a2,"main_att")] "  RELGAP=" %14.10f (a1[1,colnumb(a1,"main_att")]-a2[1,colnumb(a2,"main_att")])/a2[1,colnumb(a2,"main_att")]

dis "DONE_L1E"
