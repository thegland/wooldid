clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

set obs 20
gen i = _n
gen region = mod(_n,4)+1
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
gen post = t>5
gen zone = region*10 + post
gen zfe = 2.7*mod(zone,7) - 1.8*region
gen tfe = 0.4*t
gen dose = 1 + 2*runiform()
gen trtpost = (!missing(ttre) & t>=ttre)
gen eps = 0.5*rnormal()
gen y = sfe + cfe + tfe + zfe + 0.35*dose + 1.3*dose*trtpost + eps

* diagnose the fe() failure
cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(dydx) robust fe(i.zone)
dis "DIAG_FEi_RC=" _rc
cap noisily wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(dydx) robust fe(zone)
dis "DIAG_FEbare_RC=" _rc

*======================================================================
* B (retry). PLACEMENT INVARIANCE, reghdfe both times.
*    fe(zone) -> absorb();  controls(i.zone) -> main equation.
*======================================================================
foreach et in dydx dyex eydx eyex {
  qui cap wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust fe(zone)
  dis "B_FE_`et'_RC=" _rc
  mat bF = e(b)
  local r2F = e(r2)
  local nF  = e(N)
  qui cap wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust controls(i.zone)
  dis "B_CT_`et'_RC=" _rc
  mat bT = e(b)
  local r2T = e(r2)
  local nT  = e(N)
  dis "B_`et'_R2  FE=" %18.14f `r2F' "  CT=" %18.14f `r2T' "   N: " `nF' " / " `nT'
  local cF = colnumb(bF,"main_ame")
  local cT = colnumb(bT,"main_ame")
  local aF = colnumb(bF,"main_att")
  local aT = colnumb(bT,"main_att")
  dis "B_`et'_ATT FE=" %21.12f bF[1,`aF'] "  CT=" %21.12f bT[1,`aT'] "  RELGAP=" %14.10f (bF[1,`aF']-bT[1,`aT'])/bT[1,`aT']
  dis "B_`et'_AME FE=" %21.12f bF[1,`cF'] "  CT=" %21.12f bT[1,`cT'] "  RELGAP=" %14.10f (bF[1,`cF']-bT[1,`cT'])/bT[1,`cT']
}

* binary semielasticity placement
qui cap wooldid y i t ttre, robust semielasticity fe(zone)
dis "C2_FE_RC=" _rc
mat sF = e(b)
local r2F = e(r2)
qui cap wooldid y i t ttre, robust semielasticity controls(i.zone)
dis "C2_CT_RC=" _rc
mat sT = e(b)
local r2T = e(r2)
dis "C2_R2 FE=" %18.14f `r2F' "  CT=" %18.14f `r2T'
dis "C2_SEMI FE=" %21.12f sF[1,colnumb(sF,"main_att")] "  CT=" %21.12f sT[1,colnumb(sT,"main_att")] "  RELGAP=" %14.10f (sF[1,colnumb(sF,"main_att")]-sT[1,colnumb(sT,"main_att")])/sT[1,colnumb(sT,"main_att")]
qui cap wooldid y i t ttre, robust fe(zone)
mat aF = e(b)
qui cap wooldid y i t ttre, robust controls(i.zone)
mat aT = e(b)
dis "C2_ATT  FE=" %21.12f aF[1,colnumb(aF,"main_att")] "  CT=" %21.12f aT[1,colnumb(aT,"main_att")] "  RELGAP=" %14.10f (aF[1,colnumb(aF,"main_att")]-aT[1,colnumb(aT,"main_att")])/aT[1,colnumb(aT,"main_att")]

*======================================================================
* E. second seed, engine parity only (robustness of the 12-14% gap)
*======================================================================
set seed 77712345
replace eps = 0.5*rnormal()
replace dose = 1 + 2*runiform()
replace y = sfe + cfe + tfe + zfe + 0.35*dose + 1.3*dose*trtpost + eps
foreach et in dydx eydx eyex {
  qui cap wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') robust
  mat bD = e(b)
  qui cap wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') customvce(robust)
  mat bC = e(b)
  local cD = colnumb(bD,"main_ame")
  local cC = colnumb(bC,"main_ame")
  dis "E_`et'_AME DEF=" %21.12f bD[1,`cD'] "  CVE=" %21.12f bC[1,`cC'] "  RELGAP=" %14.10f (bD[1,`cD']-bC[1,`cC'])/bC[1,`cC']
}
qui cap wooldid y i t ttre, robust semielasticity
mat sD2 = e(b)
qui cap wooldid y i t ttre, customvce(robust) semielasticity
mat sC2 = e(b)
dis "E_semi DEF=" %21.12f sD2[1,colnumb(sD2,"main_att")] "  CVE=" %21.12f sC2[1,colnumb(sC2,"main_att")] "  RELGAP=" %14.10f (sD2[1,colnumb(sD2,"main_att")]-sC2[1,colnumb(sC2,"main_att")])/sC2[1,colnumb(sC2,"main_att")]

dis "DONE_L1D"
