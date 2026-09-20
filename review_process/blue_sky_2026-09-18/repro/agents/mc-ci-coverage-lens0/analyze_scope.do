clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 40404
set sortseed 40404

* ===================== PART 1: analyse the MC =====================
use mcres, clear
count
scalar tgt_IM  = (1.00+1.20+0.80)/3
scalar tgt_att = (8*1.00 + 6*1.20 + 4*0.80)/18
dis "TARGET_IM(cluster-equal-weighted) = " tgt_IM
dis "TARGET_ATT(cell-weighted)         = " tgt_att

foreach arm in A B C {
  quietly sum b`arm'
  scalar mn`arm' = r(mean)
  gen byte cown`arm' = (l`arm'<=mn`arm' & u`arm'>=mn`arm') if !missing(l`arm')
}
foreach arm in A B {
  gen byte ctgt`arm' = (l`arm'<=tgt_IM & u`arm'>=tgt_IM) if !missing(l`arm')
}
gen byte ctgtC = (lC<=tgt_att & uC>=tgt_att) if !missing(lC)
gen double lenA = uA-lA
gen double lenB = uB-lB
gen double lenC = uC-lC

dis "================ MC RESULTS (R=120) ================"
foreach arm in A B C {
  quietly count if !missing(se`arm')
  local nn = r(N)
  quietly sum b`arm'
  local mb = r(mean)
  local sdb = r(sd)
  quietly sum se`arm'
  local ms = r(mean)
  quietly sum cown`arm'
  local co = r(mean)
  quietly sum ctgt`arm'
  local ct = r(mean)
  quietly sum len`arm'
  local ml = r(mean)
  dis "ARM `arm': n=`nn' mean(b)=" %8.5f `mb' " sd(b)=" %8.5f `sdb' ///
      " mean(se)=" %8.5f `ms' " se/sd=" %7.4f `ms'/`sdb' ///
      " meanCIlen=" %7.4f `ml' " cov_ownmean=" %6.3f `co' " cov_truth=" %6.3f `ct'
}
dis "3*binomSE at .95 with R=" _N " : " 3*sqrt(.95*.05/_N)
dis "ARM_A_is: impvalues cluster(unit)  [cluster FINER than i=cg]"
dis "ARM_B_is: impvalues cluster(cg)    [cluster == i, textbook IM]"
dis "ARM_C_is: analytic  cluster(unit)"
quietly sum lenB
local lb = r(mean)
quietly sum lenA
dis "CI_LENGTH_RATIO A/B = " r(mean)/`lb'

* ===================== PART 2: scope, one deterministic draw ==========
clear
set obs 60
gen unit = _n
gen cg = ceil(_n/15)
gen ttre = .
replace ttre = 5 if cg==1
replace ttre = 7 if cg==2
replace ttre = 9 if cg==3
gen double ui = rnormal(0,0.80)
gen double uhet = rnormal(0,0.30)
expand 12
bysort unit: gen t = _n
gen double te = 0.05*t
gen double eps = rnormal(0,0.60)
gen double effc = 0
replace effc = 1.00 if cg==1
replace effc = 1.20 if cg==2
replace effc = 0.80 if cg==3
gen byte post = !missing(ttre) & t>=ttre
gen double y = 2 + ui + te + post*(effc + uhet) + eps
* continuous dose, varying WITHIN cohort-period cells
gen double dose = 0
replace dose = runiform(0.2,2.0) if post
gen double ycont = 2 + ui + te + post*(0.45*dose) + eps

dis "=== P2a: EVENT STUDY horizons, impvalues cluster(unit) ==="
capture noisily wooldid y cg t ttre, att cluster(unit) impvalues esprelength(2) espostlength(2)
dis "P2a_RC=" _rc
capture noisily matrix list e(wdidmainresults)
matrix ESU = e(wdidmainresults)

dis "=== P2b: EVENT STUDY horizons, impvalues cluster(cg) ==="
capture noisily wooldid y cg t ttre, att cluster(cg) impvalues esprelength(2) espostlength(2)
dis "P2b_RC=" _rc
capture noisily matrix list e(wdidmainresults)
matrix ESG = e(wdidmainresults)

dis "=== ES SE RATIO per row (unit-cluster / cg-cluster) ==="
local nr = rowsof(ESU)
local nrg = rowsof(ESG)
dis "rows_unit=`nr' rows_cg=`nrg'"
forvalues j = 1/`nr' {
  local rn : word `j' of `: rownames ESU'
  if `j' <= `nrg' {
    dis "ESROW `j' `rn' : b_u=" %9.5f ESU[`j',1] " se_u=" %9.6f ESU[`j',2] ///
        " b_g=" %9.5f ESG[`j',1] " se_g=" %9.6f ESG[`j',2] ///
        " ratio=" %8.5f ESU[`j',2]/ESG[`j',2]
  }
}
dis "MULTIPLICITY_PRED_ratio(k=3,m=15,G=45) = " sqrt(3*15*2/(44*45))

dis "=== P2c: CONTINUOUS TREATMENT ame, impvalues cluster(unit) ==="
capture noisily wooldid ycont cg t ttre, att cluster(unit) impvalues contreat(dose)
dis "P2c_RC=" _rc
capture noisily matrix list e(wdidmainresults)
capture matrix CU = e(wdidmainresults)

dis "=== P2d: CONTINUOUS TREATMENT ame, impvalues cluster(cg) ==="
capture noisily wooldid ycont cg t ttre, att cluster(cg) impvalues contreat(dose)
dis "P2d_RC=" _rc
capture noisily matrix list e(wdidmainresults)
capture matrix CG = e(wdidmainresults)

capture {
  local nc = rowsof(CU)
  forvalues j = 1/`nc' {
    local rn : word `j' of `: rownames CU'
    dis "CROW `j' `rn' : se_u=" %10.7f CU[`j',2] " se_g=" %10.7f CG[`j',2] ///
        " ratio=" %8.5f CU[`j',2]/CG[`j',2]
  }
}
dis "P2_TAIL_RC=" _rc
dis "DONE_SCOPE"
