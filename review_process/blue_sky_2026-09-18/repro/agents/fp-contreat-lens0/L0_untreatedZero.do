*========================================================================
* LENS 0 / INDEPENDENT REPRODUCTION of fp-contreat-E
* Built from the CLAIM alone, on my own DGP, with my own controls.
*========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
set maxvar 20000

dis "=== BUILD ASSERT ==="
which wooldid

*---------------------------------------------------------------
* helper: print the main results table in a greppable tagged form
*---------------------------------------------------------------
cap prog drop grab
prog def grab
  args tag
  cap mat R = e(wdidmainresults)
  if _rc {
    dis "RES `tag' | NOMATRIX"
    exit
  }
  local nr = rowsof(R)
  local nc = colsof(R)
  local rn : rownames R
  local cn : colnames R
  dis "COLS `tag' | `cn'"
  forvalues r = 1/`nr' {
    local nm : word `r' of `rn'
    local line "RES `tag' | `nm'"
    forvalues c = 1/`nc' {
      local v = string(R[`r',`c'],"%20.10f")
      local line "`line' | `v'"
    }
    dis "`line'"
  }
end

*===============================================================
* DESIGN H : untreatedZero premise, HOMOGENEOUS dose effect
*   truth: ame = 0.80 under ANY cell weighting
*===============================================================
cap prog drop mkpanel
prog def mkpanel
  args nsub nper noisesd
  clear
  set obs 8
  gen grp = _n
  gen ttre = .
  replace ttre = 4 in 1
  replace ttre = 5 in 2
  replace ttre = 6 in 3
  replace ttre = 7 in 4
  replace ttre = 5 in 5
  * groups 6,7,8 never treated
  expand `nsub'
  bysort grp: gen sub = _n
  expand `nper'
  bysort grp sub: gen t = _n
  gen byte post = (ttre < . & t >= ttre)
  gen relyr = t - ttre
  gen double dose = 0
  replace dose = 0.5 + 2.5*runiform() if post
  gen double ag = 1.5*grp
  gen double ft = 0.7*t
  gen double cellint = 0
  replace cellint = 0.5*grp + 0.3*relyr if post
  gen double noise = `noisesd'*rnormal()
end

*---- H exact fit ---------------------------------------------
mkpanel 15 9 0
gen double y = ag + ft + cellint + 0.80*dose + noise
dis "=== H panel built ==="
count
count if post
sum dose if post
sum dose if !post

dis "HTRUTH | homogeneous beta = 0.8000000000 ; any weighting gives 0.8"

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) robust
dis "H1_RC=" _rc
grab H1_default_robust

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) ccnoabsorb robust
dis "H2_RC=" _rc
grab H2_ccnoabsorb

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) customvce(robust)
dis "H3_RC=" _rc
grab H3_customvce_rob

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) customvce(cluster grp)
dis "H4_RC=" _rc
grab H4_customvce_clu

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) unconditionalse robust
dis "H5_RC=" _rc
grab H5_uncondse

*---- H noisy (realism: does the wrong answer look ordinary?) ---
mkpanel 15 9 0.30
gen double y = ag + ft + cellint + 0.80*dose + noise
dis "HNOISY | truth still 0.80 under any weighting"

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) robust
dis "HN1_RC=" _rc
grab HN1_default_robust

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) customvce(robust)
dis "HN2_RC=" _rc
grab HN2_customvce_rob

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) customvce(cluster grp)
dis "HN3_RC=" _rc
grab HN3_customvce_clu

*===============================================================
* DESIGN X : heterogeneous, DISTINCT beta per (grp,relyr) cell
*   the gap default-minus-customvce should equal exactly one of them
*===============================================================
mkpanel 15 9 0
gen double bcell = 0
replace bcell = 0.50 + 0.13*grp + 0.07*relyr if post
gen double y = ag + ft + cellint + bcell*dose + noise

dis "=== X: distinct true cell slopes ==="
preserve
  keep if post
  bysort grp relyr: keep if _n==1
  sort grp relyr
  list grp relyr bcell, clean noobs
restore

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) robust
dis "X1_RC=" _rc
grab X1_default_robust
cap mat bX1 = e(b)

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) customvce(robust)
dis "X2_RC=" _rc
grab X2_customvce_rob

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) unconditionalse robust
dis "X3_RC=" _rc
grab X3_uncondse

*===============================================================
* DESIGN B : CONTROL. dose defined everywhere, contreatcontrols(basic).
*   default and customvce MUST agree here if customvce only changes the VCE
*===============================================================
clear
set obs 8
gen grp = _n
gen ttre = .
replace ttre = 4 in 1
replace ttre = 5 in 2
replace ttre = 6 in 3
replace ttre = 7 in 4
replace ttre = 5 in 5
expand 15
bysort grp: gen sub = _n
expand 9
bysort grp sub: gen t = _n
gen byte post = (ttre < . & t >= ttre)
gen relyr = t - ttre
gen double dose = 0.5 + 2.5*runiform()
gen double y = 1.5*grp + 0.7*t + 0.4*dose
replace y = y + 0.5*grp + 0.3*relyr + 0.80*dose if post

dis "=== B control design: dose nonzero everywhere, basic controls ==="
cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(basic) robust
dis "B1_RC=" _rc
grab B1_default_robust

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(basic) customvce(robust)
dis "B2_RC=" _rc
grab B2_customvce_rob

cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(basic) unconditionalse robust
dis "B3_RC=" _rc
grab B3_uncondse

*===============================================================
* DESIGN P : poisson. untreatedZero premise, count-style outcome.
*===============================================================
mkpanel 15 9 0
gen double y = exp(0.10*grp + 0.05*t + 0.20*post + 0.25*dose)
dis "=== P homogeneous: true dose semi-elasticity = 0.25 everywhere ==="
cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) poisson robust
dis "P1_RC=" _rc
grab P1_poisson_homog

mkpanel 15 9 0
gen double bcell = 0
replace bcell = 0.10 + 0.04*grp + 0.02*relyr if post
gen double y = exp(0.10*grp + 0.05*t + 0.20*post + bcell*dose)
dis "=== P heterogeneous: true semi-elasticities vary by cell ==="
preserve
  keep if post
  bysort grp relyr: keep if _n==1
  sort grp relyr
  list grp relyr bcell, clean noobs
restore
cap noi wooldid y grp t ttre, contreat(dose) contreatcontrols(untreatedZero) poisson robust
dis "P2_RC=" _rc
grab P2_poisson_het

dis "=== L0 DONE ==="
