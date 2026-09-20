* ============================================================================
* R2_E : SKEPTIC's OWN repro of seeded finding fp-contreat-E  (corrected oracle)
*
* THE MODEL wooldid is supposed to fit under contreatcontrols(untreatedZero):
*     y = FE_i + FE_t + sum_c [ kappa_c * 1{cell=c} + beta_c * 1{cell=c} * dose ]
* i.e. NO pooled baseline dose slope, because the dose is identically 0 off the
* treated cohort-periods so a pooled slope is exactly collinear with the cell
* slopes.  ado:1311-1315 says so in as many words.
*
* ORACLE A (no pooled dose)   : reghdfe y ib0.cell ib0.cell#c.dose, absorb(i.i i.t)
* ORACLE B (WITH pooled dose) : reg    y i.i i.t c.dose ib0.cell ib0.cell#c.dose
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R2_E_untreatedZero.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771102
set sortseed 771102

which wooldid

* ---------------- data: noisy, realistic-shape panel ----------------------
clear
set obs 720                        /* 6 groups x 12 units x 10 periods */
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/12)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
qui replace ttre = 8 if i == 3

gen byte post = (ttre < . & t >= ttre)
gen double dose = 0
qui replace dose = 0.5 + 3*runiform() if post == 1

gen double bit = 0
qui replace bit = 0.60 + 0.30*i - 0.09*(t-ttre) if post == 1
gen double kit = 0
qui replace kit = 0.40 + 0.15*i + 0.05*(t-ttre) if post == 1

bysort unit (t): gen double aunit = 10 + 4*runiform() if _n == 1
bysort unit (t): replace aunit = aunit[1]
gen double y = aunit + 0.35*t + kit + bit*dose + rnormal(0, 0.5)

egen cell = group(i t) if post == 1
qui replace cell = 0 if cell >= .

tempfile F
qui save `F'

qui sum bit if post == 1
local DGPTRUE = r(mean)
dis as result "DGP truth (mean beta_it over treated post rows) = " %18.10f `DGPTRUE'
levelsof cell if post == 1, local(CL)

* ---------------- ORACLE A : the intended model ---------------------------
use `F', clear
qui reghdfe y ib0.cell ib0.cell#c.dose , absorb(i.i i.t) vce(robust)
gen double bA = 0
foreach c of local CL {
  qui replace bA = _b[`c'.cell#c.dose] if cell == `c'
}
qui sum bA if post == 1
local OR_A = r(mean)
dis as result "ORACLE A (no pooled dose)  ame = " %18.10f `OR_A'

* ---------------- ORACLE B : the model the reg path actually builds -------
use `F', clear
qui reg y i.i i.t c.dose ib0.cell ib0.cell#c.dose , vce(robust)
gen double bB = 0
local DROPPED
foreach c of local CL {
  local bc = _b[`c'.cell#c.dose]
  qui replace bB = `bc' if cell == `c'
  if `bc' == 0 local DROPPED `DROPPED' `c'
}
qui sum bB if post == 1
local OR_B = r(mean)
local POOLED = _b[dose]
dis as result "ORACLE B (WITH pooled dose) ame = " %18.10f `OR_B'
dis as result "   pooled dose slope = " %18.10f `POOLED' "   dropped cell(s): `DROPPED'"

* identity: ORACLE B ame == ORACLE A ame - (ORACLE A slope of the dropped cell)
use `F', clear
qui reghdfe y ib0.cell ib0.cell#c.dose , absorb(i.i i.t) vce(robust)
local BDROP .
foreach c of local DROPPED {
  local BDROP = _b[`c'.cell#c.dose]
}
dis as result "   ORACLE A slope of dropped cell = " %18.10f `BDROP'
dis as result "   IDENTITY  A_ame - b_dropped     = " %18.10f `OR_A' - `BDROP' ///
              "   (should equal ORACLE B ame)"

* ---------------- wooldid, every engine -----------------------------------
foreach spec in "robust" "ccnoabsorb robust" "customvce(robust)" "customvce(cluster unit)" "unconditionalse robust" "cluster(unit)" {
  use `F', clear
  cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) `spec' regtol(14)
  dis "RC=" _rc " for `spec'"
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "WD | `spec' | ame = " %18.10f M[2,1] " se = " %14.8f M[2,2] ///
        " | A " %18.10f `OR_A' " | B " %18.10f `OR_B'
  }
}

* ---------------- the cautions the user actually sees (NOT quiet) ---------
dis as result "===== NOISY customvce(cluster unit) leg ====="
use `F', clear
cap noi wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) customvce(cluster unit) regtol(14)
dis "RC_NOISY=" _rc

* ---------------- POISSON legs -------------------------------------------
* homogeneous true semi-elasticity (prober predicts exactly 0 is reported)
use `F', clear
qui gen double xbh = 0.5 + 0.10*i + 0.04*t
qui replace xbh = xbh + 0.10 + 0.20*dose if post == 1
qui gen double ypoh = rpoisson(exp(xbh))
qui replace ypoh = 0 if ypoh >= .
tempfile PH
qui save `PH'

* heterogeneous true semi-elasticity (prober's mechanism predicts 0 too;
* the drop mechanism predicts mean(gamma) - gamma_dropped)
use `F', clear
qui gen double gam = 0
qui replace gam = 0.10 + 0.06*i - 0.02*(t-ttre) if post == 1
qui gen double xbx = 0.5 + 0.10*i + 0.04*t
qui replace xbx = xbx + 0.10 + gam*dose if post == 1
qui gen double ypox = rpoisson(exp(xbx))
qui replace ypox = 0 if ypox >= .
qui sum gam if post == 1
local GTRUE = r(mean)
tempfile PX
qui save `PX'
dis as result "POISSON heterogeneous DGP mean semi-elasticity = " %18.10f `GTRUE'

dis as result "===== poisson HOMOGENEOUS (truth 0.20), NOISY ====="
use `PH', clear
cap noi wooldid ypoh i t ttre , att contreat(dose) contreatcontrols(untreatedZero) poisson robust regtol(9)
dis "RC_PH=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis as result "WD-POIS-HOM | ame = " %18.10f M[2,1] " se = " %14.8f M[2,2]
}

dis as result "===== poisson HETEROGENEOUS, quiet ====="
use `PX', clear
cap noi qui wooldid ypox i t ttre , att contreat(dose) contreatcontrols(untreatedZero) poisson robust regtol(9)
dis "RC_PX=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis as result "WD-POIS-HET | ame = " %18.10f M[2,1] " se = " %14.8f M[2,2] " truth " %12.6f `GTRUE'
}

* poisson oracles, built by hand
use `PX', clear
cap noi qui ppmlhdfe ypox ib0.cell ib0.cell#c.dose , absorb(i.i i.t) vce(robust) d
dis "RC_PPMLA=" _rc
if _rc == 0 {
  gen double gA = 0
  foreach c of local CL {
    qui replace gA = _b[`c'.cell#c.dose] if cell == `c'
  }
  qui sum gA if post == 1
  dis as result "ORACLE-POIS A (no pooled dose)  ame = " %18.10f r(mean)
}
use `PX', clear
cap noi qui ppmlhdfe ypox c.dose ib0.cell ib0.cell#c.dose , absorb(i.i i.t) vce(robust) d
dis "RC_PPMLB=" _rc
if _rc == 0 {
  gen double gB = 0
  foreach c of local CL {
    qui replace gB = _b[`c'.cell#c.dose] if cell == `c'
  }
  qui sum gB if post == 1
  dis as result "ORACLE-POIS B (WITH pooled dose) ame = " %18.10f r(mean)
}

dis "DONE_R2"
