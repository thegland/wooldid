* ============================================================================
* R1_E : SKEPTIC's OWN repro of seeded finding fp-contreat-E
*   claim: contreatcontrols(untreatedZero) silently gains the pooled baseline
*   dose control on the customvce()/unconditionalse/poisson paths, so the ame
*   is wrong there while the default reghdfe path is right.
*
* Differences from the prober's repro (deliberate, to kill harness artifacts):
*   - NOISY data (their design was exact-fit, R2=1)
*   - 6 cohort groups x 12 units x 10 periods, 3 never-treated groups
*   - independent OLS oracles built by hand OUTSIDE wooldid for BOTH engines
*   - the identity  ame_customvce == ame_oracle - bhat(dropped cell)  is tested
*     by reading the omitted column out of my own regress fit
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R1_E_untreatedZero.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771102
set sortseed 771102

which wooldid

* ---------------- data: realistic-ish messy panel -------------------------
clear
set obs 720                        /* 6 groups x 12 units x 10 periods */
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/12)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
qui replace ttre = 8 if i == 3
* i = 4,5,6 never treated

gen byte post = (ttre < . & t >= ttre)
* dose is EXACTLY zero off treated cohort-periods (the untreatedZero premise)
gen double dose = 0
qui replace dose = 0.5 + 3*runiform() if post == 1

* heterogeneous cell-level dose slope and cell-level intercept shift
gen double bit = 0
qui replace bit = 0.60 + 0.30*i - 0.09*(t-ttre) if post == 1
gen double kit = 0
qui replace kit = 0.40 + 0.15*i + 0.05*(t-ttre) if post == 1

* unit and time effects carry most of the LEVEL of y (as in real panels)
bysort unit (t): gen double aunit = 10 + 4*runiform() if _n == 1
bysort unit (t): replace aunit = aunit[1]
gen double gt = 0.35*t
gen double y = aunit + gt + kit + bit*dose + rnormal(0, 0.5)

* treated (i,t) cell ids; 0 on every untreated/reference row
egen cell = group(i t) if post == 1
qui replace cell = 0 if cell >= .

tempfile F
qui save `F'

qui sum bit if post == 1
local DGPTRUE = r(mean)
dis as result "DGP truth (mean beta_it over treated post rows) = " %18.10f `DGPTRUE'

* ---------------- ORACLE 1: the absorbing fit (what the default path does) --
use `F', clear
qui reghdfe y ib0.cell##c.dose , absorb(i.i i.t) vce(robust)
* att-weighted ame = observation-weighted mean of the cell dose slopes
gen double bhat = 0
levelsof cell if post == 1, local(CL)
foreach c of local CL {
  qui replace bhat = _b[`c'.cell#c.dose] if cell == `c'
}
qui sum bhat if post == 1
local OR_HDFE = r(mean)
dis as result "ORACLE reghdfe  ame = " %18.10f `OR_HDFE'

* ---------------- ORACLE 2: the plain-regress fit WITH the pooled dose control
* this is what wooldid builds on the customvce/unconditionalse path: the pooled
* baseline dose slope is pasted in FRONT of the cell interactions
use `F', clear
qui reg y i.i i.t c.dose ib0.cell##c.dose , vce(robust)
gen double bhat2 = 0
local DROPPED
local DROPB .
foreach c of local CL {
  local bc = _b[`c'.cell#c.dose]
  qui replace bhat2 = `bc' if cell == `c'
  if `bc' == 0 local DROPPED `DROPPED' `c'
}
qui sum bhat2 if post == 1
local OR_REG = r(mean)
local POOLED = _b[dose]
dis as result "ORACLE reg+pooled ame = " %18.10f `OR_REG'
dis as result "   pooled dose slope = " %18.10f `POOLED'
dis as result "   cell(s) with a ZERO slope (dropped for collinearity): `DROPPED'"

* what the true slope of the dropped cell is, per ORACLE 1
use `F', clear
qui reghdfe y ib0.cell##c.dose , absorb(i.i i.t) vce(robust)
foreach c of local DROPPED {
  dis as result "   ORACLE1 slope of dropped cell `c' = " %18.10f _b[`c'.cell#c.dose]
  local BDROP = _b[`c'.cell#c.dose]
}
dis as result "   IDENTITY CHECK  ame_oracle1 - bdrop = " %18.10f `OR_HDFE' - `BDROP'

* ---------------- wooldid, every engine -----------------------------------
foreach spec in "robust" "ccnoabsorb robust" "customvce(robust)" "customvce(cluster unit)" "unconditionalse robust" {
  use `F', clear
  cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) `spec' regtol(14)
  dis "RC_`spec'=" _rc
  if _rc == 0 {
    mat M = e(wdidmainresults)
    local rn : rowfullnames M
    dis as result "WD | `spec' | rows: `rn'"
    dis as result "WD | `spec' | ame = " %18.10f M[2,1] " se = " %14.8f M[2,2] ///
        " | vs ORACLE1 " %18.10f `OR_HDFE' " | vs ORACLE2 " %18.10f `OR_REG'
  }
}

dis "DONE_R1"
