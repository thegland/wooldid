* ============================================================================
* R3_A : SKEPTIC's OWN repro of seeded finding fp-contreat-A
*   claim: contreatelasticitytype(eydx)/(eyex) on the default reghdfe path
*   divide by reghdfe's `xb', which EXCLUDES the absorbed fixed effects, so the
*   reported elasticity is not "% change in y" at all.
*
* Design of the test (independent of the prober's):
*   - noisy panel; unit/cohort effects carry most of the LEVEL of y
*   - dydx legs first: if the two engines agree on dydx they are fitting the
*     same model, so any eydx disagreement is purely the DENOMINATOR
*   - outside-wooldid oracles that mimic wooldid's own spellings exactly
*     (a `ringer' copy of the dose, so margins cannot see the pooled slope)
*   - closed-form target: mean over treated post rows of b_cell / yhat_full
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R3_A_elasticity.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

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

* dose defined everywhere (so the pooled `basic' control is a real control)
gen double dose = 1 + 2*runiform()
gen double bit = 0
qui replace bit = 0.60 + 0.30*i - 0.09*(t-ttre) if post == 1
gen double kit = 0
qui replace kit = 0.40 + 0.15*i + 0.05*(t-ttre) if post == 1
bysort unit (t): gen double aunit = 20 + 6*runiform() if _n == 1
bysort unit (t): replace aunit = aunit[1]
gen double y = aunit + 0.35*t + 0.9*dose + kit + bit*dose + rnormal(0, 0.5)
gen double ringer = dose            /* bit-exact copy, as wooldid builds */
gen byte nullv = 1                  /* wooldid's one-level `null' variable */
egen cell = group(i t) if post == 1
qui replace cell = 0 if cell >= .
tempfile F
qui save `F'
levelsof cell if post == 1, local(CL)

* ---------- ORACLE H : reghdfe spelling (ccnoabsorb-equivalent) ------------
use `F', clear
qui reghdfe y ib0.cell##c.dose , absorb(i.i i.t i.nullv##c.ringer) vce(robust)
qui margins , dydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "ORACLE H  dydx = " %18.12f r[1,1]
qui margins , eydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "ORACLE H  eydx = " %18.12f r[1,1]
qui margins , eyex(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "ORACLE H  eyex = " %18.12f r[1,1]
* what reghdfe's own xb prediction looks like
qui predict double xbH , xb
qui sum xbH
dis as result "ORACLE H  xb: mean " %14.6f r(mean) "  sd " %14.6f r(sd)
qui sum y
dis as result "         y : mean " %14.6f r(mean) "  sd " %14.6f r(sd)

* ---------- ORACLE R : plain regress spelling (customvce-equivalent) -------
use `F', clear
qui reg y i.i i.t i.nullv##c.ringer ib0.cell##c.dose , vce(robust)
qui margins , dydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "ORACLE R  dydx = " %18.12f r[1,1]
qui margins , eydx(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "ORACLE R  eydx = " %18.12f r[1,1]
qui margins , eyex(dose) subpop(if post == 1) noestimcheck
mat r = r(b)
dis as result "ORACLE R  eyex = " %18.12f r[1,1]
* closed form: mean of b_cell / full prediction over treated post rows
qui predict double yhat , xb
gen double bc = 0
foreach c of local CL {
  qui replace bc = _b[`c'.cell#c.dose] if cell == `c'
}
gen double ratio = bc/yhat
qui sum ratio if post == 1
dis as result "CLOSED FORM mean(b_cell/yhat_full) = " %18.12f r(mean)
gen double ratio2 = bc*dose/yhat
qui sum ratio2 if post == 1
dis as result "CLOSED FORM mean(b_cell*dose/yhat_full) = " %18.12f r(mean)

* ---------- wooldid ------------------------------------------------------
foreach et in dydx dyex eydx eyex {
  foreach spec in "robust" "ccnoabsorb robust" "customvce(robust)" {
    use `F', clear
    cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) ///
        contreatelasticitytype(`et') `spec' regtol(14)
    dis "RC=" _rc " | `et' | `spec'"
    if _rc == 0 {
      mat M = e(wdidmainresults)
      dis as result "WD | `et' | `spec' | ame = " %18.12f M[2,1] " se = " %14.8f M[2,2]
    }
  }
}

* ---------- binary-side semielasticity (advisor: do not inherit) ----------
foreach spec in "robust" "customvce(robust)" {
  use `F', clear
  cap noi qui wooldid y i t ttre , att semielasticity `spec' regtol(14)
  dis "RC=" _rc " | semielasticity | `spec'"
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "WD | semielast | `spec' | att = " %18.12f M[1,1]
  }
}
* binary oracles
use `F', clear
qui reghdfe y ib0.cell , absorb(i.i i.t) vce(robust)
dis as result "(binary oracle fits omitted - see R5)"

dis "DONE_R3"
