* ============================================================================
* LENS 4 (documented-behaviour) adjudication of fp-contreat-D
*
* Question the lens must settle: fe() under customvce()/unconditionalse enters
* the regression as a BARE token (continuous).  Is that a defect, or documented
* design?  sthlp:1052 documents it verbatim with the same example, and
* review_process/rulings_2026-09-07.md BS-1 is a documentation-only ruling.
*
* What is NOT yet shown by either the prober or the skeptic: whether the
* DOCUMENTED REMEDY actually restores the absorbed answer.  If fe(i.unit) under
* customvce() recovers what fe(unit) under reghdfe recovers, the capability is
* intact and the finding is refuted.  If it does NOT, there is a live residual.
*
* DGP: the prober's own probe08 fixture, regenerated here so this file is
* self-contained.  y is exactly linear in (unit FE, t, dose, treatment), so the
* absorbed fit is exact and the true average marginal effect of dose on the
* treated is known in closed form.
*
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do L4_fe_remedy.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242

dis "WHICH_BUILD:"
which wooldid

clear
set obs 400
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/10)
gen double ttre = .
qui replace ttre = 5 if i == 1
qui replace ttre = 7 if i == 2
bysort unit (t): gen double ufe = 30 + 20*runiform() if _n==1
bysort unit (t): replace ufe = ufe[1]
gen double dose = 1 + 2*runiform()
gen double kappatrue = 0
gen double betatrue  = 0
qui replace kappatrue = 1.0 + 0.2*(t-ttre) if ttre < . & t >= ttre
qui replace betatrue  = 2.0 - 0.3*(t-ttre) if ttre < . & t >= ttre
gen double y = ufe + 0.5*t + 0.9*dose + kappatrue + betatrue*dose
tempfile F
qui save `F'

* closed-form truth: average of (0.9 + betatrue) over treated post rows
qui gen byte post = (ttre < . & t >= ttre)
qui gen double dydose = 0.9 + betatrue
qui sum dydose if post == 1
dis "TRUTH_AME=" %20.12f r(mean)

* ---------------------------------------------------------------- reference
dis "===== REF: fe(unit) + cluster(unit)  [reghdfe, genuinely absorbed] ====="
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) fe(unit) cluster(unit) regtol(14)
dis "RC_REF=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis "REF_att=" %20.12f M[1,1] "  REF_ame=" %20.12f M[2,1] "  REF_r2=" %12.8f e(r2)
}

* ------------------------------------------------- the bare token (finding D)
dis "===== BARE: fe(unit) + customvce(cluster unit)  [documented continuous] ====="
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) fe(unit) customvce(cluster unit) regtol(14)
dis "RC_BARE=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis "BARE_att=" %20.12f M[1,1] "  BARE_ame=" %20.12f M[2,1] "  BARE_r2=" %12.8f e(r2)
}

* -------------------------------------------- the HELP's remedy: fe(i.unit)
dis "===== HELPFIX: fe(i.unit) + customvce(cluster unit) ====="
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) fe(i.unit) customvce(cluster unit) regtol(14)
dis "RC_HELPFIX=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis "HELPFIX_att=" %20.12f M[1,1] "  HELPFIX_ame=" %20.12f M[2,1] "  HELPFIX_r2=" %12.8f e(r2)
}

* ------------- the RUNTIME CAUTION's remedy: move it to controls with i.
dis "===== CTRLFIX: controls(i.unit) + customvce(cluster unit), no fe() ====="
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) controls(i.unit) customvce(cluster unit) regtol(14)
dis "RC_CTRLFIX=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis "CTRLFIX_att=" %20.12f M[1,1] "  CTRLFIX_ame=" %20.12f M[2,1] "  CTRLFIX_r2=" %12.8f e(r2)
}

* ----------------------------- same question on the unconditionalse spelling
dis "===== HELPFIX_UNC: fe(i.unit) + unconditionalse cluster(unit) ====="
use `F', clear
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) fe(i.unit) unconditionalse cluster(unit) regtol(14)
dis "RC_HELPFIX_UNC=" _rc
if _rc == 0 {
  mat M = e(wdidmainresults)
  dis "HELPFIXUNC_att=" %20.12f M[1,1] "  HELPFIXUNC_ame=" %20.12f M[2,1] "  HELPFIXUNC_r2=" %12.8f e(r2)
}

* --- does the caution still fire (and is it still accurate) with i. already written?
dis "===== CAUTION TEXT CHECK: fe(i.unit) run above should still print the caution ====="
dis "DONE"
