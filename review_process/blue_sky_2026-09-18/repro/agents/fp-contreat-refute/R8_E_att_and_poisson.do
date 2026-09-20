* ============================================================================
* R8 : two follow-ups
*   (i)  does finding E's inserted pooled dose control also move the ATT row
*        (not just the ame)?  -> severity of E
*   (ii) does the ppmlhdfe rc-3200 fallback fire on EVERY poisson + contreat()
*        run, whatever contreatcontrols()?  The absorbed term wooldid builds is
*        i.<constant>##c.(ringer) - a one-level factor - which is exactly the
*        condition the fallback's own caution names.  If it always fires, the
*        caution blames the user's fe() for wooldid's own construction.
*   (iii) does fe(i.unit) work on the customvce path (finding D's documented
*        remedy)?
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R8_E_att_and_poisson.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771102
set sortseed 771102

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

* true ATT on the treated = mean over treated post rows of (kit + bit*dose)
qui gen double teff = kit + bit*dose
qui sum teff if post == 1
dis as result "DGP true ATT = " %16.10f r(mean)

dis as result "----- (i) untreatedZero: does the ATT row move too? -----"
foreach spec in "robust" "customvce(robust)" "unconditionalse robust" {
  use `F', clear
  cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) `spec' regtol(14)
  dis "RC=" _rc " `spec'"
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "  UZ `spec' | att = " %16.10f M[1,1] " (se " %10.6f M[1,2] ") | ame = " %16.10f M[2,1]
  }
}

dis as result "----- (ii) poisson + contreat: does the rc-3200 fallback always fire? -----"
use `F', clear
qui gen double xbp = 0.5 + 0.10*i + 0.04*t
qui replace xbp = xbp + 0.10 + 0.20*dose if post == 1
qui gen double ypo = rpoisson(exp(xbp))
qui replace ypo = 0 if ypo >= .
* a dose defined everywhere, so `basic' is a legitimate control
qui gen double dose2 = 1 + 2*runiform()
qui gen double xbq = 0.5 + 0.10*i + 0.04*t + 0.15*dose2
qui replace xbq = xbq + 0.10 + 0.20*dose2 if post == 1
qui gen double ypo2 = rpoisson(exp(xbq))
qui replace ypo2 = 0 if ypo2 >= .
tempfile P
qui save `P'

dis as result "  ..... poisson + contreatcontrols(basic) (dose defined everywhere) ....."
use `P', clear
cap noi wooldid ypo2 i t ttre , att contreat(dose2) contreatcontrols(basic) poisson robust regtol(9)
dis "RC_PB=" _rc

dis as result "  ..... poisson + contreatcontrols(byTime) ....."
use `P', clear
cap noi qui wooldid ypo2 i t ttre , att contreat(dose2) contreatcontrols(byTime) poisson robust regtol(9)
dis "RC_PT=" _rc

dis as result "  ..... plain binary poisson (no contreat) ....."
use `P', clear
cap noi qui wooldid ypo2 i t ttre , att poisson robust regtol(9)
dis "RC_PBIN=" _rc

dis as result "----- (iii) fe(i.unit) on the customvce path (D's documented remedy) -----"
use `F', clear
cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) fe(i.unit) customvce(cluster unit) regtol(14)
dis "RC_FEI=" _rc
if _rc == 0 {
  dis as result "  fe(i.unit) customvce: R2 = " %10.6f e(r2)
  mat M = e(wdidmainresults)
  dis as result "  fe(i.unit) customvce: ame = " %16.10f M[2,1]
}
use `F', clear
cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) fe(unit) customvce(cluster unit) regtol(14)
dis "RC_FEBARE=" _rc
if _rc == 0 {
  dis as result "  fe(unit) BARE customvce: R2 = " %10.6f e(r2)
  mat M = e(wdidmainresults)
  dis as result "  fe(unit) BARE customvce: ame = " %16.10f M[2,1]
}

dis "DONE_R8"
