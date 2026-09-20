* ============================================================================
* R5_A : how big is the fp-contreat-A elasticity error on a REALISTIC panel?
*  The gap is driven by how much of y's LEVEL the absorbed effects carry.
*  Here unit fixed effects are absorbed through fe(i.unit) - the ordinary
*  applied choice - and the unit effects carry most of the dispersion of y.
*  fe(i.unit) (not fe(unit)) is used on BOTH legs, as sthlp:1051-1053 instructs.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do R5_A_magnitude.do
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606060
set sortseed 606060

* unitsd = dispersion of the absorbed unit effects; sweep it
foreach usd in 0 1 4 8 {
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
  gen double dose = 1 + 2*runiform()
  gen double bit = 0
  qui replace bit = 0.60 + 0.30*i - 0.09*(t-ttre) if post == 1
  gen double kit = 0
  qui replace kit = 0.40 + 0.15*i + 0.05*(t-ttre) if post == 1
  bysort unit (t): gen double aunit = `usd'*runiform() if _n == 1
  bysort unit (t): replace aunit = aunit[1]
  gen double y = 25 + aunit + 0.35*t + 0.9*dose + kit + bit*dose + rnormal(0, 0.5)
  egen cell = group(i t) if post == 1
  qui replace cell = 0 if cell >= .
  tempfile F
  qui save `F'

  local D_ .
  local C_ .
  local DY .
  local CY .
  foreach spec in "robust" "customvce(robust)" {
    use `F', clear
    cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) ///
        contreatelasticitytype(eydx) fe(i.unit) `spec' regtol(14)
    if _rc == 0 {
      mat M = e(wdidmainresults)
      if "`spec'" == "robust" local D_ = M[2,1]
      else local C_ = M[2,1]
    }
    else dis as error "RC=" _rc " eydx `spec' usd=`usd'"
    use `F', clear
    cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) ///
        contreatelasticitytype(dydx) fe(i.unit) `spec' regtol(14)
    if _rc == 0 {
      mat M = e(wdidmainresults)
      if "`spec'" == "robust" local DY = M[2,1]
      else local CY = M[2,1]
    }
  }
  dis as result "usd=`usd' | dydx: reghdfe " %16.10f `DY' "  regress " %16.10f `CY' ///
     "  reldif " %10.3e reldif(`DY',`CY')
  dis as result "usd=`usd' | eydx: reghdfe " %16.10f `D_' "  regress " %16.10f `C_' ///
     "  reldif " %10.3e reldif(`D_',`C_')
}
dis "DONE_R5"
