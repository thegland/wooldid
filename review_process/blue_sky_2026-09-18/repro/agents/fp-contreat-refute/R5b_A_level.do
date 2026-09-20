* R5b: the fp-contreat-A gap as a function of how big E[y] is relative to the
* dispersion the absorbed effects carry.  Same panel, only the base level of y
* moves.  Realistic analogue: a rate or a per-capita outcome with a small mean
* and large between-unit dispersion.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909090
set sortseed 909090

foreach base in 40 20 10 6 {
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
  bysort unit (t): gen double aunit = 6*runiform() if _n == 1
  bysort unit (t): replace aunit = aunit[1]
  gen double y = `base' + aunit + 0.35*t + 0.9*dose + kit + bit*dose + rnormal(0, 0.5)
  tempfile F
  qui save `F'
  qui sum y
  local YM = r(mean)
  local YS = r(sd)

  local D_ .
  local C_ .
  foreach spec in "robust" "customvce(robust)" {
    use `F', clear
    cap noi qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) ///
        contreatelasticitytype(eydx) fe(i.unit) `spec' regtol(14)
    if _rc == 0 {
      mat M = e(wdidmainresults)
      if "`spec'" == "robust" local D_ = M[2,1]
      else local C_ = M[2,1]
    }
    else dis as error "RCFAIL=" _rc " base=`base' `spec'"
  }
  local RD = reldif(`D_',`C_')
  dis as result "base=`base' ymean=" %8.3f `YM' " ysd=" %8.3f `YS' ///
     " | eydx reghdfe=" %16.10f `D_' " regress=" %16.10f `C_' " REL=" %9.4f `RD'
}
dis "DONE_R5b"
