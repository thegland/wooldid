clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"

* ---------- realistic NOISY linear design (same seed as run 2 part D) ----------
clear
set seed 918273
set sortseed 918273
set obs 2160
gen long unit = ceil(_n/9)
gen int  t    = _n - (unit-1)*9
gen byte i    = ceil(unit/40)
gen double ttre = .
qui replace ttre = 3 if i == 1
qui replace ttre = 5 if i == 2
qui replace ttre = 7 if i == 3
gen byte tc = (ttre < . & t >= ttre)
gen double r = t - ttre if tc == 1
gen double dose = 0
qui replace dose = exp(rnormal(0.2,0.55)) if tc == 1
gen double beta = 0
qui replace beta = 1.45 - 0.31*i + 0.22*r - 0.06*r*r if tc == 1
gen double kappa = 0
qui replace kappa = 0.9 - 0.15*i + 0.13*r if tc == 1
by unit, sort: gen double ue = rnormal(0,0.6) if _n == 1
by unit: replace ue = ue[1]
gen double y = (1.1 + 0.47*i) + (-0.3 + 0.21*t - 0.012*t*t) + kappa + beta*dose + ue + rnormal(0,0.4)
qui sum beta if tc == 1
local TRUE2 = r(mean)

qui gen double _tt = t - ttre + 1 if ttre < . & t >= ttre
qui replace _tt = 0 if _tt == .
qui gen double _ii = i if _tt > 0
qui replace _ii = 0 if _ii == .
qui gen byte tf = _tt > 0
qui gen double ringer = dose
qui gen byte nullv = 1
tempfile G
qui save `G'

dis as result "TRUE (population) ame = " %14.8f `TRUE2'

* ORACLE A: the CORRECT untreatedZero spec, no ringer anywhere
qui reghdfe y i._ii#i._tt#i.tf#c.dose i._ii#i._tt#i.tf , absorb(i.i i.t) vce(cluster unit) tol(1e-12)
local sum = 0
local cnt = 0
qui levelsof _ii if _ii > 0, local(il)
foreach a of local il {
  forval b = 1/7 {
    cap local cf = _b[`a'._ii#`b'._tt#1.tf#c.dose]
    if _rc == 0 & abs(`cf') > 0 {
      local sum = `sum' + `cf'
      local cnt = `cnt' + 1
      if `a' == 3 & `b' == 3 local BHATDROP = `cf'
    }
  }
}
dis as result "ORACLE A (correct spec): mean of `cnt' cell slopes = " %14.8f `sum'/`cnt'
dis as result "   estimated slope of the cell wooldid drops (i=3,relyr=3) = " %14.8f `BHATDROP'
dis as result "   => predicted broken value = " %14.8f `sum'/`cnt' - `BHATDROP'

* ORACLE B: the SAME spec but WITH the pooled ringer (what customvce/uncond do)
qui reg y i.i i.t i.nullv##c.(ringer) i._ii#i._tt#i.tf#c.dose i._ii#i._tt#i.tf , vce(cluster unit)
local sum2 = 0
local cnt2 = 0
foreach a of local il {
  forval b = 1/7 {
    cap local cf = _b[`a'._ii#`b'._tt#1.tf#c.dose]
    if _rc == 0 {
      qui count if _ii == `a' & _tt == `b'
      if r(N) > 0 {
        local sum2 = `sum2' + `cf'
        local cnt2 = `cnt2' + 1
      }
    }
  }
}
dis as result "ORACLE B (ringer inserted): mean of `cnt2' cell slopes = " %14.8f `sum2'/`cnt2' "   ringer coef = " %14.8f _b[ringer]
dis "SENT_ORACLE_DONE"

foreach spec in "att cluster(unit)" "att customvce(cluster unit)" "att unconditionalse cluster(unit)" {
  use `G', clear
  cap qui wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(12)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "REAL | `spec' | ame = " %14.8f M[2,1] "  se = " %14.8f M[2,2] "  ci = [" %10.6f M[2,5] "," %10.6f M[2,6] "]  p = " %10.6f M[2,4]
  }
  else dis as error "REAL | `spec' | RC=" _rc
}
dis "SENT_LIN_DONE"

* ---------- realistic COUNT outcome, poisson, dose zero before treatment ----------
use `G', clear
set seed 55221
set sortseed 55221
qui gen double xbc = (0.9 + 0.20*i) + (0.04*t) + 0.25*tc + 0.25*dose + ue*0.3
qui gen long ycount = rpoisson(exp(xbc))
qui sum ycount, detail
dis as result "count outcome: mean = " %10.4f r(mean) "  p50 = " r(p50) "  max = " r(max) "  zeros = "
qui count if ycount == 0
dis as result "   n zero = " r(N) " of " _N
qui save `G', replace
dis as result "TRUE poisson dose semi-elasticity = 0.25 everywhere"
foreach spec in "att poisson cluster(unit)" "att poisson robust" "att poisson customvce(cluster unit)" {
  use `G', clear
  cap qui wooldid ycount i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(10)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "COUNT | `spec' | ame = " %14.8f M[2,1] "  se = " %14.8f M[2,2] "  ci = [" %10.6f M[2,5] "," %10.6f M[2,6] "]   (TRUE = 0.25)"
  }
  else dis as error "COUNT | `spec' | RC=" _rc
}
* control: the same count data with contreatcontrols(basic) is a DIFFERENT, identified spec
dis "SENT_COUNT_DONE"
dis "DONE_OOB_UZ3"
