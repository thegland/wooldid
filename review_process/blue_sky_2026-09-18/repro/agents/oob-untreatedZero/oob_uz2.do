clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 77042
set sortseed 77042

* ---- rebuild the SAME exact-fit design as run 1 ----
clear
set obs 378
gen long unit = ceil(_n/9)
gen int  t    = _n - (unit-1)*9
gen byte i    = ceil(unit/7)
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
gen double y = (1.1 + 0.47*i) + (-0.3 + 0.21*t - 0.012*t*t) + kappa + beta*dose
tempfile F
qui save `F'

* ---- A. TRUE att and TRUE per-horizon ame ----
gen double effect = kappa + beta*dose
qui sum effect if tc == 1
dis as result "TRUE att (mean treatment effect on treated rows) = " %18.12f r(mean)
qui sum beta if tc == 1
local TRUE = r(mean)
dis as result "TRUE ame = " %18.12f `TRUE'
dis as text "TRUE per-horizon ame (relyr = t-ttre+1):"
forval h = 1/7 {
  qui sum beta if tc == 1 & (t - ttre + 1) == `h'
  if r(N) > 0 dis "   relyr `h' : n=" r(N) "  true ame = " %14.8f r(mean) "   minus 0.72 = " %14.8f r(mean)-0.72
}
dis "SENT_A_DONE"

* ---- B. what does reghdfe do with the absorbed i.null##c.ringer? ----
use `F', clear
qui gen double _tt = t - ttre + 1 if ttre < . & t >= ttre
qui replace _tt = 0 if _tt == .
qui gen double _ii = i if _tt > 0
qui replace _ii = 0 if _ii == .
qui gen byte tf = _tt > 0
qui gen double ringer = dose
qui gen byte nullv = 1
dis as text "==== B1: reghdfe WITH the absorbed ringer set (what wooldid's default path does) ===="
cap noi reghdfe y i._ii#i._tt#i.tf#c.dose i._ii#i._tt#i.tf , absorb(i.i i.t i.nullv##c.ringer) tol(1e-14)
dis "   e(df_a) = " e(df_a) "   e(r2) = " %18.14f e(r2)
qui levelsof _ii if _ii > 0, local(il)
foreach a of local il {
  forval b = 1/7 {
    cap local cf = _b[`a'._ii#`b'._tt#1.tf#c.dose]
    if _rc == 0 dis "   B1 cell i=`a' relyr=`b'  coef = " %14.8f `cf'
  }
}
dis as text "==== B2: reghdfe WITHOUT the absorbed ringer set (control) ===="
cap noi reghdfe y i._ii#i._tt#i.tf#c.dose i._ii#i._tt#i.tf , absorb(i.i i.t) tol(1e-14)
dis "   e(df_a) = " e(df_a) "   e(r2) = " %18.14f e(r2)
foreach a of local il {
  forval b = 1/7 {
    cap local cf = _b[`a'._ii#`b'._tt#1.tf#c.dose]
    if _rc == 0 dis "   B2 cell i=`a' relyr=`b'  coef = " %14.8f `cf'
  }
}
dis "SENT_B_DONE"

* ---- C. what does verbose tell the user was dropped? ----
dis as text "@@@@ VERBOSE customvce(robust) BEGIN @@@@"
use `F', clear
cap noi wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) customvce(robust) verbose regtol(14)
dis as text "@@@@ VERBOSE customvce(robust) END @@@@"
dis "SENT_C_DONE"

* ---- G. unconditionalse alone: what does the user see? ----
dis as text "@@@@ unconditionalse ALONE BEGIN @@@@"
use `F', clear
cap noi wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) unconditionalse regtol(14)
dis as text "@@@@ unconditionalse ALONE END (rc=" _rc ") @@@@"
dis "SENT_G_DONE"

* ---- F. event-study horizons ----
foreach spec in "att espost(3)" "att espost(3) customvce(robust)" {
  use `F', clear
  cap qui wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(14)
  if _rc == 0 {
    dis as result "==== ES | `spec' ===="
    mat M = e(wdidmainresults)
    mat li M, format(%14.8f)
  }
  else dis as error "ES | `spec' | RC=" _rc
}
dis "SENT_F_DONE"

* ---- E. lattice variant (shape of the help's own example 8) ----
foreach spec in "att" "att customvce(robust)" "att unconditionalse robust" {
  use `F', clear
  cap qui wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) lattice(baseTreatedPost 3) regtol(14)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "LATTICE3 | `spec' | att = " %16.10f M[1,1] "  ame = " %16.10f M[2,1] "  (no-lattice truth ame = " %16.10f `TRUE' ")"
  }
  else dis as error "LATTICE3 | `spec' | RC=" _rc
}
dis "SENT_E_DONE"

* ---- D. REALISTIC noisy design, clustered SEs ----
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
qui sum beta if i == 3 & t == 9
local BD2 = r(mean)
dis as result "REALISTIC: TRUE ame = " %16.10f `TRUE2' "   beta(3,9) = " %16.10f `BD2' "   predicted broken = " %16.10f `TRUE2'-`BD2'
tempfile G
qui save `G'
foreach spec in "att cluster(unit)" "att customvce(cluster unit)" "att unconditionalse cluster(unit)" "att poisson cluster(unit)" {
  use `G', clear
  cap qui wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(12)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "REAL | `spec' | ame = " %14.8f M[2,1] "  se = " %14.8f M[2,2] "  ci = [" %10.6f M[2,5] "," %10.6f M[2,6] "]  truth = " %14.8f `TRUE2'
  }
  else dis as error "REAL | `spec' | RC=" _rc
}
dis "SENT_D_DONE"
dis "DONE_OOB_UZ2"
