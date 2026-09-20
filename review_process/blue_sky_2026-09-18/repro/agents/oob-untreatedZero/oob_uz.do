* ============================================================================
* oob-untreatedZero : INDEPENDENT design (not the prober's) testing whether
* contreatcontrols(untreatedZero) silently acquires the pooled "basic" dose
* control on the unconditionalse / customvce() / poisson engine paths.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 77042
set sortseed 77042

* ---------------- MY design: 6 groups, 7 units each, 9 periods --------------
clear
set obs 378
gen long unit = ceil(_n/9)
gen int  t    = _n - (unit-1)*9
gen byte i    = ceil(unit/7)
gen double ttre = .
qui replace ttre = 3 if i == 1
qui replace ttre = 5 if i == 2
qui replace ttre = 7 if i == 3
* i = 4,5,6 never treated

gen byte tc = (ttre < . & t >= ttre)          /* treated cohort-period       */
gen double r = t - ttre if tc == 1            /* 0-based relative period     */

gen double dose = 0
qui replace dose = exp(rnormal(0.2,0.55)) if tc == 1     /* EXACTLY 0 elsewhere */

gen double beta = 0
qui replace beta = 1.45 - 0.31*i + 0.22*r - 0.06*r*r if tc == 1
gen double kappa = 0
qui replace kappa = 0.9 - 0.15*i + 0.13*r if tc == 1

gen double y = (1.1 + 0.47*i) + (-0.3 + 0.21*t - 0.012*t*t) + kappa + beta*dose

tempfile F
qui save `F'

qui sum beta if tc == 1
local TRUE = r(mean)
dis as result "== TRUE ame (obs mean of beta over treated cells) = " %18.12f `TRUE'
dis as result "   TRUE hex: " %21x `TRUE'

* predicted dropped cell = last in factor order = i=3, relyr 3 (t=9, r=2)
qui sum beta if i == 3 & t == 9
local BPRED = r(mean)
dis as result "== beta at (i=3,t=9)  = " %18.12f `BPRED' "   hex: " %21x `BPRED'
dis as result "== predicted broken ame = TRUE - beta(3,9) = " %18.12f `TRUE'-`BPRED' "  hex: " %21x (`TRUE'-`BPRED')

dis as text "== cell table (i, t, relyr[1-based], beta, n) =="
preserve
  keep if tc == 1
  gen double relyr = t - ttre + 1
  collapse (mean) beta (count) n=y, by(i t relyr)
  list i t relyr beta n, sepby(i) noobs
  mkmat i t relyr beta, matrix(CB)
restore
local nCB = rowsof(CB)
dis "SENT_SETUP_DONE"

* helper: given an implied dropped slope, name the closest cell
prog def scancell
  args BD
  local n = rowsof(CB)
  local best = .
  local bi = .
  forval q = 1/`n' {
    local d = abs(CB[`q',4] - (`BD'))
    if `d' < `best' {
      local best = `d'
      local bi = `q'
    }
  }
  if `bi' < . {
    dis as text "      closest cell: i=" CB[`bi',1] " t=" CB[`bi',2] " relyr=" CB[`bi',3] ///
      "  beta=" %16.10f CB[`bi',4] "  |BD-beta| = " %10.3e `best' ///
      "  rel = " %10.3e `best'/abs(CB[`bi',4])
  }
end

* ---------------- which cell does the sweep drop? my own replication --------
use `F', clear
qui gen double _tt = t - ttre + 1 if ttre < . & t >= ttre
qui replace _tt = 0 if _tt == .
qui gen double _ii = i if _tt > 0
qui replace _ii = 0 if _ii == .
qui gen byte tf = _tt > 0
qui gen double ringer = dose
qui gen byte nullv = 1
dis as text "==== MY OWN REPLICATION OF wooldid's ado:1485 SPELLING ===="
cap noi reg y i.i i.t i.nullv##c.(ringer) i._ii#i._tt#i.tf#c.dose i._ii#i._tt#i.tf
dis as text "---- colnames of e(b), marking omissions ----"
local cn : colnames e(b)
foreach c of local cn {
  if strpos("`c'","o.") != 0 dis as error "OMITTED COLUMN: `c'"
}
dis as text "---- cell dose coefficients from my replication ----"
forval q = 1/`nCB' {
  local a = CB[`q',1]
  local b = CB[`q',3]
  cap local cf = _b[`a'._ii#`b'._tt#1.tf#c.dose]
  if _rc == 0 dis "  cell i=`a' relyr=`b'  reg coef = " %16.10f `cf' "   true beta = " %16.10f CB[`q',4] "   diff = " %16.10f `cf'-CB[`q',4]
  else dis "  cell i=`a' relyr=`b'  NOT IN e(b)"
}
dis "  ringer coef = " %16.10f _b[ringer]
dis "SENT_REPL_DONE"

* ---------------- the wooldid runs -----------------------------------------
local speclist `" "att" "att robust" "att cluster(unit)" "att ccnoabsorb" "att ccnoabsorb robust" "att customvce(robust)" "att customvce(cluster unit)" "att unconditionalse" "att unconditionalse robust" "'
local k = 0
foreach spec of local speclist {
  local ++k
  use `F', clear
  cap qui wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(14)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    if `k' == 1 {
      dis as text "---- e(wdidmainresults) shape on spec 1 ----"
      mat li M
    }
    local A = M[2,1]
    local BD = `TRUE' - `A'
    dis as result "UZ | `spec' | ame = " %18.12f `A' "  relerr_vs_truth = " %10.3e abs(`A'-`TRUE')/abs(`TRUE') ///
      "  relerr_vs_pred = " %10.3e abs(`A'-(`TRUE'-`BPRED'))/abs(`TRUE'-`BPRED') "  r2 = " %18.14f e(r2)
    dis as result "      ame hex = " %21x `A' "   implied beta_dropped = " %16.10f `BD' "  hex " %21x `BD'
    scancell `BD'
  }
  else dis as error "UZ | `spec' | RC=" _rc
}
dis "SENT_MAIN_DONE"

* ---------------- comparison: contreatcontrols(basic) ----------------------
foreach spec in "att" "att customvce(robust)" {
  use `F', clear
  cap qui wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(basic) regtol(14)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "BASIC | `spec' | ame = " %18.12f M[2,1] "  hex " %21x M[2,1] "  (UZ truth = " %18.12f `TRUE' ")"
  }
  else dis as error "BASIC | `spec' | RC=" _rc
}
dis "SENT_BASIC_DONE"

* ---------------- poisson twins (exact-fit log-linear) ---------------------
use `F', clear
qui gen double xbh = (0.4 + 0.11*i) + (0.05*t) + 0.30*tc + 0.25*dose
qui gen double yh  = exp(xbh)
qui gen double gam = beta/4
qui gen double xbz = (0.4 + 0.11*i) + (0.05*t) + 0.30*tc + gam*dose
qui gen double yz  = exp(xbz)
qui sum gam if tc == 1
local TRUEP = r(mean)
qui sum gam if i == 3 & t == 9
local GPRED = r(mean)
dis as result "== TRUE poisson het semi-elasticity ame = " %18.12f `TRUEP' "  gam(3,9) = " %18.12f `GPRED' ///
  "  predicted broken = " %18.12f `TRUEP'-`GPRED'
qui save `F', replace

foreach spec in "att poisson" "att poisson robust" "att poisson ccnoabsorb" "att poisson customvce(robust)" "att poisson unconditionalse robust" {
  use `F', clear
  cap qui wooldid yh i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(12)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "POIS-HOM | `spec' | ame = " %18.12f M[2,1] "  hex " %21x M[2,1] "   (TRUE = 0.25 everywhere; prediction if broken = 0)"
  }
  else dis as error "POIS-HOM | `spec' | RC=" _rc
}
foreach spec in "att poisson" "att poisson robust" "att poisson customvce(robust)" "att poisson unconditionalse robust" {
  use `F', clear
  cap qui wooldid yz i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(12)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    local A = M[2,1]
    dis as result "POIS-HET | `spec' | ame = " %18.12f `A' "  truth = " %18.12f `TRUEP' ///
      "  relerr_vs_pred = " %10.3e abs(`A'-(`TRUEP'-`GPRED'))/abs(`TRUEP'-`GPRED') "  implied gam_dropped = " %16.10f `TRUEP'-`A'
  }
  else dis as error "POIS-HET | `spec' | RC=" _rc
}
dis "SENT_POIS_DONE"

* ---------------- WHAT THE USER SEES (verbatim, non-quiet) ------------------
foreach spec in "att" "att ccnoabsorb" "att customvce(robust)" "att unconditionalse robust" {
  dis as text "@@@@@@@@@@@@@@@@ VISIBLE OUTPUT BEGIN : `spec' @@@@@@@@@@@@@@@@"
  use `F', clear
  cap noi wooldid y i t ttre , `spec' contreat(dose) contreatcontrols(untreatedZero) regtol(14)
  dis as text "@@@@@@@@@@@@@@@@ VISIBLE OUTPUT END : `spec' (rc=" _rc ") @@@@@@@@@@@@@@@@"
}
dis as text "@@@@@@@@@@@@@@@@ VISIBLE OUTPUT BEGIN : att poisson @@@@@@@@@@@@@@@@"
use `F', clear
cap noi wooldid yz i t ttre , att poisson contreat(dose) contreatcontrols(untreatedZero) regtol(12)
dis as text "@@@@@@@@@@@@@@@@ VISIBLE OUTPUT END : att poisson (rc=" _rc ") @@@@@@@@@@@@@@@@"

dis "DONE_OOB_UZ"
