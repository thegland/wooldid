* ===========================================================================
* FINDING E (self-contained repro; run: stata-mp -e do <thisfile>)
*
* CLAIM: contreatcontrols(untreatedZero) is defined as the one control type
* that adds NO baseline dose control (ado:1311-1315: "if always zero when not
* treated, the only thing to do is not control for anything!"; sthlp:995 says
* it "cannot be paired with any other type"). But the block that builds the
* `basic' pooled dose control at ado:1234 is gated on
*     strpos(contreatcontrols,"basic") | unconditionalse | customvce() | poisson
* so on ANY of those three engine paths wooldid silently inserts the pooled
* i.null##c.(ringer) baseline dose control into an untreatedZero run. Because
* the dose is exactly 0 outside treated cohort-periods, that control is not
* redundant - it absorbs part of the treatment dose response - and the reported
* average marginal effect is wrong.
*
* ORACLE: exact-fit design. dose == 0 on every untreated observation, and
*     y = a_i + g_t + kappa_it + beta_it*dose  on treated post rows.
* The truth is mean(beta_it) over the treated post rows.
* ===========================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31415
set sortseed 31415

clear
set obs 320                       /* 4 cohort groups x 10 units x 8 periods */
gen long unit = ceil(_n/8)
gen int  t    = _n - (unit-1)*8
gen byte i    = ceil(unit/10)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
* i = 3,4 never treated
gen double dose = 0
qui replace dose = 1 + 2*runiform() if ttre<. & t>=ttre     /* ZERO when untreated */
gen double b1 = 0
gen double k0 = 0
qui replace k0 = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre<. & t>=ttre
qui replace b1 = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre<. & t>=ttre
gen double y = 2 + 0.6*i + 0.35*t + k0 + b1*dose
tempfile F
qui save `F'

qui sum b1 if ttre<. & t>=ttre
local TRUE = r(mean)
dis as result "TRUE ame = " %16.10f `TRUE'

foreach spec in "robust" "ccnoabsorb robust" "customvce(robust)" "customvce(cluster unit)" "unconditionalse robust" {
  use `F', clear
  cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(untreatedZero) `spec' regtol(14)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "UZ | `spec' | ame = " %16.10f M[2,1] "  truth = " %16.10f `TRUE' ///
        "  rel.err = " %10.3e abs(M[2,1]-`TRUE')/abs(`TRUE') "  R2 = " %14.12f e(r2)
  }
  else dis as error "UZ | `spec' | RC=" _rc
}

* poisson twin: same gate (ado:1234 includes poisson), exact-fit log-linear
use `F', clear
qui gen double xbi = 0.5 + 0.10*i + 0.04*t
qui replace xbi = xbi + 0.10 + 0.20*dose if ttre<. & t>=ttre
qui gen double ypo = exp(xbi + rnormal(0,0.15))
* on the xb reporting scale the true dose semi-elasticity is exactly 0.20
foreach spec in "robust" "ccnoabsorb robust" {
  cap qui wooldid ypo i t ttre , att contreat(dose) contreatcontrols(untreatedZero) poisson `spec' regtol(9)
  if _rc == 0 {
    mat M = e(wdidmainresults)
    dis as result "UZ-POISSON | `spec' | ame = " %16.10f M[2,1] "   (TRUE dose semi-elasticity = 0.20)"
  }
  else dis as error "UZ-POISSON | `spec' | RC=" _rc
}
dis "DONE_C"
