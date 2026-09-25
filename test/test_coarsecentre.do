* test_coarsecentre.do - the centred coarse cohort controls (author ruling 2026-09-25,
* econometric review findings F1a and F1b).
* Run from test/: stata-mp -e do test_coarsecentre.do -> expect COARSECENTRE_ALL_PASS
*
* coarsecohortcontrols() and ccc_absorb() used to build i._i#i._t#(vars) on the raw controls,
* so each treated cohort-period cell dummy read the effect at a control value of zero (the
* base category for an i. control), and the reported att moved when a constant was added to
* the control (-0.53 for c.x, -49.4 for c.x+100, against a true att of 2.00 on the fixture
* below). They now enter each control once as a pooled regressor and give each treated cell a
* slope on a copy centred at the cell's weighted mean (0 on reference rows), so the cell
* dummy reads the cell's average outcome less its predicted untreated outcome at the pooled
* slope.
*
* The fixture is the reviewers' design: 180 units, 8 periods, onsets 4 and 6 for 60 units
* each, 60 never treated, tau = 1 + .5(x - 5) + 2g, true att 1.9995; i is the onset cohort, so
* every treated cell holds many units and the within-cell slopes are estimated.
* C1 location and scale: coarsecohortcontrols(c.x), (c.x+100) and (c.10x-3) give the same
*    results (the whole e(wdidmainresults), 1e-8), with the event-study window on, so the
*    pre-treatment cells are centred too
* C2 ccc_absorb(): the same three spellings, the same results, and equal to
*    coarsecohortcontrols() (the absorbed and the explicit slopes are one model)
* C3 i. controls: coarsecohortcontrols(i.g) and (i.g2) with g2 = 1 - g (the other base level)
*    agree, and so do (c.x i.g) and (c.x+100 i.g2)
* C4 did_imputation: without the window the att is 1.92997, the reviewers' hand-built centred
*    model (reghdfe on the cell dummies, x and the centred cell slopes), and it is within 5e-4
*    of did_imputation's controls(x) (measured
*    1.92997 against 1.93019; they are close, not equal, because treated rows still inform the
*    unit effects, so the help claims no equivalence). did_imputation must be installed
* C5 the cell collapse is exact: the default and nofastcollapse agree to 1e-10 on
*    coarsecohortcontrols(c.x i.g) and on ccc_absorb(x) with the event-study window
* C6 poisson (ppmlhdfe): coarsecohortcontrols(c.x) and (c.x+100), and ccc_absorb(x), agree on
*    the att to 1e-6 (the ppmlhdfe convergence floor); C6b the plain-basic poisson path under
*    contreat(), where the ccc_absorb() slopes leave absorb(), is shift invariant too
* C7 bootstrap(standard) reps(20) completes, and the x and x+100 runs agree at the same rseed:
*    a replication that reused the observed run's cell means would not be shift invariant,
*    since each resample has its own cells and weights
* C8 permute(shuffle) permreps(20) permonsetcohorts completes, and the x and x+100 runs agree
*    at the same rseed: each permutation reshuffles the cells, so this holds only if every
*    replication centres on its own cells
clear all
set more off
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

cap program drop cc_fixture
program define cc_fixture
  clear
  set seed 20260925
  qui set obs 180
  gen unit = _n
  gen ttre = cond(unit <= 60, 4, cond(unit <= 120, 6, .))
  gen ai = rnormal()
  gen g = runiform() < .5
  qui expand 8
  bys unit: gen t = _n
  gen x = 5 + rnormal() + .3*ai
  gen double xs = x + 100
  gen double x10 = 10*x - 3
  gen g2 = 1 - g
  gen D = ttre < . & t >= ttre
  gen tau = 1 + 0.5*(x - 5) + 2*g
  gen y = ai + .2*t + .8*x + D*tau + rnormal()
  gen coh = cond(ttre == ., 99, ttre)
  gen yc = rpoisson(exp(1 + y/10))
  bys unit (t): gen double dose = .5 + runiform() if _n == 1
  bys unit (t): replace dose = dose[1]
end

* ---- C1 and C2: location and scale invariance, both options, event-study window on
local es esprelength(2) espostlength(2)
foreach opt in coarse ccc {
  local k 0
  foreach v in x xs x10 {
    local ++k
    local spec coarsecohortcontrols(c.`v')
    if "`opt'" == "ccc" local spec ccc_absorb(`v')
    cc_fixture
    cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) `es' `spec'
    local rc_`opt'`k' = _rc
    cap matrix drop M_`opt'`k'
    if `rc_`opt'`k'' == 0 matrix M_`opt'`k' = e(wdidmainresults)
  }
}
local c1d .
local c2d .
local c12d .
if `rc_coarse1' == 0 & `rc_coarse2' == 0 & `rc_coarse3' == 0 local c1d = max(mreldif(M_coarse1, M_coarse2), mreldif(M_coarse1, M_coarse3))
if `rc_ccc1' == 0 & `rc_ccc2' == 0 & `rc_ccc3' == 0 local c2d = max(mreldif(M_ccc1, M_ccc2), mreldif(M_ccc1, M_ccc3))
if `rc_coarse1' == 0 & `rc_ccc1' == 0 local c12d = mreldif(M_coarse1, M_ccc1)
di "C1 coarse rc=`rc_coarse1'/`rc_coarse2'/`rc_coarse3' mreldif x vs x+100, 10x-3 = `c1d'"
di "C2 ccc_absorb rc=`rc_ccc1'/`rc_ccc2'/`rc_ccc3' mreldif = `c2d' | ccc_absorb vs coarse = `c12d'"
if `c1d' > 1e-8 {
  di as error "C1 FAIL: coarsecohortcontrols() must be invariant to the control's location and units"
  local ++fails
}
if `c2d' > 1e-8 | `c12d' > 1e-8 {
  di as error "C2 FAIL: ccc_absorb() must be invariant to the control's location and units and equal coarsecohortcontrols()"
  local ++fails
}

* ---- C3: i. controls do not depend on the base level
local c3n 0
foreach pair in "i.g|i.g2" "c.x i.g|c.xs i.g2" {
  gettoken s1 s2 : pair, parse("|")
  local s2 = subinstr("`s2'", "|", "", 1)
  cc_fixture
  cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) coarsecohortcontrols(`s1')
  local r1 = _rc
  cap matrix drop C3a
  if `r1' == 0 matrix C3a = e(wdidmainresults)
  cc_fixture
  cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) coarsecohortcontrols(`s2')
  local r2 = _rc
  local d .
  if `r1' == 0 & `r2' == 0 local d = mreldif(C3a, e(wdidmainresults))
  di "C3 coarsecohortcontrols(`s1') vs (`s2'): rc `r1'/`r2' mreldif=`d'"
  if `d' > 1e-8 local ++c3n
}
if `c3n' > 0 {
  di as error "C3 FAIL: an i. coarse control must not depend on its base level"
  local ++fails
}

* ---- C4: did_imputation
cap which did_imputation
local c4w = _rc
local c4b .
if `c4w' == 0 {
  cc_fixture
  gen Ei = ttre
  cap qui did_imputation y unit t Ei, fe(unit t) controls(x)
  if _rc == 0 local c4b = _b[tau]
}
cc_fixture
cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) coarsecohortcontrols(c.x)
local c4rc = _rc
local c4a .
if `c4rc' == 0 local c4a = el(e(wdidmainresults), 1, 1)
di "C4 wooldid att=" %12.8f `c4a' " did_imputation controls(x)=" %12.8f `c4b' " (which rc `c4w')"
if `c4w' != 0 | `c4rc' != 0 | abs(`c4a' - `c4b') > 5e-4 | abs(`c4a' - 1.92997) > 5e-5 {
  di as error "C4 FAIL: the centred coarse control must give the hand-built 1.92997 and match did_imputation's controls(x) to 5e-4 on this fixture (did_imputation must be installed)"
  local ++fails
}

* ---- C5: cell collapse exactness
local c5n 0
foreach spec in "coarsecohortcontrols(c.x i.g)" "ccc_absorb(x) esprelength(2) espostlength(2)" {
  cc_fixture
  cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) `spec'
  local r1 = _rc
  local fc = e(wdid_fcol)
  cap matrix drop C5a C5b
  if `r1' == 0 {
    matrix C5a = e(wdidmainresults)
    matrix C5b = e(V)
  }
  cc_fixture
  cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) `spec' nofastcollapse
  local r2 = _rc
  local d .
  if `r1' == 0 & `r2' == 0 local d = max(mreldif(C5a, e(wdidmainresults)), mreldif(C5b, e(V)))
  di "C5 `spec': rc `r1'/`r2' collapsed=`fc' mreldif default vs nofastcollapse=`d'"
  if `d' > 1e-10 | "`fc'" != "1" local ++c5n
}
if `c5n' > 0 {
  di as error "C5 FAIL: the cell collapse must engage and be exact under the centred coarse controls"
  local ++fails
}

* ---- C6: poisson
local c6n 0
cc_fixture
cap qui wooldid yc coh t ttre, att fe(unit) robust poisson coarsecohortcontrols(c.x)
local r1 = _rc
local a1 .
if `r1' == 0 local a1 = el(e(wdidmainresults), 1, 1)
cc_fixture
cap qui wooldid yc coh t ttre, att fe(unit) robust poisson coarsecohortcontrols(c.xs)
local r2 = _rc
local a2 .
if `r2' == 0 local a2 = el(e(wdidmainresults), 1, 1)
cc_fixture
cap qui wooldid yc coh t ttre, att fe(unit) robust poisson ccc_absorb(x)
local r3 = _rc
local a3 .
if `r3' == 0 local a3 = el(e(wdidmainresults), 1, 1)
di "C6 poisson att coarse(c.x)=`a1' coarse(c.x+100)=`a2' ccc_absorb(x)=`a3' (rc `r1'/`r2'/`r3')"
if reldif(`a1', `a2') > 1e-6 | reldif(`a1', `a3') > 1e-6 {
  di as error "C6 FAIL: the poisson path must centre the coarse controls too"
  local ++fails
}
* C6b the poisson plain-basic path, which moves the ccc_absorb() slopes out of absorb() and
* into the regression (the spelling the r(3200) retry also uses): x and x+100 agree
local c6bn 0
foreach v in x xs {
  cc_fixture
  cap qui wooldid yc coh t ttre, att fe(unit) robust poisson contreat(dose) contreatcontrols(basic) ccc_absorb(`v')
  local c6brc_`v' = _rc
  cap matrix drop C6b_`v'
  if `c6brc_`v'' == 0 matrix C6b_`v' = e(wdidmainresults)
}
local c6bd .
if `c6brc_x' == 0 & `c6brc_xs' == 0 local c6bd = mreldif(C6b_x, C6b_xs)
di "C6b poisson plain basic ccc_absorb: rc `c6brc_x'/`c6brc_xs' mreldif x vs x+100=`c6bd'"
if `c6bd' > 1e-6 {
  di as error "C6b FAIL: the poisson plain-basic path must carry the centred ccc_absorb() slopes"
  local ++fails
}

* ---- C7: bootstrap replications centre on their own cells and weights
local c7n 0
foreach v in x xs {
  cc_fixture
  cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) coarsecohortcontrols(c.`v') bootstrap(standard) reps(20) rseed(4242) dots(0)
  local c7rc_`v' = _rc
  local c7ok_`v' = e(boot_completed)
  cap matrix drop C7_`v' C7V_`v'
  if `c7rc_`v'' == 0 {
    matrix C7_`v' = e(wdidmainresults)
    matrix C7V_`v' = e(V)
  }
}
local c7d .
if `c7rc_x' == 0 & `c7rc_xs' == 0 local c7d = max(mreldif(C7_x, C7_xs), mreldif(C7V_x, C7V_xs))
di "C7 bootstrap(standard) reps(20): rc `c7rc_x'/`c7rc_xs' completed `c7ok_x'/`c7ok_xs' mreldif x vs x+100=`c7d'"
if `c7rc_x' != 0 | `c7rc_xs' != 0 | `c7ok_x' != 20 | `c7ok_xs' != 20 | `c7d' > 1e-6 {
  di as error "C7 FAIL: bootstrap replications must complete and centre on their own resample"
  local ++fails
}

* ---- C8: permutation replications centre on their own cells
foreach v in x xs {
  cc_fixture
  cap qui wooldid y coh t ttre, att fe(unit) cluster(unit) coarsecohortcontrols(c.`v') permute(shuffle) permreps(20) permonsetcohorts rseed(5151) dots(0)
  local c8rc_`v' = _rc
  cap matrix drop C8_`v' C8d_`v'
  if `c8rc_`v'' == 0 {
    matrix C8_`v' = e(wdidmainresults)
    cap matrix C8d_`v' = e(permdist)
  }
}
local c8d .
if `c8rc_x' == 0 & `c8rc_xs' == 0 {
  local c8d = mreldif(C8_x, C8_xs)
  cap local c8d = max(`c8d', mreldif(C8d_x, C8d_xs))
}
di "C8 permute(shuffle) permreps(20): rc `c8rc_x'/`c8rc_xs' mreldif x vs x+100 (results and permutation distribution)=`c8d'"
if `c8rc_x' != 0 | `c8rc_xs' != 0 | `c8d' > 1e-6 {
  di as error "C8 FAIL: permutation replications must complete and centre on their own cells"
  local ++fails
}

if `fails' == 0 di as result "COARSECENTRE_ALL_PASS"
else di as error "COARSECENTRE_FAILS=`fails'"
