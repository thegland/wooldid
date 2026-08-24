* reghdfe_oracle.do — EXTERNAL ORACLE for wooldid's event-study return values
*   e(b_es_att) / e(V_es_att) (the objects fed to honestdid and the Roth-style
*   pretrendpower diagnostics), validated against a reghdfe TWFE event study run
*   DIRECTLY on the same data — an external benchmark sharing no wooldid code path
*   beyond the data. (Both sides use reghdfe as the OLS engine, so the engine is
*   common; what this pins is wooldid's sample/cell construction, margins
*   aggregation, labeling/window logic, and b/V extraction plumbing.)
*
* Design: COMMON-ONSET panel, never-treated controls. 40 units, t = 1..12, onset 7
*   for units 1..20; y = unit effect + time effect + noise + true effect
*   (0.5 + 0.06*rel) from onset. Fixed seed; deterministic.
*
* THE FORMULATION THAT IS AN EXACT ORACLE (measured in the 2026-08-24 prototype):
*   pass a SINGLE TREATED COHORT id to wooldid (icoh = 1 for all treated units,
*   own unit id for controls) while clustering on the unit. Then wooldid's internal
*   regression  reghdfe y i._i#i._t#i.targetflag, absorb(i.icoh i.t) vce(cluster i)
*   has exactly ONE treated cell per relative period, so the margins pass returns
*   the cell coefficients themselves and the identical submatrix of the reghdfe
*   VCV. The external analogue is
*       reghdfe y D_rel..., absorb(icoh t) vce(cluster i)
*   with one dummy per treated relative period over the FULL sample support
*   (reference rel = -1 omitted; under esfixedbaseperiod wooldid does NOT trim
*   out-of-window horizons — they stay in the regression as their own cells and
*   are merely not reported, so the oracle must include them too).
*   Measured agreement: b mreldif 2.6e-16, V mreldif 0 (bitwise), att 7.0e-17.
*
* Labeling/window logic (verified empirically): with esfixedbaseperiod and
*   esrelativeto(-1), esprelength(p) reports the p pre horizons rel = -(p+1)..-2
*   labeled pre(p+1)..pre2 (the reference -1 = "pre1" is omitted), and
*   espostlength(q) reports post0..post(q-1); e(numpre_att) = p.
*
* Headline att: treated-observation-weighted average of the post-period cell
*   coefficients over the FULL post support (including unreported horizons);
*   weights = per-horizon treated obs counts (unweighted) or aweight sums ([aw]).
*   Measured: reldif <= 1.4e-16 across all variants.
*
* Variants covered: (1) baseline symmetric window 4/4; (2) analytic weights
*   [aw=w]; (3) asymmetric window esprelength(2) espostlength(4); (4) poisson via
*   ppmlhdfe (b/V on the xb scale; skipped, uncounted, if ppmlhdfe is absent).
*   Plus (5) the unit-level-i design (each treated unit its own cohort =
*   singleton cells): b and att still match the pooled reghdfe event study
*   (measured 4.1e-15 / 1.4e-16) but V is checked ONLY as an info line — with
*   singleton cells every treated observation is fit exactly by its own dummy
*   (ref obs by the absorbed unit FE), treated residuals are identically zero,
*   treated clusters' sandwich scores vanish, and wooldid's V carries
*   control-cluster noise only. That is a structurally different (Wooldridge
*   saturated-cell) V estimator, not a bug; measured gap vs the pooled-ES V:
*   mreldif 1.2e-2, per-horizon SE ratios 0.84..1.07.
*
* Tolerances: b mreldif <= 1e-8, V mreldif <= 1e-6, att reldif <= 1e-8 (all
*   measured at <= 4.2e-15, so no loosening was needed anywhere).
*
* Args: 1 = path to wooldid.ado (file or dir); 2 = unused.
* Run:  "<stata>" -e do reghdfe_oracle.do <adopath>
*    -> expect REGHDFE_ORACLE_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 20260824
set sortseed 20260824

local ADO `1'
if `"`ADO'"' == "" local ADO /Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'

global pf 0

**************************** synthetic common-onset panel ****************************
set obs 40
gen double i = _n
gen double a = rnormal()
gen double ttre = cond(i <= 20, 7, .)
expand 12
bysort i: gen double t = _n
gen double dt = 0.2*t + 0.5*cos(t)
gen double rel = t - ttre
gen double tau = cond(ttre < . & t >= ttre, 0.5 + 0.06*(t - ttre), 0)
gen double y = a + dt + tau + rnormal()*0.4
gen double icoh = cond(ttre < ., 1, i)
gen double w = 0.5 + runiform()*2
gen double mu = exp(0.3*a + 0.05*t + cond(ttre < . & t >= ttre, 0.3 + 0.02*(t - ttre), 0))
gen double ypois = rpoisson(mu)
* one dummy per treated relative period over the FULL support, reference -1 omitted
levelsof rel if ttre < . & rel != -1, local(relvals)
local dlist
foreach r of local relvals {
  if `r' < 0 local nm Dm`=abs(`r')'
  else local nm Dp`r'
  qui gen byte `nm' = (ttre < . & rel == `r')
  local dlist `dlist' `nm'
}
tempfile base
save `base'

* helper: map wooldid's ES labels to oracle dummy names and build OB/OV from the
* current e() (oracle fit), aligned 1:1 with WB's columns; posts local onames
cap program drop _oracle_align
program define _oracle_align
  local wlbl : colnames WB
  local onames
  foreach l of local wlbl {
    if substr("`l'", 1, 3) == "pre" local onames `onames' Dm`=substr("`l'", 4, .)'
    else local onames `onames' Dp`=substr("`l'", 5, .)'
  }
  local K : word count `wlbl'
  matrix OB = J(1, `K', .)
  matrix OV = J(`K', `K', .)
  forval j = 1/`K' {
    local nj : word `j' of `onames'
    matrix OB[1, `j'] = _b[`nj']
    forval k = 1/`K' {
      local nk : word `k' of `onames'
      matrix OV[`j', `k'] = e(V)["`nj'", "`nk'"]
    }
  }
end

**************************** O1-O5: baseline cohort formulation, window 4/4 ****************************
use `base', clear
cap noi wooldid y icoh t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "O1 FAIL: baseline cohort run rc=" _rc
  global pf = $pf + 1
  dis as error "REGHDFE_ORACLE_FAILS=$pf"
  exit, clear
}
matrix WB = e(b_es_att)
matrix WV = e(V_es_att)
scalar watt = _b[main_att]
local wlbl : colnames WB
if "`wlbl'" != "pre5 pre4 pre3 pre2 post0 post1 post2 post3" | e(numpre_att) != 4 {
  dis as error "O2 FAIL: labels [`wlbl'] numpre_att=" e(numpre_att) " (want pre5..pre2 post0..post3, numpre 4)"
  global pf = $pf + 1
}
use `base', clear
cap noi reghdfe y `dlist', absorb(icoh t) vce(cluster i)
if _rc != 0 {
  dis as error "O1b FAIL: oracle reghdfe rc=" _rc
  global pf = $pf + 1
  dis as error "REGHDFE_ORACLE_FAILS=$pf"
  exit, clear
}
_oracle_align
if mreldif(WB, OB) > 1e-8 {
  dis as error "O3 FAIL: baseline b mreldif = " %12.4e mreldif(WB, OB) " (tol 1e-8)"
  global pf = $pf + 1
}
else dis "O3 info: baseline b mreldif = " %12.4e mreldif(WB, OB)
if mreldif(WV, OV) > 1e-6 {
  dis as error "O4 FAIL: baseline V mreldif = " %12.4e mreldif(WV, OV) " (tol 1e-6)"
  global pf = $pf + 1
}
else dis "O4 info: baseline V mreldif = " %12.4e mreldif(WV, OV)
scalar onum = 0
scalar oden = 0
foreach r of local relvals {
  if `r' >= 0 {
    qui count if ttre < . & rel == `r'
    scalar onum = onum + r(N)*_b[Dp`r']
    scalar oden = oden + r(N)
  }
}
if reldif(watt, onum/oden) > 1e-8 {
  dis as error "O5 FAIL: att " %21.0g watt " vs count-weighted oracle post avg " %21.0g onum/oden " (reldif " %12.4e reldif(watt, onum/oden) ", tol 1e-8)"
  global pf = $pf + 1
}
else dis "O5 info: att aggregation reldif = " %12.4e reldif(watt, onum/oden)

**************************** O6-O8: analytic weights [aw=w] ****************************
use `base', clear
cap noi wooldid y icoh t ttre [aw=w], att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "O6 FAIL: [aw=w] run rc=" _rc
  global pf = $pf + 1
}
else {
  matrix WB = e(b_es_att)
  matrix WV = e(V_es_att)
  scalar watt = _b[main_att]
  use `base', clear
  cap noi reghdfe y `dlist' [aw=w], absorb(icoh t) vce(cluster i)
  if _rc != 0 {
    dis as error "O6b FAIL: [aw=w] oracle reghdfe rc=" _rc
    global pf = $pf + 1
  }
  else {
    _oracle_align
    if mreldif(WB, OB) > 1e-8 {
      dis as error "O6 FAIL: [aw] b mreldif = " %12.4e mreldif(WB, OB) " (tol 1e-8)"
      global pf = $pf + 1
    }
    else dis "O6 info: [aw] b mreldif = " %12.4e mreldif(WB, OB)
    if mreldif(WV, OV) > 1e-6 {
      dis as error "O7 FAIL: [aw] V mreldif = " %12.4e mreldif(WV, OV) " (tol 1e-6)"
      global pf = $pf + 1
    }
    else dis "O7 info: [aw] V mreldif = " %12.4e mreldif(WV, OV)
    scalar onum = 0
    scalar oden = 0
    foreach r of local relvals {
      if `r' >= 0 {
        qui sum w if ttre < . & rel == `r', meanonly
        scalar onum = onum + r(sum)*_b[Dp`r']
        scalar oden = oden + r(sum)
      }
    }
    if reldif(watt, onum/oden) > 1e-8 {
      dis as error "O8 FAIL: [aw] att " %21.0g watt " vs w-sum-weighted oracle post avg " %21.0g onum/oden " (reldif " %12.4e reldif(watt, onum/oden) ", tol 1e-8)"
      global pf = $pf + 1
    }
    else dis "O8 info: [aw] att aggregation reldif = " %12.4e reldif(watt, onum/oden)
  }
}

**************************** O9-O12: asymmetric window esprelength(2) espostlength(4) ****************************
use `base', clear
cap noi wooldid y icoh t ttre, att cluster(i) esprelength(2) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "O9 FAIL: asymmetric run rc=" _rc
  global pf = $pf + 1
}
else {
  matrix WB = e(b_es_att)
  matrix WV = e(V_es_att)
  scalar watt = _b[main_att]
  local wlbl : colnames WB
  if "`wlbl'" != "pre3 pre2 post0 post1 post2 post3" | e(numpre_att) != 2 {
    dis as error "O10 FAIL: asym labels [`wlbl'] numpre_att=" e(numpre_att) " (want pre3 pre2 post0..post3, numpre 2)"
    global pf = $pf + 1
  }
  use `base', clear
  cap noi reghdfe y `dlist', absorb(icoh t) vce(cluster i)
  if _rc != 0 {
    dis as error "O9b FAIL: asym oracle reghdfe rc=" _rc
    global pf = $pf + 1
  }
  else {
    _oracle_align
    if mreldif(WB, OB) > 1e-8 {
      dis as error "O11 FAIL: asym b mreldif = " %12.4e mreldif(WB, OB) " (tol 1e-8)"
      global pf = $pf + 1
    }
    else dis "O11 info: asym b mreldif = " %12.4e mreldif(WB, OB)
    if mreldif(WV, OV) > 1e-6 {
      dis as error "O12 FAIL: asym V mreldif = " %12.4e mreldif(WV, OV) " (tol 1e-6)"
      global pf = $pf + 1
    }
    else dis "O12 info: asym V mreldif = " %12.4e mreldif(WV, OV)
    scalar onum = 0
    scalar oden = 0
    foreach r of local relvals {
      if `r' >= 0 {
        qui count if ttre < . & rel == `r'
        scalar onum = onum + r(N)*_b[Dp`r']
        scalar oden = oden + r(N)
      }
    }
    if reldif(watt, onum/oden) > 1e-8 {
      dis as error "O12b FAIL: asym att reldif = " %12.4e reldif(watt, onum/oden) " (tol 1e-8)"
      global pf = $pf + 1
    }
    else dis "O12b info: asym att aggregation reldif = " %12.4e reldif(watt, onum/oden)
  }
}

**************************** O13-O15: poisson via ppmlhdfe (skipped, uncounted, if absent) ****************************
cap which ppmlhdfe
if _rc != 0 dis "O13 info: ppmlhdfe not installed; poisson variant skipped (uncounted)"
else {
  use `base', clear
  cap noi wooldid ypois icoh t ttre, att cluster(i) poisson esprelength(4) espostlength(4) esfixedbaseperiod
  if _rc != 0 {
    dis as error "O13 FAIL: poisson run rc=" _rc
    global pf = $pf + 1
  }
  else {
    matrix WB = e(b_es_att)
    matrix WV = e(V_es_att)
    scalar watt = _b[main_att]
    use `base', clear
    cap noi ppmlhdfe ypois `dlist', absorb(icoh t) vce(cluster i)
    if _rc != 0 {
      dis as error "O13b FAIL: oracle ppmlhdfe rc=" _rc
      global pf = $pf + 1
    }
    else {
      _oracle_align
      if mreldif(WB, OB) > 1e-8 {
        dis as error "O14 FAIL: poisson b mreldif = " %12.4e mreldif(WB, OB) " (tol 1e-8; both on the xb scale)"
        global pf = $pf + 1
      }
      else dis "O14 info: poisson b mreldif = " %12.4e mreldif(WB, OB)
      if mreldif(WV, OV) > 1e-6 {
        dis as error "O15 FAIL: poisson V mreldif = " %12.4e mreldif(WV, OV) " (tol 1e-6)"
        global pf = $pf + 1
      }
      else dis "O15 info: poisson V mreldif = " %12.4e mreldif(WV, OV)
      scalar onum = 0
      scalar oden = 0
      foreach r of local relvals {
        if `r' >= 0 {
          qui count if ttre < . & rel == `r'
          scalar onum = onum + r(N)*_b[Dp`r']
          scalar oden = oden + r(N)
        }
      }
      if reldif(watt, onum/oden) > 1e-8 {
        dis as error "O15b FAIL: poisson att (xb scale) reldif = " %12.4e reldif(watt, onum/oden) " (tol 1e-8)"
        global pf = $pf + 1
      }
      else dis "O15b info: poisson att aggregation reldif = " %12.4e reldif(watt, onum/oden)
    }
  }
}

**************************** O16-O17: unit-level i (singleton cells) — b and att ONLY ****************************
* V is NOT asserted here: with i = unit id every treated cell is a singleton, the
* saturated regression fits treated obs exactly, treated clusters' sandwich scores
* vanish, and wooldid's V carries control-cluster noise only — a different
* (Wooldridge saturated-cell) V estimator from the pooled event study's. Measured
* V mreldif ~1.2e-2 (SE ratios 0.84..1.07) is printed as an uncounted info line.
use `base', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "O16 FAIL: unit-level-i run rc=" _rc
  global pf = $pf + 1
}
else {
  matrix WB = e(b_es_att)
  matrix WV = e(V_es_att)
  scalar watt = _b[main_att]
  use `base', clear
  cap noi reghdfe y `dlist', absorb(i t) vce(cluster i)
  if _rc != 0 {
    dis as error "O16b FAIL: unit-level oracle reghdfe rc=" _rc
    global pf = $pf + 1
  }
  else {
    _oracle_align
    if mreldif(WB, OB) > 1e-8 {
      dis as error "O16 FAIL: unit-level b mreldif = " %12.4e mreldif(WB, OB) " (tol 1e-8)"
      global pf = $pf + 1
    }
    else dis "O16 info: unit-level b mreldif = " %12.4e mreldif(WB, OB)
    dis "O16 info (uncounted): unit-level V mreldif vs pooled-ES V = " %12.4e mreldif(WV, OV) " (structural difference, expected ~1e-2; see header)"
    scalar onum = 0
    scalar oden = 0
    foreach r of local relvals {
      if `r' >= 0 {
        qui count if ttre < . & rel == `r'
        scalar onum = onum + r(N)*_b[Dp`r']
        scalar oden = oden + r(N)
      }
    }
    if reldif(watt, onum/oden) > 1e-8 {
      dis as error "O17 FAIL: unit-level att reldif = " %12.4e reldif(watt, onum/oden) " (tol 1e-8)"
      global pf = $pf + 1
    }
    else dis "O17 info: unit-level att aggregation reldif = " %12.4e reldif(watt, onum/oden)
  }
}

if $pf == 0 dis "REGHDFE_ORACLE_ALL_PASS"
else dis as error "REGHDFE_ORACLE_FAILS=$pf"
exit, clear
