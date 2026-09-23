* test_poismargins.do - the poisson margins reference-scale gate.
* Run from test/: stata-mp -e do test_poismargins.do
*
* Since the 2026-09-20 removals the gate is UNCONDITIONAL: every margins call wooldid
* makes under poisson is a plain dydx or dyex of the linear index, so both injectors name
* predict(xb) whenever poisson is on and are empty off it. predict() reaches margins'
* ANALYTIC Jacobian, which is exact for a contrast that is linear in b; the
* expression(predict(xb)) branch that used to stand beside it forced a NUMERICAL Jacobian
* and existed only for the ey* operators, which logged the sign-unconstrained linear index.
* Those operators are gone (poisexpresults and semielasticity removed outright,
* contreatelasticitytype narrowed to dydx/dyex), so nothing reaches margins on a non-xb
* scale any more and the branch went with them.
*
* That collapse is numerically inert by construction, not merely by measurement: on the
* surviving path semielasticity defaulted to dydx and BS-6 remapped eydx -> dydx, so the
* expression() branch was already unreachable under poisson-without-poisexpresults. It was
* also checked empirically - pinned c1d8979 against the patched build, in separate Stata
* processes, over four poisson paths (att, att+contreat, att+contreat+dyex, poisson event
* study), b and diagonal V dumped at %21x: bit-identical on every cell.
*
* P1 the xb reporting scale is unchanged (the documented purpose of the gate)
* P2 plain dydx poisson event study runs and produces finite inference
* P3 the surviving contreat operators dydx and dyex both run under poisson on the xb scale
* P4 non-poisson runs are untouched by the gate entirely
* P5 one cluster stops under poisson (rc 459) as it does under reghdfe: the too-few-clusters
*    guard reads e(N_clust) as well as e(df_r), because ppmlhdfe never posts e(df_r)
*
* The retired option NAMES are pinned in test_si89.do (T1/T2), which is their cert host:
* test/suite/configs_ext.csv cannot carry them, because coverage_check.py treats an option
* token it cannot resolve against the syntax line as FATAL (author ruling 2026-09-11).
* Prints POISMARGINS_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

use "`REPO'/test/data/dta_binary.dta", clear
tempfile b
save `b'

* ---- P1: the xb scale. The att on the log scale must stay the log-scale value:
* exp() of it is the count-scale ratio, so a scale regression would move it hugely.
qui wooldid ycount i t ttre, att cluster(i) poisson
scalar p1b = _b[main_att]
di "P1_xb_scale_att=" p1b
if abs(p1b) > 5 {
  di as error "P1 FAIL: att is implausible for a log-scale report (scale regression?)"
  local ++fails
}

* ---- P2: plain dydx with an event study - the gated path - runs with finite inference
use `b', clear
cap qui wooldid ycount i t ttre, att cluster(i) poisson esprelength(2) espostlength(2)
local p2rc = _rc
matrix P2 = e(wdidmainresults)
di "P2_dydx_es_rc=" `p2rc' " se=" P2[1,2]
if `p2rc' != 0 | P2[1,2] >= . | P2[1,2] <= 0 {
  di as error "P2 FAIL: gated dydx path rc " `p2rc' " or non-finite SE"
  local ++fails
}

* ---- P3: both surviving contreat operators run under poisson. dyex logs the REGRESSOR,
* not the prediction, so it takes the predict(xb) form too (F032, author ruling 2026-09-06)
* and is not affected by the ey* removal. This leg is what stops a future narrowing of the
* gate from quietly taking dyex with it.
foreach v in dydx dyex {
  use "`REPO'/test/data/dta_cont.dta", clear
  cap qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) ///
      poisson contreatelasticitytype(`v') safety(off)
  local p3rc = _rc
  di "P3_contreat_`v'_rc=" `p3rc'
  if `p3rc' != 0 {
    di as error "P3 FAIL: contreatelasticitytype(`v') under poisson aborted (rc `p3rc')"
    local ++fails
  }
}

* ---- P4: non-poisson is untouched - the gate never fires
use `b', clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
scalar p4 = _b[main_att]
di "P4_nonpois_att=" p4
if p4 >= . {
  di as error "P4 FAIL: non-poisson run broken"
  local ++fails
}

* ---- P5: onecluster-guard-poisson. the worker's too-few-clusters guard keyed on e(df_r)
* alone, which ppmlhdfe never posts (missing, and `. < 1' is false), so the same one-cluster
* input stopped the reghdfe path with rc 459 and completed under poisson with every SE
* missing. the guard now also reads e(N_clust), which ppmlhdfe does post. P5a is the leg
* that flips (rc 0 on the unpatched build); P5b is the reghdfe parity anchor; P5c is an
* unclustered poisson run, where e(N_clust) is unposted and the added clause must stay false.
use `b', clear
gen byte one = 1
cap qui wooldid ycount i t ttre, att cluster(one) poisson esprelength(1) espostlength(1)
local p5arc = _rc
di "P5a_pois_onecluster_rc=" `p5arc'
if `p5arc' != 459 {
  di as error "P5a FAIL: poisson + cluster(one) rc " `p5arc' " (expected 459, as on the reghdfe path)"
  local ++fails
}
use `b', clear
gen byte one = 1
cap qui wooldid y i t ttre, att cluster(one) esprelength(1) espostlength(1)
local p5brc = _rc
di "P5b_reghdfe_onecluster_rc=" `p5brc'
if `p5brc' != 459 {
  di as error "P5b FAIL: reghdfe + cluster(one) rc " `p5brc' " (expected 459)"
  local ++fails
}
use `b', clear
cap qui wooldid ycount i t ttre, att robust poisson esprelength(1) espostlength(1)
local p5crc = _rc
local p5cse .
cap matrix P5 = e(wdidmainresults)
if _rc == 0 local p5cse = P5[1,2]
di "P5c_pois_robust_rc=" `p5crc' " se=" `p5cse'
if `p5crc' != 0 | `p5cse' >= . | `p5cse' <= 0 {
  di as error "P5c FAIL: unclustered poisson rc " `p5crc' " or non-finite SE (the cluster-count clause must not fire when e(N_clust) is unposted)"
  local ++fails
}

if `fails' == 0 di as result "POISMARGINS_ALL_PASS"
else di as error "POISMARGINS_FAILS=`fails'"
