* test_poismargins.do - the poisson margins reference-scale gate.
* Run from test/: stata-mp -e do test_poismargins.do
*
* The gate uses predict(xb) where the margins derivative operator is dydx or (since
* 2026-09-06, F032) dyex - the analytic Jacobian, exact for a linear-in-b contrast; dyex logs
* the regressor, not the prediction - and KEEPS expression(predict(xb)) for the ey*
* operators, where both code paths take ln of the prediction and the poisson prediction is
* the sign-unconstrained linear index. Which form aborts there is data-dependent and goes
* in BOTH directions, so today's rc behaviour - including its failures - is preserved
* deliberately.
*
* P1 the xb reporting scale is unchanged (the documented purpose of the option)
* P2 poisexpresults still selects the exponentiated scale and is untouched by the gate
* P3 plain dydx poisson runs and produces finite inference
* P4 semielasticity under poisson runs to completion. Since BS-6 it is computed as dydx
*    on the xb scale, so the ey* branch is no longer reached on this combination; the
*    expression() lines stay as the guard for a future operator
* P5 contreat elasticity operators run; as P4, eydx is computed as dydx since BS-6
* P6 non-poisson runs are untouched by the gate entirely
* P7 (BS-6) poisson without poisexpresults computes semielasticity as dydx on the xb
*    scale, so the reported number is the plain poisson one; the header still says
*    semi-elasticities, the caution says why, and poisexpresults is untouched
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

* ---- P2: poisexpresults selects the other scale and must differ from P1
use `b', clear
qui wooldid ycount i t ttre, att cluster(i) poisson poisexpresults
scalar p2b = _b[main_att]
di "P2_expscale_att=" p2b
if reldif(p1b, p2b) < 1e-6 {
  di as error "P2 FAIL: poisexpresults did not change the reporting scale"
  local ++fails
}

* ---- P3: plain dydx with an event study - the gated path - runs with finite inference
use `b', clear
cap qui wooldid ycount i t ttre, att cluster(i) poisson esprelength(2) espostlength(2)
local p3rc = _rc
matrix P3 = e(wdidmainresults)
di "P3_dydx_es_rc=" `p3rc' " se=" P3[1,2]
if `p3rc' != 0 | P3[1,2] >= . | P3[1,2] <= 0 {
  di as error "P3 FAIL: gated dydx path rc " `p3rc' " or non-finite SE"
  local ++fails
}

* ---- P4: semielasticity under poisson still runs to completion. Since BS-6 (author
* ruling 2026-09-07) the ey* branch is no longer REACHED on this combination - the
* operator is remapped to dydx on the xb scale, where a plain marginal effect already IS
* the semi-elasticity - so what this leg now pins is that the remap did not break the run.
* The expression() lines are kept as the guard for a future operator; P7 pins the remap.
use `b', clear
cap qui wooldid ycount i t ttre, att cluster(i) poisson semielasticity
di "P4_eydx_rc=" _rc
if _rc != 0 {
  di as error "P4 FAIL: semielasticity aborted (rc " _rc "); the ey* branch must keep expression()"
  local ++fails
}

* ---- P5: contreat elasticity operators still run. As P4: since BS-6 the eydx request is
* computed as dydx on the xb scale, and the expression() line below it stands as the guard
use "`REPO'/test/data/dta_cont.dta", clear
cap qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) poisson contreatelasticitytype(eydx) safety(off)
di "P5_contreat_eydx_rc=" _rc " (contract: whatever it was, preserved)"

* ---- P6: non-poisson is untouched - the gate never fires
use `b', clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
scalar p6 = _b[main_att]
di "P6_nonpois_att=" p6
if p6 >= . {
  di as error "P6 FAIL: non-poisson run broken"
  local ++fails
}

* ---- P7 (BS-6): under poisson WITHOUT poisexpresults, semielasticity reports the plain
* poisson marginal effect - the semi-elasticity of E[y] on the xb scale - rather than
* margins' eydx, which takes a further log of the linear index and reports a number that
* depends on the units of y. The operator is remapped; the raw request is kept for the
* header and posted in e(wdid_semielasraw). poisexpresults is untouched.
use `b', clear
qui log using pois_p7.txt, text replace name(pm)
cap noi wooldid ycount i t ttre, att cluster(i) poisson semielasticity
local p7rc = _rc
qui log close pm
matrix P7A = e(wdidmainresults)
local p7op = "`e(wdid_semielas)'"
local p7raw = "`e(wdid_semielasraw)'"
local p7c 0
qui file open p7h using pois_p7.txt, read text
file read p7h p7l
while r(eof) == 0 {
  if strpos(`"`macval(p7l)'"', "Caution: under poisson without poisexpresults") != 0 local p7c 1
  file read p7h p7l
}
file close p7h
cap erase pois_p7.txt
use `b', clear
qui wooldid ycount i t ttre, att cluster(i) poisson
matrix P7B = e(wdidmainresults)
use `b', clear
qui wooldid ycount i t ttre, att cluster(i) poisson poisexpresults semielasticity
local p7op2 = "`e(wdid_semielas)'"
local p7raw2 = "`e(wdid_semielasraw)'"
di "P7 BS-6: rc=`p7rc' operator=`p7op' raw=`p7raw' caution=`p7c' | remapped b=" %10.7f P7A[1,1] " se=" %10.7f P7A[1,2] " vs plain poisson b=" %10.7f P7B[1,1] " se=" %10.7f P7B[1,2] " | poisexp operator=`p7op2' raw=[`p7raw2']"
if `p7rc' != 0 | "`p7op'" != "dydx" | "`p7raw'" != "eydx" | `p7c' != 1 | reldif(P7A[1,1], P7B[1,1]) > 1e-12 | reldif(P7A[1,2], P7B[1,2]) > 1e-12 | "`p7op2'" != "eydx" | "`p7raw2'" != "" {
  di as error "P7 FAIL: the poisson xb-scale remap did not fire, or moved poisexpresults with it"
  local ++fails
}

if `fails' == 0 di as result "POISMARGINS_ALL_PASS"
else di as error "POISMARGINS_FAILS=`fails'"
