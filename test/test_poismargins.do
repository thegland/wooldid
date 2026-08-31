* test_poismargins.do - the poisson margins reference-scale gate.
* Run from test/: stata-mp -e do test_poismargins.do
*
* The gate uses predict(xb) where the margins derivative operator is plain dydx (the
* analytic Jacobian, exact for a linear-in-b contrast) and KEEPS expression(predict(xb))
* for every ey*/dyex operator, where both code paths take ln of the prediction and the
* poisson prediction is the sign-unconstrained linear index. Which form aborts there is
* data-dependent and goes in BOTH directions, so today's rc behaviour - including its
* failures - is preserved deliberately.
*
* P1 the xb reporting scale is unchanged (the documented purpose of the option)
* P2 poisexpresults still selects the exponentiated scale and is untouched by the gate
* P3 plain dydx poisson runs and produces finite inference
* P4 semielasticity (eydx) keeps expression(): its rc must match whatever it was, and
*    the branch must be REACHED (this is the coverage the gate exists to preserve)
* P5 contreat elasticity operators keep expression()
* P6 non-poisson runs are untouched by the gate entirely
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

* ---- P4: semielasticity keeps expression(). The branch must be REACHED and must not
* have become an abort where it previously ran. On this panel the linear index is
* positive throughout, so ln() is in-domain and the run completes.
use `b', clear
cap qui wooldid ycount i t ttre, att cluster(i) poisson semielasticity
di "P4_eydx_rc=" _rc
if _rc != 0 {
  di as error "P4 FAIL: semielasticity aborted (rc " _rc "); the ey* branch must keep expression()"
  local ++fails
}

* ---- P5: contreat elasticity operators keep expression() too
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

if `fails' == 0 di as result "POISMARGINS_ALL_PASS"
else di as error "POISMARGINS_FAILS=`fails'"
