* test_fastcollapse.do - the fastcollapse option: equivalence, gate, and plumbing.
* Run from test/: stata-mp -e do test_fastcollapse.do
*
* fastcollapse hands every margins call one row per (_i,_t,targetflag[,subgroup]) cell
* carrying the cell's SUMMED weight, keeping the cell's FIRST row. Inside its gate the
* treatment derivative is constant within a cell, so the retained row's other values
* never enter the estimand.
*
* NOT BIT-IDENTICAL BY DESIGN. The controls cancel algebraically between xb(tf=1) and
* xb(tf=0) but not bitwise, leaving a residue of order ULP(|xb|). Measured to move no
* more than permuting the input rows moves the uncollapsed path itself. The thresholds
* below are therefore set to separate ULP-class (~1e-13) from the real breakages the
* refused options produce (7e-4 to 4.5e-1) - four orders of headroom either side.
*
* NON-VACUITY IS ASSERTED, NOT ASSUMED. test/data/dta_binary.dta is DEGENERATE at the
* main margins sites - its subpop cells are singletons, so the collapse is a no-op
* exactly where it matters and equality is free. F0 pins that this suite runs on data
* with genuine aggregation; without it every equivalence check below passes vacuously.
*
* F0 the test data actually aggregates (subpop rows per cell >> 1)
* F1 fastcollapse agrees with the uncollapsed path on b and V
* F2 ... and on the event study and the joint tests
* F3 fastcollapse is compatible with fastnose
* F4 every incompatible option is refused with r(198), no safety(off) escape
* F5 the option reaches permute replications
* F6 the option reaches bootstrap replications
* F7 poisson (dydx scale) is allowed and agrees; poisexpresults is refused
* Prints FASTCOLLAPSE_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

local ADMIN "`REPO'/test/suite/data/dta_s_admin.dta"

* ---- F0: NON-VACUITY. Count rows and cells in the post-treatment treated subpop the
* main margins site averages over. A ratio near 1 would make every check below free.
use "`ADMIN'", clear
qui gen byte __tf = (year >= ttre) & ttre < .
qui count if __tf == 1
local subrows = r(N)
qui gegen long __cell = group(cohid year) if __tf == 1
qui su __cell, meanonly
local subcells = r(max)
local ratio = `subrows'/`subcells'
di "F0_subpop_rows=" `subrows' " cells=" `subcells' " ratio=" %8.2f `ratio'
if `ratio' < 10 {
  di as error "F0 FAIL: subpop ratio `ratio' - the collapse barely aggregates here, so the equivalence checks would pass vacuously"
  local ++fails
}

* ---- F1/F2: equivalence on b, V, the event study and the joint tests
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests
matrix bA = e(b)
matrix VA = e(V)
matrix eA = e(b_es_att)
scalar jA = e(joint_horiz)
local rcA = _rc

use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests fastcollapse
matrix bB = e(b)
matrix VB = e(V)
matrix eB = e(b_es_att)
scalar jB = e(joint_horiz)
local rcB = _rc

local db = mreldif(bA, bB)
local dv = mreldif(VA, VB)
local de = mreldif(eA, eB)
di "F1_rc=" `rcA' "/" `rcB' " mreldif_b=" `db' " mreldif_V=" `dv'
di "F2_mreldif_es=" `de' " joint_p=" jA " vs " jB
if `rcA' != 0 | `rcB' != 0 | `db' > 1e-8 | `dv' > 1e-8 {
  di as error "F1 FAIL: rc `rcA'/`rcB', mreldif b=`db' V=`dv' (want <1e-8)"
  local ++fails
}
if `de' > 1e-8 | reldif(jA, jB) > 1e-6 {
  di as error "F2 FAIL: event study mreldif=`de', joint p reldif=" reldif(jA, jB)
  local ++fails
}
* structural identity: same dimensions and same names
if colsof(bA) != colsof(bB) | rowsof(VA) != rowsof(VB) {
  di as error "F2 FAIL: structural change - b or V changed dimensions"
  local ++fails
}
local cnA : colnames bA
local cnB : colnames bB
if "`cnA'" != "`cnB'" {
  di as error "F2 FAIL: e(b) column names changed"
  local ++fails
}

* ---- F3: compatible with fastnose
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) fastnose
matrix bC = e(b)
local rcC = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) fastnose fastcollapse
local rcD = _rc
matrix bD = e(b)
local dcd = mreldif(bC, bD)
di "F3_fastnose_rc=" `rcC' "/" `rcD' " mreldif_b=" `dcd'
if `rcC' != 0 | `rcD' != 0 | `dcd' > 1e-8 {
  di as error "F3 FAIL: fastnose+fastcollapse rc `rcC'/`rcD', mreldif=`dcd'"
  local ++fails
}

* ---- F4: the gate. Each of these must be refused with r(198), and safety(off) must
* NOT open an escape (these are measured wrong answers, not risks).
use "`REPO'/test/data/dta_cont.dta", clear
cap qui wooldid y i t ttre, att contreat(intensity) contreatcontrols(basic) cluster(i) fastcollapse safety(off)
local g1 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust interactivecontrols(x1 x2) fastcollapse safety(off)
local g2 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att unconditionalse fastcollapse safety(off)
local g3 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust semielasticity fastcollapse safety(off)
local g4 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) impvalues fastcollapse safety(off)
local g5 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att customvce(cluster state) fastcollapse safety(off)
local g6 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust subgroups(region) fastcollapse safety(off)
local g7 = _rc
di "F4_gate_rcs contreat=" `g1' " interactive=" `g2' " uncond=" `g3' " semielas=" `g4' " impvalues=" `g5' " customvce=" `g6' " subgroups=" `g7'
forval g = 1/7 {
  if `g`g'' != 198 {
    di as error "F4 FAIL: gate `g' returned rc `g`g'' (want 198) - safety(off) must not open an escape"
    local ++fails
  }
}

* ---- F5: the option reaches permute replications and the run stays clean.
* dta_binary is used here because permute() requires clusters to be subsets of the
* treatment-onset groups, which dta_s_admin does not satisfy for any of its cluster
* candidates. This is a PLUMBING check, not an equivalence claim: dta_binary is
* degenerate for the collapse (see the F0 note), so no numeric agreement is asserted
* from it. The forwarding itself is exercised end to end - a rep that could not parse
* the option would abort the run.
use "`REPO'/test/data/dta_binary.dta", clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) fastcollapse
local p1 = _rc
scalar pv1 = e(permute_p_att)
use "`REPO'/test/data/dta_binary.dta", clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20)
local p2 = _rc
scalar pv2 = e(permute_p_att)
di "F5_permute_rc=" `p1' "/" `p2' " p=" pv1 " vs " pv2
if `p1' != 0 | `p2' != 0 {
  di as error "F5 FAIL: permute+fastcollapse rc `p1'/`p2'"
  local ++fails
}

* ---- F6: the option reaches bootstrap replications
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) bootstrap(standard) reps(15) rseed(7) fastcollapse
local b1 = _rc
di "F6_bootstrap_rc=" `b1'
if `b1' != 0 {
  di as error "F6 FAIL: bootstrap+fastcollapse rc `b1'"
  local ++fails
}

* ---- F7: poisson on the dydx (log) scale is allowed and agrees; the exponentiated
* reporting scale is not linear in the coefficients and is refused.
use "`ADMIN'", clear
qui wooldid ycount cohid year ttre, att poisson robust esprelength(2) espostlength(2)
matrix bP = e(b)
local rp1 = _rc
use "`ADMIN'", clear
cap qui wooldid ycount cohid year ttre, att poisson robust esprelength(2) espostlength(2) fastcollapse
local rp2 = _rc
matrix bQ = e(b)
local dpq = mreldif(bP, bQ)
use "`ADMIN'", clear
cap qui wooldid ycount cohid year ttre, att poisson poisexpresults robust fastcollapse safety(off)
local rp3 = _rc
di "F7_poisson_rc=" `rp1' "/" `rp2' " mreldif=" `dpq' " poisexpresults_rc=" `rp3'
if `rp1' != 0 | `rp2' != 0 | `dpq' > 1e-8 {
  di as error "F7 FAIL: poisson dydx rc `rp1'/`rp2', mreldif=`dpq'"
  local ++fails
}
if `rp3' != 198 {
  di as error "F7 FAIL: poisexpresults returned rc `rp3' (want 198)"
  local ++fails
}

* ---- F8: REGRESSION. Before subgroups() was refused (2026-09-01) this exact command
* returned a hard r(503) WITH fastcollapse and rc 0 without it: collapsed, the att and
* att_it estimands lost different numbers of horizons to the omission-marker drop, and
* the non-conformable concatenation aborted the run. It must now be a clean refusal, and
* the command without the option must still complete.
use "`REPO'/test/suite/data/survey_rcs.dta", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], att att_it subgroups(group) ///
    fe(group_time) esprelength(2) espostlength(2) jointtests cluster(state) safety(off) fastcollapse
local r8a = _rc
use "`REPO'/test/suite/data/survey_rcs.dta", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], att att_it subgroups(group) ///
    fe(group_time) esprelength(2) espostlength(2) jointtests cluster(state) safety(off)
local r8b = _rc
di "F8_r503regression with_option=" `r8a' " (want 198, was 503)  without=" `r8b' " (want 0)"
if `r8a' != 198 | `r8b' != 0 {
  di as error "F8 FAIL: rc `r8a'/`r8b' - the r(503) regression is not pinned"
  local ++fails
}

if `fails' == 0 di as result "FASTCOLLAPSE_ALL_PASS"
else di as error "FASTCOLLAPSE_FAILS=`fails'"
