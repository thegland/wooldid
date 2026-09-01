* test_fastcollapse.do - the cell collapse: default-on behaviour, its off-switch, and
* the eligibility gate. Run from test/: stata-mp -e do test_fastcollapse.do
*
* The collapse hands every margins call one row per (_i,_t,targetflag) cell carrying the
* cell's SUMMED weight, keeping the cell's FIRST row. Inside the gate the treatment
* derivative is constant within a cell, so the retained row's other values never enter
* the estimand.
*
* IT IS ON BY DEFAULT (2026-09-01). That inverts the usual test polarity: an equivalence
* check written against "the command without fastcollapse" would now compare the collapsed
* path to ITSELF and pass for free. Every comparison below therefore uses
* {fullprecisionslow} as the uncollapsed reference, and G1/G2 assert the default is really
* engaging.
*
* NOT BIT-IDENTICAL BY DESIGN. The controls cancel algebraically between xb(tf=1) and
* xb(tf=0) but not bitwise, leaving a residue of order ULP(|xb|). Thresholds below
* separate that (~1e-13) from the breakages the ineligible options produce (7e-4 to
* 4.5e-1) - four orders of headroom either side.
*
* NON-VACUITY IS ASSERTED, NOT ASSUMED. test/data/dta_binary.dta is DEGENERATE for the
* collapse (its main-site subpop cells are singletons, so the collapse is a measured no-op
* there and equality is free). F0 pins that the equivalence work runs on aggregating data.
*
* F0 the test data actually aggregates (subpop rows per cell >> 1)
* G1 the default ENGAGES: default != fullprecisionslow, and e(wdid_fcol) reports it
* G2 the default is the SAME PATH as an explicit fastcollapse (bit-identical)
* F1 default agrees with fullprecisionslow on b and V
* F2 ... and on the event study and the joint tests
* F3 compatible with fastnose
* F4 an explicit fastcollapse with an ineligible option is refused r(198), no safety(off)
*    escape; and fastcollapse + fullprecisionslow is refused as contradictory
* F5 WITHOUT an explicit request every ineligible option runs normally and silently skips
* F6 the off-switch reaches permute and bootstrap replications
* F7 poisson (dydx scale) collapses and agrees; poisexpresults does not collapse
* F8 REGRESSION: the subgroups()+fe()+jointtests command that hard-errored r(503) under the
*    collapse is refused when asked for explicitly and runs uncollapsed by default
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

* ---- the three arms: default (collapsed), fullprecisionslow (not), explicit request
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests
matrix bD = e(b)
matrix VD = e(V)
matrix eD = e(b_es_att)
scalar jD = e(joint_horiz)
local fcolD = e(wdid_fcol)
local rcD = _rc

use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests fullprecisionslow
matrix bS = e(b)
matrix VS = e(V)
matrix eS = e(b_es_att)
scalar jS = e(joint_horiz)
local fcolS = e(wdid_fcol)
local rcS = _rc

use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests fastcollapse
matrix bE = e(b)
local fcolE = e(wdid_fcol)

* ---- G1: the default really engages. If this ever reads 0 differences, the default has
* silently stopped applying and every equivalence check below is meaningless.
local ddefslow = mreldif(bD, bS)
di "G1_default_fcol=" `fcolD' " slow_fcol=" `fcolS' " explicit_fcol=" `fcolE' " default_vs_slow_mreldif=" `ddefslow'
if `fcolD' != 1 | `fcolS' != 0 | `fcolE' != 1 {
  di as error "G1 FAIL: e(wdid_fcol) should be 1/0/1 for default/fullprecisionslow/explicit, got `fcolD'/`fcolS'/`fcolE'"
  local ++fails
}
if `ddefslow' == 0 {
  di as error "G1 FAIL: the default produced results identical to fullprecisionslow - the collapse is NOT engaging, and every check below passes vacuously"
  local ++fails
}

* ---- G2: default and an explicit request are the same code path, so bit-identical
local ddefexp = mreldif(bD, bE)
di "G2_default_vs_explicit_mreldif=" `ddefexp'
if `ddefexp' != 0 {
  di as error "G2 FAIL: default and explicit fastcollapse differ (`ddefexp') - they must be the same code path"
  local ++fails
}

* ---- F1/F2: the collapsed default agrees with the uncollapsed reference
local db = mreldif(bD, bS)
local dv = mreldif(VD, VS)
local de = mreldif(eD, eS)
di "F1_rc=" `rcD' "/" `rcS' " mreldif_b=" `db' " mreldif_V=" `dv'
di "F2_mreldif_es=" `de' " joint_p=" jD " vs " jS
if `rcD' != 0 | `rcS' != 0 | `db' > 1e-8 | `dv' > 1e-8 {
  di as error "F1 FAIL: rc `rcD'/`rcS', mreldif b=`db' V=`dv' (want <1e-8)"
  local ++fails
}
if `de' > 1e-8 | reldif(jD, jS) > 1e-6 {
  di as error "F2 FAIL: event study mreldif=`de', joint p reldif=" reldif(jD, jS)
  local ++fails
}
* structural identity: same dimensions and same names
if colsof(bD) != colsof(bS) | rowsof(VD) != rowsof(VS) {
  di as error "F2 FAIL: structural change - b or V changed dimensions"
  local ++fails
}
local cnD : colnames bD
local cnS : colnames bS
if "`cnD'" != "`cnS'" {
  di as error "F2 FAIL: e(b) column names changed"
  local ++fails
}

* ---- F3: composes with fastnose
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) fastnose fullprecisionslow
matrix bC = e(b)
local rcC = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) fastnose
local rcF = _rc
matrix bF = e(b)
local fcolF = e(wdid_fcol)
local dcf = mreldif(bC, bF)
di "F3_fastnose_rc=" `rcC' "/" `rcF' " fcol=" `fcolF' " mreldif_b=" `dcf'
if `rcC' != 0 | `rcF' != 0 | `fcolF' != 1 | `dcf' > 1e-8 {
  di as error "F3 FAIL: fastnose rc `rcC'/`rcF', fcol `fcolF' (want 1), mreldif=`dcf'"
  local ++fails
}

* ---- F4: an EXPLICIT request for an ineligible configuration is refused r(198), and
* safety(off) does NOT open an escape (these are measured wrong answers, not risks).
* Plus the contradictory pair.
use "`REPO'/test/data/dta_cont.dta", clear
cap qui wooldid y i t ttre, att contreat(intensity) contreatcontrols(basic) cluster(i) fastcollapse safety(off)
local g1 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust interactivecontrols(x1 x2) fastcollapse safety(off)
local g2 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust unconditionalse fastcollapse safety(off)
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
cap qui wooldid y cohid year ttre, att cluster(state) subgroups(region) fastcollapse safety(off)
local g7 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust fastcollapse fullprecisionslow
local g8 = _rc
di "F4_gate_rcs contreat=" `g1' " interactive=" `g2' " uncond=" `g3' " semielas=" `g4' " impvalues=" `g5' " customvce=" `g6' " subgroups=" `g7' " contradictory=" `g8'
forval g = 1/8 {
  if `g`g'' != 198 {
    di as error "F4 FAIL: gate `g' returned rc `g`g'' (want 198)"
    local ++fails
  }
}

* ---- F5: THE DEFAULT-ON CONTRACT. Without an explicit fastcollapse, every one of those
* options must run NORMALLY and simply not collapse - an ordinary contreat() or
* subgroups() run must not start erroring because the collapse became the default.
* (unconditionalse needs robust/cluster and region varies within cohid, so the
* unconditionalse and subgroups(region) commands carry robust / safety(off): each
* runs rc 0 under fullprecisionslow, measured 2026-09-01, so a nonzero rc here is the
* default collapse and nothing else.)
use "`REPO'/test/data/dta_cont.dta", clear
cap qui wooldid y i t ttre, att contreat(intensity) contreatcontrols(basic) cluster(i)
local s1 = _rc
local f1 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust interactivecontrols(x1 x2)
local s2 = _rc
local f2 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust unconditionalse
local s3 = _rc
local f3 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust semielasticity
local s4 = _rc
local f4 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) impvalues
local s5 = _rc
local f5 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att customvce(cluster state)
local s6 = _rc
local f6 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) subgroups(region) safety(off)
local s7 = _rc
local f7 = e(wdid_fcol)
di "F5_silentskip rcs=" `s1' "/" `s2' "/" `s3' "/" `s4' "/" `s5' "/" `s6' "/" `s7' "  fcol=" `f1' "/" `f2' "/" `f3' "/" `f4' "/" `f5' "/" `f6' "/" `f7'
forval k = 1/7 {
  if `s`k'' != 0 {
    di as error "F5 FAIL: ineligible option `k' errored (rc `s`k'') without an explicit fastcollapse - it must run normally and skip the collapse"
    local ++fails
  }
  if `f`k'' != 0 {
    di as error "F5 FAIL: ineligible option `k' reports e(wdid_fcol)=`f`k'' - it must be 0"
    local ++fails
  }
}

* ---- F6: the off-switch reaches the replication engines (and the default runs them).
* dta_binary is used for permute because permute() requires clusters to be subsets of
* the treatment-onset groups, which dta_s_admin does not satisfy. PLUMBING checks: a rep
* that could not parse the forwarded option would abort the run.
use "`REPO'/test/data/dta_binary.dta", clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) fullprecisionslow
local p1 = _rc
use "`REPO'/test/data/dta_binary.dta", clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20)
local p2 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) bootstrap(standard) reps(15) rseed(7) fullprecisionslow
local b1 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) bootstrap(standard) reps(15) rseed(7)
local b2 = _rc
di "F6_reps permute_slow_rc=" `p1' " permute_default_rc=" `p2' " boot_slow_rc=" `b1' " boot_default_rc=" `b2'
if `p1' != 0 | `p2' != 0 | `b1' != 0 | `b2' != 0 {
  di as error "F6 FAIL: replication rcs `p1'/`p2'/`b1'/`b2'"
  local ++fails
}

* ---- F7: poisson on the dydx (log) scale collapses and agrees; the exponentiated
* reporting scale is not linear in the coefficients and does not collapse.
use "`ADMIN'", clear
qui wooldid ycount cohid year ttre, att poisson robust esprelength(2) espostlength(2) fullprecisionslow
matrix bP = e(b)
local rp1 = _rc
use "`ADMIN'", clear
cap qui wooldid ycount cohid year ttre, att poisson robust esprelength(2) espostlength(2)
local rp2 = _rc
matrix bQ = e(b)
local fcolQ = e(wdid_fcol)
local dpq = mreldif(bP, bQ)
use "`ADMIN'", clear
cap qui wooldid ycount cohid year ttre, att poisson poisexpresults robust
local rp3 = _rc
local fcolR = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid ycount cohid year ttre, att poisson poisexpresults robust fastcollapse safety(off)
local rp4 = _rc
di "F7_poisson_rc=" `rp1' "/" `rp2' " fcol=" `fcolQ' " mreldif=" `dpq' " poisexp_rc=" `rp3' " poisexp_fcol=" `fcolR' " poisexp_explicit_rc=" `rp4'
if `rp1' != 0 | `rp2' != 0 | `fcolQ' != 1 | `dpq' > 1e-8 {
  di as error "F7 FAIL: poisson dydx rc `rp1'/`rp2', fcol `fcolQ', mreldif=`dpq'"
  local ++fails
}
if `rp3' != 0 | `fcolR' != 0 | `rp4' != 198 {
  di as error "F7 FAIL: poisexpresults rc `rp3' fcol `fcolR' explicit rc `rp4' - it must run uncollapsed by default and be refused when asked for"
  local ++fails
}

* ---- F8: REGRESSION. Before subgroups() was excluded (2026-09-01) this exact command
* returned a hard r(503) WITH the collapse and rc 0 without it: collapsed, the att and
* att_it estimands lost different numbers of horizons to the omission-marker drop, and
* the non-conformable concatenation aborted the run. An explicit request must be a clean
* refusal, and the plain command must complete UNCOLLAPSED under the default.
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
local f8b = e(wdid_fcol)
di "F8_r503regression explicit=" `r8a' " (want 198, was 503)  default=" `r8b' " (want 0) default_fcol=" `f8b' " (want 0)"
if `r8a' != 198 | `r8b' != 0 | `f8b' != 0 {
  di as error "F8 FAIL: rc `r8a'/`r8b' fcol `f8b' - the r(503) regression is not pinned"
  local ++fails
}

if `fails' == 0 di as result "FASTCOLLAPSE_ALL_PASS"
else di as error "FASTCOLLAPSE_FAILS=`fails'"
