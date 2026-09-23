* test_fastcollapse.do - the cell collapse: default-on behaviour, its off-switch, and
* the eligibility gate. Run from test/: stata-mp -e do test_fastcollapse.do
*
* The collapse hands every margins call one row per (_i,_t,targetflag) cell carrying the
* cell's SUMMED weight, keeping the cell's FIRST row. Inside the gate the treatment
* derivative is constant within a cell, so the retained row's other values never enter
* the estimand.
*
* IT IS ON BY DEFAULT (2026-09-01). That inverts the usual test polarity: an equivalence
* check written against "the command without an off-switch" would now compare the collapsed
* path to ITSELF and pass for free. Every comparison below therefore uses
* {nofastcollapse} as the uncollapsed reference, and G1/G2 assert the default is really
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
* G1 the default ENGAGES: default != nofastcollapse, and e(wdid_fcol) reports it
* G2 THE CLEAN BREAK: the two retired spellings, and an abbreviation of the new one, are
*    all refused r(198)
* F1 default agrees with nofastcollapse on b and V
* F2 ... and on the event study and the joint tests
* F3 compatible with fastnose
* F4 an ineligible option runs normally, skips the collapse, and DISCLOSES the skip in a
*    caution naming the option; safety(off) does not silence it and nofastcollapse does;
*    a subgroups() run is eligible and collapses (subgroups() left the ineligible list
*    2026-09-01; contreat() left it 2026-09-02, so F4's contreat leg now carries
*    contreatpoly(2))
* F5 every ineligible option runs normally and skips the collapse, while a subgroups() run
*    and a bare contreat() run collapse
* F6 the off-switch reaches permute and bootstrap replications
* F7 poisson (dydx scale) collapses and agrees. Its poisexpresults half went with that
*    option on 2026-09-20; no poisson reporting scale is uncollapsible any more
* F8 REGRESSION: the subgroups()+fe()+jointtests command that hard-errored r(503) under the
*    collapse now completes collapsed (safety(off) => the screen prints its caution)
* F9 THE SUBGROUP SCREEN: a healthy subgroup run collapses with an empty e(wdid_fcguard); a
*    numerically singular subgroup VCE under the default safety falls back to the full data
*    (e(wdid_fcguard)=fallback) and every result read after the fallback is BIT-IDENTICAL to
*    nofastcollapse; the same run under safety(off) stays collapsed with
*    e(wdid_fcguard)=caution and rc 0
* F10 THE 2026-09-01 ADVERSARIAL PROBE'S PINS: (a) a subgroup whose pre window holds only
*    the reference period (no pre-trend joint test in either arm) does not trip the screen;
*    (b) a rank-1 subgroup-contrast VCE with a tiny largest eigenvalue is flagged at the
*    level scale; (c) an observed-run fallback carries into the bootstrap replications
* F11 THE CONTINUOUS TREATMENT: contreat() left the ineligible list 2026-09-02, so the
*    collapse now carries the cell's estimand-weighted mean of the dose on the surviving row
*    and puts the user's cohort variable in the key. (a) five eligible reads - plain att; all
*    four estimands with an event study and jointtests; contreatelasticitytype(dyex); poisson
*    on the default dydx scale; subgroups(); poisson+dyex (since 2026-09-06, F032) - collapse
*    and agree with nofastcollapse, and so do lattice() and makecfxplots, whose dose bin joins
*    the cell key, across all four estimands under aweights on every e() matrix and scalar;
*    (b) each of the two refused reads (contreatpoly(>1), the dose inside controls()) skips
*    and stays bit-identical to it;
*    (c) each of those skips announces itself in a caution naming its own token, which
*    safety(off) does not silence, while a plain contreat() run collapses quietly; (d) the two provisos the worker's form exists for -
*    a customweightmultiplier() of 0 on the pick row (the mean is written onto every row of
*    the cell, not by gcollapse (mean) [aw=w] merge) and a float-stored dose (the recast)
* F12 THE TWO GUARDS ADDED ON TOP OF F11 (2026-09-02). (a) THE GATE UNABBREVIATES: the
*    contreat-in-controls token tokenises controls(), coarsecohortcontrols(), fe() and
*    ccc_absorb() and looks for the treatment variable among the tokens, but
*    coarsecohortcontrols() is declared string and so is NOT unabbreviated by Stata before the
*    gate sees it - c.dose_var#c.x1 with a dose named dose_vary walked past the literal
*    comparison and the run COLLAPSED when it had to refuse. The gate now runs cap unab on
*    each token. Four legs: the exact name, the abbreviation (the regression), an unrelated
*    pair (no false positive), and the varlist route. (b) THE CELL MEAN IS MISSING-AWARE PER
*    VARIABLE: each averaged variable carries its own weight sum over the rows where THAT
*    variable is non-missing. Both halves are pinned - the caller's guarantee that no missing
*    dose can reach the worker, and wooldid_fcol called DIRECTLY on a constructed fixture,
*    which is the only way to exercise a defensive guard.
* F13 (BS-28, 2026-09-07) THE ONE-REALIZED-SUBGROUP LEVEL MARGINS TAKE THE SCREEN. The
* dedicated margins call that runs when only one subgroup is realized, and the pre-period
* fallback after a deferred pwcompare, now sit inside the same two-pass wooldid_fcdeg
* screen the pwcompare pass uses. All three legs assert the guard did NOT fire and the
* numbers match nofastcollapse: a spurious fire on a healthy run is the risk here.
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

* ---- the two arms: default (collapsed) and nofastcollapse (not)
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests
matrix bD = e(b)
matrix VD = e(V)
matrix eD = e(b_es_att)
scalar jD = e(joint_horiz)
local fcolD = e(wdid_fcol)
local rcD = _rc

use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) jointtests nofastcollapse
matrix bS = e(b)
matrix VS = e(V)
matrix eS = e(b_es_att)
scalar jS = e(joint_horiz)
local fcolS = e(wdid_fcol)
local rcS = _rc

* ---- G1: the default really engages. If this ever reads 0 differences, the default has
* silently stopped applying and every equivalence check below is meaningless.
local ddefslow = mreldif(bD, bS)
di "G1_default_fcol=" `fcolD' " slow_fcol=" `fcolS' " default_vs_slow_mreldif=" `ddefslow'
if `fcolD' != 1 | `fcolS' != 0 {
  di as error "G1 FAIL: e(wdid_fcol) should be 1/0 for default/nofastcollapse, got `fcolD'/`fcolS'"
  local ++fails
}
if `ddefslow' == 0 {
  di as error "G1 FAIL: the default produced results identical to nofastcollapse - the collapse is NOT engaging, and every check below passes vacuously"
  local ++fails
}

* ---- G2: THE CLEAN BREAK (option consolidation, author ruling 2026-09-16). fastcollapse is
* gone and fullprecisionslow was renamed rather than aliased, so an archived do-file naming
* either stops with r(198) instead of running something the author did not ask for. The new
* spelling takes no abbreviation, the convention both retired options used, because
* nofastcollapse shortened to nofast is read by Stata's syntax parser as the negation of
* fastnose and would be swallowed without a word (measured 2026-09-16).
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att robust fastcollapse
local gx1 = _rc
cap qui wooldid y cohid year ttre, att robust fullprecisionslow
local gx2 = _rc
cap qui wooldid y cohid year ttre, att robust nofastc
local gx3 = _rc
di "G2_retired fastcollapse=" `gx1' " fullprecisionslow=" `gx2' " abbrev_nofastc=" `gx3' " (want 198 x3)"
if `gx1' != 198 | `gx2' != 198 | `gx3' != 198 {
  di as error "G2 FAIL: a retired or abbreviated spelling was accepted (rcs `gx1'/`gx2'/`gx3') - the break must be clean"
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
qui wooldid y cohid year ttre, att robust esprelength(2) espostlength(2) fastnose nofastcollapse
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

* ---- F4: THE DISCLOSURE. Each of these configurations is ineligible for the collapse.
* Until the 2026-09-16 option consolidation the way to learn that was to type fastcollapse
* and read the r(198) it raised; that option is gone and its diagnostic is now a caution
* every such run prints, which is what this leg asserts. Each gate must (a) run to rc 0,
* (b) report e(wdid_fcol) = 0, and (c) name its own reason in the caution. safety(off) is
* carried where the run needs it for an unrelated screen and must not silence the caution.
* Gate 1 was a bare contreat() run until 2026-09-02, when contreat() itself left the
* ineligible list (F11); it now carries contreatpoly(2), one of the two contreat reads that
* ARE still gated, so this leg keeps covering the contreat family on dta_cont. The
* acceptance side of the new contract is pinned in F11c.
local F4D1 `REPO'/test/data/dta_cont.dta
local F4V1 y i t ttre
local F4O1 att contreat(intensity) contreatcontrols(basic) contreatpoly(2) cluster(i) safety(off)
local F4W1 contreatpoly(>1)
local F4D2 `ADMIN'
local F4V2 y cohid year ttre
local F4O2 att robust interactivecontrols(x1 x2) safety(off)
local F4W2 interactivecontrols()
local F4D3 `ADMIN'
local F4V3 y cohid year ttre
local F4O3 att robust unconditionalse safety(off)
local F4W3 unconditionalse
* Gate 4 is the treatment variable inside controls(): like contreatpoly(>1), it is a
* contreat read that is NOT linear in the dose inside a cell. contreatpoly stays at its
* default of 1 here so that contreat-in-controls is the FIRST token in fcbad and the caution
* test below matches on it. The path is written out because the CCD local is defined at F11.
local F4D4 `REPO'/test/suite/data/continuous_control_dense.dta
local F4V4 y_cont_vary cohort t ttre
local F4O4 att contreat(dose_vary) contreatcontrols(basic) controls(c.dose_vary#c.x1) cluster(person_id) safety(off)
local F4W4 contreat-in-controls
* 2026-09-06 (F062): impvalues left the ineligible list - the cluster variable joins the
* cell key - so an impvalues run collapses (see F15); the gate count below keeps its shape
* with customvce() as gate 5's stand-in
local F4D5 `ADMIN'
local F4V5 y cohid year ttre
local F4O5 att customvce(cluster state) safety(off)
local F4W5 customvce()
* The five gates are the five reasons on the ineligible list, one each.
forval g = 1/5 {
  use "`F4D`g''", clear
  qui log using fcol_f4.txt, text replace name(fc)
  cap noi wooldid `F4V`g'' , `F4O`g''
  local g`g' = _rc
  local q`g' = e(wdid_fcol)
  qui log close fc
  local lg = fileread("fcol_f4.txt")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  local c`g' = strpos(`"`lg'"', "Caution: the default cell collapse does not apply with: `F4W`g''") > 0
  local r`g' = strpos(`"`lg'"', "Specify nofastcollapse to silence this") > 0
}
cap erase fcol_f4.txt
di "F4_rcs contreatpoly=" `g1' " interactive=" `g2' " uncond=" `g3' " ctl=" `g4' " customvce=" `g5' " (want 0 x5)"
di "F4_fcol=" `q1' "/" `q2' "/" `q3' "/" `q4' "/" `q5' " (want 0 x5)  caution=" `c1' "/" `c2' "/" `c3' "/" `c4' "/" `c5' " remedy=" `r1' "/" `r2' "/" `r3' "/" `r4' "/" `r5' " (want 1 x5)"
forval g = 1/5 {
  if `g`g'' != 0 {
    di as error "F4 FAIL: gate `g' (`F4W`g'') returned rc `g`g'' (want 0) - an ineligible request must run, not error"
    local ++fails
  }
  if `q`g'' != 0 {
    di as error "F4 FAIL: gate `g' (`F4W`g'') reports e(wdid_fcol)=`q`g'' (want 0)"
    local ++fails
  }
  if !`c`g'' {
    di as error "F4 FAIL: gate `g' skipped the collapse without a caution naming `F4W`g''"
    local ++fails
  }
  if !`r`g'' {
    di as error "F4 FAIL: gate `g' (`F4W`g'') cautioned without naming nofastcollapse as the way to silence it"
    local ++fails
  }
}
* the off-switch is the ONE thing that silences it: the collapse is off because the user
* asked for that, so there is nothing to disclose
use "`ADMIN'", clear
qui log using fcol_f4.txt, text replace name(fc)
cap noi wooldid y cohid year ttre, att robust interactivecontrols(x1 x2) nofastcollapse
local g7 = _rc
qui log close fc
local lg = fileread("fcol_f4.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local c7 = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
cap erase fcol_f4.txt
di "F4_offswitch_silent rc=" `g7' " (want 0) caution_present=" `c7' " (want 0)"
if `g7' != 0 | `c7' {
  di as error "F4 FAIL: nofastcollapse rc `g7', caution present `c7' - the user's own switch must say nothing"
  local ++fails
}
* subgroups() is NOT on the ineligible list: the run collapses and says nothing
use "`ADMIN'", clear
qui log using fcol_f4.txt, text replace name(fc)
cap noi wooldid y cohid year ttre, att cluster(state) subgroups(region) safety(off)
local g8 = _rc
local f8 = e(wdid_fcol)
qui log close fc
local lg = fileread("fcol_f4.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local c8 = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
cap erase fcol_f4.txt
di "F4_subgroups rc=" `g8' " (want 0)  fcol=" `f8' " (want 1) caution_present=" `c8' " (want 0)"
if `g8' != 0 | `f8' != 1 | `c8' {
  di as error "F4 FAIL: subgroups() rc `g8' fcol `f8' caution `c8' - subgroups() must be eligible and quiet"
  local ++fails
}

* ---- F5: THE DEFAULT-ON CONTRACT. Every one of those options must run NORMALLY and simply
* not collapse - a run wooldid cannot collapse must not start erroring because the collapse
* became the default - and a subgroups() run must collapse (fcol 1) and complete. F4 asserts
* the caution each of these prints; this leg asserts they still RUN. (unconditionalse needs
* robust/cluster and region varies within cohid, so the unconditionalse and subgroups(region)
* commands carry robust / safety(off): each runs rc 0 under nofastcollapse, measured
* 2026-09-01, so a nonzero rc here is the default collapse and nothing else.)
* Option 1 followed gate 1 on 2026-09-02: a bare contreat() run now COLLAPSES, so the
* silent-skip leg carries contreatpoly(2) and the bare run moved to its own leg (s8) below,
* which keeps the original protection - an ordinary contreat() run must not start erroring.
use "`REPO'/test/data/dta_cont.dta", clear
cap qui wooldid y i t ttre, att contreat(intensity) contreatcontrols(basic) contreatpoly(2) cluster(i)
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
* option 4 is F4's gate 4, the treatment variable inside controls()
use "`REPO'/test/suite/data/continuous_control_dense.dta", clear
cap qui wooldid y_cont_vary cohort t ttre, att contreat(dose_vary) contreatcontrols(basic) controls(c.dose_vary#c.x1) cluster(person_id)
local s4 = _rc
local f4 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att customvce(cluster state)
local s5 = _rc
local f5 = e(wdid_fcol)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) subgroups(region) safety(off)
local s7 = _rc
local f7 = e(wdid_fcol)
di "F5_skip rcs=" `s1' "/" `s2' "/" `s3' "/" `s4' "/" `s5' "/" `s7' "  fcol=" `f1' "/" `f2' "/" `f3' "/" `f4' "/" `f5' "/" `f7' " (want 0 x5 then 1)"
forval k = 1/5 {
  if `s`k'' != 0 {
    di as error "F5 FAIL: ineligible option `k' errored (rc `s`k'') - it must run normally and skip the collapse"
    local ++fails
  }
  if `f`k'' != 0 {
    di as error "F5 FAIL: ineligible option `k' reports e(wdid_fcol)=`f`k'' - it must be 0"
    local ++fails
  }
}
if `s7' != 0 | `f7' != 1 {
  di as error "F5 FAIL: subgroups() run rc `s7' fcol `f7' - it must complete and collapse (fcol 1)"
  local ++fails
}
* s8: the bare contreat() run this section used to hold. It must still complete, and since
* 2026-09-02 it must also collapse - the same shape as the subgroups() leg above.
use "`REPO'/test/data/dta_cont.dta", clear
cap qui wooldid y i t ttre, att contreat(intensity) contreatcontrols(basic) cluster(i)
local sct8 = _rc
local fct8 = e(wdid_fcol)
di "F5_contreat_default rc=" `sct8' " (want 0) fcol=" `fct8' " (want 1)"
if `sct8' != 0 | `fct8' != 1 {
  di as error "F5 FAIL: bare contreat() run rc `sct8' fcol `fct8' - it must complete and collapse (fcol 1)"
  local ++fails
}

* ---- F6: the off-switch reaches the replication engines (and the default runs them).
* dta_binary is used for permute because permute() requires clusters to be subsets of
* the treatment-onset groups, which dta_s_admin does not satisfy. PLUMBING checks: a rep
* that could not parse the forwarded option would abort the run.
use "`REPO'/test/data/dta_binary.dta", clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) nofastcollapse
local p1 = _rc
use "`REPO'/test/data/dta_binary.dta", clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20)
local p2 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) bootstrap(standard) reps(15) rseed(7) nofastcollapse
local b1 = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) bootstrap(standard) reps(15) rseed(7)
local b2 = _rc
di "F6_reps permute_slow_rc=" `p1' " permute_default_rc=" `p2' " boot_slow_rc=" `b1' " boot_default_rc=" `b2'
if `p1' != 0 | `p2' != 0 | `b1' != 0 | `b2' != 0 {
  di as error "F6 FAIL: replication rcs `p1'/`p2'/`b1'/`b2'"
  local ++fails
}

* ---- F7: poisson on the dydx (log) scale collapses and agrees. This leg used to carry a
* second half asserting that poisexpresults - the exponentiated reporting scale, which is
* NOT linear in the coefficients - ran uncollapsed and disclosed the skip. That option was
* removed on 2026-09-20 (author ruling), and with it the only poisson reporting scale the
* collapse could not reproduce: every margins call wooldid now makes under poisson is a
* plain dydx or dyex of the linear index. What remains is the agreement check, which is the
* half that was always load-bearing.
use "`ADMIN'", clear
qui wooldid ycount cohid year ttre, att poisson robust esprelength(2) espostlength(2) nofastcollapse
matrix bP = e(b)
local rp1 = _rc
use "`ADMIN'", clear
cap qui wooldid ycount cohid year ttre, att poisson robust esprelength(2) espostlength(2)
local rp2 = _rc
matrix bQ = e(b)
local fcolQ = e(wdid_fcol)
local dpq = mreldif(bP, bQ)
di "F7_poisson_rc=" `rp1' "/" `rp2' " fcol=" `fcolQ' " mreldif=" `dpq'
if `rp1' != 0 | `rp2' != 0 | `fcolQ' != 1 | `dpq' > 1e-8 {
  di as error "F7 FAIL: poisson dydx rc `rp1'/`rp2', fcol `fcolQ', mreldif=`dpq'"
  local ++fails
}

* ---- F8: REGRESSION. Before the subgroup screen existed this exact command returned a
* hard r(503) WITH the collapse and rc 0 without it: collapsed, the att and att_it
* estimands lost different numbers of horizons to the omission-marker drop, and the
* non-conformable concatenation aborted the run. The screen catches the singular VCE
* that drives that; under safety(off) the run stays collapsed with the caution
* (e(wdid_fcguard) = caution) and must COMPLETE. The off-switch arm is the control: the
* same command with nofastcollapse never reaches the collapse, so it completes with
* e(wdid_fcol) = 0 and an empty guard.
use "`REPO'/test/suite/data/survey_rcs.dta", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], att att_it subgroups(group) ///
    fe(group_time) esprelength(2) espostlength(2) jointtests cluster(state) safety(off) nofastcollapse
local r8a = _rc
local f8a = e(wdid_fcol)
local g8a "`e(wdid_fcguard)'"
use "`REPO'/test/suite/data/survey_rcs.dta", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], att att_it subgroups(group) ///
    fe(group_time) esprelength(2) espostlength(2) jointtests cluster(state) safety(off)
local r8b = _rc
local f8b = e(wdid_fcol)
local g8b "`e(wdid_fcguard)'"
di "F8_r503regression offswitch=" `r8a' " (want 0) fcol=`f8a' (want 0) guard=`g8a' (want empty)  default=" `r8b' " (want 0, was 503) default_fcol=" `f8b' " (want 1) guard=`g8b' (want caution)"
if `r8a' != 0 | `f8a' != 0 | "`g8a'" != "" | `r8b' != 0 | `f8b' != 1 | "`g8b'" != "caution" {
  di as error "F8 FAIL: rc `r8a'/`r8b' fcol `f8a'/`f8b' guard `g8a'/`g8b' - the r(503) regression is not pinned"
  local ++fails
}

* ---- F9: THE SUBGROUP SCREEN. Three runs of one degenerate command and two healthy ones.
* The degenerate command is survey_rcs with the cell shock removed and fe(i.group##i.t):
* the clustered VCE is singular in the subgroup-contrast direction, so the collapsed and
* full-data reads at the subgroup sites disagree beyond residue (a measured divergence,
* 2026-09-01). Under the default safety the screen must fall back to the full data at the
* first flagged site, and everything read after it must be BIT-IDENTICAL to
* nofastcollapse (mreldif exactly 0 on the posted subgroup matrices), while the main
* results read before it stay within residue. Under safety(off) the run stays collapsed
* with the caution. The healthy command must pass the screen untouched: guard empty,
* collapsed, and within residue of nofastcollapse.
local RCS "`REPO'/test/suite/data/survey_rcs.dta"
local DEG "att subgroups(group) fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state)"
use "`RCS'", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], `DEG'
local r9a = _rc
local f9a = e(wdid_fcol)
local g9a "`e(wdid_fcguard)'"
cap matrix B9a = e(b)
cap matrix P9a = e(pwsgs)
cap matrix S9a1 = e(wdidsgresults_1)
cap matrix S9a2 = e(wdidsgresults_2)
use "`RCS'", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], `DEG' nofastcollapse
local r9u = _rc
local g9u "`e(wdid_fcguard)'"
cap matrix B9u = e(b)
cap matrix P9u = e(pwsgs)
cap matrix S9u1 = e(wdidsgresults_1)
cap matrix S9u2 = e(wdidsgresults_2)
local d9b = mreldif(B9a, B9u)
local d9p = mreldif(P9a, P9u)
local d9s = max(mreldif(S9a1, S9u1), mreldif(S9a2, S9u2))
use "`RCS'", clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], `DEG' safety(off)
local r9c = _rc
local f9c = e(wdid_fcol)
local g9c "`e(wdid_fcguard)'"
di "F9_screen default: rc=" `r9a' " fcol=" `f9a' " guard=`g9a' (want fallback)  slow: rc=" `r9u' " guard=`g9u' (want empty)"
di "F9_screen post-fallback mreldif pwsgs=" `d9p' " sgresults=" `d9s' " (want exactly 0)  pre-flag mreldif b=" `d9b' " (want < 1e-8)"
di "F9_screen safety(off): rc=" `r9c' " fcol=" `f9c' " guard=`g9c' (want caution)"
if `r9a' != 0 | `f9a' != 1 | "`g9a'" != "fallback" | `r9u' != 0 | "`g9u'" != "" {
  di as error "F9 FAIL: default-safety degenerate run rc `r9a' fcol `f9a' guard `g9a' (slow rc `r9u' guard `g9u') - the screen must fall back"
  local ++fails
}
if `d9p' != 0 | `d9s' != 0 | `d9b' >= 1e-8 {
  di as error "F9 FAIL: after the fallback pwsgs/sgresults must be bit-identical to nofastcollapse (mreldif `d9p'/`d9s') and b within residue (`d9b')"
  local ++fails
}
if `r9c' != 0 | `f9c' != 1 | "`g9c'" != "caution" {
  di as error "F9 FAIL: safety(off) degenerate run rc `r9c' fcol `f9c' guard `g9c' - it must stay collapsed with the caution"
  local ++fails
}
* Two healthy runs, both on aggregating data. (i) binary_balanced under the conventional
* VCE: balanced cohorts, so the subgroup contrasts share one se legitimately - the screen
* must not mistake that for degeneracy. (ii) the admin data at the cohort level with a
* subgroup that splits cohorts (share weighting under safety(off)): the cell collapse's
* headline case, where the subgroup reads are 3-5x faster collapsed (measured 2026-09-01:
* 2.1 s vs 9.7 s) and must agree to residue with an empty guard.
local BB "`REPO'/test/suite/data/binary_balanced.dta"
use "`BB'", clear
cap qui wooldid y cohort t ttre, att att_it esprelength(2) espostlength(2) subgroups(sg_three) jointtests
local r9h = _rc
local f9h = e(wdid_fcol)
local g9h "`e(wdid_fcguard)'"
cap matrix P9h = e(pwsgs)
cap matrix S9h = e(wdidsgresults_1)
use "`BB'", clear
cap qui wooldid y cohort t ttre, att att_it esprelength(2) espostlength(2) subgroups(sg_three) jointtests nofastcollapse
local r9hu = _rc
cap matrix P9hu = e(pwsgs)
cap matrix S9hu = e(wdidsgresults_1)
local d9hp = mreldif(P9h, P9hu)
local d9hs = mreldif(S9h, S9hu)
di "F9_healthy rc=" `r9h' "/" `r9hu' " fcol=" `f9h' " (want 1) guard=`g9h' (want empty)  mreldif pwsgs=" `d9hp' " sgresults=" `d9hs' " (want < 1e-8)"
if `r9h' != 0 | `r9hu' != 0 | `f9h' != 1 | "`g9h'" != "" | `d9hp' >= 1e-8 | `d9hs' >= 1e-8 {
  di as error "F9 FAIL: healthy subgroup run rc `r9h'/`r9hu' fcol `f9h' guard `g9h' mreldif `d9hp'/`d9hs' - it must collapse untouched and agree"
  local ++fails
}
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att att_it subgroups(region) esprelength(2) espostlength(2) safety(off)
local r9k = _rc
local f9k = e(wdid_fcol)
local g9k "`e(wdid_fcguard)'"
cap matrix P9k = e(pwsgs)
cap matrix S9k = e(wdidsgresults_1)
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att att_it subgroups(region) esprelength(2) espostlength(2) safety(off) nofastcollapse
local r9ku = _rc
cap matrix P9ku = e(pwsgs)
cap matrix S9ku = e(wdidsgresults_1)
local d9kp = mreldif(P9k, P9ku)
local d9ks = mreldif(S9k, S9ku)
di "F9_healthy_admin rc=" `r9k' "/" `r9ku' " fcol=" `f9k' " (want 1) guard=`g9k' (want empty)  mreldif pwsgs=" `d9kp' " sgresults=" `d9ks' " (want < 1e-8)"
if `r9k' != 0 | `r9ku' != 0 | `f9k' != 1 | "`g9k'" != "" | `d9kp' >= 1e-8 | `d9ks' >= 1e-8 {
  di as error "F9 FAIL: healthy admin subgroup run rc `r9k'/`r9ku' fcol `f9k' guard `g9k' mreldif `d9kp'/`d9ks' - it must collapse untouched and agree"
  local ++fails
}

* ---- F10: THE 2026-09-01 ADVERSARIAL PROBE'S THREE PINS (7 angles, ~500 configs).
* (a) FALSE ALARM REMOVED: an early-treated cohort used as a subgroup has a pre-treatment
*     window holding only the reference period, so its pre-trend joint test is absent in
*     both arms (b = 0, V = 0 structurally). The overjoint screen used to flag that as
*     "zero variance" and cost the speedup at default safety; it is now gated on r(p).
*     Pin: guard empty, collapsed, bit-identical to nofastcollapse.
local DB "`REPO'/test/data/dta_binary.dta"
use "`DB'", clear
qui sum t, meanonly
local t2 = r(min) + 1
qui gen byte sgd2 = 1 + (mod(i, 7) == 0)
qui replace ttre = `t2' if sgd2 == 2
cap qui wooldid y i t ttre, att robust esprelength(2) espostlength(2) subgroups(sgd2) jointtests
local r10a = _rc
local f10a = e(wdid_fcol)
local g10a "`e(wdid_fcguard)'"
local j10a "`e(joint_horiz_s1)'"
local k10a "`e(joint_horiz_s2)'"
cap matrix P10a = e(pwsgs)
cap qui wooldid y i t ttre, att robust esprelength(2) espostlength(2) subgroups(sgd2) jointtests nofastcollapse
local r10au = _rc
local j10au "`e(joint_horiz_s1)'"
local k10au "`e(joint_horiz_s2)'"
cap matrix P10au = e(pwsgs)
local d10a = mreldif(P10a, P10au)
di "F10a rc=" `r10a' "/" `r10au' " fcol=" `f10a' " (want 1) guard=`g10a' (want empty) joint_horiz_s1 `j10a'/`j10au' s2 [`k10a']/[`k10au'] (want both empty) mreldif pwsgs=" `d10a'
if `r10a' != 0 | `r10au' != 0 | `f10a' != 1 | "`g10a'" != "" | "`j10a'" != "`j10au'" | "`k10a'`k10au'" != "" | `d10a' >= 1e-8 {
  di as error "F10a FAIL: reference-only subgroup pre window rc `r10a'/`r10au' fcol `f10a' guard `g10a' joint `j10a'/`j10au' [`k10a']/[`k10au'] mreldif `d10a' - the absent joint test must not trip the screen"
  local ++fails
}
* (b) SLIP CLOSED: a rank-1 subgroup-contrast VCE whose largest eigenvalue is itself ~1e-8
*     of the level variance passed the contrast-relative singular-direction test and kept a
*     collapsed joint test that differed from the full data by an order of magnitude
*     (1.5e-82 vs 3.3e-83; a known-degenerate dataset). The test is now anchored to the
*     level scale too. Pin: the att read flags (guard=fallback) and the all-subgroups-equal
*     joint tests match nofastcollapse exactly.
use "`BB'", clear
qui gen double pw1 = 1
qui replace pw1 = 1 + 1e-3 if person_id == 101
qui replace pw1 = 1 + 2e-3 if person_id == 201
cap qui wooldid y cohort t ttre [pw=pw1], att att_it esprelength(2) espostlength(2) subgroups(sg_three) cluster(cohort)
local r10b = _rc
local g10b "`e(wdid_fcguard)'"
local ja10b = e(joint_att_allsgssame)
local jp10b = e(joint_pre_att_allsgssame)
cap qui wooldid y cohort t ttre [pw=pw1], att att_it esprelength(2) espostlength(2) subgroups(sg_three) cluster(cohort) nofastcollapse
local r10bu = _rc
local ja10bu = e(joint_att_allsgssame)
local jp10bu = e(joint_pre_att_allsgssame)
local da10b = reldif(`ja10b', `ja10bu')
local dp10b = reldif(`jp10b', `jp10bu')
di "F10b rc=" `r10b' "/" `r10bu' " guard=`g10b' (want fallback) joint_att reldif=" `da10b' " joint_pre_att reldif=" `dp10b' " (want < 1e-8)"
if `r10b' != 0 | `r10bu' != 0 | "`g10b'" != "fallback" | `da10b' >= 1e-8 | `dp10b' >= 1e-8 {
  di as error "F10b FAIL: rank-1 contrast VCE rc `r10b'/`r10bu' guard `g10b' reldif `da10b'/`dp10b' - the att read must flag and the joint tests must match the full data"
  local ++fails
}
* (c) FALLBACK REACHES THE BOOTSTRAP REPS: reps run under fastnose, where the screen has no
*     r(V) to read, so after an observed-run fallback they are sent to the full data
*     outright. Before, on this exactly degenerate contrast the collapsed replicates were
*     bit-constant and the bootstrapped subgroup-equality test came back unavailable while
*     nofastcollapse reported one. Pin: guard=fallback and pwsgs / joint test / e(V)
*     match nofastcollapse exactly under the same seed.
local RCS "`REPO'/test/suite/data/survey_rcs.dta"
use "`RCS'", clear
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) fe(group_time) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) bootstrap(standard) reps(5) rseed(123) dots(0)
local r10c = _rc
local g10c "`e(wdid_fcguard)'"
local ja10c "`e(joint_att_allsgssame)'"
cap matrix P10c = e(pwsgs)
cap matrix V10c = e(V)
use "`RCS'", clear
cap qui wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) fe(group_time) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) bootstrap(standard) reps(5) rseed(123) dots(0) nofastcollapse
local r10cu = _rc
local ja10cu "`e(joint_att_allsgssame)'"
cap matrix P10cu = e(pwsgs)
cap matrix V10cu = e(V)
local dp10c = mreldif(P10c, P10cu)
local dv10c = mreldif(V10c, V10cu)
di "F10c rc=" `r10c' "/" `r10cu' " guard=`g10c' (want fallback) joint_att `ja10c'/`ja10cu' (want equal) mreldif pwsgs=" `dp10c' " V=" `dv10c' " (want < 1e-8)"
if `r10c' != 0 | `r10cu' != 0 | "`g10c'" != "fallback" | "`ja10c'" != "`ja10cu'" | `dp10c' >= 1e-8 | `dv10c' >= 1e-8 {
  di as error "F10c FAIL: bootstrap after fallback rc `r10c'/`r10cu' guard `g10c' joint `ja10c'/`ja10cu' mreldif `dp10c'/`dv10c' - the reps must follow the observed run to the full data"
  local ++fails
}

* ---- F11: THE CONTINUOUS TREATMENT (2026-09-02). contreat() used to sit on the ineligible
* list outright. The collapse now engages for it: the user's cohort variable joins the cell
* key (without it the reference-period rows pool across every cohort sharing an onset, and a
* read weighted by w_att_it consumes a mean built under w_att - measured: a pre-trend joint p
* moving .0956 -> .1112), and the surviving row carries the cell's ESTIMAND-WEIGHTED MEAN of
* the dose and of wooldid's internal ringer copy, written onto every row of the cell before
* the first-row pick. What stays refused is the reads that are not linear in the dose. A dose
* bin - the lattice() variable, and the plot bins makecfxplots reads - joins the cell key.
* F11.0 non-vacuity; F11a the eligible reads agree with nofastcollapse; F11b every gated
* contreat token skips and stays bit-identical; F11c the explicit-request refusals and
* the plain-contreat acceptance; F11d the two provisos (a zero-weight pick row, a float dose).
local CCD "`REPO'/test/suite/data/continuous_control_dense.dta"
local DC  "`REPO'/test/data/dta_cont.dta"

* ---- F11.0: NON-VACUITY, twice over. (i) the collapse must actually aggregate under the
* contreat key - rows per (cohort, t) cell far above 1; (ii) the DOSE must vary inside those
* cells, which is what makes the carried mean load-bearing: with a constant-within-cell dose
* the surviving row's own draw already equals the mean and F11a would pass without the new
* code (measured on the pre-mean form, which this fixture breaks by 5.6% on the att se).
use "`CCD'", clear
qui count
local ctrows = r(N)
qui gegen long __ccell = group(cohort t)
qui su __ccell, meanonly
local ctcells = r(max)
local ctratio = `ctrows'/`ctcells'
qui gegen double __dsd = sd(dose_vary), by(cohort t)
qui su __dsd, meanonly
local ctdsd = r(min)
di "F11_0_rows=" `ctrows' " cells=" `ctcells' " rows_per_cell=" %6.2f `ctratio' " min_within_cell_dose_sd=" `ctdsd'
if `ctratio' < 10 {
  di as error "F11.0 FAIL: rows per contreat cell `ctratio' - the collapse barely aggregates here, so F11a would pass vacuously"
  local ++fails
}
if `ctdsd' <= 0 | `ctdsd' >= . {
  di as error "F11.0 FAIL: min within-cell dose sd `ctdsd' - the dose must vary inside a cell or the carried mean is not exercised"
  local ++fails
}
* the float-dose proviso pin (F11d3) is only a pin while dta_cont's intensity is stored float
use "`DC'", clear
local ctftype : type intensity
di "F11_0_dta_cont_intensity_type=`ctftype' (want float)"
if "`ctftype'" != "float" {
  di as error "F11.0 FAIL: dta_cont intensity is `ctftype', not float - F11d3 no longer tests the recast-to-double proviso"
  local ++fails
}

* ---- F11a: THE COLLAPSE ENGAGES WITH contreat() AND AGREES WITH nofastcollapse.
* Five reads, each compared on e(b), e(V) and the reported table e(wdidmainresults):
* (1) a plain att run, (2) all four estimands with an event study and the joint tests,
* (3) contreatelasticitytype(dyex), (4) poisson on the default dydx scale, (5) subgroups(),
* (6) poisson with contreatelasticitytype(dyex) (eligible since 2026-09-06, F032: the analytic
* Jacobian form; the ame rows agree to 6.3e-16, the pretreatment CI bound to 4.5e-12),
* (7)-(9) lattice() and makecfxplots, whose dose bin joins the cell key (compared on every
* e() matrix and scalar at 1e-10, see the block before F11b).
* Residues measured 2026-09-02 sit at 0..9.1e-14, four orders inside the 1e-8 threshold, and
* the gate-off breakages of F11b's forms run 1e-5 to 0.8 relative - the same headroom F1 has.

* (1) plain att
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, att contreat(dose_vary) contreatcontrols(basic) robust
local ra1 = _rc
local fa1 = e(wdid_fcol)
matrix b11a1 = e(b)
matrix V11a1 = e(V)
matrix M11a1 = e(wdidmainresults)
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, att contreat(dose_vary) contreatcontrols(basic) robust nofastcollapse
local ra1u = _rc
local fa1u = e(wdid_fcol)
matrix b11a1u = e(b)
matrix V11a1u = e(V)
matrix M11a1u = e(wdidmainresults)
local da1b = mreldif(b11a1, b11a1u)
local da1v = mreldif(V11a1, V11a1u)
local da1m = mreldif(M11a1, M11a1u)
di "F11a1_plain rc=" `ra1' "/" `ra1u' " fcol=" `fa1' "/" `fa1u' " (want 1/0) mreldif b=" `da1b' " V=" `da1v' " table=" `da1m'
if `ra1' != 0 | `ra1u' != 0 | `fa1' != 1 | `fa1u' != 0 | `da1b' > 1e-8 | `da1v' > 1e-8 | `da1m' > 1e-8 {
  di as error "F11a1 FAIL: plain contreat rc `ra1'/`ra1u' fcol `fa1'/`fa1u' mreldif `da1b'/`da1v'/`da1m' - contreat() must collapse and agree"
  local ++fails
}

* (2) all four estimands, event study, joint tests. The joint-test p-values are the reads
* proviso (iii) protects, so they are compared as scalars (full precision) rather than macros.
local CT2 att att_it att_i att_itime contreat(dose_vary) contreatcontrols(basic) cluster(person_id) esprelength(2) espostlength(2) jointtests
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, `CT2'
local ra2 = _rc
local fa2 = e(wdid_fcol)
matrix b11a2 = e(b)
matrix V11a2 = e(V)
matrix M11a2 = e(wdidmainresults)
matrix E11a2 = e(b_es_att)
scalar jh2 = e(joint_horiz)
scalar jc2 = e(joint_Choriz_i)
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, `CT2' nofastcollapse
local ra2u = _rc
local fa2u = e(wdid_fcol)
matrix b11a2u = e(b)
matrix V11a2u = e(V)
matrix M11a2u = e(wdidmainresults)
matrix E11a2u = e(b_es_att)
scalar jh2u = e(joint_horiz)
scalar jc2u = e(joint_Choriz_i)
local da2b = mreldif(b11a2, b11a2u)
local da2v = mreldif(V11a2, V11a2u)
local da2m = mreldif(M11a2, M11a2u)
local da2e = mreldif(E11a2, E11a2u)
local da2jh = reldif(jh2, jh2u)
local da2jc = reldif(jc2, jc2u)
di "F11a2_all4es rc=" `ra2' "/" `ra2u' " fcol=" `fa2' "/" `fa2u' " (want 1/0) mreldif b=" `da2b' " V=" `da2v' " table=" `da2m' " es=" `da2e'
di "F11a2_joint reldif joint_horiz=" `da2jh' " joint_Choriz_i=" `da2jc' " (want < 1e-6; values " jc2 "/" jc2u ")"
if `ra2' != 0 | `ra2u' != 0 | `fa2' != 1 | `fa2u' != 0 | `da2b' > 1e-8 | `da2v' > 1e-8 | `da2m' > 1e-8 | `da2e' > 1e-8 {
  di as error "F11a2 FAIL: all-four-estimand contreat rc `ra2'/`ra2u' fcol `fa2'/`fa2u' mreldif `da2b'/`da2v'/`da2m'/`da2e'"
  local ++fails
}
if `da2jh' > 1e-6 | `da2jc' > 1e-6 {
  di as error "F11a2 FAIL: joint p reldif `da2jh'/`da2jc' - proviso (iii) (the user cohort in the key) is what keeps these exact"
  local ++fails
}
* structural identity, as F2 asserts for the binary path
if colsof(b11a2) != colsof(b11a2u) | rowsof(V11a2) != rowsof(V11a2u) {
  di as error "F11a2 FAIL: structural change - b or V changed dimensions under the collapse"
  local ++fails
}
local cn2 : colnames b11a2
local cn2u : colnames b11a2u
if "`cn2'" != "`cn2u'" {
  di as error "F11a2 FAIL: e(b) column names changed under the collapse"
  local ++fails
}

* (3) contreatelasticitytype(dyex) - linear in the dose, so it stays eligible
local CT3 att att_it contreat(dose_vary) contreatcontrols(basic) contreatelasticitytype(dyex) robust esprelength(2) espostlength(2)
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, `CT3'
local ra3 = _rc
local fa3 = e(wdid_fcol)
matrix b11a3 = e(b)
matrix V11a3 = e(V)
matrix M11a3 = e(wdidmainresults)
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, `CT3' nofastcollapse
local ra3u = _rc
matrix b11a3u = e(b)
matrix V11a3u = e(V)
matrix M11a3u = e(wdidmainresults)
local da3b = mreldif(b11a3, b11a3u)
local da3v = mreldif(V11a3, V11a3u)
local da3m = mreldif(M11a3, M11a3u)
di "F11a3_dyex rc=" `ra3' "/" `ra3u' " fcol=" `fa3' " (want 1) mreldif b=" `da3b' " V=" `da3v' " table=" `da3m'
if `ra3' != 0 | `ra3u' != 0 | `fa3' != 1 | `da3b' > 1e-8 | `da3v' > 1e-8 | `da3m' > 1e-8 {
  di as error "F11a3 FAIL: dyex rc `ra3'/`ra3u' fcol `fa3' mreldif `da3b'/`da3v'/`da3m' - dyex is linear in the dose and must collapse"
  local ++fails
}

* (6) poisson with contreatelasticitytype(dyex) - eligible since 2026-09-06 (F032): the
* analytic-Jacobian form under poisson makes the collapsed ame rows exact, so the run
* collapses without a word
local CT6 att poisson contreat(dose_vary) contreatcontrols(basic) contreatelasticitytype(dyex) robust esprelength(2) espostlength(2)
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, `CT6'
local ra6 = _rc
local fa6 = e(wdid_fcol)
cap matrix b11a6 = e(b)
cap matrix V11a6 = e(V)
cap matrix M11a6 = e(wdidmainresults)
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, `CT6' nofastcollapse
local ra6u = _rc
cap matrix b11a6u = e(b)
cap matrix V11a6u = e(V)
cap matrix M11a6u = e(wdidmainresults)
local da6b = mreldif(b11a6, b11a6u)
local da6v = mreldif(V11a6, V11a6u)
local da6m = mreldif(M11a6, M11a6u)
use "`CCD'", clear
qui log using fcol_x.txt, text replace name(fc)
cap noi wooldid y_count cohort t ttre, `CT6'
local ra6x = _rc
qui log close fc
local lg = fileread("fcol_x.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local ca6 = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
cap erase fcol_x.txt
di "F11a6_poisson_dyex rc=" `ra6' "/" `ra6u' "/" `ra6x' " fcol=" `fa6' " (want 1) caution=" `ca6' " (want 0) mreldif b=" `da6b' " V=" `da6v' " table=" `da6m'
if `ra6' != 0 | `ra6u' != 0 | `ra6x' != 0 | `ca6' | `fa6' != 1 | `da6b' > 1e-8 | `da6v' > 1e-8 | `da6m' > 1e-8 {
  di as error "F11a6 FAIL: poisson+dyex rc `ra6'/`ra6u'/`ra6x' caution `ca6' fcol `fa6' mreldif `da6b'/`da6v'/`da6m' - eligible since F032, so it must collapse exactly and say nothing"
  local ++fails
}

* (4) poisson on the default dydx scale (predict(xb), differentiated analytically)
local CT4 att poisson contreat(dose_vary) contreatcontrols(basic) robust esprelength(2) espostlength(2)
use "`CCD'", clear
qui wooldid y_count cohort t ttre, `CT4'
local ra4 = _rc
local fa4 = e(wdid_fcol)
matrix b11a4 = e(b)
matrix V11a4 = e(V)
matrix M11a4 = e(wdidmainresults)
use "`CCD'", clear
qui wooldid y_count cohort t ttre, `CT4' nofastcollapse
local ra4u = _rc
matrix b11a4u = e(b)
matrix V11a4u = e(V)
matrix M11a4u = e(wdidmainresults)
local da4b = mreldif(b11a4, b11a4u)
local da4v = mreldif(V11a4, V11a4u)
local da4m = mreldif(M11a4, M11a4u)
di "F11a4_poisson rc=" `ra4' "/" `ra4u' " fcol=" `fa4' " (want 1) mreldif b=" `da4b' " V=" `da4v' " table=" `da4m'
if `ra4' != 0 | `ra4u' != 0 | `fa4' != 1 | `da4b' > 1e-8 | `da4v' > 1e-8 | `da4m' > 1e-8 {
  di as error "F11a4 FAIL: poisson+contreat rc `ra4'/`ra4u' fcol `fa4' mreldif `da4b'/`da4v'/`da4m'"
  local ++fails
}

* (5) subgroups(). sg_three is constant within a cohort on this fixture, so it does not split
* a cell; the run must pass the subgroup screen untouched (guard empty) rather than fall back,
* which is what makes the comparison a test of the collapse and not of the fallback.
local CT5 att att_it contreat(dose_vary) contreatcontrols(basic) cluster(person_id) esprelength(2) espostlength(2) subgroups(sg_three)
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, `CT5'
local ra5 = _rc
local fa5 = e(wdid_fcol)
local ga5 "`e(wdid_fcguard)'"
matrix b11a5 = e(b)
matrix V11a5 = e(V)
matrix M11a5 = e(wdidmainresults)
matrix P11a5 = e(pwsgs)
matrix S11a5 = e(wdidsgresults_1)
use "`CCD'", clear
qui wooldid y_cont_vary cohort t ttre, `CT5' nofastcollapse
local ra5u = _rc
matrix b11a5u = e(b)
matrix V11a5u = e(V)
matrix M11a5u = e(wdidmainresults)
matrix P11a5u = e(pwsgs)
matrix S11a5u = e(wdidsgresults_1)
local da5b = mreldif(b11a5, b11a5u)
local da5v = mreldif(V11a5, V11a5u)
local da5m = mreldif(M11a5, M11a5u)
local da5p = mreldif(P11a5, P11a5u)
local da5s = mreldif(S11a5, S11a5u)
di "F11a5_subgroups rc=" `ra5' "/" `ra5u' " fcol=" `fa5' " (want 1) guard=[`ga5'] (want empty) mreldif b=" `da5b' " V=" `da5v' " table=" `da5m'
di "F11a5_subgroups mreldif pwsgs=" `da5p' " sgresults_1=" `da5s'
if `ra5' != 0 | `ra5u' != 0 | `fa5' != 1 | "`ga5'" != "" | `da5b' > 1e-8 | `da5v' > 1e-8 | `da5m' > 1e-8 | `da5p' > 1e-8 | `da5s' > 1e-8 {
  di as error "F11a5 FAIL: contreat+subgroups rc `ra5'/`ra5u' fcol `fa5' guard `ga5' mreldif `da5b'/`da5v'/`da5m'/`da5p'/`da5s'"
  local ++fails
}

* (7)-(9) THE DOSE BIN IN THE CELL KEY. lattice() gives every bin its own intercept and dose
* slope, and makecfxplots reads margins over() its plot bins, so both split a cohort-by-period
* cell. The worker puts the bin in the cell key - the lattice variable, and cfxlattice() under
* makecfxplots - so every read stays linear in the dose inside a cell-by-bin group. Each leg
* runs all four estimands under aweights and must (a) collapse (e(wdid_fcol) 1, e(wdid_fcskip)
* empty), (b) print no collapse caution, and (c) agree with nofastcollapse on EVERY posted e()
* matrix (mreldif, row and column names equal) and e() scalar (reldif) at 1e-10, both runs
* posting the same names. e(wdid_fcncell) is the collapse's own record and is left out.
* (7) an auto lattice built within cohort, with an event study and the joint tests;
* (8) a user-supplied lattice variable read by makecfxplots, safety(off) so that the cfx
*     matrices of all four estimands, not only att's, enter the comparison;
* (9) makecfxplots with cfxlattice() and no lattice(), so the plot bin alone splits the cell;
*     the bin is missing on some treated rows, whose rows must stay in the main reads (their
*     own cells) while the cfx reads omit them, exactly as on the full data.
cap program drop fc_esnap
program define fc_esnap
  * copies every posted e() matrix to <stub>m<k> and e() scalar to <stub>s<k>, k indexing the
  * SORTED name lists, which come back in the caller's locals <stub>_mn and <stub>_sn
  args stub
  local mn : e(matrices)
  local mn : list sort mn
  local sn : e(scalars)
  local fcx wdid_fcncell
  local sn : list sn - fcx
  local sn : list sort sn
  local k 0
  foreach m of local mn {
    local ++k
    matrix `stub'm`k' = e(`m')
  }
  local k 0
  foreach m of local sn {
    local ++k
    scalar `stub's`k' = e(`m')
  }
  c_local `stub'_mn `mn'
  c_local `stub'_sn `sn'
end
cap program drop fc_ecmp
program define fc_ecmp, rclass
  * r(worst): the largest mreldif/reldif over the matched names, missing when the two name
  * sets differ or a pair does not conform or differs in its row/column names; r(where) names
  * the worst item; r(nm)/r(ns) count the matrices and scalars compared
  syntax , a(name) b(name) amn(string) bmn(string) asn(string) bsn(string)
  local worst 0
  local where
  if "`amn'" != "`bmn'" | "`asn'" != "`bsn'" {
    return scalar worst = .
    return local where "e() name sets differ"
    return scalar nm = 0
    return scalar ns = 0
    exit
  }
  local k 0
  foreach m of local amn {
    local ++k
    local rna : rowfullnames `a'm`k'
    local rnb : rowfullnames `b'm`k'
    local cna : colfullnames `a'm`k'
    local cnb : colfullnames `b'm`k'
    local d .
    if `"`rna'"' == `"`rnb'"' & `"`cna'"' == `"`cnb'"' {
      cap local d = mreldif(`a'm`k', `b'm`k')
      if _rc local d .
    }
    if missing(`d') | `d' > `worst' {
      local worst = `d'
      local where `m'
    }
  }
  return scalar nm = `k'
  local k 0
  foreach m of local asn {
    local ++k
    local d = reldif(scalar(`a's`k'), scalar(`b's`k'))
    if missing(`d') | `d' > `worst' {
      local worst = `d'
      local where `m'
    }
  }
  return scalar ns = `k'
  return scalar worst = `worst'
  return local where `where'
end
local CT7 att att_it att_i att_itime contreat(dose_vary) contreatcontrols(basic) lattice(byCohort 3) cluster(person_id) esprelength(2) espostlength(2) jointtests
local CT8 att att_it att_i att_itime contreat(dose_vary) contreatcontrols(basic) lattice(custom_lattice) makecfxplots cfxplottypes(post pre 1) cluster(person_id) esprelength(2) espostlength(2) safety(off)
local CT9 att att_it att_i att_itime contreat(dose_vary) contreatcontrols(basic) makecfxplots cfxlattice(cfxb) cluster(person_id) esprelength(2) espostlength(2) safety(off)
local CN7 lattice(byCohort 3)
local CN8 lattice(custom_lattice)+makecfxplots
local CN9 makecfxplots+cfxlattice()
use "`CCD'", clear
qui gen cfxb = mod(unit_in_cohort, 3)
qui replace cfxb = . if ttre < . & mod(person_id, 11) == 0
tempfile f11bin
qui save `f11bin'
forval k = 7/9 {
  cap matrix drop fcdm* fcum*
  cap scalar drop fcds* fcus*
  local fcd_mn
  local fcu_mn
  local fcd_sn
  local fcu_sn
  local ncfx`k' = 0
  use `f11bin', clear
  qui log using fcol_x.txt, text replace name(fc)
  cap noi wooldid y_cont_vary cohort t ttre [aw=analytic_weight], `CT`k''
  local rb`k' = _rc
  local fb`k' = e(wdid_fcol)
  local sb`k' "`e(wdid_fcskip)'"
  local nb`k' = e(wdid_fcncell)
  if `rb`k'' == 0 fc_esnap fcd
  qui log close fc
  local lg = fileread("fcol_x.txt")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  local cb`k' = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
  use `f11bin', clear
  cap qui wooldid y_cont_vary cohort t ttre [aw=analytic_weight], `CT`k'' nofastcollapse
  local rbu`k' = _rc
  local fbu`k' = e(wdid_fcol)
  if `rbu`k'' == 0 fc_esnap fcu
  local db`k' = .
  local wb`k' "a run failed"
  local mb`k' = 0
  if `rb`k'' == 0 & `rbu`k'' == 0 {
    fc_ecmp , a(fcd) b(fcu) amn(`fcd_mn') bmn(`fcu_mn') asn(`fcd_sn') bsn(`fcu_sn')
    local db`k' = r(worst)
    local wb`k' "`r(where)'"
    local mb`k' = r(nm)
    foreach m of local fcd_mn {
      if substr("`m'", 1, 3) == "cfx" local ++ncfx`k'
    }
  }
  di "F11a`k' `CN`k'' rc=" `rb`k'' "/" `rbu`k'' " fcol=" `fb`k'' "/" `fbu`k'' " (want 1/0) fcskip=[`sb`k''] cells=" `nb`k'' " caution=" `cb`k'' " (want 0) matrices=" `mb`k'' " cfx=" `ncfx`k'' " worst=" `db`k'' " at [`wb`k''] (want <= 1e-10)"
  if `rb`k'' != 0 | `rbu`k'' != 0 | `fb`k'' != 1 | `fbu`k'' != 0 | "`sb`k''" != "" | `cb`k'' | missing(`db`k'') | `db`k'' > 1e-10 {
    di as error "F11a`k' FAIL: `CN`k'' rc `rb`k''/`rbu`k'' fcol `fb`k''/`fbu`k'' fcskip [`sb`k''] caution `cb`k'' worst `db`k'' at [`wb`k''] - the dose bin is in the cell key, so the run must collapse quietly and agree with nofastcollapse"
    local ++fails
  }
}
cap erase fcol_x.txt
* the cfx legs must actually carry cfx matrices for all four estimands (safety(off)), or (c)
* compared nothing the bin changes
forval k = 8/9 {
  if `ncfx`k'' < 4 {
    di as error "F11a`k' FAIL: only `ncfx`k'' cfx matrices posted - the leg does not exercise the plot bins"
    local ++fails
  }
}

* ---- F11b: EVERY GATED CONTREAT READ SKIPS, UNDER THE CAUTION F4 ASSERTS. Each of the gated
* reads must run normally, report e(wdid_fcol)=0, and - because it then takes the SAME
* code path as nofastcollapse - reproduce it BIT-IDENTICALLY (mreldif exactly 0, as measured
* 2026-09-02). Anything else means the gate leaked. TWO contreat tokens are gated:
* contreatpoly(>1) and the contreat variable inside controls(). lattice() and makecfxplots
* are eligible (the dose bin is in the cell key) and are asserted in F11a (7)-(9).
local BOPT1 att contreat(dose_vary) contreatcontrols(basic) contreatpoly(2) robust esprelength(2) espostlength(2)
local BOPT2 att att_it contreat(dose_vary) contreatcontrols(basic) controls(c.dose_vary#c.x1) cluster(person_id) esprelength(2) espostlength(2)
local BY1 y_cont_vary
local BY2 y_cont_vary
local BNM1 contreatpoly(2)
local BNM2 contreat-in-controls
forval k = 1/2 {
  * dropped first, so a run that posted nothing is a loud mreldif error rather than a
  * silent comparison of the previous iteration's matrices
  cap matrix drop b11b V11b M11b b11bu V11bu M11bu
  use "`CCD'", clear
  cap qui wooldid `BY`k'' cohort t ttre, `BOPT`k''
  local bq`k' = _rc
  local bf`k' = e(wdid_fcol)
  cap matrix b11b = e(b)
  cap matrix V11b = e(V)
  cap matrix M11b = e(wdidmainresults)
  use "`CCD'", clear
  cap qui wooldid `BY`k'' cohort t ttre, `BOPT`k'' nofastcollapse
  local bqu`k' = _rc
  cap matrix b11bu = e(b)
  cap matrix V11bu = e(V)
  cap matrix M11bu = e(wdidmainresults)
  local bd`k' = max(mreldif(b11b, b11bu), mreldif(V11b, V11bu), mreldif(M11b, M11bu))
  * every surviving leg is a gated read: it must skip the collapse (fcol 0) and therefore
  * take the same code path as nofastcollapse, reproducing it BIT-IDENTICALLY. The BS-6
  * exception that used to sit here went with contreatelasticitytype(eydx) on 2026-09-20.
  local wantfcol`k' 0
  local wanttol`k' 0
  di "F11b`k' `BNM`k'' rc=" `bq`k'' "/" `bqu`k'' " fcol=" `bf`k'' " (want `wantfcol`k'') worst mreldif b/V/table=" `bd`k'' " (want <= `wanttol`k'')"
  if `bq`k'' != 0 | `bqu`k'' != 0 {
    di as error "F11b FAIL: read `BNM`k'' errored (rc `bq`k''/`bqu`k'')"
    local ++fails
  }
  if `bf`k'' != `wantfcol`k'' {
    di as error "F11b FAIL: read `BNM`k'' reports e(wdid_fcol)=`bf`k'' (want `wantfcol`k'')"
    local ++fails
  }
  if `bd`k'' > `wanttol`k'' {
    di as error "F11b FAIL: read `BNM`k'' differs from nofastcollapse by `bd`k'' (want <= `wanttol`k'')"
    local ++fails
  }
}

* ---- F11c: each of those two gated reads DISCLOSES its skip in the caution, naming its own
* token. Until 2026-09-20 a sixth leg (poisson + eydx, eligible since BS-6) collapsed and said
* nothing, which is why this loop used to carry a silent case; that leg is gone. This is
* F4's contract on the continuous-treatment gate tokens: the old blanket contreat() refusal is
* gone, only the reads above are gated, and each gated one is now announced rather than either
* silent or fatal. safety(off) is carried as before and must not silence it.
local CW1 contreatpoly(>1)
local CW2 contreat-in-controls
forval k = 1/2 {
  use "`CCD'", clear
  qui log using fcol_f11c.txt, text replace name(fc)
  cap noi wooldid `BY`k'' cohort t ttre, `BOPT`k'' safety(off)
  local cx`k' = _rc
  local cf`k' = e(wdid_fcol)
  qui log close fc
  local lg = fileread("fcol_f11c.txt")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  local cc`k' = strpos(`"`lg'"', "Caution: the default cell collapse does not apply with: `CW`k''") > 0
}
cap erase fcol_f11c.txt
di "F11c_rcs poly2=" `cx1' " ctl=" `cx2' " (want 0 x2)"
di "F11c_caution=" `cc1' "/" `cc2' " (want 1 x2)"
forval k = 1/2 {
  if `cx`k'' != 0 {
    di as error "F11c FAIL: `BNM`k'' returned rc `cx`k'' (want 0)"
    local ++fails
  }
  if !`cc`k'' {
    di as error "F11c FAIL: `BNM`k'' skipped the collapse without a caution naming its token"
    local ++fails
  }
}
use "`CCD'", clear
qui log using fcol_f11c.txt, text replace name(fc)
cap noi wooldid y_cont_vary cohort t ttre, att contreat(dose_vary) contreatcontrols(basic) robust
local cx7 = _rc
local cf7 = e(wdid_fcol)
qui log close fc
local lg = fileread("fcol_f11c.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local cc7 = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
cap erase fcol_f11c.txt
di "F11c_plain rc=" `cx7' " (want 0) fcol=" `cf7' " (want 1) caution_present=" `cc7' " (want 0)"
if `cx7' != 0 | `cf7' != 1 | `cc7' {
  di as error "F11c FAIL: plain contreat() rc `cx7' fcol `cf7' caution `cc7' - contreat() left the ineligible list 2026-09-02 and must collapse quietly"
  local ++fails
}

* ---- F11d: THE TWO PROVISOS THE WORKER'S FORM EXISTS FOR.
* (1) and (2) pin the mean-on-EVERY-ROW form. The pick row is the cell's first in-sample row,
* chosen without reference to any weight, and customweightmultiplier() zeroes an estimand
* weight WITHOUT removing the row from e(sample) (a zero sampling weight cannot: Stata drops
* those rows). The earlier gcollapse (mean) [aw=w] merge form wrote MISSING onto zero-weight
* rows, so a zeroed pick row kept its raw dose: measured, a pre-trend horizon moving
* p .167 -> .0022 with its CI crossing zero. (1) is the brief's named case, cwm = 0 on the
* first row of every cell; (2) is its non-vacuity twin - cwm = 0 on ELEVEN of every twelve
* rows, so a zero-weight row is picked in 143 of the 144 cells whatever order wooldid leaves
* the data in, and the cell mean must come back as the surviving row's dose exactly.
* (3) pins the recast to double: dta_cont's intensity is stored float, and a float store of
* the cell mean injects ~2.4e-8 relative into every dose-scaled read.
local CTD att att_it att_i att_itime contreat(dose_vary) contreatcontrols(basic) cluster(person_id) esprelength(2) espostlength(2) jointtests
* the cwm fixtures address cells by row block, so first pin that 12-row blocks ARE the cells
use "`CCD'", clear
qui gegen long __ccell = group(cohort t)
qui gen long __blk = ceil(_n/12)
qui gegen long __bmin = min(__ccell), by(__blk)
qui gegen long __bmax = max(__ccell), by(__blk)
qui count if __bmin != __bmax
local dblkbad = r(N)
qui su __blk, meanonly
local dnblk = r(max)
di "F11d_blocks nblocks=" `dnblk' " (want 144) rows_spanning_two_cells=" `dblkbad' " (want 0)"
if `dnblk' != `ctcells' | `dblkbad' != 0 {
  di as error "F11d FAIL: the 12-row blocks the cwm fixtures use are not the collapse cells (`dnblk' blocks vs `ctcells' cells, `dblkbad' straddling rows)"
  local ++fails
}

* (1) cwm = 0 on the first row of every cell, with pw
use "`CCD'", clear
qui gen double w0f = cond(mod(_n,12)==1, 0, 1)
qui count if w0f == 0
local dnz1 = r(N)
qui wooldid y_cont_vary cohort t ttre [pw=probability_weight], `CTD' customweightmultiplier(w0f)
local rd1 = _rc
local fd1 = e(wdid_fcol)
matrix b11d1 = e(b)
matrix V11d1 = e(V)
matrix M11d1 = e(wdidmainresults)
use "`CCD'", clear
qui gen double w0f = cond(mod(_n,12)==1, 0, 1)
qui wooldid y_cont_vary cohort t ttre [pw=probability_weight], `CTD' customweightmultiplier(w0f) nofastcollapse
local rd1u = _rc
matrix b11d1u = e(b)
matrix V11d1u = e(V)
matrix M11d1u = e(wdidmainresults)
local dd1b = mreldif(b11d1, b11d1u)
local dd1v = mreldif(V11d1, V11d1u)
local dd1m = mreldif(M11d1, M11d1u)
di "F11d1_cwm0_firstrow zeroed=" `dnz1' " (want 144) rc=" `rd1' "/" `rd1u' " fcol=" `fd1' " (want 1) mreldif b=" `dd1b' " V=" `dd1v' " table=" `dd1m'
if `dnz1' != `ctcells' {
  di as error "F11d1 FAIL: the fixture zeroed `dnz1' rows, not one per cell (`ctcells')"
  local ++fails
}
if `rd1' != 0 | `rd1u' != 0 | `fd1' != 1 | `dd1b' > 1e-8 | `dd1v' > 1e-8 | `dd1m' > 1e-8 {
  di as error "F11d1 FAIL: cwm 0 on the first row of every cell rc `rd1'/`rd1u' fcol `fd1' mreldif `dd1b'/`dd1v'/`dd1m' - the mean must be written onto the zero-weight pick row too"
  local ++fails
}

* (2) cwm = 0 on all but the 5th row of every cell
use "`CCD'", clear
qui gen double w01 = cond(mod(_n,12)==5, 1, 0)
qui count if w01 == 0
local dnz2 = r(N)
qui wooldid y_cont_vary cohort t ttre, `CTD' customweightmultiplier(w01)
local rd2 = _rc
local fd2 = e(wdid_fcol)
matrix b11d2 = e(b)
matrix V11d2 = e(V)
matrix M11d2 = e(wdidmainresults)
use "`CCD'", clear
qui gen double w01 = cond(mod(_n,12)==5, 1, 0)
qui wooldid y_cont_vary cohort t ttre, `CTD' customweightmultiplier(w01) nofastcollapse
local rd2u = _rc
matrix b11d2u = e(b)
matrix V11d2u = e(V)
matrix M11d2u = e(wdidmainresults)
local dd2b = mreldif(b11d2, b11d2u)
local dd2v = mreldif(V11d2, V11d2u)
local dd2m = mreldif(M11d2, M11d2u)
di "F11d2_cwm0_allbut1 zeroed=" `dnz2' " (want " `ctrows'-`ctcells' ") rc=" `rd2' "/" `rd2u' " fcol=" `fd2' " (want 1) mreldif b=" `dd2b' " V=" `dd2v' " table=" `dd2m'
if `dnz2' != `ctrows' - `ctcells' {
  di as error "F11d2 FAIL: the fixture zeroed `dnz2' rows, not all but one per cell (`=`ctrows'-`ctcells'')"
  local ++fails
}
if `rd2' != 0 | `rd2u' != 0 | `fd2' != 1 | `dd2b' > 1e-8 | `dd2v' > 1e-8 | `dd2m' > 1e-8 {
  di as error "F11d2 FAIL: cwm 0 on eleven of every twelve rows rc `rd2'/`rd2u' fcol `fd2' mreldif `dd2b'/`dd2v'/`dd2m' - the picked zero-weight row must still carry the cell mean"
  local ++fails
}

* (3) a float-stored dose
use "`DC'", clear
qui wooldid y i t ttre, att att_it contreat(intensity) contreatcontrols(basic)
local rd3 = _rc
local fd3 = e(wdid_fcol)
matrix b11d3 = e(b)
matrix V11d3 = e(V)
matrix M11d3 = e(wdidmainresults)
use "`DC'", clear
qui wooldid y i t ttre, att att_it contreat(intensity) contreatcontrols(basic) nofastcollapse
local rd3u = _rc
matrix b11d3u = e(b)
matrix V11d3u = e(V)
matrix M11d3u = e(wdidmainresults)
local dd3b = mreldif(b11d3, b11d3u)
local dd3v = mreldif(V11d3, V11d3u)
local dd3m = mreldif(M11d3, M11d3u)
di "F11d3_float_dose rc=" `rd3' "/" `rd3u' " fcol=" `fd3' " (want 1) mreldif b=" `dd3b' " V=" `dd3v' " table=" `dd3m'
if `rd3' != 0 | `rd3u' != 0 | `fd3' != 1 | `dd3b' > 1e-8 | `dd3v' > 1e-8 | `dd3m' > 1e-8 {
  di as error "F11d3 FAIL: float dose rc `rd3'/`rd3u' fcol `fd3' mreldif `dd3b'/`dd3v'/`dd3m' - without the recast to double a float dose injects ~2.4e-8"
  local ++fails
}

* ---- F12: THE TWO GUARDS ADDED ON TOP OF THE CONTINUOUS-TREATMENT COLLAPSE (2026-09-02).
* F12a THE GATE UNABBREVIATES. The contreat-in-controls token refuses the collapse when the
* treatment variable ALSO appears inside controls(), coarsecohortcontrols(), fe() or
* ccc_absorb(): the dose derivative then varies inside a cell on a variable the collapse does
* not average. The gate tokenises those option strings and compares each token to the
* treatment variable's name - but controls(), fe() and ccc_absorb() are declared varlists, so
* Stata unabbreviates them before the gate sees them, while coarsecohortcontrols() is declared
* string and does NOT get unabbreviated. So c.dose_var#c.x1 with a dose named dose_vary walked
* past the literal comparison and the run COLLAPSED when it had to refuse. The gate now runs
* cap unab on each token and looks for the treatment variable in the expansion, falling back
* to the literal comparison when unab fails.
* F12b THE CELL MEAN IS MISSING-AWARE, PER VARIABLE. The worker used to build ONE shared
* weight sum and divide every averaged variable's sum(w*x) by it; each averaged variable now
* carries its own weight sum, restricted to the rows where THAT variable is non-missing. With
* a non-missing dose the two forms are the same number bit for bit, so this is a DEFENSIVE
* guard: it only matters if a missing value ever reaches the worker, where the shared form
* would divide a numerator taken over the non-missing rows by a weight total taken over all of
* them, and would silently rewrite an all-missing cell's dose to 0. A missing dose cannot
* reach the worker through wooldid - the caller zeroes the if-marker on any row with a missing
* contreat() value - so both halves are pinned: that the caller's guarantee holds, and then
* wooldid_fcol called DIRECTLY, which is the only way to exercise the guard itself.
local CCD "`REPO'/test/suite/data/continuous_control_dense.dta"

* ---- F12a.0: NON-VACUITY. F12a2 and F12a4 only test the unabbreviation while dose_var is an
* ABBREVIATION of dose_vary and not a variable in its own right. If the fixture ever grows a
* literal dose_var, both legs would be testing the plain literal comparison instead.
use "`CCD'", clear
cap confirm variable dose_var, exact
local gaex = _rc
cap unab gaunab : dose_var
local gaunrc = _rc
di "F12a_0_abbrev exact_rc=" `gaex' " (want nonzero) unab_rc=" `gaunrc' " (want 0) expands_to=[`gaunab'] (want dose_vary)"
if `gaex' == 0 | `gaunrc' != 0 | "`gaunab'" != "dose_vary" {
  di as error "F12a.0 FAIL: dose_var is not a bare abbreviation of dose_vary here (exact rc `gaex', unab rc `gaunrc', expansion [`gaunab']) - F12a2/F12a4 would not be exercising the unabbreviation at all"
  local ++fails
}

* ---- F12a: the four legs. Each is run three ways: the default (does the gate fire?),
* nofastcollapse (the reference), and an explicit fastcollapse (is the refusal loud?).
* Legs 1/2/4 must refuse, and with the collapse off they take the SAME code path as
* nofastcollapse, so they must be bit-identical to it (measured exactly 0, 2026-09-02).
* Leg 3 must NOT fire - coarsecohortcontrols(c.x1#c.x2) never mentions the dose - so it
* collapses and agrees to residue (measured 9.5e-16), and its explicit request is accepted.
* safety(off) rides on the three refusing legs to pin that it opens no escape.
local GA1 att contreat(dose_vary) contreatcontrols(basic) robust coarsecohortcontrols(c.dose_vary#c.x1)
local GA2 att contreat(dose_vary) contreatcontrols(basic) robust coarsecohortcontrols(c.dose_var#c.x1)
local GA3 att contreat(dose_vary) contreatcontrols(basic) robust coarsecohortcontrols(c.x1#c.x2)
local GA4 att contreat(dose_vary) contreatcontrols(basic) robust controls(c.dose_var#c.x1)
local GAN1 ccc-c.dose_vary#c.x1
local GAN2 ccc-c.dose_var#c.x1-ABBREV
local GAN3 ccc-c.x1#c.x2-unrelated
local GAN4 controls-c.dose_var#c.x1-ABBREV
local GAF1 0
local GAF2 0
local GAF3 1
local GAF4 0
* GAX: does the gated read announce itself? 1 = caution naming the token, 0 = silence
* (leg 3's term does not touch the dose, so the collapse engages and nothing is said)
local GAX1 1
local GAX2 1
local GAX3 0
local GAX4 1
local GAS1 safety(off)
local GAS2 safety(off)
local GAS3
local GAS4 safety(off)
local GAD1 0
local GAD2 0
local GAD3 1e-8
local GAD4 0
local GAR1 the gate is off, so both arms are the same code path and must be bit-identical
local GAR2 the gate is off, so both arms are the same code path and must be bit-identical
local GAR3 the collapse engages here and must agree with the full data to residue
local GAR4 the gate is off, so both arms are the same code path and must be bit-identical
forval k = 1/4 {
  * dropped first, so a run that posted nothing is a loud mreldif error rather than a
  * silent comparison of the previous iteration's matrices
  cap matrix drop b12a V12a M12a b12au V12au M12au
  use "`CCD'", clear
  cap qui wooldid y_cont_vary cohort t ttre, `GA`k''
  local gr`k' = _rc
  local gf`k' = e(wdid_fcol)
  cap matrix b12a = e(b)
  cap matrix V12a = e(V)
  cap matrix M12a = e(wdidmainresults)
  use "`CCD'", clear
  cap qui wooldid y_cont_vary cohort t ttre, `GA`k'' nofastcollapse
  local gru`k' = _rc
  cap matrix b12au = e(b)
  cap matrix V12au = e(V)
  cap matrix M12au = e(wdidmainresults)
  local gd`k' = max(mreldif(b12a, b12au), mreldif(V12a, V12au), mreldif(M12a, M12au))
  use "`CCD'", clear
  qui log using fcol_x.txt, text replace name(fc)
  cap noi wooldid y_cont_vary cohort t ttre, `GA`k'' `GAS`k''
  local gy`k' = _rc
  qui log close fc
  local lg = fileread("fcol_x.txt")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  local gx`k' = strpos(`"`lg'"', "Caution: the default cell collapse does not apply with: contreat-in-controls") > 0
  cap erase fcol_x.txt
  di "F12a`k' `GAN`k'' rc=" `gr`k'' "/" `gru`k'' "/" `gy`k'' " fcol=" `gf`k'' " (want `GAF`k'') worst mreldif b/V/table=" `gd`k'' " (want <= `GAD`k'') caution=" `gx`k'' " (want `GAX`k'')"
  if `gr`k'' != 0 | `gru`k'' != 0 {
    di as error "F12a FAIL: `GAN`k'' errored (rc `gr`k''/`gru`k'') - a gated read must run normally on the full data, never start erroring"
    local ++fails
  }
  if `gf`k'' != `GAF`k'' {
    di as error "F12a FAIL: `GAN`k'' reports e(wdid_fcol)=`gf`k'', want `GAF`k'' - the contreat-in-controls gate must see through the abbreviation and must not fire on an unrelated term"
    local ++fails
  }
  if `gd`k'' > `GAD`k'' {
    di as error "F12a FAIL: `GAN`k'' differs from nofastcollapse by `gd`k'' (want <= `GAD`k'') - `GAR`k''"
    local ++fails
  }
  if `gy`k'' != 0 | `gx`k'' != `GAX`k'' {
    di as error "F12a FAIL: `GAN`k'' rc `gy`k'' caution `gx`k'', want rc 0 and caution `GAX`k'' (safety(off) must not silence the disclosure)"
    local ++fails
  }
}

* ---- F12b1: THE CALLER'S GUARANTEE. wooldid zeroes the if-marker on any row whose contreat()
* value is missing, so no missing dose can reach the worker. Blank the dose on every 7th row
* and check that every one of those rows leaves e(sample) - under the default (collapsed) path
* and under nofastcollapse alike. Measured 2026-09-02: 1728 rows, 246 blanked, e(N) falling
* 1728 -> 1482 in both arms with no missing-dose row left in e(sample).
use "`CCD'", clear
qui count
local gbn = r(N)
qui wooldid y_cont_vary cohort t ttre, att contreat(dose_vary) contreatcontrols(basic) robust
local gbbase = e(N)
di "F12b1_baseline rows=" `gbn' " e(N)=" `gbbase' " (want equal)"
if `gbbase' != `gbn' {
  di as error "F12b1 FAIL: the unblanked run keeps `gbbase' of `gbn' rows - the e(N) delta below no longer isolates the missing-dose drop"
  local ++fails
}
local GBW att contreat(dose_vary) contreatcontrols(basic) robust
local GBN1 default
local GBN2 nofastcollapse
local GBF1 1
local GBF2 0
local GBO1
local GBO2 nofastcollapse
forval k = 1/2 {
  use "`CCD'", clear
  qui replace dose_vary = . if mod(_n,7) == 0
  qui count if missing(dose_vary)
  local gbm`k' = r(N)
  cap qui wooldid y_cont_vary cohort t ttre, `GBW' `GBO`k''
  local gbr`k' = _rc
  local gbf`k' = e(wdid_fcol)
  local gbN`k' = e(N)
  qui count if e(sample)
  local gbes`k' = r(N)
  qui count if e(sample) & missing(dose_vary)
  local gbbad`k' = r(N)
  local gbwant`k' = `gbn' - `gbm`k''
  di "F12b1_`GBN`k'' rc=" `gbr`k'' " fcol=" `gbf`k'' " (want `GBF`k'') blanked=" `gbm`k'' " e(N)=" `gbN`k'' " esample_rows=" `gbes`k'' " (want `gbwant`k'') esample_and_missing_dose=" `gbbad`k'' " (want 0)"
  if `gbm`k'' <= 0 {
    di as error "F12b1 FAIL: the blanking fixture blanked `gbm`k'' rows - with nothing blanked the caller's guarantee is not exercised"
    local ++fails
  }
  if `gbr`k'' != 0 | `gbf`k'' != `GBF`k'' {
    di as error "F12b1 FAIL: `GBN`k'' arm rc `gbr`k'' fcol `gbf`k'' (want 0 / `GBF`k'') - a blanked dose must not change whether the run completes or collapses"
    local ++fails
  }
  if `gbN`k'' != `gbwant`k'' | `gbes`k'' != `gbwant`k'' | `gbbad`k'' != 0 {
    di as error "F12b1 FAIL: `GBN`k'' arm e(N) `gbN`k'' / e(sample) rows `gbes`k'' (want `gbwant`k'') with `gbbad`k'' missing-dose rows still in e(sample) - the caller must zero the if-marker on every one of them, or a missing dose reaches wooldid_fcol"
    local ++fails
  }
}

* ---- F12b2: THE WORKER, CALLED DIRECTLY. The only way to exercise a guard the caller makes
* unreachable. Thirteen rows, three cells of the key, one out-of-sample row (marker 0, weight
* 1000, dose 9999) that must never enter anything:
*   cell 1  weights 1 2 3 4, doses 10 20 30 40      -> weighted mean 30
*   cell 2  weights 1 2 3 4, doses 10 20 . .        -> weighted mean over the rows that HAVE
*                                                      a value, (10 + 40)/3 = 50/3 = 16.667
*   cell 3  weights 1 2 3 4, doses . . . .          -> stays MISSING
* The shared weight sum the guard replaced gives 50/10 = 5 for cell 2 and rewrites cell 3 to
* 0/10 = 0 (measured on the pre-guard form's arithmetic, 2026-09-02). The doses are read into
* scalars rather than macros so the comparison is at full precision, as F11a2 does for the
* joint p-values. Cell weight sums are 10 each, which also pins that the out-of-sample row's
* weight of 1000 never entered.
clear
set obs 13
gen byte __es = (_n > 1)
gen long __k = cond(_n<=5, 1, cond(_n<=9, 2, 3))
gen double __w = cond(_n==1, 1000, mod(_n-2, 4) + 1)
gen double __d = .
qui replace __d = 9999 in 1
qui replace __d = 10*(_n-1) if inrange(_n,2,5)
qui replace __d = 10 in 6
qui replace __d = 20 in 7
tempfile F12FIX
qui save "`F12FIX'"

* (i) the weighted path
use "`F12FIX'", clear
cap wooldid_fcol __es , key(__k) wvars(__w) contreat(__d)
local ubr = _rc
local ubn = r(ncell)
local ubN = _N
scalar f12wd1 = __d[1]
scalar f12wd2 = __d[2]
scalar f12wd3 = __d[3]
scalar f12ww1 = __w[1]
scalar f12ww2 = __w[2]
scalar f12ww3 = __w[3]
local ubd2 = reldif(f12wd2, 50/3)
di "F12b2_weighted rc=" `ubr' " ncell=" `ubn' " N=" `ubN' " (want 0/3/3)"
di "F12b2_weighted dose nomissing=" %20.15g f12wd1 " (want 30)  partmissing=" %20.15g f12wd2 " (want 50/3, reldif " `ubd2' "; a shared weight sum gives 5)  allmissing=" f12wd3 " (want missing; a shared weight sum gives 0)"
di "F12b2_weighted cell weight sums=" f12ww1 "/" f12ww2 "/" f12ww3 " (want 10 each)"
if `ubr' != 0 | `ubn' != 3 | `ubN' != 3 {
  di as error "F12b2 FAIL: the weighted wooldid_fcol call returned rc `ubr' ncell `ubn' N `ubN' - it must leave exactly one row per cell"
  local ++fails
}
if f12wd1 != 30 {
  di as error "F12b2 FAIL: the no-missing cell's dose came back " f12wd1 ", want the weighted mean 30"
  local ++fails
}
if `ubd2' > 1e-15 {
  di as error "F12b2 FAIL: the part-missing cell's dose came back " f12wd2 " (reldif `ubd2' from 50/3) - it must be the weighted mean over the rows that HAVE a value; a weight sum taken over the whole cell gives 50/10 = 5"
  local ++fails
}
if f12wd3 < . {
  di as error "F12b2 FAIL: the all-missing cell's dose came back " f12wd3 ", want missing - a weight sum taken over the whole cell divides a zero numerator by a positive total and silently rewrites the dose to 0"
  local ++fails
}
if f12ww1 != 10 | f12ww2 != 10 | f12ww3 != 10 {
  di as error "F12b2 FAIL: cell weight sums " f12ww1 "/" f12ww2 "/" f12ww3 ", want 10 each (cell 1's 10 is also what pins that the out-of-sample row's weight of 1000 stayed out)"
  local ++fails
}

* (ii) the no-weight-variable path: with no weight variable each non-missing row counts 1, so
* the same three cells must give the PLAIN mean over the rows that have a value.
use "`F12FIX'", clear
cap wooldid_fcol __es , key(__k) wvars() contreat(__d)
local unr = _rc
local unn = r(ncell)
local unN = _N
scalar f12ud1 = __d[1]
scalar f12ud2 = __d[2]
scalar f12ud3 = __d[3]
di "F12b2_noweight rc=" `unr' " ncell=" `unn' " N=" `unN' " (want 0/3/3) dose nomissing=" f12ud1 " (want 25) partmissing=" f12ud2 " (want 15) allmissing=" f12ud3 " (want missing)"
if `unr' != 0 | `unn' != 3 | `unN' != 3 {
  di as error "F12b2 FAIL: the wvars()-empty wooldid_fcol call returned rc `unr' ncell `unn' N `unN' - it must leave exactly one row per cell"
  local ++fails
}
if f12ud1 != 25 | f12ud2 != 15 | f12ud3 < . {
  di as error "F12b2 FAIL: the wvars()-empty path gave doses " f12ud1 "/" f12ud2 "/" f12ud3 ", want the plain mean over the non-missing rows: 25 / 15 / missing"
  local ++fails
}

* ---- F13 (F061, 2026-09-06): e(wdid_fcncell) is the realised cell count - positive and
* below the row count on aggregating data, MISSING when the collapse did not run (the
* off-switch, or a silent skip) - and e(wdid_fcskip) names the skip reason(s)
use "`ADMIN'", clear
local n13 = _N
qui wooldid y cohid year ttre, att robust
local nc13 = e(wdid_fcncell)
local fs13 "`e(wdid_fcskip)'"
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att robust nofastcollapse
local nc13s = e(wdid_fcncell)
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att customvce(cluster state)
local nc13i = e(wdid_fcncell)
local fs13i "`e(wdid_fcskip)'"
di "F13 fcncell default=`nc13' slow=`nc13s' customvce=`nc13i' (rows `n13') fcskip default=[`fs13'] customvce=[`fs13i']"
if `nc13' >= . | `nc13' < 1 | `nc13' >= `n13' | `nc13s' < . | `nc13i' < . | "`fs13'" != "" | "`fs13i'" != "customvce()" {
  di as error "F13 FAIL: e(wdid_fcncell)/e(wdid_fcskip) do not record the collapse outcome"
  local ++fails
}

* ---- F14 (F037, 2026-09-06): the verbose collapse-decision line names nofastcollapse as
* the reason, alone and beside an ineligible option (which used to be blamed instead)
use "`ADMIN'", clear
qui log using fcol_f14.txt, text replace name(fc)
cap noi wooldid y cohid year ttre, att robust nofastcollapse verbose
local rc14a = _rc
qui log close fc
local lg = fileread("fcol_f14.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local ok14a = strpos(`"`lg'"', "on the full data (nofastcollapse).") > 0
use "`ADMIN'", clear
qui log using fcol_f14.txt, text replace name(fc)
cap noi wooldid y cohid year ttre, att customvce(cluster state) nofastcollapse verbose
local rc14b = _rc
qui log close fc
local lg = fileread("fcol_f14.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local ok14b = strpos(`"`lg'"', "(nofastcollapse; the default cell collapse would not apply with customvce() either)") > 0
local ok14c = strpos(`"`lg'"', "does not apply with customvce() (see nofastcollapse") == 0
cap erase fcol_f14.txt
di "F14 rcs=`rc14a'/`rc14b' alone=`ok14a' with_impvalues=`ok14b' old_text_absent=`ok14c'"
if `rc14a' != 0 | `rc14b' != 0 | !`ok14a' | !`ok14b' | !`ok14c' {
  di as error "F14 FAIL: the verbose collapse line does not name nofastcollapse as the reason"
  local ++fails
}

* ---- F15 (F062, 2026-09-06): impvalues collapses. With the cluster variable in the cell
* key the per-cluster reads are exact: (a) a cluster that differs from i and splits cells -
* i coded cohort-by-region on treated rows with every never-treated unit pooled under one
* code, and cluster(region). The treated levels of i nest in the cluster, so the IM screen
* passes, while the never-treated pool spans the four regions, so the key change splits its
* cells (240 cells against 195 (i, year) pairs) - and (b) a cluster NESTING i (cohid itself,
* 60 cells) both agree with nofastcollapse on b, V and the main table at 1e-12; the
* collapse engages (e(wdid_fcol) = 1) and the run says nothing.
*
* (c) A cluster CROSSED with i - state x cohid on this panel, 3 treated cohort levels across
* 48 states, whose 48 per-cluster reads are re-weightings of 3 cohorts' coefficients - is
* refused r(459) by the IM screen under every safety() setting, collapsed or not, so the
* crossed geometry is pinned as a refusal with safety(off) given.
use "`ADMIN'", clear
gen cohreg = cond(ttre == ., 1, (ttre - 2000)*10 + region)
qui wooldid y cohreg year ttre, att cluster(region) impvalues
matrix bF15a = e(b)
matrix VF15a = e(V)
matrix MF15a = e(wdidmainresults)
local fc15a = e(wdid_fcol)
local nc15a = e(wdid_fcncell)
qui gegen long __f15it = group(cohreg year)
qui sum __f15it
local np15a = r(max)
use "`ADMIN'", clear
gen cohreg = cond(ttre == ., 1, (ttre - 2000)*10 + region)
qui wooldid y cohreg year ttre, att cluster(region) impvalues nofastcollapse
local d15a = max(mreldif(bF15a, e(b)), mreldif(VF15a, e(V)), mreldif(MF15a, e(wdidmainresults)))
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att cluster(cohid) impvalues
matrix bF15b = e(b)
matrix VF15b = e(V)
matrix MF15b = e(wdidmainresults)
local fc15b = e(wdid_fcol)
local nc15b = e(wdid_fcncell)
use "`ADMIN'", clear
qui wooldid y cohid year ttre, att cluster(cohid) impvalues nofastcollapse
local d15b = max(mreldif(bF15b, e(b)), mreldif(VF15b, e(V)), mreldif(MF15b, e(wdidmainresults)))
use "`ADMIN'", clear
gen cohreg = cond(ttre == ., 1, (ttre - 2000)*10 + region)
qui log using fcol_x.txt, text replace name(fc)
cap noi wooldid y cohreg year ttre, att cluster(region) impvalues
local rc15c = _rc
qui log close fc
local lg = fileread("fcol_x.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local cc15 = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
cap erase fcol_x.txt
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) impvalues safety(off)
local rx15a = _rc
use "`ADMIN'", clear
cap qui wooldid y cohid year ttre, att cluster(state) impvalues nofastcollapse safety(off)
local rx15b = _rc
di "F15 pooled: fcol=`fc15a' cells=`nc15a' (i,year)=`np15a' mreldif=`d15a' | nesting: fcol=`fc15b' cells=`nc15b' mreldif=`d15b' | rc=`rc15c' caution=`cc15' (want 0) | crossed safety(off) rc=`rx15a'/`rx15b' (want 459)"
if `fc15a' != 1 | `fc15b' != 1 | `d15a' > 1e-12 | `d15b' > 1e-12 | `nc15a' <= `np15a' | `rc15c' != 0 | `cc15' {
  di as error "F15 FAIL: impvalues under the collapse is not exact or did not engage"
  local ++fails
}
if `rx15a' != 459 | `rx15b' != 459 {
  di as error "F15c FAIL: impvalues with a cluster crossing i ran under safety(off) (rc `rx15a'/`rx15b', want 459)"
  local ++fails
}


* ---- F16 (bug sweep BS-27, 2026-09-06): cfxlattice() without makecfxplots is ignored, as
* announced, and does not decide the collapse, so that run collapses and says nothing about
* it. The real bin-readers (makecfxplots with cfxlattice(), a user lattice() variable) put
* their bin in the cell key, so they too collapse without a collapse caution
use "`REPO'/test/data/dta_cont.dta", clear
gen zbin = mod(i, 3)
tempfile f16
save `f16'
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) att esprelength(1) espostlength(1) cfxlattice(zbin)
local f16rc = _rc
local f16col = e(wdid_fcol)
local f16skip `e(wdid_fcskip)'
use `f16', clear
qui log using fcol_x.txt, text replace name(fc)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) att esprelength(1) espostlength(1) cfxlattice(zbin)
local f16rcb = _rc
qui log close fc
local lg = fileread("fcol_x.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local f16cb = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
use `f16', clear
qui log using fcol_x.txt, text replace name(fc)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) att esprelength(1) espostlength(1) makecfxplots cfxlattice(zbin)
local f16rcc = _rc
local f16colc = e(wdid_fcol)
qui log close fc
local lg = fileread("fcol_x.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local f16cc = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
use `f16', clear
qui log using fcol_x.txt, text replace name(fc)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) att esprelength(1) espostlength(1) lattice(zbin)
local f16rcd = _rc
local f16cold = e(wdid_fcol)
qui log close fc
local lg = fileread("fcol_x.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local f16cd = strpos(`"`lg'"', "Caution: the default cell collapse does not apply") > 0
cap erase fcol_x.txt
di "F16 cfxlattice alone: rc=`f16rc' fcol=`f16col' fcskip=[`f16skip'] caution=`f16cb' (want 0) | makecfxplots+cfxlattice rc=`f16rcc' fcol=`f16colc' caution=`f16cc' | lattice rc=`f16rcd' fcol=`f16cold' caution=`f16cd' (want fcol 1, caution 0 each)"
if `f16rc' != 0 | `f16col' != 1 | "`f16skip'" != "" | `f16rcb' != 0 | `f16cb' | `f16rcc' != 0 | `f16colc' != 1 | `f16cc' | `f16rcd' != 0 | `f16cold' != 1 | `f16cd' {
  di as error "F16 FAIL: cfxlattice() alone must not decide the collapse, and the real bin-readers must collapse with their bin in the cell key, without a collapse caution"
  local ++fails
}

* ---- F13 (BS-28): the one-realized-subgroup level margins take the collapse screen ----
* the dominant risk of BS-28 is a SPURIOUS fire on a healthy run, so all three legs assert
* the guard stayed empty and the numbers did not move against nofastcollapse

* F13a SITE A/B inertness, no missing subgroup values
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg1 = 7
qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg1)
local rc13a = _rc
local fcol13a = e(wdid_fcol)
local g13a = "`e(wdid_fcguard)'"
matrix F13AD = e(wdidsgresults_7)
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg1 = 7
qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg1) nofastcollapse
local rc13as = _rc
local g13as = "`e(wdid_fcguard)'"
matrix F13AS = e(wdidsgresults_7)
local d13a = mreldif(F13AD, F13AS)
if `rc13a' != 0 | `rc13as' != 0 | `fcol13a' != 1 | "`g13a'" != "" | "`g13as'" != "" | `d13a' >= 1e-8 {
  di as error "F13a FAIL: rc `rc13a'/`rc13as' fcol `fcol13a' guard [`g13a']/[`g13as'] mreldif `d13a'"
  local ++fails
}

* F13b SITE A/B on the BS-28 exposure shape: one realized subgroup, treated rows with a
* missing subgroup value, so the level margin is a strict subset of the main ATT
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sgm1 = 1
qui replace sgm1 = . if mod(i,5) == 0
qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sgm1)
local rc13b = _rc
local nsg13b = wordcount("`e(wdid_sglist)'")
local g13b = "`e(wdid_fcguard)'"
matrix F13BD = e(wdidsgresults_1)
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sgm1 = 1
qui replace sgm1 = . if mod(i,5) == 0
qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sgm1) nofastcollapse
local rc13bs = _rc
local g13bs = "`e(wdid_fcguard)'"
matrix F13BS = e(wdidsgresults_1)
local d13b = 0
if `rc13b' == 0 & `rc13bs' == 0 {
  forval cc = 1/4 {
    matrix F13BDC = F13BD[1..., `cc']
    matrix F13BSC = F13BS[1..., `cc']
    local dd = mreldif(F13BDC, F13BSC)
    if `dd' > `d13b' local d13b = `dd'
  }
}
if `rc13b' != 0 | `rc13bs' != 0 | `nsg13b' != 1 | "`g13b'" != "" | "`g13bs'" != "" | `d13b' >= 1e-8 {
  di as error "F13b FAIL: rc `rc13b'/`rc13bs' sglist n=`nsg13b' guard [`g13b']/[`g13bs'] max col mreldif `d13b'"
  local ++fails
}

* F13c SITE C: the pre-period fallback after a deferred multi-subgroup pwcompare
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg2 = 1
qui replace sg2 = 2 if mod(i,7) == 0
qui drop if sg2 == 2 & t < ttre - 1 & ttre != .
qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) esfixedbaseperiod esrelativeto(-1) subgroups(sg2)
local rc13c = _rc
local nsg13c = wordcount("`e(wdid_sglist)'")
local g13c = "`e(wdid_fcguard)'"
matrix F13C1 = e(wdidsgresults_1)
matrix F13C2 = e(wdidsgresults_2)
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg2 = 1
qui replace sg2 = 2 if mod(i,7) == 0
qui drop if sg2 == 2 & t < ttre - 1 & ttre != .
qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) esfixedbaseperiod esrelativeto(-1) subgroups(sg2) nofastcollapse
local rc13cs = _rc
local g13cs = "`e(wdid_fcguard)'"
matrix F13C1S = e(wdidsgresults_1)
matrix F13C2S = e(wdidsgresults_2)
local d13c1 = mreldif(F13C1, F13C1S)
local d13c2 = mreldif(F13C2, F13C2S)
di "F13 BS-28: a rc=`rc13a'/`rc13as' fcol=`fcol13a' d=`d13a' | b rc=`rc13b'/`rc13bs' nsg=`nsg13b' d=`d13b' | c rc=`rc13c'/`rc13cs' nsg=`nsg13c' d=`d13c1'/`d13c2' guards [`g13a'][`g13b'][`g13c']"
if `rc13c' != 0 | `rc13cs' != 0 | `nsg13c' != 2 | "`g13c'" != "" | "`g13cs'" != "" | `d13c1' >= 1e-8 | `d13c2' >= 1e-8 {
  di as error "F13c FAIL: rc `rc13c'/`rc13cs' sglist n=`nsg13c' guard [`g13c']/[`g13cs'] mreldif `d13c1'/`d13c2'"
  local ++fails
}

if `fails' == 0 di as result "FASTCOLLAPSE_ALL_PASS"
else di as error "FASTCOLLAPSE_FAILS=`fails'"
