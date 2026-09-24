* test_bootstrap_gates.do - commit-B bootstrap option surface + wooldid_bootctl gates.
* Run from test/: stata-mp -e do test_bootstrap_gates.do
* Commit-B contract: every gate fires informatively; a fully valid request still exits
* r(198) at the loud "engine not yet enabled" tail. Prints BOOTGATES_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

* deterministic staggered panel: 12 units, unit-level rows; s = coarse grouping (3 states)
clear
set obs 12
gen unit = _n
gen s = ceil(unit/4)
gen ttre = cond(unit<=4, 4, cond(unit<=8, 6, .))
expand 8
bysort unit: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
tempfile panel
save `panel'

cap prog drop chk
prog def chk
  args num expected_rc got_rc
  if `got_rc' != `expected_rc' {
    di as error "G`num' FAIL: rc " `got_rc' " != expected `expected_rc'"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
end
global BOOTG_FAILS 0

cap prog drop usegate
prog def usegate
  * a saver's failed write is a caution since 2026-09-21, so a run's rc 0 no longer proves
  * the file exists: name the miss as a FAIL and load an empty stand-in carrying the columns
  * the checks that follow touch, so they fail on content rather than abort the suite r(601)
  args file label
  cap confirm file `file'.dta
  if _rc {
    di as error "`label' FAIL: bootsave() wrote no `file'.dta (e(boot_save_rc)=`e(boot_save_rc)')"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
    clear
    qui set obs 1
    gen _rep = 0
    gen w_0 = .
    gen advW = .
    gen jointW_0 = .
    exit
  }
  use `file', clear
end

* G1 stray companions without bootstrap()
cap qui wooldid y unit t ttre, att cluster(unit) reps(50)
chk 1 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) rseed(42)
chk 1b 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootalpha(equal)
chk 1c 198 _rc

* G2 bad master token
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(wild) reps(50)
chk 2 198 _rc

* G3 reps mandatory
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard)
chk 3 198 _rc

* G4 companion validation
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) level(5)
chk 4a 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) dots(-2)
chk 4b 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) rseed(abc)
chk 4c 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) rseed(2.5)
chk 4d 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(50) bootalpha(nonsense)
chk 4e 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) bootalpha(equal)
chk 4f 198 _rc
* G4g (change audit round 2, U08/U07 dots-negative): dots() is shared with permute(), and its
* range check sat below wooldid_bootctl's permute() exit, so a negative dots() under permute()
* ran silently as dots(0) while bootstrap() refused it (4b). Both engines refuse it now, with
* the same message
qui log using boot_gates_g4g.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(7) dots(-2)
local g4grc = _rc
qui log close bg
chk 4g 198 `g4grc'
if strpos(fileread("boot_gates_g4g.txt"), "dots() must be a nonnegative integer") == 0 {
  di as error "G4g FAIL: a negative dots() under permute() is not refused by the dots() range check"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g4g.txt

* G5 interaction refusals
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) permute(shuffle) permreps(20)
chk 5a 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) impvalues
chk 5b 198 _rc
* 5c ruled 2026-09-21: fastnose alongside bootstrap() cautions and is dropped rather than
* refused (the replications already run fast unless boottstat asks otherwise). The boottstat
* pairing still refuses - test_bootstrap_tstat T8 pins that leg.
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) fastnose
chk 5c 0 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) unconditionalse
chk 5d 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) customvce(robust)
chk 5e 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) oldsyntax
chk 5f 198 _rc
cap qui wooldid y unit t ttre, att att_it cluster(unit) bootstrap(standard) reps(50) suppresssgprimaryeffects subgroups(s)
chk 5g 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit s) bootstrap(standard) reps(50)
chk 5h 198 _rc

* G6 bayes + cell-weighted: refused under default safety and under safety(off) alike;
* safety(off) is not a way past it (6b; message-checked below)
cap qui wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(50)
chk 6a 198 _rc
cap noi wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(30) safety(off) rseed(70101)
chk 6b 198 _rc

* G7 standard + cell-weighted + nesting cluster (cluster == i) -> RUNS (re-id leg)
cap qui wooldid y unit t ttre, att_it cluster(unit) bootstrap(standard) reps(20) rseed(70102)
chk 7 0 _rc
* G7b coarser nesting cluster (s nests unit? NO - s is COARSER than unit: each unit
* wholly inside one s, so cluster(s) nests i=unit) -> gates pass -> tail
cap qui wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) rseed(70103)
chk 7b 0 _rc

* G8 CONVERTED 2026-09-16 (T2) from a refusal test to a positive one, case number kept
* so the history reads. It used to pin r(198), "observation-level resampling cannot
* re-identify duplicated treatment cells". The leg now keeps i and resamples
* observations, and duplicate (i,t) rows are a cell with more rows and more weight
* rather than a contradiction - so the case now pins the branch it used to exclude:
* the run completes, NO replication fails, and the point estimates are bitwise the
* analytic run's, since a bootstrap replaces the inference and must never move the
* estimate. A branch that fragmented units into single-period units would fail
* replications on inestimable horizons long before it moved a coefficient.
* rseed() is load-bearing: without it the draw depends on whatever RNG state the cases
* above happened to leave, so the completion count is not a property of this case at all.
* 50/50 is the MEASURED value at rseed(9090) on this fixture, pinned deliberately.
qui wooldid y unit t ttre, att_it robust
tempname G8B
matrix `G8B' = e(b)
cap qui wooldid y unit t ttre, att_it robust bootstrap(standard) reps(50) rseed(9090)
chk 8 0 _rc
if _rc == 0 {
  if e(boot_completed) != 50 | e(boot_failed) != 0 {
    di as error "G8b FAIL: " e(boot_completed) " of 50 replications completed, " e(boot_failed) " failed (measured 50/0 at rseed(9090))"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
  if mreldif(e(b), `G8B') != 0 {
    di as error "G8c FAIL: bootstrap(standard) moved the point estimates (mreldif " mreldif(e(b), `G8B') ")"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}

* G9 CONVERTED 2026-09-16 (T2) the same way, case number kept. i = s (3 coarse groups),
* cluster = unit (12), so the cluster is STRICTLY FINER than i - unit varies within s.
* It used to pin r(459), "strict subsets of your cluster() variable". i is an
* estimand-definition choice and cluster() a claim about dependence and assignment, so
* the pairing is legitimate and the leg keeps i here rather than splitting a straddling
* unit in two. The fixture's own shape is asserted first: without it a later change to
* the panel could turn this back into a nesting case and the test would still "pass"
* while no longer being about the feature.
qui wooldid y s t ttre, att_it cluster(unit)
tempname G9B
matrix `G9B' = e(b)
local g9nc = e(N_clust)
qui levelsof s, local(g9s)
qui levelsof unit, local(g9u)
if `: word count `g9s'' != 3 | `: word count `g9u'' != 12 | `g9nc' != 12 {
  di as error "G9a FAIL: fixture is not finer-than-i (s levels `: word count `g9s'', unit levels `: word count `g9u'', e(N_clust) `g9nc'; expected 3, 12, 12)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap qui wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(50) rseed(9090)
chk 9 0 _rc
local g9done = e(boot_completed)
local g9fail = e(boot_failed)
if _rc == 0 {
  if `g9done' != 50 | `g9fail' != 0 {
    di as error "G9b FAIL: " `g9done' " of 50 replications completed, " `g9fail' " failed (measured 50/0 at rseed(9090))"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
  if mreldif(e(b), `G9B') != 0 {
    di as error "G9c FAIL: bootstrap(standard) moved the point estimates (mreldif " mreldif(e(b), `G9B') ")"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}
* G9d the PAIRED CONTROL, and the reason G9b can be pinned at all. Measured 2026-09-16:
* on this 12-cluster fixture a replication can lose a cohort entirely and fail, and that
* happens at the SAME seeds and the SAME rate under keep-i (i = s, cluster = unit) as
* under the group() branch (i = unit, cluster = unit) - across seeds 1-5, both arms gave
* 50, 49, 50, 50, 50 completions, while dropping to cluster(s) (3 clusters) cost 15-20 of
* 50 in BOTH regimes. So a lost replication here is a property of the cluster count, not
* of the branch, and the discriminating assertion is that the two regimes agree at one
* seed rather than that either reaches 50. Pinning only G9b would make a fixture change
* look like a branch defect; pinning only G9d would miss a branch that failed outright.
cap qui wooldid y unit t ttre, att_it cluster(unit) bootstrap(standard) reps(50) rseed(9090)
if _rc == 0 & e(boot_completed) != `g9done' {
  di as error "G9d FAIL: keep-i completed `g9done' of 50 where the group() branch on the same fixture and seed completed " e(boot_completed) "; the two regimes should lose replications at the same rate here"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G10 standard + att only + cluster(unit): runs (count path)
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(70104)
chk 10 0 _rc

* G11 single cluster -> r459
cap qui wooldid y unit t ttre if s==1, att cluster(s) bootstrap(standard) reps(50)
chk 11 459 _rc

* G12 update outranks (pre-existing semantics unchanged)
cap qui wooldid, update bootstrap(standard)
local g12rc = _rc
if `g12rc' != 0 & `g12rc' != 198 {
  * update may legitimately fail rc!=0 in offline test envs; accept 0 or the
  * documented offline-update rc but never a bootctl-shaped refusal message
}

* G14 jointtests + bootstrap: RUNS (boot Wald; exempt from the SI-A1 screen)
cap qui wooldid y unit t ttre, att esprelength(2) espostlength(1) cluster(unit) bootstrap(bayes) reps(20) rseed(2) jointtests
chk 14 0 _rc
* G15 pretrendpower + bootstrap: RUNS as of 2026-09-04. The plan-A13 refusal was lifted
* on the author's ruling: the Roth diagnostics now run from the shell AFTER the bootstrap
* rewiring, so Roth's own formulas consume the BOOTSTRAP event-study covariance. The
* refusal had rested on an ordering constraint, not an impossibility.
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) esfixedbaseperiod bootstrap(bayes) reps(20) pretrendpower(80) rseed(70105)
chk 15 0 _rc
* G15b it must actually CONSUME the bootstrap covariance, not merely run: the power
* figures have to move against the analytic fit on the same data. A run that silently
* kept reading the analytic V_es would pass G15 and be wrong, so pin the movement.
cap matrix G15b_boot = e(roth_power_att)
local g15rc = _rc
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) esfixedbaseperiod pretrendpower(80)
cap matrix G15b_anl = e(roth_power_att)
local g15rc2 = _rc
local g15moved 0
local g15cells 0
if `g15rc' == 0 & `g15rc2' == 0 {
  forvalues r = 1/`=rowsof(G15b_anl)' {
    forvalues c = 1/`=colsof(G15b_anl)' {
      if G15b_anl[`r',`c'] < . & G15b_boot[`r',`c'] < . {
        local ++g15cells
        if reldif(G15b_anl[`r',`c'], G15b_boot[`r',`c']) > 1e-10 local ++g15moved
      }
    }
  }
}
di "G15b_cells=`g15cells' moved=`g15moved' (want cells>0 and moved>0)"
if `g15cells' == 0 | `g15moved' == 0 {
  di as error "G15b FAIL: pretrendpower under bootstrap did not consume the bootstrap covariance"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* G15c hbar, the weighted mean post-treatment relative period, is a property of the data and not
* of the covariance: the worker computes it once and the shell hands it back to the power algebra
* after the rewiring, so a bootstrap run must report the analytic run's hbar. Not read off the
* G15/G15b pair, whose hbar is 1.625: that is exactly representable, so it survives the hand-back
* whatever the precision. On dta_binary hbar is 27/11. The hand-back went through string(), which
* keeps 7 significant digits, so the two differed at about 1.3e-7 relative (change audit
* 2026-09-23, U09 ptp-boot-hbar-string-precision); the author ruled full precision, the shell now
* writes it in %21x, and this leg, which allowed 1e-6 until then, asserts bitwise equality. It
* also fails a bootstrap run that loses hbar or takes it from the draws
use "data/dta_binary", clear
cap qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) pretrendpower(80)
local g15crc1 = _rc
scalar G15c_anl = e(roth_hbar_att)
use "data/dta_binary", clear
cap qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) pretrendpower(80) bootstrap(standard) reps(10) rseed(1) dots(0)
local g15crc2 = _rc
scalar G15c_boot = e(roth_hbar_att)
di "G15c rc=`g15crc1'/`g15crc2' hbar analytic=" %21x G15c_anl " bootstrap=" %21x G15c_boot " reldif=" reldif(G15c_anl, G15c_boot)
if `g15crc1' != 0 | `g15crc2' != 0 | G15c_anl >= . | G15c_boot >= . | G15c_anl != G15c_boot {
  di as error "G15c FAIL: the bootstrap run's hbar (" G15c_boot ") is not the analytic run's (" G15c_anl "), rc=`g15crc1'/`g15crc2'"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
* G13 plain run untouched
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1)
chk 13 0 _rc

* ---- message grading: rcs alone can mask wrong-path failures (the oldsyntax
* normalization bug rode under matching 198s). capture output and assert the
* EXPECTED message printed for the tail-reaching and caution shapes.
cap prog drop logh
prog def logh
  * logh <marker> <needle...>: assert the log section written since `marker' log
  * was opened contains the needle string; closes the log
  * The {com} skip is load-bearing, not cosmetic (2026-09-17). Stata echoes every
  * command into EVERY open log, so the `logh mN <needle>' line itself lands in the
  * log a moment before logh closes and reads it - and the echo carries the needle.
  * Without this filter the needle always matched its own echo and every assertion
  * below passed vacuously; the sibling lognh always failed, which is how it
  * surfaced. SMCL marks echoed commands `{com}. ' and continuations `{com}> ';
  * genuine output lines are unprefixed or carry {txt}/{res}/{err}
  gettoken lg 0 : 0
  qui log close `lg'
  * gettoken leaves the separating blank on the remainder; an untrimmed needle
  * would never match a message that starts at column 0
  local needle = strtrim(`"`0'"')
  tempname fh
  local found 0
  file open `fh' using boot_gates_`lg'.smcl, read
  file read `fh' line
  while r(eof) == 0 {
    if substr(`"`macval(line)'"', 1, 5) != "{com}" {
      if strpos(`"`macval(line)'"', `"`needle'"') local found 1
    }
    file read `fh' line
  }
  file close `fh'
  if `found' == 0 {
    di as error "MSG-`lg' FAIL: expected message not found: `needle'"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
end

* lognh is logh's other half. A message assertion that only ever checks for
* PRESENCE cannot tell a correctly suppressed message from a message that was
* never written at all, so a gate with two branches needs both directions
cap prog drop lognh
prog def lognh
  * lognh <marker> <needle...>: assert the log section written since `marker' log
  * was opened does NOT contain the needle string; closes the log
  * skips echoed command lines for the reason given in logh above
  gettoken lg 0 : 0
  qui log close `lg'
  * gettoken leaves the separating blank on the remainder; an untrimmed needle
  * would never match a message that starts at column 0
  local needle = strtrim(`"`0'"')
  tempname fh
  local found 0
  file open `fh' using boot_gates_`lg'.smcl, read
  file read `fh' line
  while r(eof) == 0 {
    if substr(`"`macval(line)'"', 1, 5) != "{com}" {
      if strpos(`"`macval(line)'"', `"`needle'"') local found 1
    }
    file read `fh' line
  }
  file close `fh'
  if `found' != 0 {
    di as error "MSG-`lg' FAIL: message present that should be suppressed: `needle'"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
end

use `panel', clear
qui log using boot_gates_m1.smcl, replace name(m1) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(8)
logh m1 Bootstrap Inference
qui log using boot_gates_m2.smcl, replace name(m2) smcl
cap noi wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(30) safety(off) rseed(70106)
logh m2 Error: bootstrap(bayes) is not compatible with non-att estimands
qui log using boot_gates_m3.smcl, replace name(m3) smcl
* m3 used to pin the nesting refusal's message; T2 (2026-09-16) removed the refusal,
* so what this leg now pins is that the same request RUNS and prints the bootstrap
* footer rather than a message about nesting
cap noi wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(50) rseed(70107)
logh m3 Bootstrap Inference
qui log using boot_gates_m4.smcl, replace name(m4) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) oldsyntax rseed(70108)
logh m4 not available under oldsyntax
qui log using boot_gates_m5.smcl, replace name(m5) smcl
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) rseed(8)
logh m5 Bootstrap Inference
qui log using boot_gates_m6.smcl, replace name(m6) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(20) rseed(3)
logh m6 Bootstrap Inference

* ---- A10 (open_defects.md): the re-id leg's covariate caution, pinned on BOTH
* branches of its gate. The caution says that when clusters are duplicated, a
* covariate keyed to unit identity keeps its original codes across the copies.
* T2/T3 (2026-09-16) gated it on `wdb_nest' because under a cluster finer than or
* cross-cutting i the leg keeps i whole and re-identifies nothing, so the caution
* would be false there. Neither branch had ever run: no config wrote fe(i).
* Three legs, because the gate has two independent conjuncts:
*   m7  nesting regime + fe(i)     -> MUST fire      (the claim is true)
*   m8  finer regime  + fe(i)      -> MUST NOT fire  (the T2/T3 suppression)
*   m9  nesting regime, no fe(i)   -> MUST NOT fire  (the wdb_ikeyed conjunct)
* m8 and m9 together rule out the two ways a vacuous gate could pass m7.
* s nests unit (s = ceil(unit/4)), so cluster(s) with i=unit is the nesting
* regime and cluster(unit) with i=s is the finer one. att_it + bootstrap(standard)
* is what sets `wdb_reid' in the first place (wooldid.ado:378)
qui log using boot_gates_m7.smcl, replace name(m7) smcl
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) fe(unit) rseed(70109)
chk 7m 0 _rc
* m7/m8 pinned the re-id covariate caution, DELETED 2026-09-21 on the author's ruling: the
* right answer is to substitute the re-id'd i into the covariate payload rather than warn, and
* the scan the caution used could not see i.unit anyway. See
* review_process/open_ticket_bootreid_covariates_2026-09-21.md.
*
* The two message assertions are removed rather than inverted. lognh cannot make a negative
* assertion about a needle that appears in its OWN echoed command line: Stata writes that echo
* into the open log, and it is not always prefixed {com}, which is the same vacuity trap the
* logh comment above documents. The two RUNS and their rc assertions stay, so the re-id leg is
* still exercised; the logs are closed by hand since logh/lognh are what used to close them.
qui log close m7
qui log using boot_gates_m8.smcl, replace name(m8) smcl
cap noi wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(20) fe(s) rseed(70110)
chk 8m 0 _rc
qui log close m8

* ---- m10-m13 (open_ticket_bootreid_covariates, built 2026-09-22): a covariate that names
* i reaches the re-id leg's replications with the fresh unit code in its place.
*   m10 nesting regime + boottstat + a per-unit slope keyed to i: the same seed with and
*       without the rewrite (the WOOLDID_BOOTNOCOVSUB negative control). The point
*       estimates and the draw distribution are the same to the last bit, since drawn
*       copies are exact duplicates; the studentized inference moves, because each
*       replication's cluster-robust covariance treats the slope per copy rather than
*       pooled across copies.
*       e(boot_covsub) counts the rewritten terms and is absent under the control.
*   m11 names that merely contain the i name are left alone: a rewrite would make every
*       replication fail to find the variable.
*   m12 no rewrite outside the nesting regime (cluster finer than i keeps i whole).
*   m13 the matcher itself, including the unabbreviation a string-typed option needs.
* the covariate is a per-unit slope entered as a regressor, controls(c.x#i.unit): measured
* 2026-09-22, that form moves the studentized results on this panel and on a 24-unit one,
* whereas an absorbed unit trend fe(unit#c.t) moved nothing here and 0.07 on the other
use `panel', clear
set seed 71
gen x = rnormal()
cap noi wooldid y unit t ttre, att cluster(s) bootstrap(standard) boottstat reps(20) controls(c.x#i.unit) rseed(70112) dots(0)
chk 10m 0 _rc
local m10sub = e(boot_covsub)
matrix M10B = e(b)
matrix M10R = e(wdidmainresults)
matrix M10D = e(bootdist)
global WOOLDID_BOOTNOCOVSUB 1
use `panel', clear
set seed 71
gen x = rnormal()
cap noi wooldid y unit t ttre, att cluster(s) bootstrap(standard) boottstat reps(20) controls(c.x#i.unit) rseed(70112) dots(0)
local m10rc = _rc
global WOOLDID_BOOTNOCOVSUB
chk 10n 0 `m10rc'
local m10nsub = e(boot_covsub)
di "m10 covsub live=`m10sub' control=`m10nsub' b_reldif=" mreldif(M10B, e(b)) " draws_reldif=" mreldif(M10D, e(bootdist)) " results_reldif=" mreldif(M10R, e(wdidmainresults))
if "`m10sub'" != "1" | "`m10nsub'" != "." {
  di as error "m10 FAIL: e(boot_covsub) is `m10sub' with the rewrite and `m10nsub' under the control (want 1 and absent)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
if mreldif(M10B, e(b)) != 0 | mreldif(M10D, e(bootdist)) > 1e-10 {
  di as error "m10 FAIL: the rewrite moved the point estimates or the draws, which the exact-duplicate argument says it cannot"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
if mreldif(M10R, e(wdidmainresults)) < 1e-6 {
  di as error "m10 FAIL: the studentized results are the same with and without the rewrite, so the rewrite is not reaching the replications"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* on this 3-cluster design a draw with no treated or no control cluster is discarded, so
* the completed count is compared with the same seed run without the controls, not with reps
use `panel', clear
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) rseed(70113) dots(0)
chk 11l 0 _rc
local m11ref = e(boot_completed)
use `panel', clear
set seed 72
gen unitx = rnormal()
gen myunit = rnormal()
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) controls(unitx myunit) rseed(70113) dots(0)
chk 11m 0 _rc
if "`e(boot_covsub)'" != "" | e(boot_completed) != `m11ref' | `m11ref' == 0 {
  di as error "m11 FAIL: e(boot_covsub)=`e(boot_covsub)' completed=" e(boot_completed) " vs `m11ref' without the controls (want absent and equal): a name containing the i name was rewritten"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
cap noi wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(20) fe(s#c.t) rseed(70114) dots(0)
chk 12m 0 _rc
if "`e(boot_covsub)'" != "" | e(boot_reid) != 0 {
  di as error "m12 FAIL: e(boot_covsub)=`e(boot_covsub)' reid=" e(boot_reid) " (want absent and 0): a rewrite outside the nesting regime"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
gen unitx = 1
gen x = 1
wooldid_isub , ivar(unit) new(__wdbnewi) payload(i.unit#c.t unitx myunit unit c.x#i.unit i.unit##c.t)
if r(n) != 4 | `"`r(payload)'"' != "i.__wdbnewi#c.t unitx myunit __wdbnewi c.x#i.__wdbnewi i.__wdbnewi##c.t" {
  di as error `"m13 FAIL: matcher gave n=`r(n)' payload=`r(payload)'"'
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
preserve
clear
qui set obs 3
gen cohort = _n
gen x = _n
gen t = _n
wooldid_isub , ivar(cohort) new(__wdbnewi) payload(i.coh#c.x c.x i.cohort##c.t) unab
if r(n) != 2 | `"`r(payload)'"' != "i.__wdbnewi#c.x c.x i.__wdbnewi##c.t" {
  di as error `"m13b FAIL: unabbreviating matcher gave n=`r(n)' payload=`r(payload)'"'
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
restore
qui log using boot_gates_m9.smcl, replace name(m9) smcl
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) rseed(70111)
chk 9m 0 _rc
lognh m9 Caution: a covariate specification names the i variable

* G16 bootsave(): stray without bootstrap() -> 198; writes <name>.dta to the cwd with the draws
* bitwise equal to e(bootdist) and an observed _rep = 0 row equal to e(b); an existing file is
* refused without replace (602) and overwritten with it
cap erase "boot_gates_g16.dta"
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) bootsave(boot_gates_g16)
chk 16 198 _rc
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(15) rseed(7) bootsave(boot_gates_g16)
chk 16b 0 _rc
matrix G16_BD = e(bootdist)
matrix G16_OB = e(b)
local g16bn : colnames G16_OB
cap confirm file "boot_gates_g16.dta"
chk 16c 0 _rc
if _rc == 0 {
  preserve
  use "boot_gates_g16.dta", clear
  mata:
    BD = st_matrix("G16_BD"); OB = st_matrix("G16_OB"); nm = tokens(st_local("g16bn"))
    X = J(st_nobs(), 0, .)
    for (k=1;k<=cols(BD);k++) X = (X, st_data(., strtoname(nm[k])))
    rep = st_data(., "_rep")
    obs = X[selectindex(rep:==0), .]
    draws = select(X[selectindex(rep:>0), .], rowmissing(X[selectindex(rep:>0), .]):==0)
    // Mata's & does not short-circuit, so guard the conformability before comparing
    ok = (rows(obs)==1 & cols(obs)==cols(OB) & rows(draws)==rows(BD) & cols(draws)==cols(BD) & rows(BD) > 0)
    if (ok) ok = (max(abs(obs :- OB))==0 & max(abs(draws :- BD))==0)
    st_numscalar("g16ok", ok)
  end
  restore
  if g16ok != 1 {
    di as error "G16 FAIL: the saved draws or the observed row do not match e(bootdist) / e(b)"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(15) rseed(7) bootsave(boot_gates_g16)
chk 16d 602 _rc
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(15) rseed(7) bootsave(boot_gates_g16) replace
chk 16e 0 _rc
cap erase "boot_gates_g16.dta"


* G17 bootsave() name rules are checked BEFORE the replication loop (review 2026-09-05):
* a clash or a bad name is refused with nothing posted (e(boot_reps) absent), where it
* used to surface only after every replication had run
use `panel', clear
qui save "boot_gates_g17.dta", replace
ereturn clear
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(15) rseed(7) bootsave(boot_gates_g17)
chk 17 602 _rc
cap confirm scalar e(boot_reps)
if _rc == 0 {
  di as error "G17b FAIL: e(boot_reps) is posted, so the replication loop ran before the bootsave() clash was refused"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
ereturn clear
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(15) rseed(7) bootsave(sub/boot_gates_g17)
chk 17c 198 _rc
cap confirm scalar e(boot_reps)
if _rc == 0 {
  di as error "G17d FAIL: e(boot_reps) is posted after a bad bootsave() name"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase "boot_gates_g17.dta"


* G18 bootsave() with a long subgroup variable name: the pwse<j>_ column names are
* truncated to Stata's 32-character limit (the index keeps them unique); r(3300) before
* 2026-09-05. boottstat + cluster(unit) nests i, so the pairwise-se accumulator exists
use `panel', clear
gen sg_a_very_long_subgroup_name_xx = s
cap erase "boot_gates_g18.dta"
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) boottstat subgroups(sg_a_very_long_subgroup_name_xx) bootsave(boot_gates_g18)
chk 18 0 _rc
cap confirm file "boot_gates_g18.dta"
chk 18b 0 _rc
if _rc == 0 {
  preserve
  use "boot_gates_g18.dta", clear
  cap ds pwse1_*
  if _rc != 0 {
    di as error "G18c FAIL: no pwse1_ column in the saved draws"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
  restore
}
cap erase "boot_gates_g18.dta"


* G19 verbose + bootstrap(): the observed run's retained ___wooldidfullmodel survives the
* replications (each replication stored ITS model under that name until 2026-09-05, so
* est restore handed back the last replication's reweighted fit); replications store
* under ___wooldidrepmodel and drop it at exit
use `panel', clear
cap est drop _all
qui wooldid y unit t ttre, att cluster(unit) verbose
qui est restore ___wooldidfullmodel
matrix G19_B = e(b)
cap est drop _all
cap qui wooldid y unit t ttre, att cluster(unit) verbose bootstrap(standard) reps(10) rseed(7)
chk 19 0 _rc
cap qui est restore ___wooldidfullmodel
chk 19b 0 _rc
if _rc == 0 {
  if mreldif(e(b), G19_B) != 0 {
    di as error "G19c FAIL: ___wooldidfullmodel after a verbose bootstrap run is not the observed run's model"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}
cap qui est restore ___wooldidrepmodel
if _rc == 0 {
  di as error "G19d FAIL: a replication's ___wooldidrepmodel leaked"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap est drop _all


* G20 a failure after the replication loop is an error of the call (review 2026-09-05):
* the event-study plot save fails, the rc propagates (it was swallowed: rc 0 with a stale
* saved file), the caller's frame is restored, and the replay branch refuses until a
* wooldid call completes. The failure is forced with a read-only .gph under replace
* (r(608)): a file that merely exists without replace is refused before the loop since
* 2026-09-23 (G41), so it no longer reaches the plot stage
use `panel', clear
foreach f in boot_gates_g20_ES_att.gph boot_gates_g20_ES.dta {
  cap erase "`f'"
}
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) makeplots saveplots(boot_gates_g20)
chk 20 0 _rc
! chmod 444 boot_gates_g20_ES_att.gph
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) makeplots saveplots(boot_gates_g20) replace
chk 20b 608 _rc
qui frame
if "`r(currentframe)'" != "default" {
  di as error "G20c FAIL: current frame is `r(currentframe)' after the plot failure"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  frame change default
}
cap qui wooldid
chk 20d 198 _rc
! chmod 644 boot_gates_g20_ES_att.gph
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) makeplots saveplots(boot_gates_g20) replace
chk 20e 0 _rc
cap qui wooldid
chk 20f 0 _rc
foreach f in boot_gates_g20_ES_att.gph boot_gates_g20_ES.dta {
  cap erase "`f'"
}


* G21 dots(#) prints one mark every # replications, as documented (every replication
* printed a mark whatever # was, until 2026-09-05)
use `panel', clear
qui log using boot_gates_g21.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(7) dots(5)
local g21rc = _rc
qui log close bg
chk 21 0 `g21rc'
local lg = fileread("boot_gates_g21.txt")
local p0 = strpos(`"`lg'"', "Bootstrap replications (20): ")
if `p0' == 0 {
  di as error "G21b FAIL: no progress line in the log"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
else {
  local seg = substr(`"`lg'"', `p0' + 29, 300)
  local seg = substr(`"`seg'"', 1, strpos(`"`seg'"', "done") - 1)
  local nd = length(`"`seg'"') - length(subinstr(`"`seg'"', ".", "", .))
  if `nd' != 4 {
    di as error "G21c FAIL: dots(5) with reps(20) printed `nd' marks (expected 4)"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}
cap erase boot_gates_g21.txt


* G22 the boottstat pairwise-contrast caution names clustering at or above the i level;
* on an UNCLUSTERED run it must not print (nclus was the observation count there)
use `panel', clear
qui log using boot_gates_g22.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(7) boottstat subgroups(s)
local g22rc = _rc
qui log close bg
chk 22 0 `g22rc'
if `g22rc' == 0 {
  local lg = fileread("boot_gates_g22.txt")
  if strpos(`"`lg'"', "under clustering at or above the i level") > 0 {
    di as error "G22b FAIL: the clustered pairwise-se caution printed on an unclustered run"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}
cap erase boot_gates_g22.txt


* G23 gteffects under the physical re-id leg (cell-weighted estimand + nesting cluster):
* the replication's cells are keyed on its fresh unit codes, so they must be translated
* back to the original units before harvesting (review 2026-09-05: with codes 1..N the
* wrong cells were harvested, with any other coding nothing was). Pin: relabelling the
* units by +100 must leave the bootstrapped cell columns bitwise unchanged, and they
* must be populated
use `panel', clear
cap qui wooldid y unit t ttre, att_it cluster(unit) esprelength(1) espostlength(1) gteffects bootstrap(standard) reps(12) rseed(7)
chk 23 0 _rc
if _rc == 0 {
  matrix G23_A = e(gteffects)
  use `panel', clear
  gen unit2 = unit + 100
  cap qui wooldid y unit2 t ttre, att_it cluster(unit2) esprelength(1) espostlength(1) gteffects bootstrap(standard) reps(12) rseed(7)
  chk 23b 0 _rc
  if _rc == 0 {
    matrix G23_B = e(gteffects)
    * matrix == is missing-aware (a missing cell equals only a missing cell), where a
    * max(abs(.)) comparison would skip the blanks that the pristine build produces
    mata: A = st_matrix("G23_A"); B = st_matrix("G23_B"); ok = (rows(A) == rows(B) & cols(A) == cols(B)); if (ok) ok = (A[., 2..cols(A)] == B[., 2..cols(B)]) & (B[., 1] == A[., 1] :+ 100); nse = min((sum(A[., 6] :< .), sum(B[., 6] :< .))); st_numscalar("g23ok", ok); st_numscalar("g23nse", nse)
    if g23ok != 1 | g23nse == 0 {
      di as error "G23c FAIL: gteffects bootstrap columns are not relabelling-invariant under the re-id leg, or no cell carries a bootstrap se (ok=" g23ok ", cells with se=" g23nse ")"
      global BOOTG_FAILS = ${BOOTG_FAILS} + 1
    }
  }
}


* G24 the bias-corrected endpoints follow bstat's order-statistic rule (integer n*p averages
* the adjacent order statistics; the same p round-trip through normal(2*z0 + invnormal(p))
* that bstat performs, floating-point noise included): with z0 = 0 the BC interval equals
* bstat's e(ci_bc) exactly. Synthetic n = 40 draws 1..40, comparator 20.5 (z0 = 0), level
* 95: bstat gives BC [2, 39.5] and percentile [1.5, 39.5]; the pristine v[ceil(n*p)] rule
* gave the BC upper endpoint as 39
preserve
clear
qui set obs 40
gen double b = _n
qui bstat b, stat(20.5) n(40) level(95)
matrix G24_BS_BC = e(ci_bc)
matrix G24_BS_PC = e(ci_percentile)
restore
mata: A = (1::40); st_matrix("G24_A", A); st_matrix("G24_OB", 20.5); st_matrix("G24_CMP", 20.5)
tempname g24n g24p g24b g24pp g24bp g24z
mata: wooldid_bootci("G24_A", "G24_OB", "G24_CMP", 95, 0, "`g24n'", "`g24p'", "`g24b'", "`g24pp'", "`g24bp'", "`g24z'", .)
if `g24z'[1,1] != 0 | `g24b'[1,1] != G24_BS_BC[1,1] | `g24b'[2,1] != G24_BS_BC[2,1] | `g24p'[1,1] != G24_BS_PC[1,1] | `g24p'[2,1] != G24_BS_PC[2,1] | G24_BS_BC[2,1] != 39.5 {
  di as error "G24 FAIL: z0=" `g24z'[1,1] " BC=[" `g24b'[1,1] ", " `g24b'[2,1] "] vs bstat [" G24_BS_BC[1,1] ", " G24_BS_BC[2,1] "]; percentile=[" `g24p'[1,1] ", " `g24p'[2,1] "] vs bstat [" G24_BS_PC[1,1] ", " G24_BS_PC[2,1] "]"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G25 footer wording: the method label says observation-level on an unclustered run, and the
* contrast line under boottstat says studentized (both lines were unconditional)
use `panel', clear
qui log using boot_gates_g25.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(7)
local g25rc = _rc
qui log close bg
chk 25 0 `g25rc'
local lg = fileread("boot_gates_g25.txt")
if strpos(`"`lg'"', "observation-level bootstrap (resampling with replacement)") == 0 {
  di as error "G25b FAIL: the footer does not label the unclustered run as an observation-level bootstrap"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g25.txt
use `panel', clear
qui log using boot_gates_g25b.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) boottstat subgroups(s)
local g25rc = _rc
qui log close bg
chk 25c 0 `g25rc'
local lg = fileread("boot_gates_g25b.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "are bootstrapped as coefficient differences") > 0 | strpos(`"`lg'"', "studentized by their own analytic pairwise standard error") == 0 {
  di as error "G25d FAIL: the boottstat footer still describes the contrasts as coefficient differences"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* the positive half of the label check: this run IS clustered, so the label must say so
* (the label is keyed on e(wdid_vcespec); wooldid posts no e(clustvar))
if strpos(`"`lg'"', "cluster bootstrap (resampling with replacement)") == 0 {
  di as error "G25e FAIL: the footer does not label the clustered run as a cluster bootstrap"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g25b.txt

* G26 pretrendpower() under bootstrap on an att-only run: no 'bootstrap event-study
* covariance is unavailable' caution for families that were never estimated
use `panel', clear
qui log using boot_gates_g26.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) esfixedbaseperiod pretrendpower(50) bootstrap(standard) reps(10) rseed(7)
local g26rc = _rc
qui log close bg
chk 26 0 `g26rc'
local lg = fileread("boot_gates_g26.txt")
if strpos(`"`lg'"', "bootstrap event-study covariance is unavailable") > 0 {
  di as error "G26b FAIL: spurious covariance-unavailable caution for a family that was never estimated"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g26.txt

* G27 (F040, 2026-09-06) bootdf(clusters) without cluster() is refused (198) instead of
* resolving silently to t(N-1); bootdf(normal) and an integer stay legal unclustered
use `panel', clear
qui log using boot_gates_g27.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(7) bootdf(clusters)
local g27rc = _rc
qui log close bg
chk 27 198 `g27rc'
local lg = fileread("boot_gates_g27.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "bootdf(clusters) requires cluster()") == 0 {
  di as error "G27b FAIL: the refusal message does not name the requirement"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g27.txt
use `panel', clear
cap qui wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(7) bootdf(normal)
chk 27c 0 _rc
use `panel', clear
cap qui wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(7) bootdf(11)
chk 27d 0 _rc
if e(boot_df) != 11 {
  di as error "G27e FAIL: an integer bootdf() unclustered must post e(boot_df) = 11, got " e(boot_df)
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G28 (F050, 2026-09-06) the treated-cluster caution: fewer than 6 treated resampling
* clusters -> a caution naming permute() and boottest (pairs wording under standard, the
* moments wording under bayes); silent with 8 treated clusters; not keyed on the total
use `panel', clear
qui replace ttre = . if unit > 2
qui log using boot_gates_g28.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7)
local g28rc = _rc
qui log close bg
chk 28 0 `g28rc'
local lg = fileread("boot_gates_g28.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "only 2 treated cluster(s) are available for resampling") == 0 | strpos(`"`lg'"', "pairs cluster bootstrap, which keeps only the replications") == 0 | strpos(`"`lg'"', "permute()") == 0 | strpos(`"`lg'"', "boottest") == 0 {
  di as error "G28b FAIL: no treated-cluster caution with 2 treated clusters under bootstrap(standard)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
qui replace ttre = . if unit > 2
qui log using boot_gates_g28.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(10) rseed(7)
local g28rc = _rc
qui log close bg
chk 28c 0 `g28rc'
local lg = fileread("boot_gates_g28.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "only 2 treated cluster(s) are available for resampling") == 0 | strpos(`"`lg'"', "Bayesian bootstrap, which shares the pairs cluster bootstrap's first two moments") == 0 {
  di as error "G28d FAIL: no treated-cluster caution with 2 treated clusters under bootstrap(bayes)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
qui log using boot_gates_g28.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7)
local g28rc = _rc
qui log close bg
chk 28e 0 `g28rc'
local lg = fileread("boot_gates_g28.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "treated cluster(s) are available for resampling") > 0 {
  di as error "G28f FAIL: the treated-cluster caution fired with 8 treated clusters"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g28.txt

* G29 (F053, 2026-09-06) the statistic column is labelled t_wald / z_wald under the
* rank-based displays (and in the pctl/bc twins), t / z under bootdisplay(wald) and in the
* wald twin, t_pt under boottstat; e(boot_stat) matches the displayed label
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) subgroups(s) bootstrap(standard) reps(10) rseed(7)
local c29a : colnames e(wdidmainresults)
local c29w : colnames e(boot_wald_mainresults)
local c29p : colnames e(boot_pctl_mainresults)
local c29s : colnames e(wdidsgresults_1)
local c29pw : colnames e(pwsgs)
local c29pww : colnames e(boot_wald_pwsgs)
local b29a `e(boot_stat)'
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) bootdisplay(wald)
local c29b : colnames e(wdidmainresults)
local b29b `e(boot_stat)'
use `panel', clear
qui wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(7)
local c29c : colnames e(wdidmainresults)
local b29c `e(boot_stat)'
di "G29 default=[`c29a'] stat=`b29a' | wald=[`c29b'] stat=`b29b' | unclustered=[`c29c'] stat=`b29c'"
di "G29 twins wald=[`c29w'] pctl=[`c29p'] sg=[`c29s'] pwsgs=[`c29pw'] pwsgs_wald=[`c29pww']"
if !`: list posof "t_wald" in c29a' | `: list posof "t" in c29a' | "`b29a'" != "t_wald" ///
   | !`: list posof "t" in c29b' | `: list posof "t_wald" in c29b' | "`b29b'" != "t" ///
   | !`: list posof "z_wald" in c29c' | "`b29c'" != "z_wald" ///
   | !`: list posof "t" in c29w' | !`: list posof "t_wald" in c29p' | !`: list posof "t_wald" in c29s' ///
   | !`: list posof "t_wald" in c29pw' | !`: list posof "t" in c29pww' {
  di as error "G29 FAIL: the statistic column label does not follow the displayed family"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G30 (F038, 2026-09-06) a rank-based p of exactly 0 gets a line under the table and a
* footer sentence saying it means p < 2/n; the cell and e() keep the 0; nothing fires
* under bootdisplay(wald). The resolution prints unpadded in %9.2g - .2 at 2/10, .0004 at
* 2/4999, never 0.000 - and the boottstat footer's 2/(n + 1) prints the same way (G30e/f)
use `panel', clear
qui log using boot_gates_g30.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7)
local g30rc = _rc
qui log close bg
chk 30 0 `g30rc'
local p30 = e(wdidmainresults)[rowsof(e(wdidmainresults)), 4]
local lg = fileread("boot_gates_g30.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if `p30' != 0 | strpos(`"`lg'"', "p = 0 above means p < .2 (2/10 completed replications") == 0 | strpos(`"`lg'"', "The rank-based coefficient p-values have resolution 2/10 = .2:") == 0 {
  di as error "G30b FAIL: the p = 0 resolution line/footer did not print (att p = `p30')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
qui log using boot_gates_g30.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) bootdisplay(wald)
local g30rc = _rc
qui log close bg
chk 30c 0 `g30rc'
local lg = fileread("boot_gates_g30.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "p = 0 above means") > 0 | strpos(`"`lg'"', "have resolution 2/") > 0 {
  di as error "G30d FAIL: the resolution line fired under bootdisplay(wald)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
qui log using boot_gates_g30.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) boottstat
local g30rc = _rc
qui log close bg
chk 30e 0 `g30rc'
local n30 = e(boot_completed)
local lg = fileread("boot_gates_g30.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if `n30' != 10 | strpos(`"`lg'"', "smallest attainable value are both 2/(n + 1) = .18: a studentized p") == 0 {
  di as error "G30f FAIL: the boottstat footer does not print 2/(n + 1) as .18 at 10 completed replications (n = `n30')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* G30g (change audit round 2, U10 pfloor-bc-rank-wording; the author's ruling 2026-09-23): the
* BC p is the inverted BC interval and is off the 2/n grid, so under bootdisplay(bc) the p = 0
* line and the footer's resolution clause are omitted, here on a run whose att BC p is 0 (on
* 8d8b324 both printed); the footer keeps its joint-test sentence
use `panel', clear
qui log using boot_gates_g30.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) bootdisplay(bc)
local g30rc = _rc
qui log close bg
chk 30g 0 `g30rc'
local p30g = e(wdidmainresults)[rowsof(e(wdidmainresults)), 4]
local lg = fileread("boot_gates_g30.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local g30line = strpos(`"`lg'"', "p = 0 above means") > 0
local g30res = strpos(`"`lg'"', "have resolution 2/") > 0
local g30jt = strpos(`"`lg'"', "The joint, equality and adversarial p-values are add-one") > 0
di "G30g bc: att p = `p30g', p = 0 line `g30line', resolution clause `g30res', joint-test sentence `g30jt'"
if `p30g' != 0 | `g30line' | `g30res' | !`g30jt' {
  di as error "G30g FAIL: under bootdisplay(bc) the p = 0 line and the resolution clause must be omitted and the joint-test sentence kept (att p `p30g', line `g30line', clause `g30res', joint `g30jt')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g30.txt

* G31 (F059/F063/F066, 2026-09-06) the reproducibility record and the resample columns:
* e(boot_rngstate) is posted with or without rseed() and replays the draws exactly; the
* bootsave() file carries the engine/seed/rngstate characteristics, w_<cluster> columns
* (counts summing to G, observed row 1) under standard, Dirichlet multipliers summing to
* 1 (observed row 1/G) under bayes, no w_ columns unclustered, jointW/advW under boottstat;
* the projection line prints after the fifth replication
use `panel', clear
cap erase boot_gates_g31.dta
qui log using boot_gates_g31.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(8) rseed(7) bootsave(boot_gates_g31) dots(1)
local g31rc = _rc
qui log close bg
chk 31 0 `g31rc'
local g31rng `e(boot_rngstate)'
local lg = fileread("boot_gates_g31.txt")
if strlen("`g31rng'") < 100 | strpos(`"`lg'"', "projected for 8 replications)") == 0 {
  di as error "G31b FAIL: e(boot_rngstate) not posted or no runtime projection (rng length " strlen("`g31rng'") ")"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
preserve
usegate boot_gates_g31 G31r1
qui ds w_*
local g31w : word count `r(varlist)'
egen double g31s = rowtotal(w_*)
qui sum g31s if _rep > 0
local g31min = r(min)
local g31max = r(max)
qui sum g31s if _rep == 0
local g31obs = r(mean)
local g31eng `_dta[wooldid_engine]'
local g31seed `_dta[wooldid_seed]'
local g31rl = strlen("`_dta[wooldid_rngstate]'")
restore
di "G31 w-cols=`g31w' rowsums " `g31min' "-" `g31max' " obs=" `g31obs' " chars engine=`g31eng' seed=`g31seed' rnglen=`g31rl'"
if `g31w' != 12 | `g31min' != 12 | `g31max' != 12 | `g31obs' != 12 | "`g31eng'" != "bootstrap" | "`g31seed'" != "7" | `g31rl' < 100 {
  di as error "G31c FAIL: the bootsave() resample columns or characteristics are wrong"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* replay from the posted state: no rseed, then set rngstate and repeat
use `panel', clear
set seed 4242
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(8)
matrix G31D1 = e(bootdist)
local g31st `e(boot_rngstate)'
use `panel', clear
set rngstate `g31st'
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(8)
matrix G31D2 = e(bootdist)
if mreldif(G31D1, G31D2) != 0 {
  di as error "G31d FAIL: set rngstate e(boot_rngstate) does not replay the draws (mreldif " mreldif(G31D1, G31D2) ")"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* bayes: multipliers sum to 1, observed row 1/G
use `panel', clear
cap erase boot_gates_g31.dta
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(8) rseed(7) bootsave(boot_gates_g31)
chk 31e 0 _rc
preserve
usegate boot_gates_g31 G31r2
egen double g31s = rowtotal(w_*)
qui sum g31s if _rep > 0
local g31bmin = r(min)
local g31bmax = r(max)
qui sum g31s if _rep == 0
local g31bobs = r(mean)
restore
di "G31 bayes rowsums " %10.8f `g31bmin' "-" %10.8f `g31bmax' " obs=" %10.8f `g31bobs'
if abs(`g31bmin' - 1) > 1e-10 | abs(`g31bmax' - 1) > 1e-10 | abs(`g31bobs' - 1) > 1e-10 {
  di as error "G31f FAIL: the bayes resample columns do not sum to 1 (observed row 1/G)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* unclustered: no w_ columns; boottstat: jointW and advW columns with the observed row
use `panel', clear
cap erase boot_gates_g31.dta
cap qui wooldid y unit t ttre, att bootstrap(standard) reps(8) rseed(7) bootsave(boot_gates_g31)
chk 31g 0 _rc
preserve
usegate boot_gates_g31 G31r3
cap ds w_*
local g31u = _rc
restore
use `panel', clear
cap erase boot_gates_g31.dta
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) esfixedbaseperiod adversarial jointtests bootstrap(standard) reps(8) rseed(7) boottstat bootsave(boot_gates_g31)
chk 31h 0 _rc
preserve
usegate boot_gates_g31 G31r4
cap ds jointW*
local g31j : word count `r(varlist)'
cap confirm variable advW
local g31a = _rc
qui count if _rep == 0 & advW < .
local g31a0 = r(N)
restore
di "G31 unclustered w_ rc=`g31u' (want 111) jointW cols=`g31j' advW rc=`g31a' observed advW present=`g31a0'"
if `g31u' != 111 | `g31j' < 1 | `g31a' != 0 | `g31a0' != 1 {
  di as error "G31i FAIL: unclustered/boottstat file layout"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g31.dta
cap erase boot_gates_g31.txt

* G32 (F067, 2026-09-06) bootdisplay(bca): the acceleration equals the hand delete-one-
* cluster jackknife (G fastnose fits dropping each cluster), the BCa endpoints and p
* reconstruct from e(bootdist), e(boot_z0) and e(boot_accel) with the order-statistic rule,
* the displayed table carries them, the cells fall back to bc with a caution, and bca is
* refused with boottstat
use `panel', clear
qui log using boot_gates_g32.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(40) rseed(7) bootdisplay(bca) gteffects esprelength(1) espostlength(1)
local g32rc = _rc
qui log close bg
chk 32 0 `g32rc'
matrix G32M = e(wdidmainresults)
matrix G32D = e(bootdist)
matrix G32Z = e(boot_z0)
matrix G32A = e(boot_accel)
matrix G32CI = e(boot_ci_bca)
matrix G32P = e(boot_p_bca)
local g32jk = e(boot_jk_completed)
local g32fam `e(boot_gteffects_family)'
local g32cn : colnames G32A
local g32k : list posof "main_att" in g32cn
local lg = fileread("boot_gates_g32.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local g32caut = strpos(`"`lg'"', "cells carry the bias-corrected (bc) family under bootdisplay(bca)") > 0
* hand jackknife: drop each unit in turn, fastnose
local g32hand
forval g = 1/12 {
  use `panel', clear
  qui drop if unit == `g'
  * the SAME specification as the bootstrap run (the event-study window changes the
  * regression and with it the pooled att; gteffects does not)
  qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) fastnose
  local g32hand `g32hand' `=string(e(b)[1, colnumb(e(b), "main_att")], "%21x")'
}
mata: th = strtoreal(tokens(st_local("g32hand")))'; d2 = sum((mean(th) :- th):^2); d3 = sum((mean(th) :- th):^3); st_numscalar("g32_a", d3/(6*d2^1.5))
* hand endpoints and p from the draws with the same rule as the shipped bootqs
local g32dn : colnames G32D
local g32kd : list posof "main_att" in g32dn
mata: v = sort(select(st_matrix("G32D")[., `g32kd'], st_matrix("G32D")[., `g32kd'] :< .), 1); n = rows(v); z0 = st_matrix("G32Z")[1, `g32k']; a = st_matrix("G32A")[1, `g32k']; al = .025; w = z0 + invnormal(al); p1 = normal(z0 + w/(1 - a*w)); w = z0 + invnormal(1 - al); p2 = normal(z0 + w/(1 - a*w)); x = n*p1; lo = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))]); x = n*p2; hi = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))]); r0 = sum(v :<= 0)/n; q = (r0 <= 0 ? 0 : (r0 >= 1 ? 1 : normal((invnormal(r0) - z0)/(1 + a*(invnormal(r0) - z0)) - z0))); st_numscalar("g32_lo", lo); st_numscalar("g32_hi", hi); st_numscalar("g32_p", min((1, 2*min((q, 1 - q)))))
* the pooled att row of the main table (the last row is an event-study horizon here)
local r32 0
local g32rn : rowfullnames G32M
local g32ri 0
foreach rn of local g32rn {
  local ++g32ri
  if "`rn'" == "main:att" local r32 `g32ri'
}
if `r32' == 0 local r32 1
di "G32 a=" %12.8f G32A[1,`g32k'] " hand=" %12.8f g32_a " ci=[" %12.8f G32CI[1,`g32k'] "," %12.8f G32CI[2,`g32k'] "] hand=[" %12.8f g32_lo "," %12.8f g32_hi "] p=" %10.8f G32P[1,`g32k'] " hand=" %10.8f g32_p " table lb/ub/p=" %12.8f G32M[`r32',5] "/" %12.8f G32M[`r32',6] "/" %10.8f G32M[`r32',4] " jk=`g32jk' cells=`g32fam' caution=`g32caut'"
if reldif(G32A[1,`g32k'], g32_a) > 1e-10 | reldif(G32CI[1,`g32k'], g32_lo) > 1e-12 | reldif(G32CI[2,`g32k'], g32_hi) > 1e-12 | reldif(G32P[1,`g32k'], g32_p) > 1e-12 ///
   | G32M[`r32',5] != G32CI[1,`g32k'] | G32M[`r32',6] != G32CI[2,`g32k'] | G32M[`r32',4] != G32P[1,`g32k'] | `g32jk' != 12 | "`g32fam'" != "bc" | !`g32caut' {
  di as error "G32b FAIL: the BCa family does not reconstruct (acceleration, endpoints, p, table, cells)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) bootdisplay(bca) boottstat
chk 32c 198 _rc
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7)
cap confirm matrix e(boot_ci_bca)
chk 32d 111 _rc
cap erase boot_gates_g32.txt

* G33 (CS-87, 2026-09-06) makeplots with no event-study window under bootstrap() used to
* die r(2000) in the deferred plot call (the analytic path skips the plot); now rc 0 with
* the results intact, and makeplots WITH a window still draws
use `panel', clear
set graphics off
cap qui wooldid y unit t ttre, att cluster(unit) makeplots bootstrap(standard) reps(5) rseed(7)
chk 33 0 _rc
cap confirm matrix e(wdidmainresults)
chk 33b 0 _rc
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) makeplots bootstrap(standard) reps(5) rseed(7)
chk 33c 0 _rc

set linesize 255
cap prog drop unwrap
prog def unwrap
  args fname
  local lg = fileread("`fname'")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  c_local lg `"`lg'"'
end
* G34 (bug sweep BS-16, 2026-09-06): bootdisplay(bca) without cluster() is refused (the
* jackknife would be one fit per observation); with cluster() the jackknife runs once per
* cluster, the acceleration is posted, and the runtime projection counts the extra fits
use `panel', clear
qui log using boot_gates_g34.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(8) rseed(7) bootdisplay(bca)
local g34rc = _rc
qui log close bg
chk 34 198 `g34rc'
unwrap boot_gates_g34.txt
if strpos(`"`lg'"', "bootdisplay(bca) requires cluster()") == 0 {
  di as error "G34b FAIL: the bca-without-cluster refusal does not name the reason"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g34.txt
use `panel', clear
qui log using boot_gates_g34.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(8) rseed(7) bootdisplay(bca) dots(1)
local g34rc = _rc
qui log close bg
chk 34c 0 `g34rc'
local g34jk = e(boot_jk_completed)
mata: st_local("g34acc", strofreal(sum(st_matrix("e(boot_accel)") :< .)))
unwrap boot_gates_g34.txt
di "G34 jk_completed=`g34jk' accel non-missing=`g34acc' projection has jackknife=" strpos(`"`lg'"', "+ 12 jackknife fits")
if `g34jk' != 12 | `g34acc' == 0 | strpos(`"`lg'"', "+ 12 jackknife fits") == 0 {
  di as error "G34d FAIL: the clustered bca run should jackknife 12 clusters, post an acceleration and project the extra fits"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g34.txt

* G35 (BS-17): with two clusters the jackknife has at most two fits, so no acceleration exists;
* the BCa p-value and interval are blank, but the bootstrap SD (and b/sd) survive the
* display switch, the caution names the blank family, and the footer stops claiming an
* acceleration
use `panel', clear
gen s2 = (unit <= 6) + 1
qui log using boot_gates_g35.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(s2) bootstrap(standard) reps(8) rseed(7) bootdisplay(bca)
local g35rc = _rc
qui log close bg
chk 35 0 `g35rc'
local g35jk = e(boot_jk_completed)
mata: st_local("g35cimiss", strofreal(allof(st_matrix("e(boot_ci_bca)"), .)))
matrix G35M = e(wdidmainresults)
matrix G35B = e(boot_bc_mainresults)
local g35sd = G35M[1,2]
local g35sdbc = G35B[1,2]
mata: st_local("g35blank", strofreal(allof(st_matrix("G35M")[1, 4..6], .)))
local lg = fileread("boot_gates_g35.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local g35c = strpos(`"`lg'"', "had no usable BCa family") > 0
local g35f = strpos(`"`lg'"', "No BCa acceleration could be computed") > 0
local g35g = strpos(`"`lg'"', "The BCa acceleration comes from") > 0
di "G35 jk=`g35jk' ci_bca all missing=`g35cimiss' sd=`g35sd' (bc sd `g35sdbc') cols4-6 blank=`g35blank' caution=`g35c' footer-none=`g35f' footer-claims=`g35g'"
if `g35jk' > 2 | `g35cimiss' != 1 | `g35sd' >= . | `g35sd' != `g35sdbc' | `g35blank' != 1 | !`g35c' | !`g35f' | `g35g' {
  di as error "G35b FAIL: an unavailable BCa family should blank only the p/interval, with the caution and the corrected footer"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g35.txt

* G36 (BS-18): bootsave() under bootalpha(size) writes row 0 on the replication scale -
* the invariant sum_c n_c w_c = 1 (n_c = 8 rows per unit here, so the row total of the
* w_ columns is 1/8) holds on every row including _rep = 0. n_c, the cluster's ROW COUNT,
* is the divisor since the 2026-09-13 ruling; unweighted it equals the base-weight total
* this leg was written against, so the expected numbers are unchanged
use `panel', clear
cap erase boot_gates_g36.dta
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) bootalpha(size) reps(8) rseed(7) bootsave(boot_gates_g36)
local g36rc = _rc
chk 36 0 `g36rc'
preserve
usegate boot_gates_g36 G36r
egen double g36s = rowtotal(w_*)
qui sum g36s
local g36min = r(min)
local g36max = r(max)
qui sum g36s if _rep == 0
local g36obs = r(mean)
restore
di "G36 size: w row totals " `g36min' "-" `g36max' " obs row " `g36obs' " (expect 1/8 = .125 throughout)"
if abs(`g36min' - 1/8) > 1e-9 | abs(`g36max' - 1/8) > 1e-9 | abs(`g36obs' - 1/8) > 1e-9 {
  di as error "G36b FAIL: bootsave() row 0 under bootalpha(size) is not on the replication scale"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g36.dta

* G37 (BS-31): the bayes header names the resampling unit, as the standard one does
use `panel', clear
qui log using boot_gates_g37.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(bayes) reps(8) rseed(7)
local g37rc = _rc
qui log close bg
chk 37 0 `g37rc'
unwrap boot_gates_g37.txt
local g37a = strpos(`"`lg'"', "Bayesian bootstrap (Dirichlet reweighting over observations") > 0
local g37b = strpos(`"`lg'"', "Bayesian cluster bootstrap") > 0
cap erase boot_gates_g37.txt
use `panel', clear
qui log using boot_gates_g37.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(8) rseed(7)
local g37rc = _rc
qui log close bg
chk 37c 0 `g37rc'
unwrap boot_gates_g37.txt
local g37c = strpos(`"`lg'"', "Bayesian cluster bootstrap") > 0
local g37d = strpos(`"`lg'"', "over observations") > 0
cap erase boot_gates_g37.txt
di "G37 unclustered: observations=`g37a' cluster-label=`g37b' | clustered: cluster-label=`g37c' observations=`g37d'"
if !`g37a' | `g37b' | !`g37c' | `g37d' {
  di as error "G37e FAIL: the bayes header does not follow the resampling unit"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G38 (BS-30): adversarial under bootstrap() without boottstat has no bootstrap p-value;
* the display must not print a bare "." line, and the bare replay discloses the absence.
* G38e (2026-09-23, author ruling: no line-only results, text approved): the fitted line is
* not posted either, and the queued caution says the test failed
use `panel', clear
qui log using boot_gates_g38.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(2) espostlength(1) adversarial bootstrap(standard) reps(8) rseed(7)
local g38rc = _rc
qui log close bg
chk 38 0 `g38rc'
unwrap boot_gates_g38.txt
local g38a = strpos(`"`lg'"', "single line: .") > 0
local g38adv `e(boot_adv_absent)'
cap confirm matrix e(adversarial_line_att)
local g38l = _rc == 0
local g38t = strpos(`"`lg'"', "Caution: the adversarial pre-trend test failed: under bootstrap() it needs boottstat.") > 0
di "G38e line=`g38l' failed-caution=`g38t'"
if `g38l' | !`g38t' {
  di as error "G38e FAIL: under bootstrap() without boottstat the adversarial line must not be posted and the caution must say the test failed"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
qui log using boot_gates_g38b.txt, text replace name(bg)
cap noi wooldid
local g38rc2 = _rc
qui log close bg
unwrap boot_gates_g38b.txt
local lg2 `"`lg'"'
local g38b = strpos(`"`lg2'"', "single line: .") > 0
local g38c = strpos(`"`lg2'"', "adversarial line tests are NOT reported") > 0
cap erase boot_gates_g38.txt
cap erase boot_gates_g38b.txt
di "G38 live bare-dot=`g38a' e(boot_adv_absent)=`g38adv' | replay rc=`g38rc2' bare-dot=`g38b' disclosure=`g38c'"
if `g38a' | "`g38adv'" != "1" | `g38rc2' != 0 | `g38b' | !`g38c' {
  di as error "G38d FAIL: the absent adversarial bootstrap p prints as a bare dot or is not disclosed on replay"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G39 (CS-20, author ruling 2026-09-07 "caution, no floor"): reps() gains permreps()'s
* resolution caution and keeps its absence of a minimum. Which number binds depends on the
* family the table carries, so all three branches are checked, together with the silence
* above the threshold and the fact that a low-rep run still RUNS.
use `panel', clear
qui log using boot_gates_g39.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(3)
local g39rc = _rc
qui log close bg
chk 39 0 `g39rc'
unwrap boot_gates_g39.txt
local g39a = strpos(`"`lg'"', "the coefficient p-values fall on a grid of 0.100 (2/reps)") > 0
local g39a2 = strpos(`"`lg'"', "the only value below .05 this run can produce is exactly 0") > 0
use `panel', clear
qui log using boot_gates_g39.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(3) boottstat
local g39rcb = _rc
qui log close bg
chk 39b 0 `g39rcb'
unwrap boot_gates_g39.txt
local g39b = strpos(`"`lg'"', "the smallest studentized coefficient p-value this run can produce is 0.095") > 0
use `panel', clear
qui log using boot_gates_g39.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(3) bootdisplay(wald)
local g39rcc = _rc
qui log close bg
chk 39c 0 `g39rcc'
unwrap boot_gates_g39.txt
local g39c = strpos(`"`lg'"', "the joint, equality and adversarial p-values - one-sided add-one shares") > 0
* above the threshold nothing is said: 2/60 = .033 for the coefficient grid
use `panel', clear
qui log using boot_gates_g39.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(60) rseed(3)
local g39rcd = _rc
qui log close bg
chk 39d 0 `g39rcd'
unwrap boot_gates_g39.txt
local g39d = strpos(`"`lg'"', "Caution: with reps(") > 0
cap erase boot_gates_g39.txt
di "G39 grid=`g39a'/`g39a2' studentized=`g39b' wald-joint=`g39c' quiet-at-60=`g39d'"
if !`g39a' | !`g39a2' | !`g39b' | !`g39c' | `g39d' {
  di as error "G39e FAIL: the reps() resolution caution is wrong somewhere (grid=`g39a'/`g39a2' studentized=`g39b' wald=`g39c' fired-at-60=`g39d')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G39f (change audit 2026-09-22, U08_bootstrap / bootctl-respcaution-clusteronly): the same
* caution on runs WITHOUT cluster(). The block sat inside the cluster() branch of
* wooldid_bootctl, so an unclustered bootstrap at low reps() ran silently on every path. The
* grid branch is checked and required to print exactly once, and so is the bootdisplay(wald)
* branch, the one path whose floor no post-table footer discloses.
use `panel', clear
qui log using boot_gates_g39.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(20) rseed(3)
local g39rcf = _rc
qui log close bg
chk 39f 0 `g39rcf'
unwrap boot_gates_g39.txt
local g39f = strpos(`"`lg'"', "the coefficient p-values fall on a grid of 0.100 (2/reps)") > 0
local g39fn = (length(`"`lg'"') - length(subinstr(`"`lg'"', "Caution: with reps(", "", .))) / length("Caution: with reps(")
use `panel', clear
qui log using boot_gates_g39.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att bootstrap(standard) reps(10) rseed(3) bootdisplay(wald)
local g39rcg = _rc
qui log close bg
chk 39g 0 `g39rcg'
unwrap boot_gates_g39.txt
local g39g = strpos(`"`lg'"', "the joint, equality and adversarial p-values - one-sided add-one shares") > 0
cap erase boot_gates_g39.txt
di "G39f unclustered grid=`g39f' (printed `g39fn' times) wald-joint=`g39g'"
if !`g39f' | `g39fn' != 1 | !`g39g' {
  di as error "G39f FAIL: the reps() resolution caution is missing or repeated on an unclustered bootstrap (grid=`g39f' printed=`g39fn' wald=`g39g')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G40 (permboot_tracking audit 2026-09-13): the three boottstat refusals that test_bootstrap_joint
* :140-148 relies on but nothing pins - poisson (its replications are estimated on sampling
* weights, so a replication's own standard error is not the same functional as the observed
* one), bootstrap(bayes) (continuous Dirichlet mass is not a duplication count to studentize
* over), bootdisplay(wald) (a studentized run reports the analytic se; no bootstrap sd).
use `panel', clear
set seed 40
qui gen ycount = rpoisson(1.5 + (t>=ttre & ttre<.))
cap qui wooldid ycount unit t ttre, att cluster(unit) poisson bootstrap(standard) reps(50) boottstat
chk 40a 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(50) boottstat
chk 40b 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) boottstat bootdisplay(wald)
chk 40c 198 _rc
* and the neighbour that must STAY open: bc/bca are refused with boottstat too (T11/G32c
* already pin them); percentile is the studentized display and must run
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) boottstat bootdisplay(percentile) rseed(7) dots(0)
chk 40d 0 _rc

* ---- G38 (2026-09-21 ruling, built ed6a0da; gate added 2026-09-22): a bootsave() write that
* fails after the loop is a caution, not a refusal. Forced with a read-only target: the run
* completes rc 0, the results stand, the caution names the file, e(boot_save_rc) carries the
* code, the previous file is untouched, and a bare wooldid replays
use `panel', clear
cap erase boot_gates_g38.dta
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(8) rseed(7) bootsave(boot_gates_g38)
chk 38 0 _rc
cap confirm file boot_gates_g38.dta
chk 38a 0 _rc
! chmod 444 boot_gates_g38.dta
use `panel', clear
qui log using boot_gates_g38.txt, text replace name(g38)
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(8) rseed(8) bootsave(boot_gates_g38) replace dots(0)
local g38rc = _rc
qui log close g38
chk 38b 0 `g38rc'
local g38src `e(boot_save_rc)'
local g38lg = fileread("boot_gates_g38.txt")
if strpos(`"`g38lg'"', "Caution: bootsave(boot_gates_g38) could not produce the replication dataset") == 0 | "`g38src'" == "" | "`g38src'" == "0" {
  di as error "G38c FAIL: the failed write did not caution, or did not post e(boot_save_rc) (got `g38src')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* G38g the caution prints below the result tables, after the part(boot) footer, with the
* run's other end-of-run cautions rather than above the deferred display
local g38pf = strpos(`"`g38lg'"', "Bootstrap Inference (")
local g38pc = strpos(`"`g38lg'"', "Caution: bootsave(boot_gates_g38) could not produce the replication dataset")
if `g38pf' == 0 | `g38pc' < `g38pf' {
  di as error "G38g FAIL: the failed-write caution (at `g38pc') does not print below the bootstrap footer (at `g38pf')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap confirm matrix e(bootdist)
if _rc != 0 | e(boot_completed) != 8 {
  di as error "G38d FAIL: the results are not intact after the failed write"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap noi wooldid
chk 38e 0 _rc
preserve
usegate boot_gates_g38 G38r
qui count if _rep > 0
local g38n = r(N)
local g38seed `_dta[wooldid_seed]'
restore
if `g38n' != 8 | "`g38seed'" != "7" {
  di as error "G38f FAIL: the read-only file was altered (rows=`g38n' seed=`g38seed'; want 8 and 7)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
! chmod 644 boot_gates_g38.dta
cap erase boot_gates_g38.dta
cap erase boot_gates_g38.txt

* ---- G41 saveplots-no-precheck (U12 plots review, 2026-09-23): under bootstrap() and permute()
* the plots are drawn after the replication loop, from the rewired e(), so a saveplots() file
* that already existed without replace was found only once every replication had run - r(602)
* after the loop, the run lost - where bootsave() and permsave() refuse a clash before it. It is
* refused before the loop now (wooldid_plotpre). A seed run lists the files a call writes; then,
* for each one alone, a stand-in of that name must make the rerun without replace fail r(602)
* with no replication run (e(boot_method)/e(perm_method) empty: e() holds the observed results)
* and nothing else written, and the observed results must replay, as after a failed loop. With
* the files gone the same call without replace must still run and write the same files (a check
* that refused it would be a regression). The three calls cover every naming rule the check
* mirrors: the cfx .gph files and _cfx.dta, the event-study families and _ES.dta, and the
* subgroup _sg<case>.gph and _ES_<sgvar>_<case>.dta twins
cap prog drop g41leg
prog def g41leg
  args lab stem dataset key
  local cmd $G41CMD
  local W : dir . files "`stem'_*"
  local W : subinstr local W `"""' "", all
  foreach f of local W {
    cap erase "`f'"
  }
  use "`dataset'", clear
  cap qui wooldid `cmd' saveplots(`stem') replace
  local rc = _rc
  local W : dir . files "`stem'_*"
  local W : subinstr local W `"""' "", all
  local W : list sort W
  if `rc' != 0 | "`W'" == "" {
    di as error "G41`lab' FAIL: the seed run failed (rc `rc') or wrote no `stem'_* file"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
    exit
  }
  foreach f of local W {
    foreach g of local W {
      cap erase "`g'"
    }
    * a stand-in for the stale file: the check is on existence
    tempname fh
    file open `fh' using "`f'", write replace
    file close `fh'
    use "`dataset'", clear
    cap qui wooldid `cmd' saveplots(`stem')
    local rc = _rc
    local eng `e(`key')'
    local left : dir . files "`stem'_*"
    local left : subinstr local left `"""' "", all
    if `rc' != 602 | "`eng'" != "" | "`left'" != "`f'" {
      di as error "G41`lab' FAIL: with `f' present, rc `rc' (want 602), e(`key')=[`eng'] (want empty: refused before the loop), `stem'_* files [`left'] (want only `f')"
      global BOOTG_FAILS = ${BOOTG_FAILS} + 1
    }
    cap qui wooldid
    if _rc != 0 {
      di as error "G41`lab' FAIL: the observed results do not replay after the refusal of `f'"
      global BOOTG_FAILS = ${BOOTG_FAILS} + 1
    }
  }
  foreach g of local W {
    cap erase "`g'"
  }
  use "`dataset'", clear
  cap qui wooldid `cmd' saveplots(`stem')
  local rc = _rc
  local eng `e(`key')'
  local W2 : dir . files "`stem'_*"
  local W2 : subinstr local W2 `"""' "", all
  local W2 : list sort W2
  if `rc' != 0 | "`eng'" == "" | "`W2'" != "`W'" {
    di as error "G41`lab' FAIL: the clean run without replace: rc `rc' (want 0), e(`key')=[`eng'], files [`W2'] (want [`W'])"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
  foreach g of local W2 {
    cap erase "`g'"
  }
end
global G41CMD y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots makeplots bootstrap(standard) reps(10) rseed(7) dots(0)
g41leg a boot_gates_g41a "data/dta_cont" boot_method
global G41CMD y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots makeplots permute(shuffle) permreps(20) rseed(7) dots(0)
g41leg b boot_gates_g41b "data/dta_cont" perm_method
use `panel', clear
gen sgv = mod(unit, 2)
tempfile panelsg
save `panelsg'
global G41CMD y unit t ttre, att att_it cluster(unit) subgroups(sgv) esprelength(1) espostlength(1) makeplots bootstrap(standard) reps(10) rseed(7) dots(0)
g41leg c boot_gates_g41c "`panelsg'" boot_method
global G41CMD

* ---- G42 (author ruling 2026-09-23, save first and caution last): bootsave() writes its file
* before the deferred display and plots, so a plot that fails after the replication loop does
* not lose it, while the saver's line still prints below the result tables. A first run writes
* the draws file and the event-study .gph (G42a: rc 0, and the success line prints below the
* part(boot) footer). The draws file is then erased and the .gph made read-only, as in G20, so
* the rerun under replace fails r(608) at the plot save (G42b). The draws file must exist
* regardless, holding the rerun's replications and seed (G42c)
use `panel', clear
foreach f in boot_gates_g42.dta boot_gates_g42p_ES_att.gph boot_gates_g42p_ES.dta boot_gates_g42.txt {
  cap erase "`f'"
}
qui log using boot_gates_g42.txt, text replace name(g42)
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) dots(0) bootsave(boot_gates_g42) makeplots saveplots(boot_gates_g42p)
local g42rc = _rc
qui log close g42
local g42lg = fileread("boot_gates_g42.txt")
local g42pf = strpos(`"`g42lg'"', "Bootstrap Inference (")
local g42ps = strpos(`"`g42lg'"', "bootsave(): the replication draws were written to boot_gates_g42.dta")
if `g42rc' != 0 | `g42pf' == 0 | `g42ps' < `g42pf' {
  di as error "G42a FAIL: the seed run (rc `g42rc') did not print the bootsave() success line (at `g42ps') below the bootstrap footer (at `g42pf')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g42.dta
! chmod 444 boot_gates_g42p_ES_att.gph
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(8) dots(0) bootsave(boot_gates_g42) makeplots saveplots(boot_gates_g42p) replace
chk 42b 608 _rc
cap confirm file boot_gates_g42.dta
if _rc {
  di as error "G42c FAIL: the plot failure lost the bootsave() file"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
else {
  preserve
  use boot_gates_g42, clear
  qui count if _rep > 0
  local g42n = r(N)
  local g42seed `_dta[wooldid_seed]'
  restore
  if `g42n' != 10 | "`g42seed'" != "8" {
    di as error "G42c FAIL: the bootsave() file is not the rerun's (rows=`g42n' seed=`g42seed'; want 10 and 8)"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
}
! chmod 644 boot_gates_g42p_ES_att.gph
foreach f in boot_gates_g42.dta boot_gates_g42p_ES_att.gph boot_gates_g42p_ES.dta boot_gates_g42.txt {
  cap erase "`f'"
}

* G43 (change audit round 2, U10 replay-extra-blank-line-engines / display-test-gaps): the bare
* replay of a bootstrap() run is line-identical to the live run from the "Wooldid Estimation"
* banner to EOF, by test_fastnose T5's method. The replay branch printed a blank line before
* the part(boot) footer that the live run does not. The run has no bootsave() and stays below
* the 10,000-draw cap: the lines those print below the footer are, by design, not replayed
use `panel', clear
qui log using boot_gates_g43a.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) dots(0)
qui log close bg
local g43rc1 = _rc
qui log using boot_gates_g43b.txt, text replace name(bg)
cap noi wooldid
qui log close bg
local g43rc2 = _rc
chk 43a 0 `g43rc1'
chk 43b 0 `g43rc2'
local g43same 0
local g43n 0
mata: a_ = cat("boot_gates_g43a.txt") \ ""; b_ = cat("boot_gates_g43b.txt") \ ""
mata: ia_ = selectindex(strpos(a_, "Wooldid Estimation") :== 1); ib_ = selectindex(strpos(b_, "Wooldid Estimation") :== 1)
mata: hb_ = (rows(ia_) > 0) & (rows(ib_) > 0)
mata: sa_ = a_[|(ia_ \ 1)[1] \ rows(a_)|]; sb_ = b_[|(ib_ \ 1)[1] \ rows(b_)|]
mata: st_local("g43same", strofreal(hb_ & (rows(sa_) > 20) & (sa_ == sb_))); st_local("g43n", strofreal(rows(sa_)))
di "G43 replay line-identical from the banner=`g43same' (`g43n' live lines)"
if `g43same' != 1 {
  di as error "G43c FAIL: the bootstrap() replay differs from the live run from the banner onward (`g43n' live lines)"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
cap erase boot_gates_g43a.txt
cap erase boot_gates_g43b.txt

* G44 (change audit round 2, U08 test-gap-straddle-unclustered): a PARTIAL straddle. The
* baselines pin e(boot_nstraddle) only at 0 (re-id runs) and at every unit; nothing pinned the
* case in between, 0 < e(boot_nstraddle) < units with e(boot_reid) 0, or its caution. 24 units
* in 6 states; unit 4 moves to state 2 after t 4, so under cluster(cs) it alone straddles two
* clusters, the run keeps every unit's identity (no re-id), and the straddle caution prints
* once. The twin at cluster(s), which nests unit, re-identifies and stays silent
clear
set seed 2026
set sortseed 2026
qui set obs 24
gen unit = _n
gen s = ceil(unit/4)
gen ttre = cond(unit<=8, 4, cond(unit<=16, 6, .))
qui expand 8
bysort unit: gen t = _n
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + 0.5*(t>=ttre & ttre<.)*mod(unit,4) + rnormal()
gen cs = s
qui replace cs = 2 if unit == 4 & t > 4
tempfile strad
save `strad'
qui log using boot_gates_g44.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att_i cluster(cs) esprelength(1) espostlength(1) bootstrap(standard) reps(8) rseed(5) dots(0)
local g44rc = _rc
qui log close bg
chk 44a 0 `g44rc'
local g44ns = e(boot_nstraddle)
local g44ri = e(boot_reid)
unwrap boot_gates_g44.txt
local g44n = (length(`"`lg'"') - length(subinstr(`"`lg'"', "straddle two or more", "", .))) / length("straddle two or more")
use `strad', clear
qui log using boot_gates_g44.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att_i cluster(s) esprelength(1) espostlength(1) bootstrap(standard) reps(8) rseed(5) dots(0)
local g44rcb = _rc
qui log close bg
chk 44b 0 `g44rcb'
local g44nsb = e(boot_nstraddle)
local g44rib = e(boot_reid)
unwrap boot_gates_g44.txt
local g44nb = (length(`"`lg'"') - length(subinstr(`"`lg'"', "straddle two or more", "", .))) / length("straddle two or more")
cap erase boot_gates_g44.txt
di "G44 partial: nstraddle=`g44ns' reid=`g44ri' caution printed `g44n' time(s) | nesting twin: nstraddle=`g44nsb' reid=`g44rib' printed `g44nb'"
if "`g44ns'" != "1" | "`g44ri'" != "0" | `g44n' != 1 {
  di as error "G44c FAIL: a one-unit straddle should post e(boot_nstraddle) 1 and e(boot_reid) 0 and caution once (got `g44ns', `g44ri', `g44n')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
if "`g44nsb'" != "0" | "`g44rib'" != "1" | `g44nb' != 0 {
  di as error "G44d FAIL: the nesting twin should post e(boot_nstraddle) 0 and e(boot_reid) 1 and stay silent (got `g44nsb', `g44rib', `g44nb')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
* G44e (change audit round 2, U08 bootrun-straddle-runlevel; the author's ruling 2026-09-23,
* card m277): boottstat under cluster() takes the same leg with att alone, so the straddle is
* counted and the two scalars are posted as before, but the caution names att_i, att_it and
* att_itime, none of which is in the run, and prints only when one of them is
use `strad', clear
qui log using boot_gates_g44.txt, text replace name(bg)
cap noi wooldid y unit t ttre, att cluster(cs) esprelength(1) espostlength(1) bootstrap(standard) boottstat reps(8) rseed(5) dots(0)
local g44rcc = _rc
qui log close bg
chk 44e 0 `g44rcc'
local g44nsc = e(boot_nstraddle)
local g44ric = e(boot_reid)
unwrap boot_gates_g44.txt
local g44nc = (length(`"`lg'"') - length(subinstr(`"`lg'"', "straddle two or more", "", .))) / length("straddle two or more")
cap erase boot_gates_g44.txt
di "G44e att-only boottstat: nstraddle=`g44nsc' reid=`g44ric' caution printed `g44nc' time(s)"
if "`g44nsc'" != "1" | "`g44ric'" != "0" | `g44nc' != 0 {
  di as error "G44e FAIL: an att-only boottstat run should post e(boot_nstraddle) 1 and e(boot_reid) 0 and print no straddle caution (got `g44nsc', `g44ric', `g44nc')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}

* G45 (change audit round 2, S5 accum-writeback-cubic-cost): the per-replication write-back of
* the R-row accumulators must not scale with R. Each replication wrote its row back with
* st_matrix() + st_matrixcolstripe(), and st_matrix() rebuilds the default r1..rR row stripe on
* every call, so one call cost about 0.2 s at R = 8,000 and a run about R^3. The kernels are
* called here directly (as test_permute calls wooldid_permx), 20 times each against 8,000-row
* accumulators: bootgrab (every replication of both engines), bootwrow (the resample counts
* under bootsave()), bootcells (the gteffects cells) and bootjwr (the boottstat joint W).
* Each must take well under a second in all, write its row, and keep the stripes it was
* allocated with. The 1-second bar sits about 4x below the old cost and far above the new
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(1) gteffects
chk 45a 0 _rc
local g45R 8000
matrix G45O = e(b)
local g45cn : colfullnames e(b)
matrix G45A = J(`g45R', colsof(G45O), .)
matrix colnames G45A = `g45cn'
mata: wooldid_bootcellkeys("e(gteffects)", "G45C_gt", `g45R', 1)
mata: wooldid_bootjobs("att", "G45J")
matrix G45JW = J(`g45R', 1, .)
matrix colnames G45JW = att
local g45rn0 : rownames G45A
local g45gcn0 : colfullnames G45C_gt
timer clear 45
timer clear 46
timer clear 47
timer clear 48
forval r = 1/20 {
  timer on 45
  mata: wooldid_bootgrab("G45A", `r', "G45O")
  timer off 45
  timer on 47
  mata: wooldid_bootcells("G45C", `r', 0, 0)
  timer off 47
  timer on 48
  mata: wooldid_bootjwr("G45JW", `r', "G45J")
  timer off 48
}
preserve
clear
qui set obs 12
gen d = _n
matrix G45W = J(`g45R', 12, .)
local g45wn
forval c = 1/12 {
  local g45wn `g45wn' w_`c'
}
matrix colnames G45W = `g45wn'
forval r = 1/20 {
  timer on 46
  mata: wooldid_bootwrow("G45W", `r', "d")
  timer off 46
}
restore
qui timer list
local g45t1 = r(t45)
local g45t2 = r(t46)
local g45t3 = r(t47)
local g45t4 = r(t48)
* the rows landed and the stripes survived
mata: st_local("g45ok1", strofreal(mreldif(st_matrix("G45A")[20, .], st_matrix("G45O")) == 0 & missing(st_matrix("G45A")[21, .]) == cols(st_matrix("G45O"))))
mata: st_local("g45ok2", strofreal(st_matrix("G45W")[20, .] == (1..12) & missing(st_matrix("G45W")[21, .]) == 12))
mata: st_local("g45ok3", strofreal(!missing(st_matrix("G45C_gt")[20, 1]) & missing(st_matrix("G45C_gt")[21, 1])))
mata: st_local("g45ok4", strofreal(st_matrix("G45JW")[20, 1] < . & st_matrix("G45JW")[21, 1] >= .))
local g45cn1 : colfullnames G45A
local g45rn1 : rownames G45A
local g45wn1 : colfullnames G45W
local g45gcn1 : colfullnames G45C_gt
local g45jn1 : colfullnames G45JW
local g45str = ("`g45cn1'" == "`g45cn'") & ("`g45rn1'" == "`g45rn0'") & ("`: list retokenize g45wn1'" == "`: list retokenize g45wn'") & ("`g45gcn1'" == "`g45gcn0'") & ("`g45jn1'" == "att")
di "G45 20 calls at R=`g45R': bootgrab `g45t1' s, bootwrow `g45t2' s, bootcells `g45t3' s, bootjwr `g45t4' s | rows ok `g45ok1'`g45ok2'`g45ok3'`g45ok4', stripes kept `g45str'"
if !(`g45t1' < 1) | !(`g45t2' < 1) | !(`g45t3' < 1) | !(`g45t4' < 1) {
  di as error "G45b FAIL: a per-replication accumulator write-back still scales with the replication count (`g45t1' / `g45t2' / `g45t3' / `g45t4' s for 20 calls at R = `g45R')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
if `g45ok1' != 1 | `g45ok2' != 1 | `g45ok3' != 1 | `g45ok4' != 1 | `g45str' != 1 {
  di as error "G45c FAIL: a write-back lost its row or its stripes (rows `g45ok1'`g45ok2'`g45ok3'`g45ok4', stripes `g45str')"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
}
matrix drop G45O G45A G45W G45JW G45J G45C_gt

if ${BOOTG_FAILS} == 0 di as result "BOOTGATES_ALL_PASS"
else di as error "BOOTGATES_FAILS=${BOOTG_FAILS}"
