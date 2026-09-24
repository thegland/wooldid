* Targeted suite for F6 permutation inference: validation gates (P6*), the engine's
* exact/sampled/strata inference (P1-P4), its e() surface and additivity (P7, P11),
* cautions (P8), replay (P9), and determinism (P10); onsets after the last period (P13);
* failing replications (P14); a failed loop keeps the caller's e() (P15).
* Run from test/:  "<stata>" -e do test_permute.do   -> expect PERMUTE_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 13579
set sortseed 13579
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* helper: expect a given rc
cap program drop xrc
program define xrc
  args label want got
  if `got' != `want' {
    dis as error "`label' FAIL: rc=`got' (expected `want')"
    global pf = $pf + 1
  }
end
global pf 0

use "`DATA'/dta_binary", clear

* ---- P6a: perm suboptions without permute() are refused ----
cap wooldid y i t ttre, att cluster(i) permreps(50)
xrc P6a1 198 _rc
cap wooldid y i t ttre, att cluster(i) permexact
xrc P6a2 198 _rc
cap wooldid y i t ttre, att cluster(i) permtstat
xrc P6a3 198 _rc

* ---- P6b: incompatible inference modes ----
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50) impvalues
xrc P6b1 198 _rc
cap wooldid y i t ttre, att cluster(i) robust permute(shuffle) permreps(50) unconditionalse
xrc P6b2 198 _rc

* ---- P6c: reps/exact arithmetic ----
cap wooldid y i t ttre, att cluster(i) permute(shuffle)
xrc P6c1 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50) permexact
xrc P6c2 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(50)
xrc P6c3 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50) permprobs((ttre 5 .) (1 .5 .5))
xrc P6c4 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(nonsense) permreps(50)
xrc P6c5 198 _rc

* ---- P6d: randomization unit must be a single cluster variable ----
cap wooldid y i t ttre, att permute(shuffle) permreps(50)
xrc P6d1 198 _rc
cap wooldid y i t ttre, att cluster(i t) permute(shuffle) permreps(50)
xrc P6d2 198 _rc

* ---- P6e: cluster nesting in onset groups (mixing clusters are refused) ----
qui gen badclus = mod(i, 3)
cap wooldid y i t ttre, att cluster(badclus) permute(shuffle) permreps(50)
xrc P6e1 459 _rc
qui drop badclus

* ---- P6f: strata must contain whole clusters ----
qui gen badstr = mod(t, 2)
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50) permstrata(badstr)
xrc P6f1 459 _rc
qui drop badstr

* ---- P6k: probabilities needs a numeric cluster; P6l: the permreps floor ----
preserve
qui gen str8 sclus = "c" + string(i)
cap wooldid y i t ttre, att cluster(sclus) permute(probabilities) permreps(50) permprobs((ttre 5 .) (1 .5 .5))
xrc P6k1 109 _rc
restore
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(10)
xrc P6l1 198 _rc

* ---- P6j: extended missing codes in the timing variable are refused ----
preserve
qui replace ttre = .a if ttre >= . & mod(i, 2) == 0
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50)
xrc P6j1 459 _rc
restore

* ---- P6g: nothing to randomize ----
preserve
qui replace ttre = 5 if ttre < .
qui keep if ttre < .
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50)
xrc P6g1 459 _rc
restore

* ---- P6i: methods/options awaiting their engine commits refuse loudly ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) permute(probabilities) permexact permprobs((ttre 5 .) (1 .5 .5))
xrc P6i1 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50) permtstat fastnose
xrc P6i2 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50) permtci fastnose
xrc P6i3 198 _rc

* ---- helper panel: 6 clusters, 3 treated at t=4, strong effect ----
cap program drop mkperm
program define mkperm
  args nclus ntreat eff
  clear
  set seed 4242
  set sortseed 4242
  qui set obs `nclus'
  qui gen i = _n
  qui gen ttre = cond(_n <= `ntreat', 4, .)
  qui expand 8
  qui bysort i: gen t = _n
  qui gen y = `eff'*(t >= ttre & ttre < .) + 0.1*i + rnormal()
end

* ---- P2: exact enumeration, M = C(6,3) = 20; strong effect -> p = 2/20 (the
* observed assignment and its mirror complement) ----
mkperm 6 3 3
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
if _rc != 0 {
  dis as error "P2 FAIL: rc=" _rc
  global pf = $pf + 1
}
else {
  if e(perm_npossible) != 20 | e(perm_completed) != 20 {
    dis as error "P2 FAIL: npossible=" e(perm_npossible) " completed=" e(perm_completed)
    global pf = $pf + 1
  }
  if abs(e(perm_p_att) - 0.10) > 1e-12 {
    dis as error "P2 FAIL: exact p=" e(perm_p_att) " (expected 2/20)"
    global pf = $pf + 1
  }
  if abs(e(perm_p_att)*20 - round(e(perm_p_att)*20)) > 1e-9 {
    dis as error "P2 FAIL: exact p not a multiple of 1/M"
    global pf = $pf + 1
  }
}

* ---- P1: placebo (no effect) -> exact p large; the p-value machinery has power both ways ----
mkperm 6 3 0
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
if _rc != 0 | e(perm_p_att) < 0.2 {
  dis as error "P1 FAIL: rc=" _rc " placebo exact p=" e(perm_p_att)
  global pf = $pf + 1
}

* ---- P3: sampled p is Monte-Carlo consistent with the exact p on the same panel ----
mkperm 6 3 3
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
local pex = e(perm_p_att)
set seed 31415
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(200)
if abs(e(perm_p_att) - `pex') > 0.08 {
  dis as error "P3 FAIL: sampled p=" e(perm_p_att) " vs exact `pex'"
  global pf = $pf + 1
}

* ---- P4: strata - M is the within-strata product; every level keeps its per-stratum count ----
clear
set seed 55
set sortseed 55
qui set obs 8
qui gen i = _n
qui gen st = cond(_n <= 4, 1, 2)
qui gen ttre = cond(mod(_n,2) == 1, 5, .)
qui expand 8
qui bysort i: gen t = _n
qui gen y = 2*(t >= ttre & ttre < .) + 0.3*st + rnormal()
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permstrata(st)
if _rc != 0 | e(perm_npossible) != 36 | e(perm_completed) != 36 {
  dis as error "P4 FAIL: rc=" _rc " npossible=" e(perm_npossible) " completed=" e(perm_completed) " (expected 6x6)"
  global pf = $pf + 1
}

* ---- P5: permute(probabilities) - grammar, draws, counts, cautions ----
mkperm 6 3 3
set seed 31415
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(60) permprobs((ttre 4 .) (1 .5 .5) (2 .5 .5) (3 .5 .5) (4 .5 .5) (5 .5 .5) (6 .5 .5))
if _rc != 0 | e(perm_npossible) != 64 | !(e(perm_p_att) >= 0 & e(perm_p_att) <= 1) {
  dis as error "P5a FAIL: rc=" _rc " npossible=" e(perm_npossible) " p=" e(perm_p_att)
  global pf = $pf + 1
}
mkperm 6 3 3
cap wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(20) permprobs((ttre 4 .) (1 .6 .5) (2 .5 .5) (3 .5 .5) (4 .5 .5) (5 .5 .5) (6 .5 .5))
xrc P5b 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(20) permprobs((ttre 4 .) (1 .5 .5) (2 .5 .5) (3 .5 .5) (4 .5 .5) (5 .5 .5) (99 .5 .5))
xrc P5c 198 _rc
cap wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(20) permprobs((wrongvar 4 .) (1 .5 .5) (2 .5 .5) (3 .5 .5) (4 .5 .5) (5 .5 .5) (6 .5 .5))
xrc P5d 198 _rc
qui log using perm_p5.txt, text replace name(pm)
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(30) permprobs((ttre 4 .) (1 0 1) (2 .5 .5) (3 .5 .5) (4 .5 .5) (5 .5 .5) (6 .5 .5))
qui log close pm
if _rc != 0 | strpos(fileread("perm_p5.txt"), "zero probability under the specified design") == 0 {
  dis as error "P5e FAIL: rc=" _rc " or zero-probability caution missing"
  global pf = $pf + 1
}

* ---- P7: the full e() surface on a sampled run ----
mkperm 10 5 3
set seed 2718
cap noi wooldid y i t ttre, att att_it cluster(i) permute(shuffle) permreps(30)
if _rc != 0 {
  dis as error "P7 FAIL: rc=" _rc
  global pf = $pf + 1
}
else {
  local bad 0
  if "`e(perm_method)'" != "shuffle" local ++bad
  if e(perm_reps) != 30 local ++bad
  if e(perm_completed) + e(perm_failed) != 30 local ++bad
  if e(perm_npossible) != 252 local ++bad
  foreach ee in att att_it {
    if !(e(perm_p_`ee') >= 0 & e(perm_p_`ee') <= 1) local ++bad
    if !(e(perm_lb_`ee') <= e(perm_ub_`ee')) local ++bad
    if !(e(perm_sd_`ee') > 0 & e(perm_sd_`ee') < .) local ++bad
  }
  cap confirm matrix e(permdist)
  if _rc != 0 local ++bad
  else if rowsof(e(permdist)) != e(perm_completed) | colsof(e(permdist)) != 2 local ++bad
  if `bad' != 0 {
    dis as error "P7 FAIL: `bad' e() surface check(s)"
    global pf = $pf + 1
  }
}

* ---- P8: cautions - enumeration size and reps-vs-possible ----
qui log using perm_p8.txt, text replace name(pm)
mkperm 6 3 3
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50)
local rc8 = _rc
qui log close pm
if `rc8' != 0 | strpos(fileread("perm_p8.txt"), "exceeds 10% of the") == 0 {
  dis as error "P8 FAIL: reps>10%-of-possible caution missing"
  global pf = $pf + 1
}

* ---- P9: replay reprints the permutation block ----
mkperm 6 3 3
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
qui log using perm_p9.txt, text replace name(pm)
cap noi wooldid
qui log close pm
if _rc != 0 | strpos(fileread("perm_p9.txt"), "Randomization Inference") == 0 {
  dis as error "P9 FAIL: rc=" _rc " or replay block missing"
  global pf = $pf + 1
}

* ---- P9b: the replay is line-identical to the live run from the banner onward ----
* (change audit round 2, U10 replay-extra-blank-line-engines / display-test-gaps) test_fastnose
* T5's method on a permute() run. The replay branch printed a blank line before the part(perm)
* footer that the live run does not. No permsave() and below the 10,000-draw cap: the lines
* those print below the footer are, by design, not replayed
mkperm 6 3 3
qui log using perm_p9b_live.txt, text replace name(pm)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) rseed(99) dots(0)
qui log close pm
local rc9b1 = _rc
qui log using perm_p9b_rep.txt, text replace name(pm)
cap noi wooldid
qui log close pm
local rc9b2 = _rc
local p9bsame 0
local p9bn 0
mata: a_ = cat("perm_p9b_live.txt") \ ""; b_ = cat("perm_p9b_rep.txt") \ ""
mata: ia_ = selectindex(strpos(a_, "Wooldid Estimation") :== 1); ib_ = selectindex(strpos(b_, "Wooldid Estimation") :== 1)
mata: hb_ = (rows(ia_) > 0) & (rows(ib_) > 0)
mata: sa_ = a_[|(ia_ \ 1)[1] \ rows(a_)|]; sb_ = b_[|(ib_ \ 1)[1] \ rows(b_)|]
mata: st_local("p9bsame", strofreal(hb_ & (rows(sa_) > 20) & (sa_ == sb_))); st_local("p9bn", strofreal(rows(sa_)))
if `rc9b1' != 0 | `rc9b2' != 0 | `p9bsame' != 1 {
  dis as error "P9b FAIL: rc est=`rc9b1' rep=`rc9b2', or the replay differs from the live run from the banner onward (`p9bn' live lines)"
  global pf = $pf + 1
}
else dis "P9b info: replay line-identical from the banner onward (`p9bn' lines)"
cap erase perm_p9b_live.txt
cap erase perm_p9b_rep.txt

* ---- P10: determinism under a set seed ----
mkperm 6 3 3
set seed 777
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(25)
* held in a SCALAR: a local's decimal round-trip re-parses one ulp off (the equal-tailed
* p's k/13-type fractions do not survive it), which is not a determinism failure
tempname pd1
scalar `pd1' = e(perm_p_att)
set seed 777
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(25)
if `pd1' != e(perm_p_att) {
  dis as error "P10 FAIL: same seed gave different p"
  global pf = $pf + 1
}

* ---- P12: permtstat - t-ratio randomization; exact mirror symmetry pins p = 2/20 ----
mkperm 6 3 3
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permtstat
if _rc != 0 {
  dis as error "P12 FAIL: rc=" _rc
  global pf = $pf + 1
}
else {
  if "`e(perm_stat)'" != "t" {
    dis as error "P12 FAIL: perm_stat=[`e(perm_stat)']"
    global pf = $pf + 1
  }
  * equal-tailed convention (2026-09-05): the observed assignment is a member of its own
  * reference set, so an exact two-sided p has a floor of 2/M in both modes; require p a
  * positive multiple of 1/M at or above 2/20
  if !(e(perm_p_att) >= 0.10 - 1e-12 & abs(e(perm_p_att)*20 - round(e(perm_p_att)*20)) < 1e-9) {
    dis as error "P12 FAIL: exact t-stat p=" e(perm_p_att)
    global pf = $pf + 1
  }
}

* ---- P12b (audit 2026-09-13): the exact permtstat p RECONSTRUCTED from the posted draws on a
* placebo, where the observed t sits inside its 20-member reference set (P12 above only
* pins the 2/20 floor; no other leg reconstructs a permtstat p in exact mode) ----
mkperm 6 3 0
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permtstat
if _rc != 0 {
  dis as error "P12b FAIL: rc=" _rc
  global pf = $pf + 1
}
else {
  mata:
    A_ = st_matrix("e(permdist_b)"); S_ = st_matrix("e(permdist_se)")
    k_ = selectindex(st_matrixcolstripe("e(permdist_b)")[., 2] :== "main_att")
    cb_ = st_matrix("e(perm_comparator_b)")[1, k_]; cs_ = st_matrix("e(perm_comparator_se)")[1, k_]
    ok_ = (A_[., k_] :< .) :& (S_[., k_] :< .) :& (S_[., k_] :> 0)
    tt_ = select(A_[., k_] :/ S_[., k_], ok_)
    n_ = rows(tt_)
    to_ = cb_/cs_
    lo_ = sum(tt_ :<= to_); hi_ = sum(tt_ :>= to_)
    // exact mode: plain shares (the observed assignment is a member of its own reference
    // set), equal-tailed - the same operation order as the engine (test_permute3 wdp3_reconx)
    p_ = min((1, 2*min((lo_/n_, hi_/n_))))
    st_numscalar("__p12bn", n_); st_numscalar("__p12btail", min((lo_, hi_))); st_numscalar("__p12bp", p_)
  end
  if "`e(perm_stat)'" != "t" | e(perm_exact) != 1 | __p12bn != 20 | e(perm_npossible) != 20 | e(perm_completed) != 20 {
    dis as error "P12b FAIL: stat=[`e(perm_stat)'] exact=" e(perm_exact) " complete t draws=" __p12bn " npossible=" e(perm_npossible) " completed=" e(perm_completed) " (want t/1/20/20/20)"
    global pf = $pf + 1
  }
  else if __p12btail < 2 {
    dis as error "P12b FAIL (VACUOUS): smaller tail = " __p12btail " of 20 - the placebo's observed t is at the extreme, so this leg would only re-derive the floor"
    global pf = $pf + 1
  }
  else if __p12bp != e(perm_p_att) {
    dis as error "P12b FAIL: exact permtstat p reconstructed from e(permdist_b)/e(permdist_se) = " %21x __p12bp " vs e(perm_p_att) = " %21x e(perm_p_att)
    global pf = $pf + 1
  }
  else dis as result "P12b PASS: exact permtstat p " e(perm_p_att) " reconstructs bitwise (smaller tail " __p12btail " of 20)"
}

* ---- P13: permtci implies permtstat; the inverted-t CI exists and covers the estimate ----
mkperm 6 3 3
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permtci
if _rc != 0 {
  dis as error "P13 FAIL: rc=" _rc
  global pf = $pf + 1
}
else {
  local bad 0
  if "`e(perm_stat)'" != "t" local ++bad
  if !(e(perm_tlb_att) < . & e(perm_tub_att) < .) local ++bad
  if !(e(perm_tlb_att) <= _b[main_att] & _b[main_att] <= e(perm_tub_att)) local ++bad
  if !(e(perm_tlb_att) < e(perm_tub_att)) local ++bad
  if `bad' != 0 {
    dis as error "P13 FAIL: `bad' inverted-t CI check(s): [" e(perm_tlb_att) ", " e(perm_tub_att) "] vs " _b[main_att]
    global pf = $pf + 1
  }
}

* ---- P11: additivity under the 2026-09-05 contract - the point estimates are those of
* the fast-mode run (the permute() default), e(V) is absent by default, and under
* permanalyticse the observed run equals the full-inference run bitwise ----
mkperm 6 3 3
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fastnose
tempname bb
matrix `bb' = e(b)
local nn = e(N)
mkperm 6 3 3
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
tempname bf vf
matrix `bf' = e(b)
matrix `vf' = e(V)
mkperm 6 3 3
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20)
cap confirm matrix e(V)
if mreldif(`bb', e(b)) != 0 | e(N) != `nn' | _rc == 0 {
  dis as error "P11 FAIL: fast-mode point estimates moved under permute (b " mreldif(`bb', e(b)) ") or e(V) posted"
  global pf = $pf + 1
}
mkperm 6 3 3
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20) permanalyticse
if mreldif(`bf', e(b)) != 0 | mreldif(`vf', e(V)) != 0 {
  dis as error "P11b FAIL: observed results moved under permute + permanalyticse (b " mreldif(`bf', e(b)) ", V " mreldif(`vf', e(V)) ")"
  global pf = $pf + 1
}
if "`e(cmd)'" != "wooldid" {
  dis as error "P11 FAIL: e(cmd)=[`e(cmd)']"
  global pf = $pf + 1
}

* ---- P6h: gates leave no residue - a normal run afterwards is clean ----
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "P6h FAIL: clean run after gate refusals rc=" _rc
  global pf = $pf + 1
}

* ---- P12: the permexact enumeration ceiling is the Stata matrix row limit (2026-09-22) ----
* Ruled by the author: the cap is no longer a round number chosen for the user's patience but
* 65,534, the largest number of rows a Stata matrix holds, because the replication accumulators
* (`res' and `accname' in wooldid_permrun) are Stata matrices with one row per replication. One
* more and the run dies r(915) after paying for the whole enumeration.
*
* These legs call wooldid_permx DIRECTLY. That is deliberate: the band this ruling opened runs
* from 10,001 to 65,534 assignments, and going through the command would mean actually running
* 40,320 replications. The routine reads one row per cluster from a frame, so a synthetic frame
* exercises the cap in milliseconds while still enumerating for real.
*
* capP/capL are passed as unset TEMPNAMES, not as "", because the routine's probabilities
* branch reads `pname != "" & rows(st_matrix(pname)) > 0` and Mata does not short-circuit `&`:
* an empty pname reaches st_matrix("") and dies r(3300). Production always passes tempnames, so
* that arm is unreachable there, but a caller must not hand it an empty string.
tempname capfr capA capP capL
frame create `capfr'
frame `capfr' {
  qui set obs 8
  qui gen double ons = _n
  qui gen double str = 1
}
* 8 clusters, 8 distinct onsets = 8! = 40,320 assignments: ABOVE the old 10,000 cap and below
* the new one, so this is the case the ruling admits. It must enumerate, not refuse.
local wdp_err
local wdp_R
mata: wooldid_permx("`capfr'", "ons", "str", 0, 1, "`capA'", "`capP'", "`capL'")
if "`wdp_err'" != "" {
  dis as error "P12a FAIL: 40,320 assignments refused: `wdp_err'"
  global pf = $pf + 1
}
else if "`wdp_R'" != "40320" {
  dis as error "P12a FAIL: wdp_R=`wdp_R' (expected 40320)"
  global pf = $pf + 1
}
else if colsof(`capA') != 40320 | rowsof(`capA') != 8 {
  dis as error "P12a FAIL: assignment matrix is " rowsof(`capA') " x " colsof(`capA') " (expected 8 x 40320)"
  global pf = $pf + 1
}
* 12 clusters, 12 distinct onsets = 12! assignments, far above the ceiling: refused, and the
* message must name the NEW ceiling rather than the retired 10,000
frame `capfr' {
  qui set obs 12
  qui replace ons = _n
  qui replace str = 1
}
local wdp_err
mata: wooldid_permx("`capfr'", "ons", "str", 0, 1, "`capA'", "`capP'", "`capL'")
if strpos("`wdp_err'", "65534") == 0 {
  dis as error "P12b FAIL: over-ceiling refusal did not name 65534: [`wdp_err']"
  global pf = $pf + 1
}
frame drop `capfr'

* ---- P12c: end to end, an over-ceiling permexact run is refused by the command ----
* cheap because the count is computed and tested before anything is enumerated
clear
qui set obs 12
qui gen i = _n
* 8 treated with distinct onsets and 4 never-treated: 12!/4! = 19,958,400 assignments, far
* above the ceiling, and the never-treated give the design a control group
qui gen double ttre = cond(_n <= 8, 3 + _n, .)
qui expand 18
qui bysort i: gen t = _n
qui gen y = 2*(t >= ttre & ttre < .) + 0.1*i + rnormal()
cap noi qui wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
if _rc != 198 {
  dis as error "P12c FAIL: over-ceiling permexact rc=" _rc " (expected 198)"
  global pf = $pf + 1
}

* ---- P13: an onset after the last period is never-treated in the permutation design ----
* The worker makes an onset after max(t) + oostreatedcontrols(), net of anticipation(),
* never-treated for the estimation, and the permutation design does the same: the cluster
* map, the level accounting and nesting screen, the ever-treated count, the freedom check, the
* permprobs() level set and the assignment count all read that design onset, and a caution
* counts the clusters recoded. Fixture: 12 units x t = 1..10, units 1-2 first treated at 6,
* units 3-6 never treated, units 7-12 at the onset(s) passed. The outcome is drawn before the
* onsets of units 7-12 are set, so a late fixture and its hand-recoded twin (onset set to .)
* differ only in ttre.
cap program drop mklate
program define mklate
  args late
  clear
  set seed 1313
  set sortseed 1313
  qui set obs 12
  qui gen i = _n
  qui gen double ttre = cond(_n <= 2, 6, .)
  qui expand 10
  qui bysort i: gen t = _n
  qui gen y = 1.5*(t >= ttre & ttre < .) + 0.1*i + 0.05*t + rnormal()
  if "`late'" == "distinct" qui replace ttre = 4 + i if i >= 7
  else if "`late'" == "oos" qui replace ttre = cond(i <= 9, 12, 13) if i >= 7
  else if "`late'" != "" qui replace ttre = `late' if i >= 7
end

* P13a: six late adopters at onset 12 under shuffle. The design has 2 treated clusters and 10
* never treated, so C(12,2) = 66 assignments (13,860 when the late label counted as a level);
* the recode caution, the few-treated caution and the 10% caution all print. The p-values are
* those of the hand-recoded twin, and bitwise those of the build before the design onset
* existed (the %21x pin), because shuffle permutes positions rather than labels
tempname p13p p13b p13c
mklate 12
qui log using perm_p13a.txt, text replace name(p13)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20) rseed(1313) dots(0)
local rc13 = _rc
qui log close p13
xrc P13a1 0 `rc13'
if `rc13' == 0 {
  matrix `p13p' = e(perm_p)
  matrix `p13b' = e(b)
  matrix `p13c' = e(perm_ci)
  if e(perm_npossible) != 66 {
    dis as error "P13a2 FAIL: e(perm_npossible)=" e(perm_npossible) " (expected 66)"
    global pf = $pf + 1
  }
  if strpos(fileread("perm_p13a.txt"), "Caution: 6 cluster(s) with onsets after the last period are permuted as never treated.") == 0 {
    dis as error "P13a3 FAIL: the recode caution did not print"
    global pf = $pf + 1
  }
  if strpos(fileread("perm_p13a.txt"), "Caution: only 2 of the resampling clusters are ever treated") == 0 {
    dis as error "P13a4 FAIL: the few-treated caution did not count the design's 2 treated clusters"
    global pf = $pf + 1
  }
  if strpos(fileread("perm_p13a.txt"), "exceeds 10% of the 66 possible assignments") == 0 {
    dis as error "P13a5 FAIL: the 10% caution was not judged on the design's 66 assignments"
    global pf = $pf + 1
  }
  local p13hex
  forval k = 1/`=colsof(`p13p')' {
    local p13h : di %21x el(`p13p', 1, `k')
    local p13hex `p13hex' `p13h'
  }
  dis "P13a p hex: `p13hex'"
  if "`p13hex'" != "+1.8618618618618X-004 +1.8618618618618X-001 +1.e79e79e79e79eX-002 +1.e79e79e79e79eX-001 +1.5555555555555X-001 +1.e79e79e79e79eX-002" {
    dis as error "P13a6 FAIL: e(perm_p) moved against the pre-change build"
    global pf = $pf + 1
  }
}
mklate 12
qui replace ttre = . if ttre > 10
qui log using perm_p13a2.txt, text replace name(p13)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20) rseed(1313) dots(0)
local rc13b = _rc
qui log close p13
xrc P13a7 0 `rc13b'
if `rc13' == 0 & `rc13b' == 0 {
  if mreldif(`p13p', e(perm_p)) != 0 | mreldif(`p13c', e(perm_ci)) != 0 | mreldif(`p13b', e(b)) != 0 | e(perm_npossible) != 66 {
    dis as error "P13a8 FAIL: the late run differs from its hand-recoded twin (p " mreldif(`p13p', e(perm_p)) ", ci " mreldif(`p13c', e(perm_ci)) ", b " mreldif(`p13b', e(b)) ", M " e(perm_npossible) ")"
    global pf = $pf + 1
  }
  if strpos(fileread("perm_p13a2.txt"), "with onsets after the last period") != 0 {
    dis as error "P13a9 FAIL: the recode caution printed on data with no late onset"
    global pf = $pf + 1
  }
}

* P13b: the six late adopters at distinct onsets 11..16. As labels they made 12!/(2!4!) =
* 9,979,200 assignments and permexact was refused above its 65,534 ceiling; as the
* never-treated they are, the design has P13a's 66 and enumerates them, with the exact
* p-values of the hand-recoded twin
mklate distinct
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact dots(0)
local rc13 = _rc
xrc P13b1 0 `rc13'
if `rc13' == 0 {
  matrix `p13p' = e(perm_p)
  if e(perm_reps) != 66 | e(perm_completed) != 66 {
    dis as error "P13b2 FAIL: permexact enumerated " e(perm_reps) " assignments, completed " e(perm_completed) " (expected 66)"
    global pf = $pf + 1
  }
  mklate distinct
  qui replace ttre = . if ttre > 10
  cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact dots(0)
  local rc13b = _rc
  xrc P13b3 0 `rc13b'
  if `rc13b' == 0 {
    if mreldif(`p13p', e(perm_p)) != 0 {
      dis as error "P13b4 FAIL: exact p-values differ from the hand-recoded twin by " mreldif(`p13p', e(perm_p))
      global pf = $pf + 1
    }
  }
}

* P13c: anticipation() moves the window with the onset. Onset 11 under anticipation(1) is an
* effective onset of 10, inside the window, so the six units stay a treated level: no recode
* caution and 12!/(2!6!4!) = 13,860 assignments
mklate 11
qui log using perm_p13c.txt, text replace name(p13)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) anticipation(1) rseed(1313) dots(0)
local rc13 = _rc
qui log close p13
xrc P13c1 0 `rc13'
if `rc13' == 0 {
  if e(perm_npossible) != 13860 {
    dis as error "P13c2 FAIL: e(perm_npossible)=" e(perm_npossible) " (expected 13860)"
    global pf = $pf + 1
  }
  if strpos(fileread("perm_p13c.txt"), "with onsets after the last period") != 0 {
    dis as error "P13c3 FAIL: an onset inside the anticipation window was recoded"
    global pf = $pf + 1
  }
}

* P13d: oostreatedcontrols(2) keeps onset 12 (units 7-9) a treated level and recodes onset 13
* (units 10-12): 12!/(2!3!7!) = 7,920 assignments where the labels gave 12!/(2!4!3!3!) =
* 277,200, and the caution counts 3 clusters
mklate oos
qui log using perm_p13d.txt, text replace name(p13)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) oostreatedcontrols(2) rseed(1313) dots(0)
local rc13 = _rc
qui log close p13
xrc P13d1 0 `rc13'
if `rc13' == 0 {
  if e(perm_npossible) != 7920 {
    dis as error "P13d2 FAIL: e(perm_npossible)=" e(perm_npossible) " (expected 7920)"
    global pf = $pf + 1
  }
  if strpos(fileread("perm_p13d.txt"), "Caution: 3 cluster(s) with onsets after the last period are permuted as never treated.") == 0 {
    dis as error "P13d3 FAIL: the recode caution did not count the 3 clusters past oostreatedcontrols(2)"
    global pf = $pf + 1
  }
}

* P13e: the permprobs() level set is the design's. With the late adopters at 12 the levels are
* 6 and never-treated, so (ttre 6 .) validates and (ttre 6 12 .) is refused r(198)
mklate 12
local PP2 (ttre 6 .)
local PP3 (ttre 6 12 .)
forval c = 1/12 {
  local PP2 `PP2' (`c' .5 .5)
  local PP3 `PP3' (`c' .4 .3 .3)
}
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permprobs(`PP2') permreps(30) rseed(1313) dots(0)
xrc P13e1 0 `=_rc'
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permprobs(`PP3') permreps(30) rseed(1313) dots(0)
xrc P13e2 198 `=_rc'
foreach f in perm_p13a perm_p13a2 perm_p13c perm_p13d {
  cap erase `f'.txt
}

* ---- P14: replications that fail to estimate (change audit round 2, U07 perm-test-gaps) ----
* Every e(perm_failed) assertion elsewhere is == 0, so the badrep path was never driven. Six
* clusters, one treated at t 4, and clusters 5 and 6 carry no outcome: they stay in the design,
* so a draw that hands them the onset leaves the replication with no treated observation and it
* fails. The run must complete (30 of 40 succeed here, above the max(20, half) floor), count the
* failures in e(perm_failed), print one x per failure on the dots line, and print the
* above-10% caution once
clear
set seed 4242
set sortseed 4242
qui set obs 6
qui gen i = _n
qui gen ttre = cond(_n <= 1, 4, .)
qui expand 8
qui bysort i: gen t = _n
qui gen y = 3*(t >= ttre & ttre < .) + 0.1*i + rnormal()
qui replace y = . if i >= 5
qui log using perm_p14.txt, text replace name(pm)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(1414)
local rc14 = _rc
qui log close pm
local p14f = e(perm_failed)
local p14c = e(perm_completed)
mata: l_ = cat("perm_p14.txt"); k_ = selectindex(strpos(l_, "Permutation replications (") :== 1)
mata: d_ = (rows(k_) ? l_[k_[1]] : ""); st_local("p14x", strofreal(strlen(d_) - strlen(subinstr(d_, "x", ""))))
local lg14 = fileread("perm_p14.txt")
local p14n = (length(`"`lg14'"') - length(subinstr(`"`lg14'"', "replications failed to estimate and were excluded", "", .))) / length("replications failed to estimate and were excluded")
dis "P14 info: rc=`rc14' completed=`p14c' failed=`p14f' x-marks=`p14x' caution printed `p14n' time(s)"
if `rc14' != 0 | !(`p14f' > 4 & `p14f' < .) | `p14c' + `p14f' != 40 | `p14x' != `p14f' | `p14n' != 1 {
  dis as error "P14 FAIL: failed replications are not counted, marked or cautioned as documented"
  global pf = $pf + 1
}
cap erase perm_p14.txt

* ---- P15: a replication loop that fails keeps the caller's e() (the author's round-2 ruling
* 2026-09-23, "keep the promise": a call that fails, at any stage, leaves e() as it was before
* the call). P14's design with the outcome missing on clusters 4-8 of 8: a draw that hands the
* one onset to any of them fails, so too few replications complete for the max(20, half) floor
* and the loop errors. The engine used to leave the failed call's own observed run in e(); e()
* must now be the plain run before it, and it must still replay (the loop failure leaves the
* replay marker clear) ----
clear
set seed 4243
set sortseed 4243
qui set obs 8
qui gen i = _n
qui gen ttre = cond(_n <= 1, 4, .)
qui expand 8
qui bysort i: gen t = _n
qui gen y = 3*(t >= ttre & ttre < .) + 0.1*i + rnormal()
qui replace y = . if i >= 4
cap noi wooldid y i t ttre, att cluster(i)
local rc15a = _rc
local cl15 `"`e(cmdline)'"'
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(1515) dots(0)
local rc15 = _rc
local cl15b `"`e(cmdline)'"'
cap qui wooldid
local rc15r = _rc
dis `"P15 info: plain rc=`rc15a' | permute rc=`rc15' e(cmdline)=[`cl15b'] | replay rc=`rc15r'"'
if `rc15a' != 0 | `rc15' != 459 | `"`cl15b'"' != `"`cl15'"' | `rc15r' != 0 {
  dis as error "P15 FAIL: after the failed replication loop e() is not the earlier run's, or it does not replay"
  global pf = $pf + 1
}

if $pf == 0 dis "PERMUTE_ALL_PASS"
else dis as error "PERMUTE_FAILS=$pf"
exit, clear
