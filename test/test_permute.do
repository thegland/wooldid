* Targeted suite for F6 permutation inference: validation gates (P6*), the engine's
* exact/sampled/strata inference (P1-P4), its e() surface and additivity (P7, P11),
* cautions (P8), replay (P9), and determinism (P10).
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

if $pf == 0 dis "PERMUTE_ALL_PASS"
else dis as error "PERMUTE_FAILS=$pf"
exit, clear
