* Targeted suite for F6 permutation inference. COMMIT-A SCOPE: option plumbing and
* validation gates only - the engine refuses loudly until it lands (P7 pins that
* refusal and will be rewritten when the engine commit arrives).
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

* ---- P6g: nothing to randomize ----
preserve
qui replace ttre = 5 if ttre < .
qui keep if ttre < .
cap wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50)
xrc P6g1 459 _rc
restore

* ---- P7 (commit-A pin, REWRITTEN when the engine lands): a fully valid spec passes
* every screen and hits the loud engine-not-yet-enabled refusal ----
qui log using perm_p7.txt, text replace name(pm)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(50)
local rc7 = _rc
qui log close pm
if `rc7' != 198 | strpos(fileread("perm_p7.txt"), "engine is not yet") == 0 {
  dis as error "P7 FAIL: rc=`rc7' or the engine-pending refusal is missing"
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
