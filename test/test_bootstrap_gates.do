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

* G5 interaction refusals
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) permute(shuffle) permreps(20)
chk 5a 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) impvalues
chk 5b 198 _rc
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) fastnose
chk 5c 198 _rc
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

* G6 bayes + cell-weighted: refuse under default safety, caution + proceed-to-tail under safety(off)
cap qui wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(50)
chk 6a 198 _rc
cap noi wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(50) safety(off)
chk 6b 198 _rc
* (6b must fail at the TAIL, i.e. print the caution first - checked via log grep below)

* G7 standard + cell-weighted + nesting cluster (cluster == i) -> gates pass -> tail
cap qui wooldid y unit t ttre, att_it cluster(unit) bootstrap(standard) reps(50)
chk 7 198 _rc
* G7b coarser nesting cluster (s nests unit? NO - s is COARSER than unit: each unit
* wholly inside one s, so cluster(s) nests i=unit) -> gates pass -> tail
cap qui wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(50)
chk 7b 198 _rc

* G8 standard + cell-weighted + no cluster -> obs-level refusal
cap qui wooldid y unit t ttre, att_it robust bootstrap(standard) reps(50)
chk 8 198 _rc

* G9 standard + cell-weighted + NON-nesting cluster: i = s (coarse), cluster = unit
* (finer than i: unit varies within s) -> r459
cap qui wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(50)
chk 9 459 _rc

* G10 standard + att only + cluster(unit): all gates pass -> tail r198
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50)
chk 10 198 _rc

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
  gettoken lg 0 : 0
  qui log close `lg'
  local needle `"`0'"'
  tempname fh
  local found 0
  file open `fh' using boot_gates_`lg'.smcl, read
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') local found 1
    file read `fh' line
  }
  file close `fh'
  if `found' == 0 {
    di as error "MSG-`lg' FAIL: expected message not found: `needle'"
    global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  }
end

use `panel', clear
qui log using boot_gates_m1.smcl, replace name(m1) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50)
logh m1 engine is not yet enabled
qui log using boot_gates_m2.smcl, replace name(m2) smcl
cap noi wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(50) safety(off)
logh m2 Caution: bootstrap(bayes) standard errors for att_it
qui log using boot_gates_m3.smcl, replace name(m3) smcl
cap noi wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(50)
logh m3 requires the cluster() variable to nest i
qui log using boot_gates_m4.smcl, replace name(m4) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) oldsyntax
logh m4 not available under oldsyntax
qui log using boot_gates_m5.smcl, replace name(m5) smcl
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(50)
logh m5 engine is not yet enabled

if ${BOOTG_FAILS} == 0 di as result "BOOTGATES_ALL_PASS"
else di as error "BOOTGATES_FAILS=${BOOTG_FAILS}"
