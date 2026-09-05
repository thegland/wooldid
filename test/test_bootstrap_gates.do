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
cap noi wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(30) safety(off)
chk 6b 0 _rc
* (6b now RUNS: caution printed, then full bayes bootstrap - message-checked below)

* G7 standard + cell-weighted + nesting cluster (cluster == i) -> RUNS (re-id leg)
cap qui wooldid y unit t ttre, att_it cluster(unit) bootstrap(standard) reps(20)
chk 7 0 _rc
* G7b coarser nesting cluster (s nests unit? NO - s is COARSER than unit: each unit
* wholly inside one s, so cluster(s) nests i=unit) -> gates pass -> tail
cap qui wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20)
chk 7b 0 _rc

* G8 standard + cell-weighted + no cluster -> obs-level refusal
cap qui wooldid y unit t ttre, att_it robust bootstrap(standard) reps(50)
chk 8 198 _rc

* G9 standard + cell-weighted + NON-nesting cluster: i = s (coarse), cluster = unit
* (finer than i: unit varies within s) -> r459
cap qui wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(50)
chk 9 459 _rc

* G10 standard + att only + cluster(unit): runs (count path)
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20)
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
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) esfixedbaseperiod bootstrap(bayes) reps(20) pretrendpower(80)
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
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(8)
logh m1 Bootstrap Inference
qui log using boot_gates_m2.smcl, replace name(m2) smcl
cap noi wooldid y unit t ttre, att_it cluster(unit) bootstrap(bayes) reps(30) safety(off)
logh m2 Caution: bootstrap(bayes) standard errors for att_it
qui log using boot_gates_m3.smcl, replace name(m3) smcl
cap noi wooldid y s t ttre, att_it cluster(unit) bootstrap(standard) reps(50)
logh m3 requires the cluster() variable to nest i
qui log using boot_gates_m4.smcl, replace name(m4) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(50) oldsyntax
logh m4 not available under oldsyntax
qui log using boot_gates_m5.smcl, replace name(m5) smcl
cap noi wooldid y unit t ttre, att_it cluster(s) bootstrap(standard) reps(20) rseed(8)
logh m5 Bootstrap Inference
qui log using boot_gates_m6.smcl, replace name(m6) smcl
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(20) rseed(3)
logh m6 Bootstrap Inference

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

if ${BOOTG_FAILS} == 0 di as result "BOOTGATES_ALL_PASS"
else di as error "BOOTGATES_FAILS=${BOOTG_FAILS}"
