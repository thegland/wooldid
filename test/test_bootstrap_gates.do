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
* the event-study plot save is refused without replace (602), the rc propagates (it was
* swallowed: rc 0 with a stale saved file), the caller's frame is restored, and the
* replay branch refuses until a wooldid call completes
use `panel', clear
foreach f in boot_gates_g20_ES_att.gph boot_gates_g20_ES.dta {
  cap erase "`f'"
}
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) makeplots saveplots(boot_gates_g20)
chk 20 0 _rc
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) makeplots saveplots(boot_gates_g20)
chk 20b 602 _rc
qui frame
if "`r(currentframe)'" != "default" {
  di as error "G20c FAIL: current frame is `r(currentframe)' after the plot failure"
  global BOOTG_FAILS = ${BOOTG_FAILS} + 1
  frame change default
}
cap qui wooldid
chk 20d 198 _rc
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
* under bootdisplay(wald)
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
if `p30' != 0 | strpos(`"`lg'"', "p = 0 above means p < 0.200 (2/10 completed replications") == 0 | strpos(`"`lg'"', "The rank-based p-values have resolution 2/10 = 0.200") == 0 {
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
cap erase boot_gates_g30.txt

if ${BOOTG_FAILS} == 0 di as result "BOOTGATES_ALL_PASS"
else di as error "BOOTGATES_FAILS=${BOOTG_FAILS}"
