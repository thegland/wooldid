* wooldid regression battery. Run from a fresh subdirectory of test/ (see test/README.md):
*   mkdir test/run1 && cd test/run1 && "<stata>" -e do ../runbattery.do
* Loads the CURRENT repo wooldid.ado, runs 48 configs, dumps all e() content
* (exact %21x hex) + rc + post-run varlist + saved-file names to <cfg>.csv in cwd.
clear all
set more off
set maxvar 20000
set graphics off
set seed 98765
set sortseed 98765

local REPO "../.."
local DATA "../data"

do "`REPO'/wooldidinstall/wooldid.ado"

cap program drop dumpres
program define dumpres
  gettoken cfg 0 : 0
  gettoken rc  0 : 0
  * every remaining token is a file prefix to pin (see the files| note below)
  local svpre `0'
  tempname fh
  file open `fh' using "`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    if _rc != 0 continue
    local R = rowsof(`M')
    local C = colsof(`M')
    local rn : rowfullnames `M'
    local cn : colfullnames `M'
    file write `fh' `"matnames|`m'|`R'x`C'|`rn'|`cn'"' _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  * files| (2026-09-13, suite audit F-8; mirrors mkrun.py's line): every file in cwd whose name
  * starts with one of the prefixes given after cfg and rc. A config that writes files -
  * permsave(pfx), bootsave(pfx), saveplots(pfx) - names its prefixes on its dumpres line
  * (`dumpres perm_save `rc' pfx'); the sorted names are pinned and the files are erased so the
  * next config's glob is clean. Default prefix svp (mkrun's). compare.py reads a dump without
  * this line as an empty list, so the pre-2026-09-13 baseline files stay comparable as they are.
  if "`svpre'" == "" local svpre svp
  local svfiles
  foreach sp of local svpre {
    local sp1 : dir . files "`sp'*"
    local svfiles `svfiles' `sp1'
  }
  local svfiles : list uniq svfiles
  local svfiles : list sort svfiles
  file write `fh' `"files|`svfiles'"' _n
  foreach svf of local svfiles {
    cap erase `svf'
  }
  file close `fh'
end

* helper to reset state between configs
cap program drop freshstart
program define freshstart
  args dtafile
  ereturn clear
  use "`dtafile'", clear
end

* ---------------- configs ----------------

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats jointtests safety(off)
local rc = _rc
dumpres basic `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres allest_aw `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2) adversarial
local rc = _rc
dumpres esf_adv `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests safety(off)
local rc = _rc
dumpres sg_attit `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att_itime cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests safety(off)
local rc = _rc
dumpres sg_attitime `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) impvalues
local rc = _rc
dumpres impv `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid ycount i t ttre, att cluster(i) poisson
local rc = _rc
dumpres pois `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [pw=w], att cluster(i) esprelength(2) espostlength(2) summarystats
local rc = _rc
dumpres pw_ss `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [aw=w], att robust unconditionalse esprelength(2) espostlength(2)
local rc = _rc
dumpres aw_uncond `rc'

freshstart "`DATA'/dta_tiny"
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
local rc = _rc
dumpres tiny `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(-999)
local rc = _rc
dumpres esptrip `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) esfixedbaseperiod esprelength(3) espostlength(3) cluster(i)
local rc = _rc
dumpres cont_basic `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres cont_uz `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2) summarystats
local rc = _rc
dumpres cont_uz_ss `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre [aw=w], contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres cont_lat `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) summarystats
local rc = _rc
dumpres cont_ss `rc'

freshstart "`DATA'/dta_simple1"
cap noi wooldid y i t ttre, att cluster(i) subgroup(ttre)
local rc = _rc
dumpres simple1 `rc'

freshstart "`DATA'/dta_simple1"
cap noi wooldid y i t ttre, att cluster(i) subgroup(ttre) esprelength(3) espostlength(3) esfixedbaseperiod
local rc = _rc
dumpres simple_esf `rc'

* E1 coverage (added 2026-07-26, blessed pre-E1): the contreat SG-comparison block had
* zero battery coverage; sg_missif pins missing-sg + bare-flag-if through the SG block;
* sg_degpre pins the degenerate no-pre-period state (absent pre joints, post-only pwsgs)
freshstart "`DATA'/dta_cont"
qui gen sg = mod(i,2)+1
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests safety(off)
local rc = _rc
dumpres cont_sg `rc'

freshstart "`DATA'/dta_binary"
qui gen sgm = mod(i,3)
qui replace sgm = . if mod(i,11) == 0
cap noi wooldid y i t ttre if w > 0, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sgm) jointtests safety(off)
local rc = _rc
dumpres sg_missif `rc'

freshstart "`DATA'/dta_binary"
qui sum t, meanonly
local t2 = r(min) + 1
qui gen sgd2 = 1
qui replace sgd2 = 2 if mod(i,7) == 0
qui replace ttre = `t2' if mod(i,7) == 0
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sgd2) jointtests safety(off)
local rc = _rc
dumpres sg_degpre `rc'

* B16 coverage (added 2026-07-26, blessed pre-fix): the esf adversarial reference-row
* planting and the altmethod_adversarial pre-trend cut had no battery coverage; these
* four pin the affected cells (adv_altb: below-window reference regime with altmethod;
* adv_r2001: degenerate 1-row fit, rc 2001 at bless time; adv_gapa/adv_gapb: gap-year
* data where the reference-gap hunt aborts (rc 198) or plants at the wrong horizon)
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-3) esprelength(2) espostlength(2) adversarial altmethod_adversarial
local rc = _rc
dumpres adv_altb `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(1) espostlength(2) adversarial altmethod_adversarial
local rc = _rc
dumpres adv_r2001 `rc'

freshstart "`DATA'/dta_path_relyrgap"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(4) espostlength(4) adversarial
local rc = _rc
dumpres adv_gapa `rc'

freshstart "`DATA'/dta_path_relyrgap"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-4) esprelength(3) espostlength(3) adversarial
local rc = _rc
dumpres adv_gapb `rc'

* dupcol coverage (added 2026-07-26, blessed pre-fix): continuous + subgroups where the
* time-invariant dose makes the per-sg continuous pre-treatment margins come back
* all-omitted, tripping the sg rowname filter (silently shifted/duplicated colnames in
* e(b)/e(V)/e(wdidsgresults_*) on the pre-fix build: 45 columns, 42 unique names)
clear
set seed 447102
set obs 50
gen i = _n
gen ttre = cond(_n<=13, 5, cond(_n<=26, 8, cond(_n<=38, 10, .)))
expand 13
bysort i: gen t = _n
gen dose = 0.5 + mod(i,6) + 0.25*mod(i,7)/7
gen post = (t >= ttre & ttre < .)
gen y = 0.4*dose*post + 0.02*i + 0.08*t + sin(3*i + 2*t)
gen sg = mod(i,3)
drop post
ereturn clear
cap noi wooldid y i t ttre, contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg)
local rc = _rc
dumpres cont_sg_tidose `rc'

freshstart "`DATA'/dta_path_relyrgap"
gen dose2 = 0.5 + mod(i,5) + 0.2*mod(i,3)
gen sgr = mod(i,3)
cap noi wooldid y i t ttre, contreat(dose2) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sgr)
local rc = _rc
dumpres cont_sg_gap `rc'

* F3 coverage (added 2026-07-27, blessed on the feature build): pins e(gteffects) content
* and the off-state e() surface of a weighted two-estimand run
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) esprelength(2) espostlength(2) gteffects
local rc = _rc
dumpres gt_basic `rc'

* F6 coverage (added 2026-07-27, blessed on the feature build): pins the permutation
* e() surface; the seed is set INSIDE the config so the dump is independent of
* battery ordering (draws consume the session RNG stream)
freshstart "`DATA'/dta_binary"
set seed 246813
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20)
local rc = _rc
dumpres perm_basic `rc'

* a cell-weighted estimand under permute(): measured 2026-09-16 as having ZERO coverage - no
* att_it/att_i/att_itime config anywhere in the deck runs under permute(). This matters for the
* onset-rebuild feature specifically: permonsetcohorts rebuilds i inside every replication, so
* att_i is then the equal average over the REBUILT cohorts rather than the observed ones. Under
* permute(shuffle) the onset labels are permuted, so the level set is preserved by construction
* and the rebuilt cohort count matches the observed one; this row is the observed-side pin of
* that, taken with the option absent. Same seed convention as perm_basic (set INSIDE the config).
freshstart "`DATA'/dta_binary"
set seed 173205
cap noi wooldid y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(10)
local rc = _rc
dumpres perm_atti `rc'

* New-feature benchmark pins (added 2026-08-24, blessed in-commit; PREREG_NEWFEAT
* chunk): close the battery's coverage gaps for fastnose, pretrendpower,
* anticipation, permtstat/permtci, and permute(probabilities). Appended AFTER all
* prior configs so earlier dumps are bit-unaffected; permute configs set their own
* seed (perm_basic precedent).
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fastnose gteffects
local rc = _rc
dumpres fn_gt `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod pretrendpower(50 80)
local rc = _rc
dumpres ptp_two `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) anticipation(1)
local rc = _rc
dumpres antic_es `rc'

freshstart "`DATA'/dta_binary"
set seed 135791
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permtstat permtci
local rc = _rc
dumpres perm_tstat `rc'

* 10-cluster subset spanning all three onset blocks (2005: i 1-25, 2008: 26-50,
* never: 51-80); permprobs needs one row per in-sample cluster value
freshstart "`DATA'/dta_binary"
qui keep if inlist(i, 1, 2, 3, 26, 27, 28, 51, 52, 53, 54)
set seed 975310
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(30) permprobs((ttre 2005 2008 .) (1 .3 .3 .4) (2 .3 .3 .4) (3 .3 .3 .4) (26 .3 .3 .4) (27 .3 .3 .4) (28 .3 .3 .4) (51 .3 .3 .4) (52 .3 .3 .4) (53 .3 .3 .4) (54 .3 .3 .4))
local rc = _rc
dumpres perm_probs `rc'

* ---- 2026-09-11: bit-exact cover for the three paths the 504-config dual build cannot see.
* Every option on these three lines is NEW-ONLY, so the dual build can only ever report their
* e() objects as additions - additions are informational and are never compared, which is why
* the three rulings built on 2026-09-11 (F033 on the studentized cell p, the cfxplotnose
* carve-out under permute(), and the per-family band seed) moved nothing there. The battery is
* the tool that pins them: it compares this build against a committed baseline cell by cell.

* the ANALYTIC simultaneous band. Pins the sup-t critical value and both endpoints per family.
* The seed is now derived per family from the family's NAME (D-9, author ruling 2026-09-11),
* so these numbers are also the regression test for that derivation - and att_it is here
* precisely because it is the family whose draws used to depend on att having drawn first.
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) simultaneousbands
local rc = _rc
dumpres band_analytic `rc'

* boottstat WITH gteffects: the studentized CELL p-values, which now carry F033's add-one
* (2*(min+1)/(n+1)) like every other recentred rank-based p. Before that ruling this fixture's
* cell p-values included exact zeros; the baseline pins that they no longer can.
freshstart "`DATA'/dta_binary"
set seed 246810
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) gteffects bootstrap(standard) boottstat reps(20) rseed(24680) dots(0)
local rc = _rc
dumpres boot_ts_gt `rc'

* permute() + cfxplotnose + makecfxplots: the continuous-treatment plot cells. cfxplotnose used
* to truncate the cfx matrices to one column unconditionally, so under permute() every one of
* them was skipped by the cell rewiring and the cells stayed analytic under a permutation
* header with no message. They are permuted now, and this pins their p and shift interval.
freshstart "`DATA'/dta_cont"
set seed 864209
cap noi wooldid y i t ttre, att contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(1) espostlength(1) makecfxplots cfxplotnose permute(shuffle) permreps(20) rseed(86420) dots(0)
local rc = _rc
dumpres perm_cfxnose `rc'

* ---- 2026-09-13 (suite audit, permboot_tracking): bit-exact cover for the permute()/bootstrap()
* surfaces that no battery config, and no content-judged oldnew config, pins. Every option here is
* NEW-ONLY, so the dual build reports its e() objects as additions and never compares them; the cert
* suites either reconstruct these numbers from the posted draws (identities) or, on a few att-only
* unweighted fixtures, match the draw generators to an external oracle (DO1, O1/O2, G32,
* ri_crosscheck Sd); neither sees the weighted, poisson, stratified, re-identified, observation-level
* or subgroup/joint-wired paths bit-exactly. The battery pins those. RULES for this block: (1) every
* replication engine call sets its seed INSIDE the config AND passes rseed() - without rseed() the
* engine leaves the session stream alone and dumps c(rngstate) as e(perm_rngstate)/e(boot_rngstate),
* so the dump would pass a full battery run and fail any single-config rerun; (2) dots(0); (3) most
* configs run on an 18-cluster subset of dta_binary (6 units from each onset block, both subgroups
* in every block - the perm_probs precedent) so a replication costs ~1/5 of a full-panel one. These
* are regression pins of the code paths, NOT realistic-data evidence; the full-panel behaviour is
* the cert suites' business. Blessed additively (the 37 earlier baseline files are untouched).

* permute() + subgroups + jointtests + adversarial: the fixed-metric permutation joint tests and the
* permuted pairwise contrasts - e(perm_jointW/perm_sgW/perm_advW) and their *Wobs, e(pwsgs) with its
* perm_pwsgs_kernel, the per-subgroup tables and the joint_*_s* scalars. The largest permute surface
* with no bit pin anywhere (test_permute3 Q5 reconstructs it from the draws; nothing pins the draws).
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
set seed 112358
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) esfixedbaseperiod esrelativeto(-1) subgroups(sg) jointtests adversarial permute(shuffle) permreps(20) rseed(11235) dots(0)
local rc = _rc
dumpres perm_sgjt_adv `rc'

* permexact + permabsp + permstepdown + level(90) on the tiny panel (C(8,3) = 56 assignments): the
* exact plain-share convention (no add-one, e(perm_p_se)/e(perm_p_ci) missing), the absolute-value
* p with its symmetric dual interval, the opt-in max-|studentized| step-down (post AND pre families:
* the pre family needs two non-reference pre horizons, hence esprelength(3)), and lb_90ci columns
* under permute(). No draws are sampled, but e(perm_rngstate) is still posted, so the seed is set.
freshstart "`DATA'/dta_tiny"
set seed 314159
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) permute(shuffle) permexact permabsp permstepdown level(90) rseed(31415) dots(0)
local rc = _rc
dumpres perm_exact_absp `rc'

* weights + permstrata + permside + two estimands under permute(): the weighted replication fits
* (the only other weighted permute run anywhere is mc_ri R7, last run 2026-08-24; this is [aw] only -
* permute() under [pw] has no coverage anywhere, and a user [fw] is refused at the shell), the
* within-strata shuffle, the one-sided add-one share, and the multi-column kernel (att and att_it
* side by side).
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
set seed 271828
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) permstrata(sg) permside(right) rseed(27182) dots(0)
local rc = _rc
dumpres perm_strata_aw `rc'

* poisson under permute(): fast-mode ppmlhdfe replications, att-only (this pins the replication
* KERNEL, not poisson event-study or subgroup wiring under permute()). Nothing in the tree has ever
* run this combination (permctl does not refuse it).
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
set seed 161803
cap noi wooldid ycount i t ttre, att cluster(i) poisson permute(shuffle) permreps(20) rseed(16180) dots(0)
local rc = _rc
dumpres perm_pois `rc'

* bootstrap(bayes) + subgroups + jointtests + ties + bootdf(5) + level(90) + simultaneousbands: the
* Dirichlet reweighting path (no battery config runs bayes), the subgroup rewiring (per-sg tables,
* their family twins, e(pwsgs) from the contrast draws), the fixed-metric joint family
* (e(boot_jointW), e(boot_jt_family)), an integer reference df, the ties z0 convention and the
* non-studentized bootstrap band engine including the per-subgroup families.
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
set seed 141421
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests bootstrap(bayes) reps(10) ties bootdf(5) level(90) simultaneousbands rseed(14142) dots(0)
local rc = _rc
dumpres boot_bayes_sgjt `rc'

* boottstat + subgroups + adversarial + jointtests + esfixedbaseperiod + simultaneousbands: the
* studentized joint statistic e(boot_jointW), the studentized subgroup family (e(boot_sgW),
* e(boot_tstud_pwsgs), e(bootdist_pwse), e(boot_tpctl_pwsgs) - posted ONLY under boottstat with
* subgroups), the studentized adversarial test e(boot_advW) (exists only under boottstat,
* test_bootstrap_gates G38) and the percentile-t band engine (e(esband_engine) =
* bootstrap-studentized, test_esband B9). Every replication carries a full analytic VCE. On this
* fixture the across-subgroup equality tests are WITHHELD (joint_att_allsgssame posted as missing,
* e(boot_sgjt_reason) = nostudentizedinputs - the pairwise-se degeneracy under cluster(i) that
* test_bootstrap_joint documents); the dump pins that withheld state and its reason, not a live test.
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
set seed 173205
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) esfixedbaseperiod esrelativeto(-1) subgroups(sg) adversarial jointtests bootstrap(standard) boottstat reps(10) simultaneousbands rseed(17320) dots(0)
local rc = _rc
dumpres boot_ts_sgadv `rc'

* bootdisplay(bca) on a coarse cluster that nests i in pairs inside one onset block (9 clusters, so
* 9 jackknife fits): e(boot_accel), e(boot_z0), e(boot_ci_bca)/e(boot_p_bca), the bca mainresults
* twin and e(boot_jk_completed). Also the first battery run of the bootstrap with cluster() above i.
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
qui gen c2 = ceil(i/2)
set seed 223607
cap noi wooldid y i t ttre, att cluster(c2) esprelength(1) espostlength(1) bootstrap(standard) reps(10) bootdisplay(bca) rseed(22360) dots(0)
local rc = _rc
dumpres boot_bca_c2 `rc'

* weights + att_it + bootdf(normal) under bootstrap(standard): the physical re-identification leg
* (a cell-weighted estimand duplicates resampled clusters as whole units) under aweights, with the
* normal reference distribution overriding the clustered t(G-1) default. [aw] only: the [pw] twin,
* where the re-identified clusters carry their weights, has an rc pin (test_bootstrap_tstat T10)
* and a collapse-vs-full identity (test_fastcollapse F10c) but no bit pin.
freshstart "`DATA'/dta_binary"
qui keep if inrange(i, 1, 6) | inrange(i, 27, 32) | inrange(i, 51, 56)
set seed 244949
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) esprelength(1) espostlength(1) bootstrap(standard) reps(10) bootdf(normal) rseed(24494) dots(0)
local rc = _rc
dumpres boot_aw_reid `rc'

* att_i on the re-id leg under nesting (2026-09-16, feature tickets T2/T3): att_it was the deck's
* ONLY cell-weighted estimand under bootstrap(standard) - boot_std_reid and boot_aw_reid both carry
* att_it and nothing carried att_i. att_i averages cohort-specific effects equally, so it is the
* estimand whose weights a replication must rebuild from its own data; this pins that rebuild on
* the branch where cluster() nests i, which the finer-than-i work must leave byte-identical.
freshstart "`DATA'/dta_binary"
set seed 314159
cap noi wooldid y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) bootstrap(standard) reps(10) rseed(31415926) dots(0)
local rc = _rc
dumpres boot_reid_atti `rc'

* both gate routes at once: att_i (cell-weighted) AND boottstat under cluster(), the two distinct
* reasons wooldid_bootctl sends a run down the physical re-id leg. Nothing in the deck reached the
* leg by both routes simultaneously. Under boottstat the replication also substitutes its own VCE
* onto the resampled cluster identity, so this row pins the substitution's output under nesting.
freshstart "`DATA'/dta_binary"
set seed 271828
cap noi wooldid y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) bootstrap(standard) reps(10) boottstat rseed(27182818) dots(0)
local rc = _rc
dumpres boot_tstat_atti `rc'

* the re-id leg with cluster() STRICTLY COARSER than i, which nothing in the deck reached: every
* other config that lands on the leg clusters at cluster(i), where the cluster IS the cohort and
* group(__wdbns, i) degenerates to the drawn id alone. With 4 cohorts inside each drawn cluster the
* re-identification genuinely mints four new units per drawn copy, which is the behaviour the
* finer-than-i work must preserve on its nesting branch (keeping i there instead would let the
* equal-cell weighting count a twice-drawn cluster once). cl4 is generated, as sg/sgm/sgd2 are.
freshstart "`DATA'/dta_binary"
qui gen cl4 = ceil(i/4)
set seed 161803
cap noi wooldid y i t ttre, att att_i cluster(cl4) esprelength(2) espostlength(2) bootstrap(standard) reps(10) rseed(16180339) dots(0)
local rc = _rc
dumpres boot_reid_coarse `rc'

* att_itime on the re-id leg: measured 2026-09-16 as having ZERO bootstrap coverage anywhere in
* the deck (13 configs use att_itime, none under bootstrap()), so the third member of the
* cell-weighted family took the physical re-id leg entirely untested. It is subject to the same
* bootctl gate and the same two-regime re-identification as att_it and att_i.
freshstart "`DATA'/dta_binary"
set seed 141421
cap noi wooldid y i t ttre, att att_itime cluster(i) esprelength(2) espostlength(2) bootstrap(standard) reps(10) rseed(14142136) dots(0)
local rc = _rc
dumpres boot_itime_reid `rc'

* poisson under bootstrap(standard) on the full panel: ppmlhdfe replications on the count path,
* att-only (the replication kernel, not poisson event-study or subgroup wiring under bootstrap).
* bootstrap + poisson has no content-judged coverage anywhere (test_bootstrap_bayes O5 is rc-only).
freshstart "`DATA'/dta_binary"
set seed 264575
cap noi wooldid ycount i t ttre, att cluster(i) poisson bootstrap(standard) reps(10) rseed(26457) dots(0)
local rc = _rc
dumpres boot_pois `rc'

* the observation-level bootstrap (no cluster()) on the full panel: the resampling unit is the row,
* and bootdf defaults to normal (z_wald). The only engine path with no cluster().
freshstart "`DATA'/dta_binary"
set seed 282842
cap noi wooldid y i t ttre, att robust esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(28284) dots(0)
local rc = _rc
dumpres boot_obs `rc'

* level(90) on the analytic path: lb_90ci/ub_90ci columns, e(level) = 90 and the t critical value
* at .05 (test_level L1 hand-checks the formula; nothing pins the numbers). No RNG is consumed.
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) level(90)
local rc = _rc
dumpres level_analytic `rc'

di "BATTERY_DONE"
exit, clear
