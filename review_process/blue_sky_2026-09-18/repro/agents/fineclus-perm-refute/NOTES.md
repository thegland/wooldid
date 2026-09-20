# fineclus-perm-REFUTE notes (skeptic on FP1/FP2)

## Harness
- build asserted with `which wooldid` in every probe -> .../bluesky/build/wooldid.ado
- OWN data: 7 states, i = coh (onset grouping, 2 treated cohorts + never), cluster(state)
  strictly finer. Ragged windows (different t0/t1 per state), state level+trend
  heterogeneity. M = 7!/(3!2!2!) = 210 enumerated assignments -> p grid = 2/210 = .00952,
  much finer than the prober's M=90 (grid .0222).
- GOTCHA found: `permsave(NAME, replace)` is NOT a sub-option; `replace` is a top-level
  wooldid option (sthlp:700). `permsave(NAME, replace)` silently treats "NAME, replace"
  as the FILE NAME. Documented, not a finding, but it cost one arm.

## probe01_smoke  (panelR.dta, TE=0.6)
BASE: b=.69998408, p=.0095238 (floor, 2*1/210), SD=.28831. e(perm_irebuilt)=1,
e(perm_failed)=0, 210/210 completed. 34 s per 210-rep run.

## probe02_arms (panelR2.dta = panA, TE=0.22)  -- PREMISE OF FP1 CONFIRMED
Observed fit ALONE (no permute), %21x of e(b)[main_att]:
  none            +1.6690b859bb834X-002
  controls(cohbase)  +1.6690b859bb834X-002   BITWISE IDENTICAL
  coarsecohortcontrols(c.cohbase) +1.6690b859bb834X-002  BITWISE IDENTICAL
  fe(wavegrp)     +1.6690b859bb838X-002   (4 ulps - reghdfe convergence noise)
  fe(coh)         +1.6690b859bb82cX-002   (8 ulps)
  fe(state)       +1.b90008cd4f840X-004   (a genuinely different model - contrast arm)
=> a covariate constant within the observed i is INERT for the point estimate. TRUE.

Permutation arms, cluster(state) (rebuild fires, e(perm_irebuilt)=1, 0 failures):
  CTRL    p=.0285714  SD=.2207672
  FEWAVE  p=.0285714  SD=.2196989
  FECOH   p=.0380952  SD=.2181833
  CCC     p=.8952381  SD=1.2655240   <-- SD 5.7x the others, no failures reported
  FESTATE p=.7142857  SD=.1968140    (different model, not part of the claim)
CONTROL ARM cluster(coh) == i (no rebuild, e(perm_irebuilt)=0, M=6):
  CBASE p=.3333333 SD=.2566696171
  CCTRL p=.3333333 SD=.2566696171   <-- BIT-IDENTICAL. The option is inert with no rebuild.
BASE arm failed r(602) on the permsave gotcha; rerun in probe03.

## probe03_paired  -- running
(probe02 BASE arm rerun in probe03 as ABASE.)

## probe03_paired + probe05_pairing  -- MECHANISM CONFIRMED ON MY OWN DATA
Two independent panels (panA seed 771131, panB seed 408822), TE=0.22, M=210.
Assignment enumeration is IDENTICAL across arms (S_ABASE vs S_ACTRL onset_* columns:
ASG_MISMATCH_ROWS = 0 of 211), so the draws pair row by row.

PANEL A (b_obs = 0.35016143844303, bitwise identical in every inert arm):
  arm      p           nullSD      SD vs BASE   #assignments whose statistic changes
  BASE     .02857143   .22768479      --         --
  CTRL     .02857143   .22076719    -3.0%        208 of 210   (0 at the identity, exactly)
  FEWAVE   .02857143   .21969891    -3.5%        198 of 210
  FECOH    .03809524   .21818330    -4.2%        208 of 210
  CCC      .89523810  1.26552397   +456%         208 of 210
PANEL B (b_obs = 0.33501463824096):
  BASE     .10476190   .24825485      --
  CTRL     .11428571   .22942360    -7.6%
  FEWAVE   .08571429   .23534421    -5.2%
  FECOH    .11428571   .22928515    -7.6%
  CCC      .80000000  1.37229962   +453%
=> the reference distribution really is rewritten, and it is EXACTLY unchanged at the
   identity assignment (|diff| = +0.0000000000000X-3ff, bitwise zero).
=> but the p MOVES BY AT MOST ONE GRID STEP (2/210 = .00952) in the controls()/fe() arms,
   in BOTH directions (panel B CTRL made p LARGER). The "2x" in the original finding is
   one grid step at M=90. The null SD shrinks 3-8%, not 11-24%.
=> the LARGE, reproducible arm is coarsecohortcontrols(): null SD x5.5 on BOTH panels,
   p .029 -> .895 and .105 -> .800. The original prober DEMOTED this arm.

## probe06_ccc  -- the blow-up is the REBUILD, and only for the inert class
  cluster(state) (rebuild):  BASE p=.0286 SD=.2277
                             CCC(c.cohbase)  p=.8952 SD=1.2655  b BITWISE = BASE
                             CCC(c.a_s)      p=.4762 SD=.2169   b DIFFERENT, no blow-up
  cluster(coh)==i (no rebuild): BASE p=.3333 SD=.2566696171
                             CCC(c.cohbase)  p=.3333 SD=.2566696171  BIT-IDENTICAL
  a_s is a STATE-level control (varies within the cohort), so it is live in the observed
  fit too -> no asymmetry -> no blow-up. It is the inert-in-observed/live-in-replication
  asymmetry, not coarsecohortcontrols() as such.

## probe06 (2) -- NEW FINDING I FOUND MYSELF
interactivecontrols(cohbase) with NO permute():  prints
  "Caution: every treated cohort-period slope term for the interactive control coh..."
  (the documented sthlp:1068 reduction caution)
Same spec WITH permute(shuffle) permreps(20) permonsetcohorts: that caution is ABSENT.
Same spec WITH permute + permanalyticse: the caution is BACK.
=> permute()'s 2026-09-05 fast-mode default silently disables the one diagnostic wooldid
   owns for "this control has no within-cohort-period variation" - on the very path where
   such a control goes live inside every replication. The user never typed fastnose.

## probe07_cautions -- FP2's facts, verified
  A permute + permonsetcohorts + fe(coh)                : NO covariate caution
  B permute + permonsetcohorts + controls(cohbase)      : NO covariate caution
  C bootstrap(standard) cluster(coh)==i att_i fe(coh)   : "Caution: a covariate
    specification names the i variable coh; under resampling, ..."  PRINTED
  D same but controls(cohbase) (constant within i, NOT named i) : NOT printed
  E bootstrap cluster(state) finer than i, fe(coh)      : NOT printed (wdb_nest gate)
=> FP2's asymmetry is REAL but NARROW: the bootstrap caution is a NAME-EQUALITY tripwire
   (ado:410 `if "`wdbcv'" == "`i'"`), so the CLASS FP1 identifies (constant within i) is
   uncautioned on BOTH engines. FP2 overstates the bootstrap side's coverage.

## probe04 / probe08 / probe09 -- 24-dataset SHARP-NULL sweep, RANDOMISED assignment
6 states, cohorts redrawn at random every dataset, y with no treatment term, M=90.
SD ratio (arm nullSD / BASE nullSD):
  controls(cohbaseY)  med 0.977  [p25 .889, p75 1.016]  max 1.29   narrower in 16/24
  controls(cohidx)    med 0.981  [.945, 1.005]          max 1.23   narrower in 16/24
  fe(early)           med 0.991  [.958, 1.039]          max 1.29   narrower in 14/24
  coarsecohortcontrols(c.cohbaseY)
                      med 3.851  [2.144, 7.557] p90 17.4 max 229.6  narrower in 0/24
p-values:  BASE min .0667 med .6222 ; CCC min .4000 med .8222 ; CCC > BASE in 20/24
  CCC nullSD > 2x BASE in 19/24, > 5x in 10/24.  e(perm_failed)=0 throughout.
=> controls()/fe(): a few percent, BOTH directions -> FP1's "anti-conservative" is not a
   systematic direction, it is dataset noise around a ~2% median narrowing.
=> coarsecohortcontrols(): systematic, huge, one-directional. The test cannot produce a
   p below .40 on any of 24 datasets. THAT is the real defect inside FP1.

## Line numbers verified on the pinned build (FP1 cites 7639; actual is 7641)
  ado:7641   local repcall `y' `repi' `t'      (varlist-only substitution)
  ado:6022-6034  subgroups() refused under permonsetcohorts, with the FP1 mechanism
                 spelled out in the comment
  ado:5414-5418  the 2026-09-05 permute() fast-mode default
  ado:12706  the interactivecontrols() reduction caution (gated out under fastnose)
  ado:410    the bootstrap A10 caution -- `if "`wdbcv'" == "`i'"`, NAME EQUALITY only

## UNTESTED (honest list)
- customweightmultiplier() keyed to the observed i under the rebuild (same repopts
  pass-through; not probed for want of Stata slots)
- ccc_absorb() and interactivecontrols() as FP1-class carriers (interactivecontrols is
  NOT in the class: it changes the observed estimate, so the asymmetry does not arise)
- whether the CCC blow-up also fires under permute(probabilities) or permstrata()
