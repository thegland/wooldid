# fineclus-perm probe notes (agent: permute under finer-than-i clustering + permonsetcohorts)

## Harness
- build asserted: `which wooldid` -> .../bluesky/build/wooldid.ado (v1.4 2026-08-19 header)
- design A (`dataA.dta`, built by mkdata_A.do): 5 states, i=coh (onset grouping),
  cluster=state. coh1={s1,s2} onset 4; coh2={s3} onset 6; coh3={s4,s5} never.
  T=1..8, 40 rows. M = 5!/(2!1!2!) = 30 assignments.

## probe01_smoke  -- CLEAN
permonsetcohorts accepted, e(perm_irebuilt)=1, e(perm_npossible)=30 (hand count 30),
rep 1 of permsave == observed assignment and reproduces main_att exactly.

## probe02_identity -- CLEAN (machinery faithful)
Hand-recomputed all 31 rows of permA.dta externally, three ways:
 (A) observed-code map (onset4->1, onset6->2, missing->3)  == what the ado builds
 (B) fresh codes via `egen group(repons), missing`
 (C) relabelled observed-code map (10/20/30)
All three agree with each other and with the stored replication to ~1e-14 relative
(reghdfe convergence noise; max seen 114 ulps on one rep). NO i-code dependence.
NOTE: probe02 compared through Stata LOCALS, which round-trip through decimal, so its
equality test is only ~17 digits. The BITWISE comparison is probe09_tstat_weights.do,
which reads the saved values with %21x and reproduces b (and se under permtstat)
exactly, including under aweights and customweightmultiplier().

## probe03 (design A, BALANCED) -- instrument was BLIND
On the balanced design A every covariate option (fe(state), fe(cohcopy), controls(...))
left the reference distribution bit-identical, because finer-than-i FE are orthogonal to
the cohort-by-time design on a balanced panel. Rebuilt on a RAGGED panel (design B).

## probe04 + repro_F1  -- FINDING F1 (high)
Ragged panel, i=coh (onset grouping), cluster(state) finer, permonsetcohorts, permexact.
`controls(wave)` where wave is constant within the OBSERVED cohort (hence exactly
collinear with the i fixed effects and provably inert for the point estimate):
  point estimate  : bitwise identical with and without the option
  permutation p   : .11111111 -> .04444444   (crosses 5%)
  null-dist SD    : .77134377 -> .58882731   (narrower => ANTI-conservative)
CONTROL ARM cluster(coh) == i (no rebuild): p identical (.33333333 both). So it is the
rebuild. Mechanism reproduced by hand: in a hand-built replication the same covariate
moves the estimate by .0218 because it is no longer collinear with the REBUILT i.
No caution is printed. The bootstrap re-id leg prints the analogous A10 caution
(ado:412); the permute rebuild leg prints nothing.
Repro: repro_F1_cohortkeyed_covariate.do

## probe06_invariance -- CLEAN (mostly)
permexact: sorted reference distribution invariant to row order, to relabelling i
(100-coh, and a permutation of the codes) and to relabelling cluster ids, to <=3.6e-15.
Sampled permute at fixed rseed(12345): row order and i-relabel reproduce the draw
sequence exactly; relabelling CLUSTER ids changes the draw sequence (max per-rep diff
1.20) because the map frame is sorted by cluster value. Benign/expected, noted not filed.

## probe16_size -- KILLED AT DEADLINE, UNTESTED
80-dataset sharp-null size simulation (9 exchangeable clusters, i = 3 onset cohorts,
cluster(state) finer, permreps(99)) was killed at its deadline before printing. Not a
result either way. Exactness is nevertheless established by construction from the
identity probe + observed-is-a-member (probe14) + the hand-verified p convention.

## probes 07/08/09/10/12/13/14/15/17/18/19 -- see FINDINGS.md

## probe20_bootcmp.do -- WRITTEN BUT NOT RUN
The confirming contrast for FP-1 (bootstrap(standard) with cluster(coh) and att_i, with
and without fe(wave)/controls(cohx); the re-id leg REFINES the partition, so a
cohort-constant covariate should stay collinear there and the se/p should be identical).
Not launched: the project's 8-instance Stata cap was breached by other agents for the
whole wrap-up window (N stayed at 9-10). Run it as-is to close that leg.
