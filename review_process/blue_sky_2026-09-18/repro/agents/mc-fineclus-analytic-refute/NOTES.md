# mc-fineclus-analytic-refute — NOTES (skeptic)

Assignment: refute MC-FC-1 (finer-than-i remedy incomplete; SE shortfall is a gradient
in K = clusters per treated cell) and MC-FC-2 (sthlp:251's 0.66-0.75 band is not the size
of the at-or-above-i understatement on realistic data).

## Rule followed
Own generator (`skgen` / `ugen` / `probe_unbal.do`), own arms, own truth. The prober's
do-files were READ, never run.

## Theory worked out before any run (matters for the verdicts)
Saturated treated (i, rel) cell dummy demeans the residual across whatever sits inside
the cell. With Kc CLUSTERS inside the cell, E[sum_c score_c^2] loses exactly one of the Kc
degrees of freedom, so
        SE/SD  ~  sqrt(1 - phi/Kc),
phi = share of Var(bhat) coming from SATURATED cells (vs free pre-periods / controls).
The seeded sqrt((Kc-1)/Kc) is the phi->1 limit. cluster()==i is Kc=1.
KEY COROLLARY the seeded write-up does not state: Kc counts CLUSTERS, not units and not
noise levels, so OVER-clustering (a cluster coarser than the true dependence) makes the
reported SE SMALLER here, i.e. anti-conservative. mc2_mech.do tests exactly that on a DGP
with NO state-level dependence at all.

## Files
- pilot.do/.log        — harness check. 0.70 s per wooldid call. b bit-identical across
                         i=coh and i=state (.82521497); se .17367 vs .05666 -> 3.06x.
                         jointtests caution fires ONLY at i=state. e(joint_horiz) = the
                         pre-trend joint p.
- mc1_ladder.do        — 120 draws x K in {3,10} x arms A(i=coh,clus state,esfixedbase),
                         B(same, pooled reference), C(i=state,clus state). RUNNING.
- mc2_mech.do          — 120 draws, unit-independent DGP, arms E/F/G. QUEUED.
- probe_boot.do        — does bootstrap(standard) (remedy #2 in sthlp:251) repair K=3? QUEUED.
- probe_unbal.do       — NEW-FINDING probe: the only way to get cluster() finer than i in a
                         state-policy panel is to set i = the onset cohort, which throws
                         away the state fixed effects. Variance-only on a balanced panel;
                         what about an unbalanced one? QUEUED.

## Status
(updated as results land)

## RESULT 1 (NEW FINDING, certain) — probe_unbal2.do / probe_fe.do
ZERO-NOISE panel, homogeneous effect EXACTLY 1, 40 states, T=12, waves at t=5 and t=7,
deterministic permanent state levels s_g = 2*mod(state,7)-6.
The documented remedy geometry (i = the onset wave so that cluster(state) is finer than i)
silently BIASES THE POINT ESTIMATE on an unbalanced panel:

  case                                     i=state        i=coh (the remedy geometry)
  balanced                                 1.000000000    1.000000000
  one treated state attrits after t=8      1.000000000    1.210526316   (+21%)
  one treated state enters late at t=5     1.000000000    0.428571429   (-57%)
  ordinary ragged panel, 8% MCAR           1.000000000    1.214732454   (+21%)

probe_fe.do on the ragged case: att 1.2147, se .1296, t=9.37, p<0.001, CI [.953, 1.477],
and NO caution of any kind. fe(state) repairs it EXACTLY (1.000000000000, se 0), and so
does controls(i.state). sthlp's Standard Errors section never mentions either.
Mechanism: coarsening i replaces the state fixed effects with a single wave intercept, so
any change in the wave's state composition over time loads onto the estimate. Magnitude
scales with the dispersion of the sub-i permanent heterogeneity.

## Context checks done in the ado/help (matter for the verdicts)
- `ado:1421-1443` the jointtests caution predicate: `count distinct i` vs `count distinct
  (i, clus)` pairs. Fires ONLY when every `i` lies wholly in one cluster, i.e. Kcell = 1.
  Its closing clause claims a valid joint test "as does clustering strictly finer than the
  cohort-by-period cells" — MC1 arm A tests exactly that claim at Kcell = 3.
- `sthlp:1332`: "Where a treated (i, relative-period) cell lies wholly inside ONE cluster
  ... wooldid cautions on exactly that predicate". So the package's own model of the hazard
  is BINARY (one cluster), which is the seeded finding's point.
- NO few-treated-cluster caution exists on the analytic path. The `< 6 treated clusters`
  caution (`ado:6498`) and the `< 10 clusters` one (`ado:6479`) live inside
  `wooldid_bootctl` and fire only under bootstrap(). So the package already treats "few
  treated clusters" as a named hazard class — on the engine the help says is NOT a remedy
  for it — and says nothing on the analytic path the help calls the remedy.
- `sthlp:67-70` DOES pair a coarse `i` with `fe()`: "Choose i by the treatment heterogeneity
  ... not by the level at which you want fixed effects - those are specified separately in
  fe(), so a coarse i combines with unit-level ones." Framed as a modelling preference, not
  a correctness requirement, and absent from sthlp:255's remedy sentence. This CAPS the
  severity of SK-FC-A: it is a cross-reference + runtime-check gap, not undocumented.

## RESULT 2 — mc1_ladder.do, 120 draws per K, my own DGP (state RE + state-by-year shock
## + unit RE + iid error, 40 states x 4 units x 12 periods, waves at t=5 and t=7)
Arms: A = i=coh cluster(state) esfixedbaseperiod (cluster STRICTLY FINER than i)
      B = same, pooled reference
      C = i=state cluster(state) (cluster AT i)

              att SE/SD  att cov95 | med ES SE/SD  med ES cov95 | min ES SE/SD
  K=3   A        0.869     0.858   |    0.828         0.871     |   0.748
  K=3   B        0.884     0.900   |    0.807         0.845     |   0.737
  K=3   C        0.315     0.467   |    0.289         0.434     |   0.223
  K=10  A        1.013     0.958   |    0.986         0.938     |   0.900
  K=10  B        0.901     0.908   |    0.906         0.934     |   0.843
  K=10  C        0.634     0.800   |    0.612         0.750     |   0.483

=> MC-FC-1's K=3 rung REPRODUCES independently (0.869 / 0.858, below the 0.90/0.93 gate).
=> MC-FC-1's K=10 rung does NOT (1.013 / 0.958 here vs their 0.912 / ~0.93); mine agrees
   with the project's own 2026-09-14 figure of 0.972 / 0.934 at that same geometry.
=> sqrt((K-1)/K) is NOT a fit: it is a FLOOR implied by the stated mechanism
   (SE/SD = sqrt(1 - phi/Kc), phi<=1), and at K=10 I measure 1.013 above it while THEIR
   own K=10 measurement of 0.912 is BELOW it - which their mechanism cannot produce.
=> MC-FC-2's substance reproduces and is WORSE on the honest statistic: the ATT row (not a
   min over rows) is 0.315 SE/SD and 0.467 coverage at 6 treated states, 0.634 / 0.800 at
   20 treated states. Both outside sthlp:255's "roughly 0.66 to 0.75".

## RESULT 3 (NEW FINDING, the biggest one) — joint pre-trend test size, same run.
True pre-trend EXACTLY zero. Nominal 5%. Rejection rates:
   K=3   A (FINER, no caution)  0.617      K=3   C (AT i, caution) 0.925
   K=3   B (FINER, no caution)  0.617
   K=10  A (FINER, no caution)  0.108      K=10  C (AT i, caution) 0.550
sthlp:743 and the ado:1443 caution BOTH name "cluster strictly finer than the
cohort-by-period cells" as a FIX for the at-or-above-i oversizing. At 3 clusters per cell
that geometry rejects a true null 62% of the time and prints nothing.

## RESULT 4 — mc2_mech.do, 120 draws, DGP with NO state-level dependence AT ALL
(y = a_u + .2t + .5cos t + tau + e; units independent; `state' is only a grouping label)
Three arms, SAME draws, IDENTICAL b, so SD is shared (0.1168 on att):

 arm  geometry                   Nclust  clusters in    att       att    mean(se)
                                         treated cell  SE/SD     cov95
  E   i=coh,  cluster(state)        40         3       1.002     0.925   0.1170
  F   i=coh,  cluster(unit)        160        12       1.136     0.958   0.1327
  G   i=unit, cluster(unit)        160         1       0.362     0.567   0.0422

PAIRED ratio mean(se_E)/mean(se_F), which needs no SD at all:
  att 0.882, post0 0.867, post1 0.840, post2 0.900, post3 0.850, pre3 0.875, pre2 0.884
  mean 0.871   vs   sqrt(2/3)/sqrt(11/12) = 0.853 predicted by the seeded mechanism.
=> The mechanism MC-FC-1 states is CONFIRMED by a paired test with no SD noise, on a DGP
   with no cluster-level dependence whatever. Moving to a COARSER cluster() lowered the
   reported SE by 13%, reversing the usual conservatism of over-clustering.
=> Arm G refutes MC-FC-2's stated MECHANISM: 0.362 SE/SD and 0.567 coverage on a DGP with
   NO cluster-level time shock at all, far outside sthlp:255's 0.66-0.75. Their EXPECTED
   says "Remove it and the band reproduces"; it does not.

## RESULT 5 (NEGATIVE) — probe_ar1.do, 120 draws each arm, K=10, cluster FINER than i
   r0 iid errors  : att SE/SD 1.019, cov 0.942 ; ES range 0.869 - 1.031
   r5 AR(1) rho .5: att SE/SD 1.001, cov 0.942 ; ES range 0.843 - 1.023
The AR(1) errors in the prober's generator, the only structural difference from mine, do
NOT explain why their K=10 rung reads 0.912 where mine reads 1.013. Hypothesis refuted.

## RESULT 6 — probe_rank.do: the pre-trend covariance is NOT degenerate in the remedy geometry
   K=3   i=coh  (FINER)  rank 4/4  cond 2.23e+01  det  1.63e-06  joint p 0.037
   K=3   i=state (AT i)  rank 3/4  cond 6.20e+15  det -7.40e-25  joint p 1.0e-06
   K=10  i=coh  (FINER)  rank 4/4  cond 1.64e+01  det  2.02e-07  joint p 0.662
   K=10  i=state (AT i)  rank 3/4  cond 4.71e+15  det -2.29e-23  joint p 0.382
So SK-FC-C is not "a singular matrix inverted anyway" - the matrix is full rank and well
conditioned and the test is still 12x oversized. (The at-i arm IS rank deficient, cond 6e15,
and a p-value is printed from it; that one is the documented, cautioned case.)

## FILE INDEX (every one is self-contained: `clear all` + `adopath ++` to the pinned build)
Run any with:
  cd <this dir> && /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do <file>

  mc1_ladder.do     ~20 min  the K ladder; produces both verdicts' numbers and the joint-test
                             size table.  Sentinel MC1_DONE.
  mc2_mech.do       ~15 min  paired mechanism test on a DGP with NO state-level dependence;
                             arms E/F/G.  Sentinel MC2_DONE.
  mc2_noau.do       ~15 min  same with a_u = 0, so cluster(unit) is exactly right and E/F
                             isolates the (Kcell-1)/Kcell loss.  Sentinel MC2B_DONE.
  probe_ar1.do       ~6 min  AR(1) vs iid at K=10.  Sentinel PROBE_AR1_DONE.  NEGATIVE.
  probe_boot.do     ~60 min  does bootstrap(standard) close the K=3 gap.  PROBE_BOOT_DONE.
  repro_joint_finer.do ~25m  SK-FC-C repro with a csdid oracle and a narrow-window arm.
                             Sentinel REPRO_JF_DONE.
  repro_nocaution.do ~1 min  SK-FC-B repro: analytic silent, bootstrap cautions, same data.
  probe_rank.do      ~1 min  rank/conditioning of the pre block of e(V_es_att).
  probe_unbal2.do   ~30 sec  SK-FC-A repro, ZERO NOISE, exact arithmetic.
  probe_fe.do       ~30 sec  fe(state) / controls(i.state) repair it exactly.
  probe_unbal3.do   ~30 sec  bias scales linearly with sd of the sub-i heterogeneity.
  probe_unbal.do    ~2 min   the first, noisier version of the unbalanced probe.
  pilot.do          ~20 sec  harness check.

## STATE AT HANDOFF
FINDINGS.md is complete and is the deliverable. Verdicts: MC-FC-1 PARTIAL/medium,
MC-FC-2 PARTIAL/medium. New: SK-FC-A (high), SK-FC-C (high), SK-FC-D (medium),
SK-FC-B (low). Every number in FINDINGS.md comes from a completed run.

STILL RUNNING when this was written (none of them can change a verdict; they refine):
  repro_joint_finer.do  - a self-contained re-run of the SK-FC-C size measurement, plus a
                          csdid `estat pretrend` oracle and an esprelength(2) arm. The
                          SK-FC-C numbers already in FINDINGS.md come from mc1_ladder.do,
                          which is COMPLETE. If the esprelength(2) arm comes back near
                          nominal, SK-FC-C narrows to "ragged windows at few clusters" and
                          its confidence field should drop; if csdid also rejects 30%+,
                          severity drops to medium.
  probe_boot.do         - 26 of 30 bootstrap draws done; paired se_bootstrap/se_analytic
                          mean 1.168 (min 0.932, max 1.797). The final table adds SE/SD and
                          coverage for both arms.
  mc2_noau.do           - clean control for SK-FC-D (unit RE removed so cluster(unit) is
                          exactly right). Only affects how tightly SK-FC-D can claim the
                          0.853 constant, not the direction or the ~13%.
