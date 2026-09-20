# mc-fineclus-analytic — NOTES

TARGET: does the documented remedy (choose `i` COARSE, cluster FINER) actually deliver
correct analytic-SE coverage?  METHOD: Monte Carlo.

## FLAG RULE — WRITTEN BEFORE ANY RESULTS WERE LOOKED AT (2026-09-17)
An arm FAILS if, in the same draws where BOTH external controls (csdid, did_imputation)
sit inside the bands:
  * median over horizons of mean(reported SE) / empirical SD  < 0.90, OR
  * 95% coverage < 0.93
An arm also fails on the headline `att` row alone under the same two thresholds.
At 400-500 draws a coverage estimate carries a +/- 1.1pp standard error (sqrt(.95*.05/400)
= .0109), so a coverage of 0.93 is ~1.8 SE below 0.95 — the gate is deliberately loose.
Reported SE/SD ratios at 400 draws carry ~3.5% relative SE on the SD itself.

## DESIGN (written before results)
DGP-S ("state policy"): G=40 states x M=5 units x T=12 = 2400 rows.
  onset at the STATE level: 10 states ttre=5, 10 states ttre=8, 20 never.
  y = s_g + v_gt + a_u + 0.2 t + 0.5 cos t + tau + e
  s_g ~ N(0,.8^2) state RE; v_gt ~ N(0,.4^2) state-by-time shock (the BDM component:
  this is the variance component that a saturated state-by-relperiod dummy ABSORBS);
  a_u ~ N(0,.5^2) unit RE; e ~ N(0,.4^2).
  tau = 1.0 + .05*rel for the ttre=5 states, 0.4 + .05*rel for the ttre=8 states.
  TRUE cluster level = state.  cluster(state) is correct; cluster(unit) is too fine
  (state shocks survive) and is EXPECTED to undercover for DGP reasons.
Arms, all on the same draw:
  S1  i=unit  cluster(state)  -- cluster ABOVE i          (known-bad regime)
  S2  i=state cluster(state)  -- cluster AT i             (known-bad regime)
  S3  i=coh   cluster(state)  -- cluster FINER than i     ** THE NEW REMEDY ARM **
  S4  i=coh   cluster(unit)   -- finer still, but BELOW the DGP's dependence level
  CS  csdid long2 ivar(unit) cluster(state)               (control)
  DI  did_imputation cluster(state)                       (control)
DGP-X ("cross-cutting"): same states/units, onset assigned at the UNIT level so a state
  holds units of several cohorts -> state neither nests nor is nested in coh.
  X1 i=coh cluster(state) (CROSS-CUTTING); X2 i=unit cluster(state); controls.
TRUTH: an arm-specific zero-noise draw (all random components scaled to 1e-6), so the
  estimand each arm actually targets is recovered exactly.
INVARIANCE: b must be bit-identical across cluster() choices holding i fixed.

## STATUS LOG
(append below)

## PILOT (pilot.do / pilot.log) — 2026-09-17
Build asserted: `which wooldid` -> the pinned build path, version 1.4.
All four arms ran rc=0. Timings: i=unit 21.9s, i=state 1.5s, i=coh 0.36s.
**Point estimates are BIT-IDENTICAL across i in {unit,state,coh} and across
cluster() in {state,unit}** on the balanced panel (att .87743807 in all four).
So the MC isolates the variance path exactly. Reported att SEs on that one draw:
  i=unit  cluster(state)  se = .0853509   (V .00728474)
  i=state cluster(state)  se = .0657862   (V .00432783)
  i=coh   cluster(state)  se = .13334947  <-- the new remedy arm
  i=coh   cluster(unit)   se = .0816417   (V .00666537)
A 2x spread on identical points. wooldid's CI uses t(N_clust-1): (.87743807-.6077133)
/.13334947 = 2.0227 = t(39). So coverage is reported BOTH ways (own t, and 1.96).

## PREDICTION, LOGGED BEFORE THE MAIN RUN (so it cannot be fitted afterwards)
If the true SD of the att is ~0.11-0.13 as the hand-calc suggests, then
  S2 (i=state, cluster AT i) lands near SE/SD 0.5-0.6 — WORSE than the 0.66-0.79
     the prior MC found — and S3 (i=coh, cluster finer) lands near 1.0.
If instead S3 lands below 0.90, the documented remedy does not work and that is the
headline of this probe.

## DESIGN AMENDMENTS (advisor, before the run)
- Positive control is S2 (cluster AT i, 1.5s), not i=unit (22s). i=unit run later, fewer draws.
- Add the POOLED shape `esprelength(0)`: the prior MC's worst arm (0.661) was pooled,
  so the remedy must fix that shape too.
- Add **DGP-N**: identical code with s_g = v_gt = 0, so unit-level dependence is all there is
  and cluster(unit) is the CORRECT level. There i=coh + cluster(unit) is a real coverage claim.
- S4 (i=coh, cluster(unit)) UNDER DGP-S is EXPECTED to undercover for DGP reasons (state
  shocks live above the cluster). Logged now, before results, so it is not read as a defect.
- Record e(N_clust) and e(df_r) per arm; assert mreldif(b_S3,b_S4)==0 every draw.

## THE DOCUMENTED CLAIM UNDER TEST (build/wooldid.sthlp:251)
"The intended practice is to cluster BELOW the level of i - at a level whose units are
nested inside cohorts, so that a treated cohort-by-period cell contains several clusters.
At or above the level of i the cluster-robust covariance stops seeing the treated cells'
own variation, and two things follow that the table does not announce: the pointwise
standard errors on the att rows are understated (measured at roughly 0.66 to 0.75 of the
reference on the diagnostic fixtures) ... There is no plug-in correction ... so the
remedies are a finer cluster(), or bootstrap() or permute()".
DGP-S with i=coh + cluster(state) is EXACTLY that intended practice: each treated
cohort-by-period cell holds 10 states.

## PRELIMINARY (71 draws of DGP-S, 2026-09-17) — the headline is already clear
  arm                                   SE/SD range   cov95 range
  S2  i=state cluster(state) [AT i]     0.43 - 0.69   0.56 - 0.83   FAILS, badly
  S2p i=state cluster(state) pooled     0.51 - 0.69   0.66 - 0.89   FAILS
  S3  i=coh   cluster(state) [FINER]    0.84 - 1.16   0.92 - 0.96   PASSES  <== THE REMEDY
  S3p i=coh   cluster(state) pooled     0.89 - 1.11   0.93 - 0.99   PASSES  <== THE REMEDY
  S4  i=coh   cluster(unit)             0.46 - 0.67   0.61 - 0.80   expected (DGP: state shocks
                                                                    live ABOVE the unit) - logged
                                                                    in advance, not a defect
  CS  csdid  (control)                  0.84 - 1.26   0.90 - 0.99   in band
  DI  did_imputation (control)          0.87 - 1.11   0.90 - 0.99   in band
b invariance to cluster(): max mreldif over 71 draws = EXACTLY 0.
=> THE NEW REMEDY WORKS on the analytic path. Second result: the AT-i understatement is
   MUCH worse than sthlp:251's "roughly 0.66 to 0.75" on a realistic state-by-year-shock
   DGP - the ES rows run down to 0.43 and coverage to 0.556.

## STATUS 2026-09-17 ~01:50 — ALL MAIN MEASUREMENTS IN, FINDINGS.md WRITTEN
Waves completed (all output in out/, all analyses in RESULTS_*.txt):
  DGP-S 283 draws (mc_fineclus.do S)       — remedy PASSES (0.933/0.943), at-i FAILS (0.559/0.742)
  DGP-P 142 draws (mc_fineclus3.do)        — seeded positive control REPRODUCED (0.852/0.796
                                             vs the project's 0.838/0.794); remedy 0.951/0.957
  DGP-F 200 draws (mc_wave2.do F)          — 3 clusters/cohort: remedy 0.865/0.905, csdid 0.979/0.915
  DGP-X 200 draws (mc_wave2.do X)          — CROSS-CUTTING clean (0.999/0.950)
  DGP-U  92 draws (mc_wave2.do U)          — ragged/unequal clean (0.970/0.935)
  Ladder K=3/5/10, 150 draws each (mc_wave3.do) — wd/cs SE ratio 0.834/0.857/0.918, monotone
  probe_determ.do                          — every deterministic check CLEAN
Still running / queued: repro_remedy_works.do (self-contained repro #1),
  repro_fewclus_gap.do (self-contained repro #2), mc_fineclus.do N (redundant, lowest value).
