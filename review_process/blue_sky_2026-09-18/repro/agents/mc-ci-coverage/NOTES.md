# mc-ci-coverage — running notes

TARGET: coverage of every CI variant wooldid reports. METHOD: Monte Carlo + interval-shape sanity.
BUILD: adopath ++ .../bluesky/build  (verified `which wooldid` = pinned build in every log).

## Fixture (mkfix.doi)
Balanced staggered panel. `unit` = finest id, `cg` in 1..4 (cohorts 1,2,3 treated at 2002/2004/2005;
cg==4 never-treated), FT periods from 1999. `base = ai + gt + truetau` is NOISELESS and FIXED;
`y = base + ui + e`. So wooldid's estimator (linear in y, weights independent of y) has
E[b(y)] = b(base) exactly -> the noiseless run's point estimate IS the design truth for every
estimand. Geometry: **i = cg (4 cohorts), cluster(unit) (finer than i)** = the geometry where the
analytic path is documented sound.

## DONE
### probeA.do / probeA.log — analytic level() honouring
CLEAN. level(90/95/99/87.5): e(level) correct, columns named lb_90ci/lb_87p5ci, and the implied
critical value is EXACTLY invttail(G-1, (100-L)/200) on every row (att family and ES rows):
  90 -> 1.657759285 (= t_119), 95 -> 1.980099876, 99 -> 2.617775976, 87.5 -> 1.545006301.
No lb>=ub, no b-outside-CI, no half intervals (20 rows x 4 levels).

### probeB.do / probeB.log
- NEST: analytic 90 strictly inside 99 on all 20 rows. CLEAN.
- Analytic simultaneous bands (engine `analytic-simulated`, 10000 draws): crit 2.2895 (L=90),
  2.5645 (95), 3.1144 (99) — all strictly above the pointwise t_119 critical value; no band
  narrower than the pointwise interval; no band lb>=ub; b inside band everywhere. CLEAN.
- subgroups(sg): e(wdidsgresults_0/_1), e(pwsgs) CLEAN on all shape checks.
- e(gteffects) with level(95): 13 rows, all CLEAN.
- `unconditionalse` with a non-att estimand -> informative r(198); `impvalues` without cluster()
  -> informative r(198). (Known/by design.)

### probeC.do / probeC.log — bootstrap interval reconstruction (40 units, reps(399), bca)
percentile and bias-corrected endpoints reconstruct to EXACTLY 0.0 difference from my own
independent implementation at levels 95 and 90. `normal` differs by ~1.6e-2 = the t(G-1) vs
normal critical value (wooldid uses invt(G-1) there, my recheck passed dfden=.). No lb>=ub,
no half intervals, b inside interval in every family. (Still running for 99 / 87.5.)

## NEXT (planned)
MC arms. Decision rule PRE-REGISTERED: flag when |coverage - nominal| > 3*sqrt(p(1-p)/R).
1. analytic, R=500, + positive control (`robust` on a DGP with a unit-level shock: must undercover)
2. bootstrap(standard): normal/percentile/bc/bca at several levels, R~150, 4 parallel chunks
3. boottstat percentile-t, R~100
4. permute: perm_ci (shift) + permtci, constant-effect DGP arm
5. impvalues at a cross-cutting cluster with G~25
6. simultaneous-band SIMULTANEOUS coverage, and the boottstat "band can never be narrower"
   help claim at few clusters with an outlier.

## RESULT: analytic arm CLEAN (mc_analytic.do / mc_analytic.txt, R=500)
i=cg (4 cohorts), cluster(unit) G=40, df_r=39, balanced 40x8 panel, truth = noiseless run's e(b).
mean(b) == truth to MC error on all 16 parameters; mean(se) == sd(b).
Coverage: c95 in [0.946,0.956], c90 in [0.880,0.904], c99 in [0.986,1.000], c80 in [0.782,0.810].
Nothing exceeds 3 binomial SE. The analytic interval machinery at this geometry is SOUND.
`robust` on the same data OVER-covers (c95 ~ 0.99, mean(se) 0.266 vs sd(b) 0.202): the unit shock
cancels in the within-unit contrast, so clustering REDUCES the variance and robust overstates.
That means `robust` is NOT a usable positive control here -> mc_power.do (cluster AT i) is.

## LEAD (probeD): impvalues wipes out its whole inference surface on noisy realistic data
`wooldid y cg t ttre, att cluster(reg25) impvalues` with i=cg (4 cohorts) and reg25 = 25 regions
each holding exactly one unit per cohort: EVERY se/t/p/lb/ub is missing, with
  "Caution: ... the per-cluster estimates were numerically identical ... A collapsed IM standard
   error usually reflects deterministic or noiseless data - please verify your input".
The data is NOT noiseless (sd 0.6 idiosyncratic + 0.5 unit shock).  Diagnosis from the code
(ado:2046-2060): the per-cluster "estimates" are NOT per-cluster fits - they are
  `margins, dydx(targetflag) over(cluster)` on the ONE pooled model,
i.e. deterministic re-weightings sum_c w_{g,c} * beta_c of a single coefficient vector.  When
every cluster has the same cohort composition the weights are identical and the G "estimates"
coincide EXACTLY - which is only possible if they carry no independent sampling variation at all.
probeF.do tests the consequence: with unequal composition the spread is nonzero but should then
be a function of cell-effect heterogeneity rather than of the noise level.  Sweep sd 0.6/0.06/0.006.

## CONFIRMED: MC1 impvalues CI coverage 0.032 at nominal 0.95 (mc_imp.txt, R=187+)
mean(b)=0.9983 vs truth 1.0014 (unbiased), sd(b)=0.0885, mean reported IM se=0.00257 -> 34x
too small. Analytic clustered arm on the SAME draws: mean(se)=0.0883, sd(b)=0.0887, cov 0.967.
So the point estimate and the target are fine; only the IM standard error is wrong.
Self-contained repro: repro_MC1_impvalues_ci.do

## CLEAN so far (negative results worth reporting)
- Bootstrap interval RECONSTRUCTION: my independent implementation of normal/percentile/bc/bca
  matches e(boot_ci_*) to 0.0 difference on every draw (mcboot NATIVEDIFF column = 0 on all
  draws so far), with dfden = G-1 for the normal family. So level() is honoured exactly on all
  four bootstrap families.
- Bootstrap interval NESTING across level 80/90/95/99 on the same draws: 0 violations in 2582
  ordered pairs, including bca (nest_boot.py).
- Point estimate inside its own interval: 0 violations out of 3448 (all four families, all four
  levels) so far. (Percentile-t not yet measured.)
- Analytic simultaneous bands: crit > pointwise t critical value at 90/95/99; band never
  narrower; b always inside.

## RUNNING
mcboot_c1/c2 (bootstrap coverage, 35 draws each), mc_cheap (power check + unconditionalse +
impvalues at a NESTING cluster), probeH (permute intervals), probeE2 (band vs pointwise under
boottstat at level 99 with an asymmetric outlier cluster; + bca half-interval guard),
mc_band (simultaneous-band simultaneous coverage, analytic engine, R=400).
probeE (reps 299 version) was KILLED for cost; probeE2 replaces it.

## RESULTS (2nd batch)
### mc_cheap.txt (R=250 each)
  A_clusATi  cov95=0.544 ratio(se/sd)=0.280  FLAG  <- POSITIVE CONTROL: harness has power
  B_uncond   cov95=0.980 ratio=1.230               <- not flagged; vce(unconditional) prices in
                                                      covariate-distribution sampling my fixed
                                                      design does not have. NOT a defect.
  C_impnest  cov95=0.916 ratio=0.783               <- impvalues with the cluster NESTING i
                                                      (i=unit, cluster(reg) of 3 whole units):
                                                      mild undercoverage, not flagged. SCOPES MC1.
### mc_imp2.txt  i=cg, cluster(unit) (cluster FINER than i - the help-recommended pairing)
  im_clusfiner cov95=0.256  ratio=0.202  (R=39 so far)   <- MC1 at a geometry T4 does not touch
### mc_band.txt  analytic sup-t simultaneous band, G=40, R=400
  level 90 -> SIMULTANEOUS coverage 0.875 (3SE 0.045); level 95 -> 0.930 (3SE 0.033).
  Below nominal at BOTH levels in the same direction (1.7 and 1.8 SE) but neither flagged.
  Mechanism candidate: wooldid_esbanda simulates the critical value from rnormal() draws
  (ado:16725) while the pointwise intervals printed beside it use invttail(G-1).  At G=40 that
  is ~1%.  mc_band10.do repeats at G=12 where the gap is ~12% -> decisive either way.
### probeE2 (boottstat + simultaneousbands, 12 clusters, asymmetric outlier, level 99 and 95)
  band/pointwise width ratio 1.78 - 8.02 in every horizon of every seed so far: the band is
  always WIDER.  The help's "under boottstat the band can never come out narrower" survives.

## permute() interval anatomy (probeI.log) — for anyone reading this after me
At i=unit, cluster(unit), G=40, permreps(99), permtstat permtci, level(90):
  e(wdidmainresults) = estimate 1.093400346 | se_analytic 0.139293212 | t_perm 7.8496 |
                       p_perm 0.02 | lb_90ci 0.851974548 | ub_90ci 1.322739463
  e(perm_ci)  ll 0.851974548  ul 1.322739463         <- THE PRINTED INTERVAL
  e(perm_lb_att) = -1.6464486   e(perm_ub_att) = 1.7332201
  e(perm_tlb_att) = 0.86406123  e(perm_tub_att) = 1.3227395
DECODED: e(perm_lb/ub_*) are the level-QUANTILES of the permutation (studentized) distribution,
NOT interval endpoints.  The printed interval is the shift interval
   [ b - q_hi*se_obs , b - q_lo*se_obs ] = [1.0934 - 1.7332*0.13929, 1.0934 + 1.6464*0.13929]
and e(perm_tlb/tub_*) is permtci, symmetric about b (half width 0.22934 = q(|t|)*se).
This is why probeH's G=4 run printed lb_90ci = 0 EXACTLY: with 4 clusters the observed
assignment is the most extreme of the reference set, so q_hi = t_obs and b - t_obs*se = 0.
Structural at G=4, not a defect.
probeH: permute interval SHAPE clean at both cluster(cg) (G=4) and cluster(unit)+permonsetcohorts
(G=40, e(perm_irebuilt)=1): lb<ub, b inside, no half intervals, 90 nested in 99, 12 rows each.
probeE2: 30 horizon comparisons under boottstat+simultaneousbands at 12 clusters with an
asymmetric outlier, levels 99 and 95 -> the band was WIDER than the pointwise percentile-t in
every single one (ratio 1.78-8.02). The help claim survives. No bca half-intervals either.

## probeI ESOPT: requesting event-study horizons moves the headline att AND its SE
Identical data (fixI.dta, 40 units, i=cg, cluster(unit), N=320, df_r=39 in all four runs):
  no ES options   att = 1.093400346  se = 0.132401198
  esprelength(1)  att = 1.117099722  se = 0.158879296   (+20% se)
  esprelength(2)  att = 1.132362929  se = 0.173433656   (+31% se)
  esprelength(3)  att = 1.065211265  se = 0.211200984   (+60% se)
MECHANISM, and it IS documented - but only in the Stored Results entry for e(escomp_refN):
"Under the pooled regime that pool is the union of treated rows earlier than the requested
pre-window (relyr < -esprelength) and each treated unit's forced first in-sample observation".
So asking for pre-horizons REMOVES those periods from the reference-normalization pool.
The esprelength() option entry itself says nothing about it.  All four estimates remain
consistent for the same truth (1.0 here), so this is efficiency, not bias.  LOW severity,
doc-placement/UX.  It also explains why mc_analytic's sd(b)=0.183 (esprelength(2)) exceeds
mc_imp2's 0.125 (no ES options) on otherwise-similar designs - NOT a defect, and worth
recording so nobody else burns time on that discrepancy.

## FINAL STATE (see FINDINGS.md for the write-up)
THREE findings:
  MC1 (critical, certain)  impvalues CI coverage 0.29 at i=cg+cluster(unit) (R=200) and 0.048 at
                           a cross-cutting cluster (R=250); IM se 5.9x / 34x too small; the
                           per-cluster reads are margins over(cluster) on ONE pooled fit
                           (ado:2042/2047-2048/2054 and the three sibling sites).
                           repro: repro_MC1_impvalues_ci.do  (verified to run: R=40, cov 0.250)
  MC2 (medium, certain)    analytic sup-t simultaneousbands under-covers: 0.760 at nominal 0.90
                           and 0.848 at nominal 0.95 with G=12 (R=388); 0.875 / 0.930 at G=40
                           (R=400). crit is a NORMAL sup quantile (ado:16725) while the table's
                           pointwise intervals use invttail(G-1).
                           repro: repro_MC2_band_coverage.do (verified to run: 0.780 at nom 0.90)
  MC3 (low, certain)       esprelength()/espostlength() move the headline att and inflate its se
                           up to +60% on identical data; documented ONLY under e(escomp_refN).
                           repro: probeI.do (ESOPT block)
Everything else probed came back CLEAN - see FINDINGS.md "PROBED HARD AND CLEAN" (11 items).
Underpowered at write-up: bootstrap coverage R=36, boottstat percentile-t R~10, permute R~10.

## probeJ (bootstrap + subgroups + gteffects CI shape, levels 90/95) — CLEAN
36 bootstrap CI columns x 3 families, e(wdidsgresults_0/_1), e(pwsgs), e(gteffects) (18 cells):
zero lb>=ub, zero half intervals, zero b-outside-CI, level label follows into e(gteffects).
## VERIFIED before submitting: all 36 NATIVEDIFF values in mcboot_1.txt+mcboot_2.txt are exactly 0.
