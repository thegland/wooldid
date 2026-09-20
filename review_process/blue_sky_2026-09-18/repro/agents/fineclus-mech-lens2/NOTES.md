# FCM-4 — third lens (REALISM AND IMPACT). State as of launch.

Mechanism settled by code read, no Stata needed:
- ado:8910-8923  `smp` = (user if == 1) & !missing(clusvar).  NOTHING else is
  marked out — not y, not controls(), not weights.  Comment calls it "sample
  concept (A7): user-if sample + nonmissing cluster (permrun byte-parallel)".
- ado:8952       `wdb_C` = number of distinct clusters on that padded `smp`.
- ado:423        bootdf() DEFAULTS to `clusters` whenever cluster() is given.
- ado:10233      e(df_r) is overwritten with wdbdf (= wdb_C - 1) regardless of
                 bootdisplay().  So no opt-in beyond cluster()+bootstrap().
- ado:6460-6498  BOTH few-cluster cautions ("<10 clusters", "<6 treated
                 clusters") and the ">=2 clusters" REFUSAL are computed on
                 `smp2`, the same padded concept.
- ado:9537       wdb_C also feeds bootjnt (joint-test denominator df),
                 bootcellsx (cell df) and bootsave's observed-row weight.
- ado:7479-7486  permrun's `smp` is byte-identical — a dead cluster is also a
                 permutation unit.  DIFFERENT finding, not measured here.

Doc text that matters (sthlp bootstrap() entry, one paragraph):
  "bootstrap() requires at least two clusters IN THE ESTIMATION SAMPLE".
  The code checks the padded count, not the estimation sample.  The skeptic
  said there is no doc-vs-code divergence; that sentence says otherwise for
  the refusal, which is the same padded concept.

Probe: L2_fcm4.do — arms A1/A2 (trigger = missing merged CONTROL, not y),
B/B2 (DEFAULT invocation + lincom/test reach), D1/D2 (G_live=4, cautions),
E/E2 (G_live=1, does the >=2 refusal fire), F (dead cluster by singleton
drop, NO missing values anywhere).

Analytic magnitude grid (invttail, deterministic, no Stata):
  G_file/G_live  t(Gf-1)/t(Gl-1) = CI half-width ratio, padded vs correct
  50 / 40        2.0096/2.0227 = 0.994
  50 / 30        2.0096/2.0452 = 0.983
  24 / 12        2.0687/2.2010 = 0.940
  24 /  6        2.0687/2.5706 = 0.805
  20 /  5        2.0930/2.7764 = 0.754
  12 /  8        2.2010/2.3646 = 0.931
Material (>5% too narrow) only when G_live <~ 12 AND dead fraction >~ 30%.

## RESULTS (L2_fcm4.log, all arms completed, L2_FCM4_DONE at line 786)

| arm | trigger | printed header | printed reference | e(N_clust) | e(df_r)=e(boot_df) | e(df_r_analytic) | few-cluster cautions |
|---|---|---|---|---|---|---|---|
| A1 padded | controls(x) missing in 12 of 24 industries | N clusters = 12; N = 600 (600 obs dropped) | **t(23)** | 12 | **23** | 11 | n/a (G=12) |
| A2 dropped | same estimation sample, dead industries removed | N clusters = 12; N = 600 (0 dropped) | t(11) | 12 | 11 | 11 | n/a |
| B padded | A1 data, **DEFAULT** call (no bootdf, no bootdisplay) | as A1 | **t(23)** | 12 | **23** | 11 | n/a |
| B2 dropped | control for B | as A2 | t(11) | 12 | 11 | 11 | n/a |
| D1 padded | controls(x) missing in 20 of 24, **default call** | N clusters = 4; N = 200 (1000 obs dropped) | **t(23)** | 4 | **23** | 3 | **NONE PRINTED** |
| D2 dropped | control for D1 | N clusters = 4; N = 200 (0 dropped) | t(3) | 4 | 3 | 3 | **BOTH printed** |
| F | NO missing values anywhere: industry 1's firms each observed in one year, dropped as reghdfe singletons | N clusters = 23; N = 1150 (5 obs dropped) | **t(23)** | 23 | **23** | 22 | n/a |
| E | 1 live industry of 24 | refused r(459) | - | - | - | - | - |

Point estimates identical padded vs dropped: A1/A2/B/B2 att_i = 1.1285183513;
D1/D2 att_i = 1.0113907436.  Only the inference moves.

### What each arm settles
1. **Trigger is NOT missing-y.** A1/B: a merged control with partial coverage
   (the commonest cause of dropped rows in applied micro) does it.  **F: data
   with ZERO missing values does it too** - any internal screen that empties a
   cluster (here reghdfe's singleton drop) pads the count.  This is broader
   than either prior agent established.
2. **Zero opt-in.** Arm B is `cluster() + bootstrap(standard) + reps()` and
   nothing else.  bootdf() defaults to `clusters` (ado:423) and e(df_r) is
   overwritten regardless of bootdisplay() (ado:10233).
3. **It reaches a hand workflow.** After arm B, `test _b[main_att_i]=0` returns
   r(df_r)=23.  p at the padded df = 0.0000000000 vs 0.0000000525 at the live
   df (same b, same V).  wooldid's own footer names the consumers: "test and
   lincom run by hand ... and so do the tools that read e(V) without being
   asked: esttab, coefplot and honestdid".
4. **The printed default table is SAFE.** Under bootdisplay(percentile) the p
   and CI are rank-based and the `t_wald` column is b/se (D1: 1.0113907/
   .15146958 = 6.68).  None of them consumes df.  So the modal reader of the
   printed table loses nothing; the damage is on the e(V)+e(df_r) export path
   and under the bootdisplay(wald) opt-in.
5. **The output self-contradicts, three lines apart** (D1): "N clusters = 4;
   N = 200 (1000 obs dropped from initial sample)" ... "e(V) and e(df_r) carry
   the standard-deviation family referred to the t(23) distribution."
6. **The cautions are silenced exactly where the df error is largest** (D1 vs
   D2, identical estimation sample): the 4-cluster run gets neither "only 4
   clusters are available for the bootstrap" nor "only 4 treated cluster(s)
   are available for resampling" (MacKinnon-Webb).  Verified, as the skeptic
   claimed.
7. **The correct number is already in the same e() set**: e(df_r_analytic) is
   11 / 3 / 22 in the three padded arms - i.e. G_live-1, the value e(df_r)
   should have had.  This is an internal inconsistency, not a defensible
   sample-concept choice.

### HONEST NEGATIVE (my own hypothesis, refuted)
sthlp's bootstrap() entry says bootstrap() "requires at least two clusters IN
THE ESTIMATION SAMPLE" while that check (ado:6460-6472) runs on the padded
`smp2`.  I predicted 1 live cluster + 23 dead would slip past it.  **It does
not**: an EARLIER, sample-correct guard fires first - "Error: cluster(ind)
leaves 1 cluster(s) in the estimation sample; cluster-robust inference needs
at least 2" (arm E, r(459)).  So the refusal hole is unreachable.  Note what
this proves, though: wooldid ALREADY counts clusters on the estimation sample
elsewhere, so the padded count in the bootstrap leg is not a necessity.

### Magnitude grid (invttail, deterministic)
G_file/G_live -> padded/true CI half-width ratio -> true coverage of a "95%" CI
  3000/2800 -> 1.000 -> 95.00%     50/40 -> 0.994 -> 94.86%
    50/30   -> 0.983 -> 94.61%     51/12 -> 0.913 -> 93.02%
    24/12   -> 0.940 -> 93.71%     24/6  -> 0.805 -> 90.66%
    20/5    -> 0.754 -> 89.55%     24/4  -> 0.650 -> 86.96%
Material (>5% too narrow) only at G_live <~ 12 with dead fraction >~ 30%.

### CORRECTIONS TO THE SKEPTIC
- Their "25% too narrow" is the df effect ALONE.  Their own SE measurement
  (padded/dropped se ratios 0.989 / 1.119 / 1.176, mean 1.094) pushes the
  other way, so the NET wald interval in their own 20/5 case is ~18% too
  narrow, not 25%.  I did not re-measure the SE effect; the direction is what
  Jensen gives (E[1/m] > 1/E[m] for m ~ Binomial(G_file, G_live/G_file)) but
  the magnitude rests on 3 seeds and is not established.
- Their "no doc-vs-code divergence" is right about sthlp:532 but wrong about
  the bootstrap() entry's "at least two clusters in the estimation sample" -
  that sentence IS contradicted by the code, though as arm E shows the
  divergence is masked by an earlier guard and cannot be exploited.

### NOT MEASURED (same root cause, separate findings, not mine to file here)
- permrun's `smp` (ado:7479-7486) is byte-identical, so a dead cluster is also
  a permutation unit and a reassigned onset landing on it wastes a treated slot.
- bootdisplay(bca)'s delete-one-cluster jackknife runs over the padded set, so
  the dead clusters contribute G_dead identical pseudo-values to the
  acceleration.
- wooldid_bootjnt (ado:6738) takes the same nclus-1 as its F denominator df.
