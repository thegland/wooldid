# SKEPTIC report -- refutation pass on BT-01..BT-05
## (surface: boottstat, bootdf(), bca, simultaneousbands)

Build: the pinned bluesky build, loaded with `adopath ++`, asserted with `which wooldid`
at the top of every probe. Every dataset is mine (inline DGPs, or `gen_unbal.do` which I
wrote). The project test suite was not run and its fixtures were not read.
Every repro named below is self-contained: it sets the adopath, builds its own data and
prints a sentinel.

| id | verdict | severity as filed | severity revised |
|----|---------|-------------------|------------------|
| BT-01 | PARTIAL   | high   | medium |
| BT-02 | CONFIRMED | medium | medium (extended to 3 markout sources) |
| BT-03 | PARTIAL   | medium | low (scope is wrong: family size, not cluster count) |
| BT-04 | CONFIRMED | medium | **high** (and the advW leg is now measured, not read) |
| BT-05 | PARTIAL   | high   | low as stated -- but see SK-01, which is critical |
| SK-01 | NEW       | -      | **critical** |
| SK-02 | NEW       | -      | medium (uncertain on realism) |
| SK-03 | NEW       | -      | medium |

---

# VERDICTS ON THE SEEDED FINDINGS

## BT-01 -- band narrower than the pointwise percentile-t under boottstat: **PARTIAL**
repro: `repro_bt01_band.do` (cell A); sweep evidence `r03_bands.log`, `r08_bandsweep.log`

TRUE: the help's boottstat exemption at sthlp:528 is unsound and the band does fail to
nest the pointwise interval under `boottstat`.
* `repro_bt01_band.do` cell A (20 states, log-normal errors, 5 horizons, reps 299):
  `A_SUM m=5 inside_lo=0 inside_hi=1 narrower=1 crit=30.9196`
  horizon pre2  pointwise [-1.5106454, 66.1101377] width 67.62
                band      [-25.4339246, 25.2024119] width 50.64
  -- strictly narrower AND non-nesting on the upper side.
* `r03_bands.log` S1 (my ragged 25-state panel, log-normal, reps 399): one horizon with
  the band's upper end inside the pointwise one; band NOT narrower.
* `r03_bands.log` S2 (same data, percentile display): no violation -- so it is specific
  to the studentized engine, as the mechanism predicts.
* `r08_bandsweep.log` boottstat cells 1-6: violations in 2 of 6 (BCELL1 inside_hi +
  narrower, BCELL6 inside_hi only).

WRONG in the finding as written:
1. Its premise -- "Mathematically a simultaneous band at level() must be at least as wide
   as the pointwise interval at level() at every horizon of the family" -- is FALSE.
   Width dominance holds only between intervals of the SAME SHAPE. `wooldid_boottp`
   (ado:14332) builds the EQUAL-TAILED `[b - q_hi*se, b - q_lo*se]`; `wooldid_bootband`
   (ado:16776) builds the SYMMETRIC `b +/- q_lev(max_h |t|)*se`. Two valid level() objects
   of different shape can cross -- which is exactly what the help's OWN preceding sentence
   says about percentile/bc/bca. The defect is that the sentence EXEMPTS `boottstat` when
   percentile-t is itself a rank interval and itself asymmetric, not that a law of
   simultaneous inference was broken.
2. Their Run 1 does not demonstrate the title. Band [-6.967, 4.392] has width 11.36
   against the pointwise 10.30 -- the band is WIDER. That run shows NON-NESTING, the
   behaviour the help already documents for the other displays. Only their Run 2 (widths
   13.83 vs 14.12) is narrowness, and by 2 percent.
3. The one case where I could make the band STRICTLY narrower needed a pathological
   studentized distribution: `crit = 30.9` and a pointwise interval of [-1.5, 66.1].
   That critical value IS SK-02's inflation, so BT-01's strict narrowness and SK-02 are
   the same event, not two. On well-behaved data I get NON-NESTING only, never narrowness.
   That is honest to report, and it caps the severity.
4. Category: brief category 5 (documented-vs-actual), not category 3 (an inference path
   that does not deliver). No coverage was measured on either side; a symmetric-in-|t|
   sup band is not invalid for crossing an equal-tailed pointwise interval.
REVISED SEVERITY: **medium**, confidence certain (for the doc divergence and non-nesting).

## BT-02 -- bootdf(clusters) counts clusters outside e(sample): **CONFIRMED (extended)**
repro: `repro_bt02_bootdf.do`

Reproduced on my own 26-state panel and on THREE markout sources, not the one used:
```
A outcome missing in 4 states : e(N_clust)=22  e(boot_df)=e(df_r)=25   A_GAP=4
B same sample, states dropped : e(N_clust)=22  e(boot_df)=21           B_GAP=0
C controls(xctl) missing in 3 : e(N_clust)=23  e(boot_df)=25           C_GAP=3
D same sample, states dropped : e(N_clust)=23  e(boot_df)=22           D_GAP=0
E analysis weight 0 in 3      : e(N_clust)=23  e(boot_df)=25           E_GAP=3   <- NEW
F boottstat arm of A          : e(N_clust)=22  e(boot_df)=25           gap 4
```
A and B give the identical point estimate 1.5699826238, so it is the same estimation
sample. The printed header of A reads "N clusters = 22; N = 528 (96 obs dropped from
initial sample)"; the footer of the SAME run reads "... referred to the t(25)
distribution". invt(25,.975)=2.0595 against invt(21,.975)=2.0796.
Cause confirmed by reading: `wooldid_bootrun` ado:8911-8925 builds `smp` as
`(ifr==1) & !missing(clusv)` with no markout; ado:8954 counts `wdb_C` on it; and the
error message three lines below calls that same object "the estimation sample".
The E leg (a cluster with zero analysis weight) is new and is the most ordinary of the
three: a zero-population cell is not a data error.
SEVERITY unchanged at medium. The retracted SE-inflation sub-claim stays retracted -- I
did not re-test it and their own 3-seed measurement is consistent with a null.

## BT-03 -- the band's critical value never sees bootdf(): **PARTIAL**
repro: `repro_bt01_band.do` (cells B and C); sweep `r08_bandsweep.log` cells 7-10

It DOES fire -- but the finding's stated scope ("at modest cluster counts") is wrong.
The driver is the NUMBER OF HORIZONS IN THE FAMILY, because `q_lev(max_h |t|)` inflates
with m while `invt(G-1,.975)` does not:
```
r08 BCELL7  G=8   m=2  bandcrit 2.13844  invt(7) 2.36462  -> 2 of 2 horizons violate
r08 BCELL8  G=10  m=2  bandcrit 2.17012  invt(9) 2.26216  -> 2 of 2 horizons violate
r08 BCELL9  G=8   m=4  bandcrit 2.42834  invt(7) 2.36462  -> NO violation
r08 BCELL10 G=12  m=2  bandcrit 2.31140  invt(11)2.20099  -> NO violation
r03 S3      G=10  m=4  bandcrit 2.63837  invt(9) 2.26216  -> NO violation
r03 S5      G=10  m=4  bootdf(normal): bandcrit 2.69181 vs z 1.95996 -> NO violation
repro_bt01_band.do cell B  G=10 m=4: bandcrit 2.42613 vs invt(9) 2.26216 -> NO violation
repro_bt01_band.do cell C  G=8  m=2: bandcrit 2.06046 vs invt(7) 2.36462
                                     C_SUM inside_lo=2 inside_hi=2 narrower=2  (reproduces)
```
Corrections:
1. "the band's critical value carries no degrees-of-freedom correction" is DOCUMENTED,
   not a divergence: sthlp:554 says bootdf() "changes only the reference distribution ...
   the percentile and bias-corrected intervals are quantile-based and are unaffected".
   The finding concedes this in its own final sentence.
2. The two critical values sit within about 10 percent of one another in this corner and
   the ordering is not systematic: across my six small-G cells the band's critical value
   ranged 2.138-2.696 while invt(G-1) ranged 2.201-2.365. It is close to a coin flip,
   not a direction.
3. Their SECOND mechanism -- "the band's sd_h is computed on the LISTWISE-complete
   replications while the table's se_h is computed per column", offered as an independent
   cause of violations even when the band's critical value is the larger -- did NOT
   reproduce: in every cell of mine where bandcrit exceeded the pointwise critical value,
   no horizon violated. (The per-column/listwise split is real -- see SK-03 -- but it did
   not produce a band violation here.)
4. "the resulting band is not simultaneous at level()" is unmeasured on both sides.
REVISED SEVERITY: **low**. It is the same doc-scoping defect as BT-01, in a corner.

## BT-04 -- bootsave()'s _rep = 0 row for jointW_* (and advW*): **CONFIRMED (extended)**
repro: `repro_bt04_bootsave.do`   sentinels: `BT04_JOINT_FAIL=1`, `BT04_ADV_FAIL=1`

My own 22-state panel, cluster(state), boottstat jointtests reps(199):
```
e(joint_horiz) reported by wooldid = 0.6865671642
dump jointW_att at _rep==0         = -0.6261495363  == e(b_es_att)["pre3"] exactly
replication W_r range              = [0.4342965283, 3.56e+04]
p recomputed from the file exactly as the help documents = 1.0000000
```
The half the finding could only READ, I MEASURED. Second run with `adversarial`:
```
dump advW at _rep==0               = -0.6261495363  == e(b_es_att)["pre3"] again
replication advW range             = [0.0006373550, 4.1633e+02]
p recomputed from the file         = 1.0000000
```
So the advW leg (ado:14925-14931) fails identically. Cause confirmed by reading:
`wooldid_bootsavemk` writes `JO[1,k]` under the ACCUMULATOR's column index k, while
`jobs` (`wooldid_bootjobs`, ado:14666) has one column per (family, PRE-horizon). The
engine itself is unaffected -- it reaches the observed statistic by NAME through
`wooldid_bootjwobs` (ado:15310) -- so wooldid's own printed p is right.
REVISED SEVERITY: **high** (up from medium). The file is sold as self-sufficient ("so the
draws can be recentred from the file alone"); recomputing the documented statistic from it
returns a silently wrong number on BOTH legs, and the only tell is that a quadratic form
came out negative.

## BT-05 -- subgroup analytic SEs identical to the pooled one under cluster() >= i: **PARTIAL**
repro: `repro_bt05_balanced.do`

The exact identity is a PERFECTLY-BALANCED-PANEL artifact, and where it holds wooldid has
already detected it and withheld the affected inference.
```
BAL  24 states x 3 units x 8 years, no gaps, sg = mod(state,2), cluster(state):
     BAL_MAXRELDIF = 9.8828e-17          <- identical, as claimed
     BAL_CONTRAST_BLANKED = 1 ; e(wdid_sgdegen)=2 ; act=blank ; maxrho=2.0000
     and the run prints the named caution: "the analytic covariance of the subgroup
     effects is numerically rank-deficient in the contrast direction ... cancellation
     residue rather than an estimated variance ... use bootstrap() WITHOUT boottstat"
RAG  the same panel with 8 percent of rows deleted at random:
     RAG_MAXRELDIF = 2.9621e-02          <- the identity is GONE
```
Further controls (`r01_sgse.log`, `r06_sgse2.log`): a ragged panel with an UNBALANCED
subgroup gives ratios 0.9922 / 1.0158; a BALANCED panel with an UNBALANCED subgroup
already breaks the identity on the att row; `nofastcollapse` is bit-identical to the
default, so it is not a collapse artifact; `cluster(unit)` (below i) gives 1.096 / 1.125.
Why it is not a new defect: with a PSD covariance, "the contrast standard error is zero"
and "each subgroup's SE equals the pooled SE" are the SAME statement -- d'Vd = 0 with V
PSD forces Vd = 0, hence w'Vw constant on {sum w = 1}. wooldid documents the first (sthlp,
boottstat entry: "Under clustering at or above the i level the analytic pairwise-contrast
standard errors can be degenerate"), detects it (`wooldid_sgdeg`, ado:11092, ruled BS-15),
blanks it and names the remedy; known_bydesign section 1 carries it too. The finding
rediscovered its level-side face and did not connect the caution its own run printed.
REVISED SEVERITY: **low as stated.** But the investigation found the real defect beside
it -- SK-01.

---

# NEW FINDINGS

## SK-01 [CRITICAL, certain] The ruled subgroup-degeneracy guard misses on ordinary ragged panels, and wooldid then prints a subgroup-contrast p of 0.000 where the truth is 0.25
repro: `repro_sk01_sgcontrast.do`   sentinel: `SK01_FAIL=1`   sweep: `r07_sgdegsweep.log`
surface: `subgroups()` + `cluster()` at or above `i`, analytic path

On a perfectly balanced panel the analytic clustered contrast covariance is exactly
degenerate, `wooldid_sgdeg` sees an implied correlation of 1 (or a missing SE), blanks the
row and cautions. On a RAGGED panel the correlation lands just BELOW the ruled cut of
1 - 1e-2 and the row prints with no caution at all.
```
repro_sk01_sgcontrast.do  (20 states, 2-5 units each, 10 years, 8 percent of rows
                           dropped at random, 640 obs, i = state, cluster(state))
ANALYTIC  cluster(state) : b = -0.32158562  se = 0.0258871171  t = -12.423  p = 0.0000000
                           implied rho = 0.978609958   e(wdid_sgdegenact) = NONE
                           "p-value for joint test that all SGs have the same att: 0.0000"
BOOTSTRAP cluster(state), reps(399), no boottstat -- the remedy wooldid's own caution
          names for this design      : se = 0.310519665  t = -1.036  p = 0.250627
ROBUST                               : se = 0.273380161  t = -1.176  p = 0.240035
SK01_SE_RATIO_boot_over_analytic = 11.9951
```
Systematic, not one seed. `r07_sgdegsweep.log` runs 8 independent ragged designs
(G = 20/30, T = 8/10, 8 or 15 percent of rows dropped), 2 contrast rows each:
* 8 of 16 rows were blanked by the guard (correct behaviour)
* 8 of 16 rows LEAKED, and on EVERY one the analytic clustered contrast SE is 4.1x to
  14.5x smaller than BOTH the bootstrap and the robust SE:
```
 cell  ana_se   boot_se  rob_se   ana_p      boot_p  rob_p   rho      flagged
 1 R1  .01466   .21269   .29610   .0104      .788    .888    .98996   no
 1 R2  .04134   .22500   .32271   .0131      .438    .726    .94580   no
 3 R1  .02589   .31570   .27338   .00000000  .241    .240    .97861   no
 3 R2  .03302   .28022   .33143   .0000455   .508    .601    .98952   no
 5 R1  .02439   .31278   .25676   .00000000  .167    .111    .98425   no
 5 R2  .03622   .33925   .35519   .1225      .983    .869    .91350   no
 6 R2  .02788   .27920   .29694   .000000117 .508    .513    .97569   no
 7 R2  .04882   .19998   .37124   .1735      .722    .853    .94478   no
```
  5 of the 8 reject at the 5 percent level; 3 print p <= 5e-5.
An oracle borrowing nothing from another VCE: across the 8 cells the realized main-att
contrast estimates are {.0417, .2562, -.3216, .2640, .4102, .0263, .2734, .0929}, whose
sample SD is **0.226** -- an 8-draw Monte-Carlo estimate of the statistic's own sampling
SD. The leaked analytic SEs average 0.028; the bootstrap SEs average 0.28.
(Caveat: the 8 cells vary G and T, so 0.226 is a pooled order-of-magnitude figure.)
The DGP has a treatment effect common to every state, so the only true source of a
subgroup difference is cell composition, which is small; but the finding does not depend
on the null being exactly true -- an SE 12x below three independent estimates of the same
sampling SD is wrong whatever the true value is.

Why this is NOT the already-ruled material. known_bydesign section 1 rules analytic SEs at
cluster() >= i UNDERSTATED, measured range 0.66-0.75. On the very same runs the LEVEL SEs
sit inside that band (r06 P4: the clustered main-att SE is 0.55 of the robust one). It is
the CONTRAST that is off by an order of magnitude -- 0.09 of the robust one. If this were section 1's
understatement, the LEVEL standard error and the CONTRAST standard error would be
understated by the SAME factor. On the same run they differ by a factor of six: the level
SE is 0.55 of its robust counterpart and the contrast SE is 0.09 of its own. That is the
discriminating fact, and it is the reason this is a guard failure rather than a
rediscovery. `wooldid_sgdeg` exists precisely to withhold that row, and its own source comment anticipates this exact
failure: "A cut at 1 - 1e-3, tried first, caught the main row and left the pretreatment row
reporting a contrast standard error twenty times too small, which is the defect the ruling
names. The cut is a judgment call on a number the ruling did not fix; it lives in ONE local
below so it can be moved." The measurement says the cut that WAS chosen, 1 - 1e-2, is still
too loose once the panel is ragged: the pathology is continuous in rho and rho = 0.9786
already costs a factor of twelve.
Realistic: yes. `i` = state, `cluster(state)`, a subgroup split and a panel with holes is a
bread-and-butter applied specification, and the output is a false "these two groups differ,
p < 0.000001" with a confidence interval an order of magnitude too narrow.

## SK-02 [medium, uncertain on realism] The boottstat simultaneous critical value has no guard against a near-zero replication standard error
repro: `repro_sk02_bandcrit.do`; sweep evidence `r08_bandsweep.log` cells 1-6
surface: `simultaneousbands` + `boottstat` + `cluster()`

`wooldid_bootband`'s studentized path (ado:16776) drops a replication only when a horizon's
standard error is MISSING or NON-POSITIVE. A tiny positive one passes, and
`|b_r - b_obs| / se_r` goes to infinity. At `cluster()` == `i`, where the replication
analytic SEs are themselves unstable, the 95th percentile of the sup comes out far above
the nominal ~2.5, and the band is correspondingly vacuous, with no caution:
```
r08 BCELL1  G=20  log-normal  m=5  e(esband_crit_att) = 30.92        nband=297
r08 BCELL2  G=15  log-normal  m=5  e(esband_crit_att) = 1.17727e+14  nband=293
r08 BCELL3  G=20  GAUSSIAN    m=5  e(esband_crit_att) = 13.80        nband=295
r08 BCELL4  G=20  log-normal  m=3  e(esband_crit_att) =  4.64        nband=299
r08 BCELL5  G=12  log-normal  m=2  e(esband_crit_att) = 15.17        nband=298
r08 BCELL6  G=30  log-normal  m=5  e(esband_crit_att) =  6.77        nband=299
```
`repro_sk02_bandcrit.do` reproduces both ends of that range standalone:
```
A  G=15 log-normal m=5 : A_CRIT = 1.17727e+14  nband=293
   A_BAND pre3  [-23018976972227, 23018976972227]
   A_BAND post0 [-62463216459926, 62463216459928]      <- printed, with no caution
B  G=20 GAUSSIAN   m=5 : B_CRIT = 13.799  nband=295
   B_BAND pre3  [-4.69274, 4.78841]  halfwidth 4.74     <- about 5x the pointwise one
```
At A the 95th percentile being 1.18e14 means roughly 15 of the 293 complete
replications had a horizon standard error near machine zero. The pointwise percentile-t
rows beside it are much less affected because they use per-column 2.5/97.5 quantiles.
This is CONSERVATIVE, not anti-conservative: the band is uselessly wide rather than too
narrow. Marked `uncertain` because the 1e14 case needed log-normal errors with sigma 1.3;
but BCELL3 shows a critical value of 13.8 under plain GAUSSIAN errors at G = 20, which is
an ordinary design.

## SK-03 [medium, certain] The printed bootstrap standard error and e(V) are built on DIFFERENT completeness masks, so under bootdisplay(wald) lincom does not reproduce the printed row -- by up to 20 percent, and the gap concentrates in the subgroup block
repro: `repro_sk03_vdiag.do` (self-contained; identical to `r11_vdiag.do`); the two candidate
rules are measured in `r10_sdrule.do`
surface: `bootstrap()` + `subgroups()`, `bootdisplay(wald)`, `e(V)` / `lincom` / `test`

sthlp:554, verbatim: "Under every setting {it:e(df_r)} is aligned with the reference
distribution actually reported - the same degrees of freedom under {cmd:bootdf(clusters)}
or an explicit integer, and missing under {cmd:bootdf(normal)} so that the fallback is
chi-squared and z - so a {cmd:test} or {cmd:lincom} run by hand afterwards reproduces the
reported p-values rather than diverging from them - {bf:this holds under}
{cmd:bootdisplay(wald)}{bf:, which is what} {it:e(V)} {bf:and} {it:e(df_r)}
{bf:always carry}."  The run's footer adds that e(V) holds the bootstrap covariance
"(complete-replication listwise, e(boot_N_...))".
Measured (24 states, subgroups(sg3) with three levels, cluster(state), reps 399,
bootdisplay(wald); 368 of 399 replications are complete over all 24 coefficients):
```
WORST_PCTDIFF_PERCOL   = 0.00000   <- the PRINTED se is EXACTLY the per-column SD
WORST_PCTDIFF_LISTWISE = 24.75244  <- and is NOT the listwise SD
sqrt(diag(e(V)))       = the listwise SD, to 0.00000 percent on all 24 coefficients
worst coefficient pretreatment_SG0_att:
    printed se     0.269231486
    sqrt(V_kk)     0.215812605      -> 19.84 percent apart
    lincom _b[pretreatment_SG0_att] returns  se .2158126  t 0.81  p 0.427
                                             95% CI [-.2720803, .6208045]
second worst ES_att_SG0_pre1: printed .324726026 vs sqrt(V_kk) .293179186 (9.71 percent)
every one of the 24 differs (0.38 to 19.84 percent)
```
So under the one display the help singles out as internally consistent, the printed
standard error, statistic, p-value and interval do NOT reproduce from `e(V)`, and a
hand `lincom` diverges by up to a fifth.

WHERE the gap sits is the mechanism. Every one of the six MAIN-table coefficients is off
by 0.38-1.60 percent, which is the ordinary listwise-vs-per-column difference. The two
large gaps -- 19.84 percent on `pretreatment_SG0_att` and 9.71 percent on
`ES_att_SG0_pre1` -- are BOTH in the SG0 block. That is the block `wooldid_bootwire`'s own
comment (ado:9606-9625) flags as unguarded: "The later walks, and the `kk0`sg''/`kkP'/`kkS'
re-walks inside the subgroup loop, are unguarded". This is not a shift -- I tested for a
shift and there is none (see NEGATIVE RESULTS) -- it is column 2 of the subgroup results
matrices being filled from a different replication mask than the one `e(V)` is built on.
It bites hardest exactly where replications fail unevenly across coefficients, i.e. small
subgroups, which is also where a user is most likely to run the hand test.

---

# NEGATIVE RESULTS (probed hard, found nothing)

* **`wooldid_bootwire`'s unguarded subgroup walks.** The source comment (ado:9606-9625)
  warns that a later results-matrix walk that mis-lands "shifts every subgroup standard
  error one column with no error and no caution". I recomputed all 24 reported bootstrap
  standard errors from `e(bootdist)` name-for-name, main table and three subgroup blocks,
  and tested each against its own column AND against both neighbours (`r09_wirecheck.do`,
  `r10_sdrule.do`). No column is ever equal to a neighbour's SD; every one matches its own
  column's per-column SD to 0.00000 percent. The walk is correctly aligned on this build.
  A shift would have shown as an exact equality with the neighbour's recomputed SD.
* **`bootsave()`'s resample columns.** The help says the `_rep = 0` row of the `w_` columns
  holds 1/N under `bootalpha(size)` and 1/G under `bootalpha(equal)`. Measured
  (`r09_wirecheck.do`): N = 648, printed 0.00154320987654 = 1/648 exactly; G = 24, printed
  0.04166666666667 = 1/24 exactly. Both correct, 24 `w_` columns for 24 clusters.
* **`nofastcollapse` as an explanation for the BT-05 identity.** Bit-identical to the
  default on the same data (`r01_sgse.log` cases A and D, hex-equal SEs). Not a collapse
  artifact.
* **BT-03's second mechanism** (listwise band sd vs per-column table se producing band
  violations even when the band's critical value is larger). Did not reproduce in any of
  my cells: whenever bandcrit exceeded the pointwise critical value, no horizon violated.
