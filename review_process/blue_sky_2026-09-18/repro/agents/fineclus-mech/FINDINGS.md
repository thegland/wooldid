# FINDINGS — fineclus-mech (two-regime bootstrap unit rule)

Agent label: FCM.  Method: code-mechanism probing.
Pinned build asserted in every log (`which wooldid` -> .../bluesky/build/wooldid.ado).

## THE MECHANISM (as read on the pinned build)

wooldid_bootrun (ado:8887).

    ado:8911  tempvar smp
    ado:8912  if `"`ifr'"' != "" qui gen byte `smp' = ((`ifr') == 1)
    ado:8924  qui replace `smp' = 0 if missing(`clusv')        <- the WHOLE screen

`smp` is (a) the sample the regime predicate is evaluated on and (b) the rows the
physical re-identification leg keeps (`qui keep if smp`, ado:9345).  It applies the
user `if` and drops missing cluster ids and NOTHING ELSE — no markout for missing y,
missing controls(), missing dose, zero weights.

    ado:8986  local wdb_nest 0
    ado:8987  if "`reid'" != "" & "`clus'" != "" {
    ado:8989    qui gegen double `wdbcmin' = min(`clusv') if `smp', by(`i')
    ado:8990    qui gegen double `wdbcmax' = max(`clusv') if `smp', by(`i')
    ado:8991    cap assert `wdbcmin' == `wdbcmax' if `smp'
    ado:8992    local wdb_nest = (_rc == 0)
    ado:9348  if `wdb_nest' qui gegen double __wdbnewi = group(__wdbns `i')
    ado:9349  else          qui gen   double __wdbnewi = `i'

`reid` is set at ado:378 whenever bootstrap(standard) and (att_it|att_i|att_itime
requested, or boottstat with cluster()).

Two structural facts follow, both measured below:

1. THE FLAG IS GLOBAL, THE PROPERTY IS PER UNIT.  `assert` is one all-or-nothing
   test over every unit, so one non-nesting unit puts all units under keep-i.  The
   code's own comment (ado:9331-9334) says what that costs: "keeping i under nesting
   would make a twice-drawn cluster's cells identical duplicates that the equal-cell
   weighting counts ONCE, collapsing the composition variance of att_i toward zero."
2. e(sample) is a STRICT SUBSET of smp, so the predicate can only err one way —
   report NOT-nesting where the estimation sample does nest — and that is the
   anti-conservative direction.

The comment at ado:8975-8985 claims the `if smp` placement handles exactly this:
"Evaluated on the raw frame instead, a cluster that nests i within the sample but
straddles two i's only through rows the sample excludes would be misclassified as
non-nesting."  smp is not the sample.

## FCM-1 (mine).  One reclassified unit out of 96 shrinks the att_i bootstrap SE ~29%

repro: p08_dose.do  (pair also in p01_regime.do + p04_regime_measure.do)

Panel: 24 industries x 4 firms, T=8, i=firm, cluster(ind), industry-level treatment
heterogeneity.  Firms nest in industries except for k firms reclassified into a
neighbouring industry halfway through the panel.

  variant   straddling units   regime   BOOTSE(att_i), reps(300), same rseed
  A (k=0)   0 of 96            nest     0.14110902275043
  B (k=1)   1 of 96            keep-i   0.10057797358284

Ratio 0.713 in the SE, 0.508 in the variance.  att_i = 1.01724741722359, N = 768,
N_clust = 24 in both.  Nothing printed, nothing posted.

## FCM-2 (seeded).  The regime is decided on rows not in the estimation sample

repro: p01_regime.do then p04_regime_measure.do / p06_estimand_split.do

  variant  perturbation                                          PREDNEST  N    BOOTSE(att_i)
  C0       firm 1's y missing from 2005 on                        1(nest)  764  0.14122699826312
  C        same, AND its industry code changes on those rows      0(keep)  764  0.10043761347134

Identical estimation samples (N=764), identical point estimate
(att_i = 1.01054627838245 to every digit), identical ANALYTIC clustered SE
(0.04030370871346), opposite regimes, and a bootstrap SE ratio of 0.7112.
The only difference between the two .dta files is the value of `ind` on FOUR
rows whose y is missing -- rows the regression never sees and e(sample)
excludes.  reps(300), rseed(4242), so both runs saw the same cluster draws.

## FCM-3 (mine).  Nothing records which regime ran

No display line, no e() macro, no verbose line.  Contrast permonsetcohorts, where
the author ruled 2026-09-16 that no display line was needed BECAUSE e(perm_irebuilt)
is the record.  Here the number moves ~30% and there is no record at all.

(status of the remaining probes: NOTES.md)

## WHICH MARKOUTS VOTE ON THE REGIME AND WHICH DO NOT

The screen at ado:8924 removes exactly one thing from the detection sample:
rows with a MISSING CLUSTER VALUE.  Measured (p02_edges.do):

  T6  firm 3's INDUSTRY code is missing from 2005 on
      -> those rows leave smp, no phantom straddle, nesting regime,
         BOOTSE(att_i) = 0.13549772259762 (the nesting level)
  C   firm 3's OUTCOME is missing from 2005 on and the industry code
      changes on those same rows
      -> those rows stay in smp, phantom straddle, keep-i regime,
         BOOTSE(att_i) = 0.10043761347134

So the one markout the screen does implement is the one that could not
have produced a phantom anyway (a missing cluster id cannot straddle),
while every markout the estimator actually performs -- missing y, missing
controls(), missing dose, zero weight, the event-study window -- is left
to vote.

## BOUNDARY: cluster() == i

p02_edges.do T2: `wooldid y firm year onset, att_i cluster(firm)
bootstrap(standard) reps(100)` -> RC=0, N_clust=96,
BOOTSE(att_i)=0.09122637551236.  min==max by(i) holds trivially, so the
NESTING arm fires and `group(__wdbns, i)` degenerates to the drawn id
alone, exactly as the ado:8983-8984 comment says.  CLEAN.

## MULTI-WAY cluster()

p02_edges.do T1/T1B/T1C: `cluster(ind year)` estimates fine analytically
(att_i se .02982697) and is refused r(198) by BOTH bootstrap legs, with the
message the help promises (sthlp:532 "multi-way cluster() ... the engine
resamples one unit, so there must be exactly one"; screen at ado:6427).
The `gegen min(clusv)` in the detection block would have died on a two-
variable clusv; the screen fires first.  CLEAN.

## CONSTANT cluster()

p02_edges.do T5: r(459) "bootstrap() requires at least 2 clusters in the
estimation sample."  CLEAN.

## NOT A DEFECT (recorded so nobody re-files it)

p02_edges.do T3/T3B/T4 and most of p07_boottstat.do are INVALID TESTS of
mine, not wooldid defects: their panel assigns onset by firm, so
`wooldid y cgrp ...` has a non-constant onset inside i and is correctly
refused r(459).  Redone on a valid design in p11_finer.do.

## THE DECISIVE CONTROL (p06_estimand_split.do, reps(120), rseed(31337))

Same data pair C0 / C, one run each, all four estimands requested together:

  estimand         C0 (nesting)          C (keep-i)            ratio
  main_att         0.12587695969826      0.12587695969826      1.000000  <- BIT-IDENTICAL
  main_att_it      0.12587695969826      0.09632508066525      0.765
  main_att_i       0.12579454738696      0.09837740349041      0.782
  main_att_itime   0.12579454738696      0.09837740349041      0.782

`att` -- the observation-weighted estimand, for which a twice-drawn cluster's
copies are equivalent whether or not they are re-identified -- does not move
by one bit.  So the two runs saw exactly the same cluster draws and the
movement is entirely the unit rule.  Exactly the three CELL-WEIGHTED
estimands move, which are exactly the estimands the re-identification leg
exists to serve.  Same run also prints:

  ROWS IN e(sample) BUT NOT IN smp = 0     <- e(sample) is a subset of smp,
  ROWS IN smp BUT NOT IN e(sample) = 4        so the misdetection is
  NEST ON smp        = 0                      ONE-DIRECTIONAL: it can only
  NEST ON e(sample)  = 1                      report keep-i where the
                                              estimation sample nests, and
                                              that direction understates.

## THE PHANTOM ALSO FIRES THROUGH A MISSING CONTROL (p02_edges.do T9/T9B)

Same panel, `controls(x)`, firm 3's x missing from 2005 on:

  T9B  x missing 2005+, industry code unchanged   N=764  BOOTSE(att_i)=0.13822686104741
  T9   x missing 2005+, industry code changes on
       exactly those rows                        N=764  BOOTSE(att_i)=0.10399078061065

ratio 0.752, same point estimate (1.01945199245315) and same N.  So it is not
about the outcome specifically: EVERY markout the estimator performs is left
to vote on the regime.

## DOSE RESPONSE (p08_dose.do, reps(120), rseed(8081))

  k   straddling units   PREDNEST   se(att)              se(att_i)
  0   0 of 96            1          0.13112922202932     0.13203513225887
  1   1 of 96            0          0.13248661337362     0.10490369831254

att moves +1.0% (the data really did change on 4 rows); att_i moves -20.6%.
Opposite direction, twenty times the size.  The ANALYTIC clustered SE in the
same comparison (p04, A vs B) moves the sensible way too: 0.04087 -> 0.04499,
UP 10%.  Only the bootstrap's cell-weighted rows go down.

## INDEPENDENT ORACLE (p04_regime_measure.do)

A hand-rolled pairs-cluster bootstrap on the C0 data, 300 draws, seed 90210,
`bsample, cluster(ind) idcluster(__ns)` then `gegen __newi = group(__ns firm)`
and `wooldid y __newi year onset, att_i fastnose`:

  ORACLE_NEST  n=300  mean=1.01364576854538  SD=0.1317

i.e. the nesting neighbourhood (wooldid's C0 = 0.1412), not the keep-i one
(wooldid's C = 0.1004).

## CLEAN (probed hard, nothing found)

p11_finer.do, on a design where the onset IS constant within i
(96 firms in 12 cohort groups of 8, onset by group):

  geometry                         RC  N_clust  e(df_r)  se(att)   se(att_i)
  cluster(firm)  96, finer than i   0     96       95     0.09042   0.07575
  cluster(pair)  48, finer than i   0     48       47     0.11220   0.06916
  cluster(cgrp)  12, == i           0     12       11     0.19885   0.19581
  cluster(rowid) 768, one per obs   0    768      767     0.13166   0.15271
  no cluster() at all               0      -        -     0.14813   0.16307
  boottstat + cluster(firm/pair/cgrp)   all RC=0
  boottstat + fe(firm) + cluster(firm) + 3 one-row firms (reghdfe
      singletons): RC=0, N=744, N_clust=93 -- the P7 run-time assert
      (ado:9386-9414) did NOT fire a false positive there.

The observation-level legs giving a LARGER att_i SE than the firm-clustered
one is not a defect: with no firm fixed effect the firm-level shock
differences out within a drawn firm but not across independently drawn rows.

## THE SUBINSTR + P7 ASSERT ON THE NEW FINER-THAN-i boottstat PATH: CLEAN

p11_finer.do C_BTS_firm vs E_NOSUB (i=cgrp 12, cluster(firm) 96 -> keep-i
regime, boottstat, reps(50), rseed(77)).  e(V) under boottstat is still the
bootstrap-SD covariance, so the control is only visible in the STUDENTIZED
interval:

  run                            att_i studentized 95% CI        width
  C_BTS_firm (substitution live) [ .9158569 , 1.2092969 ]        0.29344
  E_NOSUB    (WOOLDID_BOOTNOSUB) [ .95269294, 1.1714977 ]        0.21880

width ratio 0.746, against the predicted 1/sqrt(2 - 1/G) = 0.709 at G=96
(reps(50), so the quantile is noisy).  The control reproduces the 2 - 1/G
replication-SE inflation on demand and the live path does not, so the
`subinstr' at ado:9358 IS landing on the keep-i leg.  The P7 assert did not
fire a false positive on the recommended coarse-i + fe(firm) pairing with
three reghdfe singleton firms either (D_SINGLETON_BTS RC=0, N=744,
N_clust=93).

Caveat worth recording: WOOLDID_BOOTNOSUB *skips* the P7 assert (ado:9398
gates on `wdb_sub'), so the package's own negative control cannot exercise
the check that exists to validate the substitution.  The assert is therefore
externally untestable without editing the ado.

## FCM-4 (seeded).  A cluster with no rows in e(sample) is still a resampling
## unit and still sets bootdf(clusters)

`wdb_C` (ado:8954) is counted on `smp`, which does not markout missing y.
p02_edges.do T10: industry 5's outcome is missing for all of its rows.

  plain run        N=736  e(N_clust)=23
  bootstrap run    N=736  e(N_clust)=23  e(df_r)=23   <- should be 22 on the
                                                        documented G-1 rule

So the posted e() is internally inconsistent, and ~1 draw in 24 of every
replication is spent on a cluster that contributes nothing.

## FCM-4 AT SCALE (p11_finer.do F)

i = cgrp (12 cohort groups of 8 firms), cluster(firm) = 96, and the eight
firms of cohort group 1 have no outcome at all:

  plain      N=704  e(N_clust)=88
  bootstrap  N=704  e(N_clust)=88  e(df_r)=95   (printed side by side:
                                                 "e(N_clust)-1 = 87 vs
                                                  e(df_r) = 95")

Eight of the ninety-six resampling units contribute nothing to any
replication, and bootdf(clusters) is eight degrees of freedom off the rule
the help states.

## THE DOSE RESPONSE IS A CLIFF (p08_dose.do, reps(120), rseed(8081))

  k   straddling units   PREDNEST   se(att)     se(att_i)
  0   0 of 96            1          0.13113     0.13204
  1   1 of 96            0          0.13249     0.10490
  5   5 of 96            0          0.13621     0.10528
 20  20 of 96            0          0.12888     0.10489

se(att) drifts up smoothly with k -- that is the data really changing.
se(att_i) STEPS DOWN at k = 1 and then sits FLAT to three decimal places
(0.10490 / 0.10528 / 0.10489): one straddling unit out of ninety-six buys the
whole 20%, and the next nineteen buy nothing at all.  That is the signature of
a GLOBAL flag on a PER-UNIT property.  A per-unit rule would have cost about
1/96 of the effect at k = 1.

## ORACLE, BOTH ARMS (p04_regime_measure.do) -- which number is the right one

Hand-rolled pairs-cluster bootstraps on the C0 data, 300 draws each, seed
90210, nothing of wooldid's bootstrap engine involved:

  rule                                                        SD of att_i
  ORACLE_NEST  bsample, cluster(ind) idcluster(__ns) then
               gegen __newi = group(__ns firm)                0.1317
  ORACLE_KEEP  bsample, cluster(ind) idcluster(__ns) then
               gen __newi = firm                              0.1026

against wooldid's own reps(300) numbers on the SAME estimation sample:
  C0 (nesting arm) 0.14123      C (keep-i arm) 0.10044

Both arms reproduce independently.  The textbook cluster bootstrap of an
equal-cell-weighted estimand is ORACLE_NEST (a twice-drawn cluster must
contribute its cells twice), which is also what ado:9326-9334 says the leg
exists to do.  So 0.141 is the intended answer on C's estimation sample and
0.100 is the one that is reported.

=======================================================================
## SUMMARY TABLE

  id       title                                          sev   conf     seeded
  FCM-1    global all-or-nothing regime flag on a          high  certain  no
           per-unit property; 1 straddler in 96 cuts
           se(att_i) ~20-29%
  FCM-2    the regime predicate is evaluated on `smp',     high  certain  yes
           which is NOT the estimation sample; cluster
           values on rows the regression never sees
           move se(att_i) ~25-30%
  FCM-3    nothing anywhere records which regime ran       med   certain  no
  FCM-4    a cluster with no rows in e(sample) is still    med   certain  yes
           a resampling unit and still sets
           bootdf(clusters): e(df_r)=95 vs e(N_clust)=88

FCM-1 and FCM-2 have ONE root cause and TWO independent remedies:
  (a) make the rule per unit -- the code's own comment (ado:9331-9334) says
      keep-i is wrong for a unit that nests, and the predicate is already
      computed per unit before `assert' throws the detail away;
  (b) compute the predicate on the estimation sample.  The observed run's
      e(sample) is LIVE at ado:394 (the shell has just done `qui est sto
      `wdbobs''), so the shell could pass a marker down at zero cost.  That
      alone kills FCM-2 outright even without fixing FCM-1.

## REPRO STATUS AT WRITE-UP TIME
repro_FCM2.do (the self-contained one) had printed its whole setup block
live -- 0 rows in e(sample) outside smp, 4 rows in smp outside e(sample),
NEST on smp = 0, NEST on e(sample) = 1, and the phantom=0 analytic line
(att_i 1.01054627838245, se 0.04030370871346, N=764) -- but its two
reps(200) bootstraps had not finished when this was written (14-19
competing Stata processes on a 10-core box all session).  Its lines 96-113
are the same code that produced p04_regime_measure.do's reps(300) numbers,
which DID complete: 0.14122699826312 (phantom=0) vs 0.10043761347134
(phantom=1).

NEVER RUN: p05_finer_coverage.do (size of the finer-than-i arm), and
p10_legcross.do (does att's reported SE move between the count leg and the
re-id leg depending on which other estimands were requested?).
