# SKEPTIC verdicts — "two-regime bootstrap unit rule under finer-than-i clustering"

All evidence below is from MY OWN data and MY OWN do-files.  I did not run the
original prober's scripts as evidence; I read them only to know what was claimed.
Pinned build asserted in every log (`which wooldid`).

## My DGP (used for FCM-1/2/3 and the calibration)
24 counties x 4 plants = 96 plants, T=2001-2008.  onset 2004/2006/2008 by plant
mod 4; mod4==0 never treated.  County effect heterogeneity ctau~N(1,.8^2), plant
effect ptau = ctau + N(0,.4^2).  `cluster(county)` is COARSER than `i`=plant.
**county enters the command ONLY through `cluster()`**, so moving a plant between
counties cannot move the point estimate — the cleanest possible isolation.

---

## FCM-1 — CONFIRMED (severity high, kept).  Mechanism arithmetic CORRECTED.

Repro: `repro_FCM1.do`.  Evidence: `r01_boot.log`, `r03_oracle.log`.

Dose on my own data (reps 250, rseed 606060, identical seed across k):

| k straddlers | PREDNEST | b(att_i)              | se(att_i) boot | se(att) boot |
|---|---|---|---|---|
| 0  | 1 | +1.10f479964ed18X+000 | 0.12595 | 0.12820 |
| 1  | 0 | +1.10f479964ed18X+000 | 0.09646 | 0.12844 |
| 8  | 0 | +1.10f479964ed18X+000 | 0.09760 | 0.13545 |
| 24 | 0 | +1.10f479964ed18X+000 | 0.09419 | 0.13773 |

The point estimate is **bit-identical** at every k.  One straddling plant out of 96
costs 23.4% of se(att_i); the next 23 buy nothing.  se(att) moves the other way.

**Which arm is right.**  I did not take either arm on faith.
* Data-implied answer: in this dataset the 24 county means of the plant-level
  treatment effects have SD 0.5950, so the between-county component of se(att_i)
  is 0.5950/sqrt(24) = **0.1215** (`r03_oracle.log`, EMP_SD_COUNTYMEAN_PTAU).
  The nesting arm reports 0.1257; the keep-i arm 0.0962.
* Independent hand-rolled pairs-cluster bootstrap run entirely OUTSIDE wooldid's
  engine (draw 24 counties with replacement, `joinby`, then `wooldid ... fastnose`),
  250 draws: ORACLE_NEST SD 0.12717, ORACLE_KEEP SD 0.09632.  wooldid reproduces
  BOTH rules faithfully (0.12570 / 0.09616).  The defect is *which rule it picks*,
  not the implementation of either.
* A 400-dataset Monte Carlo (`r02_mc.log`) gives SD(att_i - superpopulation ATT)
  = 0.1797 and SD(att_i - this dataset's own mean effect) = 0.0491.  The bootstrap
  is conditional on the observed sample, and this sample's between-county spread
  (0.595) happens to sit below the population value (0.833) — which is exactly why
  the nesting arm lands at 0.126 rather than 0.18.  The nesting arm is calibrated.

**CORRECTION to the claim's arithmetic.**  FCM-1 says the limiting SE ratio is
`sqrt(0.2325) = 0.482` (Bernoulli(1-1/e) variance over Poisson(1) variance).  That
drops the renormalisation: under keep-i a cluster drawn twice contributes *identical*
cells that the equal-cell weighting counts once, so the replication statistic is
`sum(I_c * thetabar_c) / sum(I_c)` — a ratio over the *realized* number of included
clusters.  Its variance is `sigma^2 (1-p)/(p G)` with p = 1-1/e, against `sigma^2/G`
under multinomial multiplicities, so the limit is `sqrt(0.3679/0.6321) = 0.763`.
Measured here: 0.09646/0.12595 = **0.766**, and 0.09616/0.12570 = **0.765**.  The
claim's own measurement (0.713) and mine bracket 0.763; 0.482 does not describe
anything either of us saw.

Implied size of a nominal 5% test on the keep-i SE, taking the nesting arm as
calibrated: 2*Phi(-1.96*0.766) = **13.3%** (the claim says ~16% off its own ratio).

**Not a doc-vs-code divergence.** sthlp:532 states the rule as a property of where
`cluster()` sits *relative to i* — "where the cluster is finer than i, **or cuts
across it**".  One straddling unit does make the two partitions cross-cutting, so
the code does what the help says.  This is a design defect (a per-unit property
collapsed to one global flag, resolving the trade-off the wrong way when 95 of 96
units nest), not a documentation mismatch.  The harm is one-directional: under a
cluster coarser than i every cell lies wholly inside one cluster, so duplication
cannot move a cell mean and keep-i can only *remove* variance.  It is always
anti-conservative.

## FCM-2 — CONFIRMED (severity high).  One sub-claim is an overread.

Repro: `repro_FCM2.do`.  Evidence: `r00_struct.log`, `r01_boot.log`, `r04_record_df.log`.

DOOR 1 (missing outcome).  Plant 1 shuts down: y missing from 2005 on.  Arm B also
records it in county 2 on exactly those rows.

| | N | b(att_i) | analytic se | NEST on smp | NEST on rows with y | boot se(att_i) | boot se(att) |
|---|---|---|---|---|---|---|---|
| A | 764 | +1.118bfbe772b97X+000 | +1.6f63d274ae0e4X-005 | 1 | 1 | 0.12570 | +1.075bb5bcac12dX-003 |
| B | 764 | +1.118bfbe772b97X+000 | +1.6f63d274ae0e4X-005 | 0 | 1 | 0.09616 | +1.075bb5bcac12aX-003 |

Identical to every hex digit except the bootstrap SEs of the three cell-weighted
estimands.  se(att) agrees to 3 ulps, so the two runs saw the same cluster draws.

DOOR 2 — **NEW, not in the original claim**: a MISSING WEIGHT opens the same door.
`[aw=wt]` with wt missing on the phantom rows: se(att_i) 0.12638 -> 0.10528
(`r04_record_df.log`, WMB0/WMB1).  So the door list is missing y, missing control
(the prober's), and missing/zero weight (mine).  A missing **i** does NOT open it:
`gegen ... by(i)` puts the missing-i rows in a group of their own (r00, D5_PRED=1).

Direction: `e(sample)` is a subset of `smp` by construction, so straddling on
e(sample) implies straddling on smp.  The misdetection is one-directional and always
toward the smaller SE.  Verified, and it is the bad direction.

**Overread to correct.**  The claim says the code comment at ado:8975-8985 "claims
the `if smp` placement handles exactly this case; it does not."  That comment says
the placement beats evaluating on the *raw frame*, which is true.  It never claims
`smp` is the estimation sample.  The genuine divergence is with sthlp:532 ("Where
the cluster nests i, each drawn copy of a cluster becomes a fresh set of units") —
on arm B's estimation sample the cluster does nest i and the copies do not become
fresh units.  That part stands.

## FCM-3 — CONFIRMED (severity medium, kept).

Repro: `repro_FCM3.do`.  Evidence: `r04_record_df.log`.
Walking every name in `e(macros)` and `e(scalars)` on both arms: 0 names containing
`nest`, `reid` or `regime` (NESTRUN_RECORD_HITS=0, KEEPRUN_RECORD_HITS=0).  A full
`verbose` run prints nothing about the regime (only the reps-resolution caution).
The two runs' printed tables differ only in the standard-error column.

## FCM-4 — CONFIRMED, severity revised UP to high; framing corrected.

Repro: `repro_FCM4.do`.  Evidence: `repro_FCM4.log`, `r04_record_df.log`, `r06_deadse.log`.

Mild case (one dead cluster of 24): `e(N_clust)=23`, `e(df_r)=e(boot_df)=23`.
Reproduced.  Consequence by itself is tiny (t(23)=2.0687 vs t(22)=2.0739).

Severe case, and this is what raises the severity.  20 counties in the file; the
outcome is collected in only 5 of them (an ordinary merged administrative file):
* printed header: `N clusters = 5; N = 160 (480 obs dropped from initial sample)`
* `e(N_clust)=5`, `e(df_r)=e(boot_df)=19`
* `bootdisplay(wald)` therefore uses t(19)=2.0930 where the live cluster count gives
  t(4)=2.7764 — the interval is **25% too narrow** (ratio 0.7539), and any hand
  `test`/`lincom` on `e(V)` inherits the wrong reference.
* control with the same estimation sample and the dead clusters dropped: `e(df_r)=4`,
  identical point estimates (att_i = 1.0836262589 both).

**Framing correction.** sthlp:532 documents bootdf(clusters) as "G being the number
of *resampling* clusters", and a dead cluster IS a resampling cluster in this code,
so help and code are self-consistent on their own terms.  The defect is that a
cluster contributing nothing to any replication is a resampling unit at all, and
that the only cluster count a user can see (`e(N_clust)`, and the printed header)
contradicts the df in the same output.

**The claim's "NOT MEASURED" item, measured**: see the r06 table below — the wasted
draws DO move the SE (mean ratio 1.094 over three seeds at 500 reps) and, more clearly,
destabilise it.  Direction and instability are established; the magnitude is indicative
only, because the design has 5 live clusters and a bootstrap SE there is noisy.

## NEW (mine) — the same dead clusters silence BOTH few-cluster cautions

Repro: `repro_FCM4.do` block A vs block B.  `wooldid_bootctl` builds its own sample
concept at ado:6459-6462 (`((if)==1) & !missing(clus)`, no outcome markout) and every
cluster-count screen hangs off it: the two-cluster refusal (ado:6464), the
"fewer than 10 clusters" caution (ado:6479) and the "fewer than 6 treated clusters"
caution (ado:6498).

In the 15-dead-of-20 run **neither caution prints**.  In the control run with the
same estimation sample and the dead rows physically removed, **both** print:
"only 5 clusters are available for the bootstrap" and "only 5 treated cluster(s)
are available for resampling ... the literature finds this bootstrap badly sized".

So the single run that most needs the few-cluster warning is the one that does not
get it, and it gets a 25%-too-narrow wald interval in the same breath.  Padding the
file with clusters that contribute nothing is enough to turn the guidance off.

---

## Measured answers to things the claims left open

### FCM-4's "NOT MEASURED": what the wasted draws do to the SE (`r06_deadse.log`)
500 replications, three seeds, the 15-dead-of-20 dataset vs the same estimation sample
with the dead clusters dropped.  se(att_i):

| seed | padded (20 resampling units, 5 live) | dropped (5 units) | ratio |
|---|---|---|---|
| 11 | 0.40684 | 0.41155 | 0.989 |
| 22 | 0.44773 | 0.39996 | 1.119 |
| 33 | 0.47193 | 0.40140 | 1.176 |
| mean | 0.44217 | 0.40430 | **1.094** |

So the dead resampling units do move the SE — about +9% on average here — and, more
telling, the padded arm's spread across seeds is 16% against the clean arm's 3%.  Drawing
20 slots when only 5 can contribute makes the live-cluster count of a replication
Binomial(20, .25) instead of a fixed 5, which fattens the replication distribution and
destabilises it.  Direction: conservative on the SE, anti-conservative on the df.
The padded arm also loses replications the clean arm does not (1 and 2 of 500 at two seeds,
draws that caught no live treated cluster).

### NEW negative result: the keep-i regime IS calibrated in its own regime (`r05_finer.log`)
§7 of known_bydesign asks for "any finer-than-`i` arm that FAILS coverage".  This one does
not.  i = cohort group (12 groups of 8 firms), cluster(firm) (96, genuinely finer than i),
treatment effects iid ACROSS FIRMS so a firm-level cluster bootstrap can capture all of the
sampling variance.  Regime predicate: PREDNEST=0 (keep-i), as intended.

* 300-dataset Monte Carlo: SD(att_i − superpopulation ATT) = **0.08290**
* mean bootstrap se(att_i) over 12 fresh datasets, reps(150) = **0.08448** (+1.9%)
* mean analytic clustered se(att_i) over the same 12 = 0.08314

The rule is right where the design intends it.  That sharpens FCM-1 rather than softening
it: the defect is not the keep-i rule, it is the global flag that applies it to a sample
where 95 of 96 units nest.

### Leads I chased and did NOT file
* **String cluster variable** (state FIPS as a string — very common).  Refused informatively
  on the plain and the bootstrap paths: "string variables not allowed in option cluster()",
  r(109).  Clean.
* **A straddle carried on rows whose `i` is MISSING.**  No phantom: `gegen ... by(i)` puts
  the missing-`i` rows in a group of their own, so the predicate still sees nesting
  (r00_struct.log, D5_PRED=1).
* **An analytic-SE lead I refuted myself.**  In `r07_anaclus.log` wooldid's analytic
  cluster(county) se(att) came back bit-identical (to 1-3 ulps) for 1, 2, 8 and 24
  straddling plants while a benchmark `reghdfe` clustered SE on the same data moved 12%.
  That looked like a VCE that could not see the row-level cluster assignment.  `r08` and
  `r09` refute it: the SE does respond to content — moving the same plant from 2003 instead
  of 2005 (1.6591a07... vs 1.659b267...), moving a different plant (1.65700d6...), moving it
  to a brand-new 25th county (1.654ba5d..., e(N_clust)=25), and half the plants straddling
  (1.4586881...) all give different numbers, and cluster(plant)/robust/region give their own.
  The k-ladder coincidence is a property of that particular symmetric ladder (same cohort,
  same within-county position, each plant moving one county up) in a saturated design, not a
  defect.  Reported here so nobody re-chases it.  Related known ruling: §1 of
  known_bydesign — analytic clustered SEs at or above `i` are understated and not to be
  re-filed.

### Repro files, all executed, all ending in their sentinel
`repro_FCM1.do` (REPRO_FCM1_DONE), `repro_FCM2.do` (REPRO_FCM2_DONE),
`repro_FCM3.do` (REPRO_FCM3_DONE), `repro_FCM4.do` (REPRO_FCM4_DONE).
Severity of the new caution-silencing finding: **medium** on its own; it co-occurs with
FCM-4's wrong df in the same run, and together they are what makes that run high.
Boundary stated honestly: `bootdf(clusters)` is the DEFAULT whenever `cluster()` is given
(ado:423, and r04's NESTRUN/KEEPRUN passed no `bootdf()` and still show `boot_df=23`), but
the DEFAULT `bootdisplay(percentile)` table is rank-based and unaffected by the df.  The
harm lands on `e(df_r)`, `e(boot_df)`, `bootdisplay(wald)`, and on every consumer the
run's own footer names as reading `e(V)` with `e(df_r)`: `test`, `lincom`, esttab,
coefplot, honestdid.
