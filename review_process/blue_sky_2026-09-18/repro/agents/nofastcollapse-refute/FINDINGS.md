# SKEPTIC report — NFC-01 (nofastcollapse vs the default cell collapse, e(V))

Everything below was measured on data I generated myself, with do-files I wrote myself.
I read the prober's script and log but never ran them as evidence.

Workdir: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse-refute

---

## VERDICT on NFC-01: PARTIAL
The phenomenon is real and I reproduced it independently, but two load-bearing parts
of the claim are wrong: the TRIGGER ("whenever a continuous, coefficient-bearing
control varies within the cell") and the ATTRIBUTION ("not solely attributable to the
collapse ... the same underlying instability" in margins). Severity stays medium.

Repro (self-contained, verified): repro/repro_skeptic_NFC01.do -> repro_skeptic_NFC01.log

### 1. CONFIRMED — the gap is real and reaches printed digits
150 units in 3 treated cohorts + never-treated, 10 periods, i = cohort, control
xpct = 50 + 8*x1 (a percentage-scale covariate varying across units inside the
cohort-year cell), att att_i att_it esprelength(3) espostlength(3):

  default vs nofastcollapse   max rel e(V) 6.59e-06   max rel diag e(V) 6.59e-06
  noise floor (arm twice)     0.00e+00
  no controls                 5.65e-16 / 5.03e-16

The move is on the DIAGONAL, i.e. on standard errors, not an off-diagonal artifact.
The printed table differs (r6.log, repro_SORTORDER.log):

  pretreatment    default collapse                  nofastcollapse
     att       se .11064425 lb -.24658272        se .11064417 lb -.24658256
     att_it    se .11036045 lb -.24568926        se .11036009 lb -.24568854
     att_i     se .11036045 ub  .18727346        se .11036009 ub  .18727275
     att_i (main)  ub .72469403                     ub .72469402

So sthlp:338 ("Standard errors and variances do not amplify at all: they were
measured at about 1e-16 relative regardless of the scale of the predictions") and
sthlp:350 ("no digit wooldid displays is affected") are both contradicted, on my data
as well as the prober's. The project's own indifference threshold is rtol 1e-10.

### 2. TRIGGER — the claim's "whenever" is FALSE
Same panel, same control, only the request changes (r4b.log PART 1):

  att                                            1.92e-16 clean
  att esprelength(1) espostlength(1)             6.47e-16 clean
  att att_i att_it (no window)                   5.41e-16 clean
  att esprelength(2) espostlength(2)             5.35e-07
  att esprelength(3) espostlength(3)             1.51e-06
  att att_i att_it esprelength(3) espostlength(3) 6.59e-06
  att esprelength(5) espostlength(5)             4.65e-07

An event-study window of two or more horizons is needed as well, and the cells that
move are the PRETREATMENT rows (pretreatment_att_it, pretreatment_att_i,
pretreatment_att, pre2). The main att_i estimand moves only 2.7e-08. jointtests is
irrelevant (6.59e-06 with and without).

The trigger is also panel-shaped, not universal. On a 6,000-row administrative micro
panel (30 states x 25 individuals x 8 years, i = state, 25 rows per collapse cell) the
same kind of control is clean at every scale I tried (r1.log, r5.log, r6.log):

  controls(z) z~N(0,1)                                   0.00e+00
  controls(age) mean 35                                  6.84e-16
  controls(inc) mean 40,000                              6.83e-16
  ladder c + 8*z, c = 0 / 50 / 5000, i = state           7.04e-16 / 7.04e-16 / 4.23e-16
  same ladder with i = 3-value onset cohort              1.07e-15 / 7.11e-16 / 7.11e-16
  3 onset cohorts, esprelength(3) espostlength(3)        6.72e-16

A within-cell-varying coefficient-bearing control is NECESSARY (a cell-constant
control is clean at 2.81e-15) but NOT SUFFICIENT.

### 3. ATTRIBUTION — the "shared instability" story is REFUTED on my data
The claim says the collapse is only the larger of two perturbations of the same
instability, citing a row permutation under nofastcollapse that moved e(V) by
2.92e-07. On my failing panel the uncollapsed path is stable:

  row permutation, collapse OFF in both runs              1.32e-15
  row permutation, collapse ON  in both runs              6.59e-06
  add a pure CONSTANT to the control, collapse OFF        nothing above 1e-9
  the same constant shift, collapse ON                    up to 3.46e-06 on se/t/ci

The constant-shift test (r4b.log PART 4) is the cleanest: adding a constant to a
regressor cannot change any reported quantity in exact arithmetic, the fixed effects
absorb it. nofastcollapse is invariant; the default is not. The DEFAULT path is the
inaccurate side, not one of two equally-perturbed paths.

Also: the prober's own repro header states arm F (shuffle, collapse off) at "~7e-16"
while their logged value is 2.92e-07 — their two numbers for the same arm are nine
orders apart, and neither reproduces here.

### 4. MECHANISM AND REMEDY — verified with a one-line patch
wooldid_fcol keeps the cell's FIRST row and its RAW controls() values; it already
replaces the contreat() dose with the cell's weighted mean (ado:1894-1901, 5645-5690).
I copied the ado to patchbuild/wooldid.ado and added one line in the collapse gate so
controls() goes through the same dose-averaging branch (armed by a global, so one
binary runs both ways). r8.log:

  P0 patch inert, default vs nofastcollapse               6.59e-06
  P1 patch ARMED, default vs nofastcollapse               4.84e-16   <- remedy works
  P2 patch armed, row shuffle, collapse ON in both        1.38e-15
  P3 patch inert, row shuffle, collapse ON in both        6.59e-06
  P4 patched default vs nofastcollapse, on e(b)           1.31e-15 (estimates unmoved)

So the retained row's raw control value IS the mechanism, and carrying the cell's
weighted mean closes the gap to machine precision without moving a point estimate.
The remedy is validated empirically, not derived: I did not isolate what inside
margins' pretreatment/horizon read (ado:2089, 2153) turns a retained covariate value
into a 1e-6 variance move, and the prober's "ULP residue of the cancelling partial
sums" is not a sufficient account of a quantity that size.
Scope for a real fix: CONTINUOUS controls only. wooldid_fcol's averaging branch takes
a numeric varlist and would choke on a factor term, and factor controls do not need it
(the prober measured controls(i.region) clean at 5e-14). coarsecohortcontrols(c.x),
which the prober measured at 1.38e-06, presumably needs the same treatment; I did not
test it.

---

## FINDING SK-01 (SEEDED: this sharpens the prober's arm G, which already said
## "under the collapse the reported variance depends on the order the user's
## observations happen to be stored in" and named the first-row mechanism.
## Mine here: the OFF arm is order-STABLE, the printed-digit evidence, the
## sthlp:357 inversion, and the validated remedy.)
Under the default collapse a plain wooldid run is not reproducible across a re-sort of
the user's own data: the printed standard errors, t, p and confidence bounds move by
up to 6.6e-06 relative.

Repro: repro/repro_SORTORDER.do -> repro_SORTORDER.log. Two identical default commands
— no option typed at all — on the same dataset, sort unit t for one and gsort -t -unit
for the other:

  RUN 1  pretreatment att    se .11064425  lb -.24658272  ub .18749341
  RUN 2  pretreatment att    se .11064417  lb -.24658256  ub .18749324
  RUN 1  pretreatment att_it se .11036045  ub .18727346
  RUN 2  pretreatment att_it se .11036009  ub .18727275
  max relative move: e(V) 6.59e-06, e(wdidmainresults) 3.81e-06
  the same two sorts under nofastcollapse: 2.45e-14

Same root cause as NFC-01, and the prober asserted it too (arm G, 2.28e-07). What is
new here is the contrast: the OFF arm under the same permutation is 1.32e-15, so this
is the collapse's own defect rather than a shared instability; plus the printed-digit
evidence and the remedy. It is the sharper user-facing statement: it needs no option
and no cross-specification comparison, and it bites the user who merges or sorts
a file before estimating. It also reverses the reassurance at sthlp:357, which tells
the reader that a row permutation of the uncollapsed path moves a near-zero coefficient
about 19 times MORE than the collapse does: for variances, under the default, the row
permutation moves things eight to nine orders of magnitude more than the same
permutation does with the collapse off (6.59e-06 vs 1.32e-15 on e(V)).

Severity medium: no decision changes at 1e-6 on a standard error, but archived-result
reproduction and byte-comparison across runs — the stated purpose of the option — are
exactly the use cases the help addresses, and the default path is the one that fails.

---

## What I attacked and did NOT find (negative results)
- The simple margins read is exact, and that is as far as it goes. The hand-built
  equivalent of the collapse outside wooldid (plain reghdfe + margins [aw=], dydx(tf)
  subpop(), one row per cell with the summed weight and the first row's covariate)
  equals the exact delta-method oracle c'Vc to 0.00e+00, with and without the control,
  at mean |xb| from 3 to 3504 (r2.do/r2.log). NOTE the limit of this: r2 exercised the
  plain post-period read, which is ALSO clean inside wooldid (bare att, 1.92e-16). The
  read that fails is the pretreatment/horizon aggregation, which r2 never reproduced
  outside wooldid — so I have NOT shown that Stata's margins is innocent in general.
- Not an off-diagonal artifact. Every failing number above is on the diagonal of e(V).
- Not driven by the control's raw magnitude alone. controls(inc) with mean 40,000 on
  the micro panel is clean at 6.83e-16; the scale ladder (which moves only the level of
  the linear prediction) rises smoothly, 3.3e-07 at c=0 to 4.7e-05 at c=5000, and
  saturates — it is an amplifier, not the trigger.
- Not the point estimates. e(b) moves at 1e-15 to 1e-13 everywhere, as the help says.
- No level-reading margins call exists in the collapsed window: every call is a
  dydx/elasticity read, so the retained row's raw covariate cancels in the point
  estimate; there is no path where it would produce a grossly wrong estimate rather
  than a 1e-6 variance move.

## Unresolved, stated plainly
Why the cohort panel fires and the 6,000-row micro panel does not is NOT settled. I
varied the control's scale (0 to 5000), the i granularity (30 states vs a 3-value
cohort id), the number of onset cohorts (2 vs 3) and the event-study window (1/1 to
5/5) on the micro panel and it stayed at 7e-16 throughout. Some further property of
the failing shape is the third ingredient; I did not find it and I am not guessing.
