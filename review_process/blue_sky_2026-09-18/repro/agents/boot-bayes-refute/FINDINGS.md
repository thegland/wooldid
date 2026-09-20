# SKEPTIC PASS on the boot-bayes findings (bb-01 / bb-02 / bb-03)
Agent: boot-bayes-refute.  Every number below is from a do-file in THIS directory, run on
data THIS directory generates.  I read the prober's scripts for design only; none of their
files is evidence for anything here.

## VERDICT SUMMARY
| id    | verdict   | severity (theirs -> mine) | why |
|-------|-----------|---------------------------|-----|
| bb-01 | PARTIAL   | high -> medium            | every number reproduces on my own data and the SCOPE IS WIDER than they state (it also fires at cluster() coarser than i, and with non-singleton cells) - but the predicate is "cluster() nests i", not "cluster() == i"; the DEFAULT ANALYTIC se is exactly as blind in the same arms, so this is not a Dirichlet-engine property; and the 4.4x uses bootstrap(standard) as truth |
| bb-02 | CONFIRMED | low -> low                | reproduces exactly; the value is internally right and only the help sentence is wrong |
| bb-03 | PARTIAL   | medium -> low             | (a) reproduces for att but is sec-1 for the analytic leg and is PRINTED AT sthlp:532 for the bootstrap legs; (b) direction reproduces, magnitude does not (3.3pp not 16.7pp, discordance 3-1 not 5-0), and the interval is NARROWER, not shifted |

New findings of my own: **N1** e(gteffects) cell CIs cover 45.0% at nominal 95% at cluster()
at-or-above i (nominal 0.950 at the recommended finer clustering) - `r11`/`r12`; **N2** the
documented 0.66-0.75 SE-understatement range is the no-unit-heterogeneity case - `r09`/`r10`;
**N3** the bayes x cell-weighted gate is regime-blind - `r02`.

---
## bb-01 — PARTIAL.  Reproduces; scope is wider; mechanism and framing need correcting.
repro: `r01_saturated.do` (own panel, their shape), `r02_scope.do`, `r05_designD.do`

### What reproduces (all on my own 28-state x 14-year panel, cluster(st), i=st, reps(200))
- att and att_it are bit-identical in b AND analytic se: b `+1.f234682b54b69X-002`,
  se `+1.2fdb96bdb85e9X-004` for both.
- bootstrap(standard) reports se(att) == se(att_it) bitwise (0.11032743 both).
- bootstrap(bayes)+safety(off): se(att_it)/se(att_it | standard) = **0.5003** with NO
  cross-cluster heterogeneity, **0.1671** with a cluster-level shock on the effect (6.0x).
- bayes se(att_it) is **bitwise identical across the two arms**: `+1.c432c1b72effdX-005`
  in both, relative move 0.0000e+00.  The exact invariance is real.

### CORRECTION 1 (the one that changes the reading): the ANALYTIC se is exactly as blind
In the same two arms the reported ANALYTIC clustered se is also **bitwise identical**:
`+1.2fdb96bdb85e9X-004` for het=0 and for het=1.5, **and for att as well as att_it**.
The bootstrap(standard) se moves by +199% between those arms; the analytic one by 0.
So the property "the standard error does not respond at all to cross-cluster
treatment-effect heterogeneity" is NOT a property of the Dirichlet engine.  It is a
property of the geometry: with i at the unit level and cluster() nesting i, the treated
cells are fitted by their own dummies, a cluster-constant shock on the post rows is
absorbed exactly, no residual moves, and every reweighting-based variance - analytic
delta-method included - is invariant to it.  bootstrap(standard) is the only one of the
three that sees it, because it is the only one that resamples the cells.
Consequence for severity: this is the sec-1 "analytic SEs are understated at cluster()
at or above i" phenomenon in its extreme form, reached by the DEFAULT inference path,
not a bayes-specific defect behind a safety gate.

### CORRECTION 2: the stated mechanism's scope clause is wrong
They wrote "SCOPE: every measurement is at cluster() == i; the magnitude at cluster()
coarser/finer than i is untested."  I tested both.
- `r02_scope.do` design C: i = state, **cluster(region) COARSER than i**: ratio
  bayes/standard on att_it = 0.6199 (homogeneous), 0.2975 (heterogeneous), and bayes
  se(att_it) is again **bitwise invariant** across arms (relmove 0.0000e+00).
  -> the effect is NOT specific to cluster() == i.
- `r02_scope.do` design B: i = onset cohort (5 groups of 8 states), **cluster(state)
  FINER than i**, 24 treated (i,t) cells of 8 rows: ratio bayes/standard on att_it =
  0.9654 (homogeneous), 0.8523 (heterogeneous), and bayes se(att_it) MOVES +53.7% with
  the heterogeneity.  -> **no collapse at all** when the cluster is finer than i.
- `r05_designD.do` design D: i = state = cluster but THREE rows per (state,year) cell
  (individual-level administrative data under a state-level policy): ratio 0.6102
  (homogeneous).  -> singleton cells are not required; what is required is that the
  Dirichlet multiplier be CONSTANT WITHIN THE CELL, i.e. that cluster() nest i.
The correct predicate is "cluster() nests i" - the identical predicate to the documented
analytic-SE weakness - not "cluster() == i", and not "singleton cells".

### CORRECTION 3: the documentation framing
They argue the package documents only a conditional ("can be too small when treatment
effects vary across clusters") and that the conditional is too narrow.  Half right.
The PRIMARY statement in both the refusal (ado:6439) and sthlp:241 is UNCONDITIONAL and
exactly describes what I measured: the standard errors "omit the cell-sampling variance
component".  Under homogeneity the omitted component is the treated cells' own
idiosyncratic draw, which is a cell-sampling component; the unconditional clause covers
it.  What is genuinely misleading is the SECOND clause of the same sentence, which a
reader takes as a sufficient condition for safety and which my het=0 arm contradicts
(ratio 0.50 with exactly zero cross-cluster effect heterogeneity).  That is a
documentation-wording defect, not the documented-vs-actual divergence they describe.

### CORRECTION 4: "4.4x" and "a 95% interval with the content of a 26% interval"
Both treat bootstrap(standard) as truth.  Section 7 of the recon note records
bootstrap(standard) at cluster() == i OVER-covering att_i (149/150), and the prober's own
Monte Carlo puts bayes att_i at 0.62 of the TRUE sampling sd, not 0.23.  The
truth-relative deficit is about 1.6x; the 4.4x (6.0x on my panel) is the deficit against
a comparator that is itself conservative in this configuration.  Report the ratio to the
truth, not to standard.

### Net
The combination is gated behind an error that names the remedy, prints a caution when
forced, and the mechanism is documented unconditionally.  What is new and worth the
author's time: the exact invariance, the fact that it is shared by the DEFAULT analytic
path, the correct predicate (cluster nests i), and the misleading second clause.
severity: **medium**, confidence certain, realistic yes.

---
## bb-02 — CONFIRMED (low).  The value is right; the help sentence is wrong.
repro: `r03_bootsave.do` (own panel, 30 states x 16 years, 15% item nonresponse on y)

    ROWS_IN_DATA=480   ROWS_WITH_Y=408   e(N)=408   e(N_clust)=30
    bootsave _rep=0 value = +1.11111111110c4X-009 = 0.0020833333
    1/nrows = 0.0020833333   MATCH_NROWS=1
    1/e(N)  = 0.0024509804   MATCH_eN   =0
    bootalpha(equal) leg: _rep=0 = 0.0333333333 = 1/G  MATCH=1

sthlp:557 says that row "holds 1/N, N the number of observations in the estimation
sample".  It holds 1/(map-frame rows) = 1/480.  Confirmed.

CORRECTION to their write-up: they call this a defect in the VALUE.  It is not.  I
checked the scale the replicate columns are on: sum over clusters of n_c*w_c on
replication 1 is 0.9999999999999999 when n_c is the MAP-FRAME row count, and 0.8488 when
it is the estimated-row count.  The saved multipliers sum to 1 over map-frame rows, so
1/480 is exactly the flat multiplier on the same scale and is the right comparator; the
ado comment at 8935-8944 says so and reasons it correctly.  Only the help's phrase
"estimation sample" is wrong.  Fix the sentence, not the divisor.
Their secondary point (the Gamma shape is also the map-frame count, so dispersion is off
by ~sqrt(n_est/n_rows)) is arithmetically right, is a footnote on an option the help
already calls unsound for panel data, and leaves the replication mean unbiased.
severity: **low**, confidence certain, realistic yes.

---
(bb-03 and my own findings below, filled in as the runs land)

---
## bb-03 — PARTIAL (severity medium -> low).  Direction right, magnitude and mechanism wrong,
## and the standard-error half of it is printed in the help.
repro: `r04_tailshape.do` (replication-distribution diagnostics), `r06_[abc].do` + `agg.py`
(my own paired coverage Monte Carlo), `r10_sesd_het.do`.

### Part (a) - "every reported SE is below the true sampling sd"
Reproduces for `att` on my own 55-panel MC of their stated design (30 clusters, 10 periods,
true att = 1, homogeneous, cluster(st), i = st, reps(79)):
mean se / true sd = standard 0.907, bayes 0.848, analytic 0.850.  All below 1, bayes lowest.
But this is (i) sec-1 for the analytic leg, and (ii) ALREADY IN THE HELP for the bootstrap
legs: sthlp:532 reports a 24-cluster staggered simulation with "standard error 0.967 of the
truth" under bootstrap(standard) and "0.877" under bootstrap(bayes).  Documented ratio
bayes/standard = 0.877/0.967 = 0.907; mine = 0.848/0.907 = 0.935.  Same result.
I did NOT re-measure their att_i figure of 0.621 - I dropped att_i from my MC for runtime.
That figure is bb-01's collapse expressed as a coverage number, not an independent fact.

### Part (b) - the percentile-interval coverage gap.  Direction reproduces, magnitude does not.
60 paired panels, same panel and same rseed for both engines (final):

    standard  percentile 56/60 = .933   mean width 0.64209   mean se/true sd .907
    bayes     percentile 54/60 = .900   mean width 0.59762   mean se/true sd .848
    analytic                                                 mean se/true sd .850
    discordant: standard covers & bayes misses 3 ; bayes covers & standard misses 1

Their 23/30 vs 28/30 (a 16.7pp gap, 5/5 one-directional discordance, sign-test p .031) does
not reproduce: mine is a 3.3pp gap with a 3-1 discordance, sign-test p = .31.  On 60 panels
a 16.7pp gap would have shown up: the standard error of a coverage difference on paired data
with 4 discordances is about 2.7pp.

### CORRECTION: "the bayes interval is shifted, not narrower" is wrong on my panels.
They could not diagnose why the two intervals sit at different locations.  On mine they do
not sit at different locations - the bayes interval is simply NARROWER:
- MC mean percentile width ratio bayes/standard = 0.931 (they reported 0.981).
- `r04_tailshape.do` at reps(999), where the endpoints are not extreme order statistics:
  P1 standard sd 0.2082 / width 0.7862 / skew -0.042 / (repmean-b)/sd -0.0023
  P1 bayes    sd 0.1958 / width 0.7547 / skew +0.015 / (repmean-b)/sd +0.0029
  P2 standard sd 0.1547 / width 0.6328 / skew +0.289
  P2 bayes    sd 0.1396 / width 0.5518 / skew +0.433
  Both replication distributions are centred on the estimate and near-symmetric; the width
  ratio is 0.96 (P1) and 0.87 (P2).
- The same panel 1 at reps(79) gives bayes skew +0.486, (repmean-b)/sd -0.188 and a width
  ratio of 0.871.  Those are reps(79) sampling noise in the draws, not engine properties:
  at reps(79) the 2.5% endpoint is essentially the second order statistic.  Part of their
  16.7pp is a replication-count choice.

### CORRECTION: the "opposite direction from sec 7" argument does not hold.
Section 7's over-coverage result is bootstrap(standard) on **att_i**; bb-03(b) is
**att**.  Different estimand, no contradiction to draw.

---
# MY OWN FINDINGS (not seeded)

## N1 (HIGH) — e(gteffects) cell standard errors, p-values and intervals are understated
## by ~4x at cluster() at-or-above i, on an ordinary state-year panel, with no caveat in
## that option's documentation.
repro: `r11_gtcell_se.do` (self-contained, 60 draws, everything computed in Stata),
corroborated by the cell columns of `r10_sesd_het.do` (a different 80 draws, different seeds).

**The headline number: the cell's own printed 95% confidence interval covers its own true
effect 45.0% of the time.**  `r11_gtcell_se.do`, 60 draws, homogeneous true effect 1:

    att   TRUE sd 0.15138   mean reported se 0.14064   se/sd 0.929
    cell  TRUE sd 0.74714   mean reported se 0.21431   se/sd 0.287
    cell reported se ranges 0.14684 to 0.28492 across draws (it is not a constant, just small)
    cell 95% interval coverage of its own true effect (= 1):  0.450   (nominal 0.95)

The att row on the same runs is fine (0.929) - it is the cells that fail.

**Scope control, and it is clean.** `r12_gtcell_finer.do` is the identical data and the
identical cell, with `i` = the onset cohort instead of the state - the help's own intended
practice, "cluster BELOW the level of i" (sthlp:255) - so each treated cell holds 7-8
states and cluster(st) is finer than i:

    att   TRUE sd 0.15138   mean reported se 0.15400   se/sd 1.017
    cell  TRUE sd 0.31685   mean reported se 0.32053   se/sd 1.012
    cell 95% interval coverage of its own true effect:  0.950   (nominal 0.95)

Exactly nominal.  So this is not "per-cell inference is intrinsically hard"; it is the
containment geometry, and the failure is confined to cluster() at-or-above i - the same
predicate as bb-01 and as sec 1.
30 states x 10 years, one row per state-year, i = st, cluster(st), true effect 1 for every
treated cell, `gteffects`.  For the cell (i = 1, relyr = 0), across 80 draws:
    true sampling sd of the cell estimate   0.8671
    mean reported analytic se of that cell  0.2166   -> se/sd = 0.250
and with a unit-level effect shock of sd 1.5 the reported se does not move at all (0.2166)
while the true sd rises to 1.5269 -> se/sd = 0.142.
The mechanism is sthlp:255's ("at or above the level of i the cluster-robust covariance
stops seeing the treated cells' own variation") - a treated (i,relyr) cell inside one
cluster is fitted by its own dummy, its residual is exactly zero, and nothing in the meat
sees it.  What is new: sthlp:255 quantifies this for the **att rows** at "roughly 0.66 to
0.75" and says nothing about e(gteffects); the gteffects help block (sthlp:1136-1151)
carries two reading notes, neither about this; and the cells' own se, p, lb and ub columns
are printed and returned at 0.25 of the right scale.  Four times worse than the number the
help gives for the rows above it.

## N2 (medium) — the documented 0.66-0.75 understatement range does not bound the att
## standard error once treatment effects vary across units: the reported SE is bitwise
## constant in that variation.
repro: `r09_analytic_blind.do` (bitwise), `r10_sesd_het.do` (the SE/SD numbers).
This SHARPENS sec 1 / sthlp:255; it does not contradict the mechanism, which that sentence
states.  What it contradicts is the only quantification on record.
`r09`: one 28-state panel, four arms identical except for the treatment-effect
heterogeneity.  Reported analytic se(att):
    arm 0 homogeneous                      +1.2e2d6dd8f2531X-003
    arm 1 unit-level effect shock sd 1.5   +1.2e2d6dd8f2531X-003   (bitwise equal)
    arm 2 unit-level effect shock sd 4.0   +1.2e2d6dd8f2531X-003   (bitwise equal)
    arm 3 per-cell effect shock sd 1.5     +1.2e2d6dd8f252dX-003   (3 ulps)
while bootstrap(standard) se(att) goes 0.152361 (arm 0) -> 0.935422 (arm 2), a factor 6.14,
and the point estimate moves 0.905 -> 1.463.
`r10`: 80 draws per arm, same design, se/sd against the TRUE sampling sd:
    unit-level effect sd 0.0   true sd 0.1515   reported se 0.1415   se/sd 0.933
    unit-level effect sd 0.5   true sd 0.1921   reported se 0.1415   se/sd 0.736
    unit-level effect sd 1.5   true sd 0.3956   reported se 0.1415   se/sd 0.358
The mean reported se is the SAME NUMBER in all three arms.  The project's own measurement
(review_process/se_coverage_2026-09-14/se_cov_mc.do) uses two DGPs whose tau is constant
within cohort - no unit-level effect heterogeneity at all - which is why it landed on
0.66-0.75; that is the middle arm here.  A reader who takes "roughly 0.66 to 0.75" as the
size of the problem is reading the homogeneous-effects case.

## N3 (low) — the bootstrap(bayes) x cell-weighted gate is regime-blind: it refuses where
## bayes is fine and lets you through where it is worst.
repro: `r02_scope.do`.
ado:6437 refuses `bootstrap(bayes)` with any cell-weighted estimand whatever the cluster
geometry.  Measured ratio bayes/standard on att_it:
    cluster FINER than i (i = cohort, cluster = state)   0.9654 (hom) / 0.8523 (het)  <- refused anyway
    cluster == i                                          0.5003 / 0.1671             <- safety(off) proceeds
    cluster COARSER than i (i = state, cluster = region)  0.6199 / 0.2975             <- safety(off) proceeds
The predicate that governs the collapse is exactly `wdb_nest`, which wooldid_bootrun already
computes at ado:8987-8998 on the resample's own sample; the gate at ado:6437 runs before it
and does not use it.  A regime-aware gate would refuse only under nesting, and its caution
could drop the misleading "when treatment effects vary across clusters" clause there.

---
# NEGATIVE RESULTS (probed hard, found nothing)

- **The two freshest centring fixes hold.** `r07_centring.do`: 26-state panel with unbalanced
  entry and aw=pop spanning p10 0.37 to p90 8.54, reps(300).  (replicate mean - b) / bootstrap
  sd: standard att -0.179, att_it -0.013, att_i -0.073; bayes(equal) att -0.054, att_it +0.027,
  att_i +0.007.  The ado comment says the pre-2026-09-13 divisor bug displaced the draws by
  9.36 bootstrap sd; nothing near that survives on weighted unbalanced data.
- **No ungated bayes collapse on the displayed rows.** `r08_gteffects.do`: with plain `att`
  (no safety gate) on the heterogeneous panel, bayes/standard se ratios are att 0.955,
  pretreatment 1.018, ES_pre2 0.907, ES_pre1 0.978, ES_post0 0.966, ES_post1 0.956,
  ES_post2 0.975.  `r08b`: over all 163 e(gteffects) cells the mean bayes/standard ratio is
  0.9414.  The cell-weighted collapse does not leak into the att-weighted surfaces.
- **The bayes cell-weighted SE is not a frozen constant.** `r05_designD.do` arm E: with the
  idiosyncratic noise sd moved 0.6 -> 1.5 and no cluster-level heterogeneity at all, bayes
  se(att_it) moves 0.05519998 -> 0.13799996, exactly the 1.5 relative change.  It is a clean
  function of the within-cell channel and exactly blind to the between-cluster one.
- **bootalpha(equal) observed row is right** (1/G), and the replicate multipliers under
  bootalpha(size) sum to 1 over the map-frame rows to 1e-16.
- **The bayes cell-weighted gate keys on the user's own safety() value**, not the working copy
  that `contreatwithin` rewrites: `safetyuser` is captured at ado:4828, one line BEFORE the
  contreatwithin relaxation at ado:4829, and bootctl is handed `si(safetyinput)`.  A
  contreatwithin run cannot slip past the refusal.  (Read, not run.)

---
## Closing note on how the three verdicts fit together
One predicate explains all of it: **does cluster() nest i?**
- Nesting (cluster == i, or coarser): every treated (i,relyr) cell sits wholly inside one
  cluster.  Its residual is exactly zero under saturation, no reweighting touches it, and
  every variance that is not built by RESAMPLING THE CELLS is blind to the cell's own
  variation.  That covers the analytic se (N2), the e(gteffects) cell se (N1), and
  bootstrap(bayes) for the cell-weighted estimands (bb-01).  bootstrap(standard) escapes
  only because its cell-weighted leg physically re-identifies the drawn copies.
- Finer than i: nothing fails.  bayes tracks standard (0.965), the analytic att se/sd is
  1.017 and the cell CI covers 0.950.
That is the help's own "cluster below the level of i" advice (sthlp:255), now with numbers
on three surfaces it does not mention.  bb-03(a)'s "bayes att_i is the worst at 0.621" is
also a nesting-regime number, not an engine property.
