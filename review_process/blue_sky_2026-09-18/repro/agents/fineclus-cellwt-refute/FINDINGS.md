# SKEPTIC PASS over FCW-01..FCW-06 (fineclus-cellwt)

Every number below comes from a do-file in this directory that I wrote from the CLAIM,
on data I generated. The prober's do-files were read for orientation and never run.

## Verdict table

| id | verdict | severity (theirs -> mine) | one line |
|---|---|---|---|
| FCW-01 | PARTIAL   | high -> high   | mechanism real, direction real; 23% is a degenerate-DGP ceiling, realistic cost is 4-16% |
| FCW-02 | CONFIRMED | high -> high   | reproduced on a realistic DGP: se(att_i) falls 13-16% while se(att), b and e(N) are bit-identical |
| FCW-03 | PARTIAL   | high -> low    | composition term measured and confirmed, but the record they "correct" had ALREADY retracted that explanation, and the attribution to the re-id leg is wrong |
| FCW-04 | PARTIAL   | medium -> low  | gate is real and wider than claimed, but the footer still states the resolution on every leg except bootdisplay(wald) |
| FCW-05 | PARTIAL   | low -> low     | no e() marker (true), but sthlp:241 documents the two-regime RULE, so "a user cannot tell" is too strong |
| FCW-06 | PARTIAL   | low -> low     | their shape really is refused, but the motivating case IS reachable when i is the onset grouping - which is what the option asserts |

---

## FCW-01 — PARTIAL.  Mechanism confirmed; "23%" is a degenerate-DGP ceiling.

**Own repro:** `r01a_straddle_real.do` (realistic arm) and `r05_emarker.do` (degenerate arm).

Realistic county x firm panel: 15 counties x 4 firms = 60 units, T = 7, 420 rows, treatment
assigned at county level, tau = 1 + a_g + b_i with a_g ~ N(0,.4^2) ACROSS counties and
b_i ~ N(0,.4^2) WITHIN them, unit RE N(0,1), idiosyncratic N(0,1).  Arm B differs from arm A
in exactly one cell of the cluster variable: firm 1 spends its last period in county 2.
reps(199), three dataset seeds x two rseeds (paired: same rseed draws the same cluster multiset
in both arms, so the comparison is regime-only).

    dseed rseed  se(att_i)B/A  se(att)B/A  b(att_i)B - A
     4101  31       0.9547       0.9971          0
     4101  32       0.9486       0.9991          0
     4102  31       0.8861       1.0134          0
     4102  32       0.8575       1.0154          0
     4103  31       0.8726       1.0142          0
     4103  32       0.8709       1.0155          0
    mean            0.8984       1.0091

So the direction, the silence and the estimand-selectivity are all real: se(att_i) falls, se(att)
does not move, and b(att_i) is bit-identical.  But the cost on a realistic DGP is ~10% (range
4.5-14%), not 23%.  The 23% is recovered only when the DGP is nearly pure composition variance:
`r05_emarker.do` (tau = 1 + .3*county exactly, noise sd .05, 15 counties x 4 units) gives
se(att_i) ratio **0.7896**, against the presence-only asymptote sqrt((1/e)/(1-1/e)) = 0.7629.
FCW-01's own noisy probe1 arm already showed 6% and one arm ABOVE 1, which is the honest number
for realistic input.

**The part of the claim that should survive verbatim, and is the reason this is still a defect:**
the predicate is GLOBAL.  In arm B, 59 of 60 units still lie wholly inside one county, and for
those 59 the project's own comment at ado:9334-9337 says keeping `i` "would make a twice-drawn
cluster's cells identical duplicates that the equal-cell weighting counts ONCE, collapsing the
composition variance of att_i toward zero."  One straddler therefore imposes the straddler's rule
on the 59 units that nest.  That is a wrong number for them, not a choice of convention.

**Severity note:** FCW-01 and FCW-02 are ONE root cause and carry ONE severity, high - not because
of the 23%, which is wrong, but because the trigger is a missing outcome or a reclassified unit,
the direction is anti-conservative, and the project's own comment says the resulting number is
wrong for the units that DO nest.  FCW-01's verdict is PARTIAL only because its headline magnitude
is a degenerate-DGP ceiling.

**Correction to the framing:** the help DOES document the two-regime rule (`sthlp:241`: "Where
cluster() nests i the copies become distinct units; where it is finer than i, or cuts across it,
the copies stay inside their cohort").  What is undocumented is that the switch is all-or-nothing
over the whole sample and what it costs.  Severity revised high -> medium.

## FCW-02 — CONFIRMED on realistic input, and it is the same root cause as FCW-01.

**Own repro:** `r02_phantom.do` (log of the completed pairs: `r02_phantom_PARTIAL_3pairs.log`).

Same realistic panel.  In BOTH arms firm 1's t = 3 outcome is missing.  Arm D additionally records
that one unestimable row under county 2.  reps(199).

    dseed rseed   se(att_i) A      se(att_i) D     ratio   se(att) A = se(att) D   e(N)
     7201  41     .24803561        .21190601       .854    .23407218  (identical)   419 / 419
     7201  42     .26059414        .22797831       .875    .24655838  (identical)   419 / 419
     7202  41     .24650980        .20780768       .843    .30014947  (identical)   419 / 419

b(att) = 1.0390337 and b(att_i) = 1.0403452 to 8 figures in both arms of pair 1 (1.1734777 /
1.1403918 in pair 3), e(N) = 419 in both, and se(att) is identical to every printed digit -
which also proves the two arms drew the same clusters and fitted the same rows.  A row that
contributes nothing to the estimation sample moves se(att_i) by 13-16%.
Source: `smp` is built from `ifr` = the RAW user `if` (ado:394 passes `ifr(`wdp_if')`,
`wdp_if` set at ado:615) plus a nonmissing cluster, with no markout.
Severity high stands.  **But FCW-01 and FCW-02 are ONE root cause** - a global predicate
evaluated on the design sample - and should be counted once.

## FCW-03 — PARTIAL.  The composition term is real and I measured it.  The "correction to the
## record" was already made by the project, and the attribution to the re-id leg is wrong.

**Own repro:** `r03_compvar.do` (mechanism), `r03b_whichleg.do` (attribution).

`r03_compvar.do`: 40 states, T = 6, one unit per state, 20 treated at t = 4, noise sd .02, i =
state, cluster(state) == i, bootstrap(standard) reps(199).  HET gives tau .4 / 1.2 (Var_s = .16),
HOM gives tau .8 for every treated state (Var_s = 0); identical noise draws, identical truth.

    het dseed rseed   se(att_i) = se(att)
     1   9001   51        .08863831
     1   9001   52        .09461377
     1   9002   51        .08991827
     1   9002   52        .09687149
     0   9001   51        .00441337
     0   9001   52        .00410414
     0   9002   51        .00565302

Predicted composition SE sqrt(.16/20) = .0894.  Measured .0886-.0969 under HET and .004-.006 under
HOM - a factor of ~20 - on data whose idiosyncratic noise is 2% of the effect.  So the mechanism
FCW-03 names is correct and quantitatively right.

**Three corrections.**
1. The record FCW-03 says it is correcting HAD ALREADY BEEN CORRECTED, in the same document that
   recorded the number.  `test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md` (F3/F4 amendment 3)
   says: "A wrong analogy, retracted here so it does not propagate. The implementation's commit
   message likened this to sthlp:255 ... That passage documents the ANALYTIC cluster-robust
   covariance ... its direction is anti-conservative ... Do not cite that passage for this
   finding."  That is FCW-03's entire "the record's explanation cannot produce over-coverage"
   argument, made by the author on 2026-09-16.
2. The record also already states the mechanism as a hypothesis: "conservatism of the percentile
   cluster bootstrap where every cohort is its own resampling cluster."  FCW-03's contribution is
   the measurement, not the explanation.
3. The attribution "the nesting re-id LEG injecting a cohort-composition component" is too narrow.
   The same record reports that `att` - which is invariant to the re-id unit rule, as FCW-01's own
   control shows (se(att) ratio 1.0029 there, 1.0091 in my r01a) - covered 149/150 on the SAME arm.
   The composition term is a property of resampling cohorts as clusters at `cluster() == i`, not of
   the unit rule.  `r03b_whichleg.do` tests this directly by asking for `att` alone, which takes the
   COUNT path and never enters the re-id leg (ado:378).
Severity high -> low: it measures a number the project recorded, explained as a hypothesis, and
ruled "recorded, not fixed" on 2026-09-16.

## FCW-05 — PARTIAL (claim true, "cannot tell" too strong).

**Own repro:** `r05_emarker.do`.  Two runs on one dataset, same rseed, differing in one cell of the
cluster variable.  I dumped every `e(scalars)` in %21x hex and every `e(macros)` for both runs:
`diff r05_A.txt r05_B.txt` is EMPTY, while se(att_i) moves from +1.0db021a8e5c64X-002 to
+1.a9dfa55549169X-003, a ratio of 0.7896.  So nothing in e() records the regime - confirmed.
**But** sthlp:241 documents the RULE, and a user can evaluate it on their own data; and the A10
caution at ado:407 is a de facto marker whenever a covariate names `i`.  "No runtime record" is
right; "a user cannot tell which variance concept produced their standard error" is too strong.

## FCW-04 — PARTIAL.  The gate is real and WIDER than claimed; the information loss is narrower.

**Own repro:** `r04_caution.do` (six runs, reps(10), one self-generated 12-county panel) and
`r04b_waldleg.do`.

    RUN1  att att_i att_it, bootstrap(standard) reps(10), no cluster()  -> rc 0, NO caution
    RUN2  same + cluster(cty)                                           -> rc 0, caution PRINTS
    RUN3  att only, no cluster()                                        -> rc 0, NO caution
    RUN4  att only, no cluster(), bootdisplay(wald)                     -> rc 0, NO caution
    RUN5  att only, no cluster(), boottstat                             -> rc 0, NO caution
    RUN6  bootstrap(bayes) reps(10), no cluster()                       -> rc 0, NO caution

So the cluster-gating is confirmed and covers all three caution variants plus bayes, not only the
percentile one FCW-04 tested.

**But the substance is not lost on most of those legs.**  The disclosure footer at ado:11784-11785
is NOT cluster-gated.  RUN1 and RUN3 still print, unprompted, "The rank-based coefficient p-values
have resolution 2/10 = 0.200: a p printed as 0 means p below that (no completed replication on the
far side of the null)", and a below-table note "p = 0 above means p < 0.200 (2/10 completed
replications, the resolution of the rank-based p)".  RUN5 prints the studentized analogue
("their resolution and their smallest attainable value are both 2/(n + 1) = 0.182").
That is the content FCW-04's "expected" paragraph says the user is denied.

The one leg where the footer is ALSO silent is `bootdisplay(wald)`, because ado:11785 opens
`else if "`wdbdsp'" != "wald"`.  RUN4's footer carries no resolution sentence at all - not even the
"joint, equality and adversarial p-values are add-one shares ... smallest value 1/(n + 1)" clause
that every other leg gets.  `r04b_waldleg.do` tests whether an unclustered wald run that actually
PRINTS joint p-values discloses their floor anywhere.
Severity medium -> low, with the wald leg carried separately below as a new finding if r04b
confirms it.

## FCW-06 — PARTIAL.  Their configuration really is refused; the motivating case is reachable.

**Own repro:** `r06_permgate.do`.  One panel, 12 states x 4 counties, T = 8, states 1-4 onset 4,
states 5-8 onset 6, 9-12 never treated.  Four arms, plus the bootstrap control:

    X1  i = state (finer than the onset grouping), cluster(county), permute(shuffle)      r459 nesting refusal
    X2  X1 + permonsetcohorts                                                             r459 bijection refusal
    Y1  i = cohort (IS the onset grouping), cluster(county) strictly finer, no option      r459 nesting refusal
    Y2  i = cohort + permonsetcohorts                    ***  rc 0, e(perm_irebuilt) = 1 ***
    Z   i = state, cluster(county), bootstrap(standard)                                    rc 0

So X1/X2 reproduce exactly as FCW-06 reports.  But Y2 is the T1 remedy working on exactly the
"coarse i, finer cluster" pairing the ticket's motivating sentence names, and FCW-06 never built
that arm.  Corrections to the claim:
- "as soon as two cohorts share an onset year" is a misstatement: a cohort IS an onset group.  The
  refusal fires when `i` is FINER than the onset grouping, i.e. when `i` is not an onset grouping -
  which is the one thing `permonsetcohorts` exists to assert and the help states twice
  (sthlp:68, sthlp:431) and again at sthlp:255.
- What survives: to reach permute() the user must define `i` as the onset grouping, and that
  changes the set `att_i` averages over.  That is a real scope limit and worth a sentence somewhere,
  but it is documented behaviour, not an uncovered case.  Severity low stands; "code defect" does
  not.

---

# NEW FINDINGS (found by me while refuting; not seeded)

## SKEP-N1 — the bootstrap's degrees of freedom come from the DESIGN sample, so a cluster whose
## outcome is never observed still buys a degree of freedom.  The same table prints both numbers.

**Own repro:** `n03_designclusters.do`.  Severity **high**, realistic **yes**.

`wooldid_bootrun` builds its cluster map on `smp` = (user `if` == 1) & nonmissing cluster
(ado:8910-8925), counts `wdb_C` there, and hands it back as `wdb_r_nclus` (ado:9537).  That count
becomes the default `bootdf(clusters)`, i.e. `e(df_r) = G - 1` (ado:423), and it keys the
"only N clusters are available for the bootstrap" and "only N treated cluster(s)" cautions
(ado:6502, ado:6520, both computed on `smp2`, the same design concept).  None of those counts is
the estimation sample.

14-county panel, T = 7, 4 units per county, outcome collected in only some counties:

    arm                                        e(N)   header          e(df_r)   analytic e(df_r)
    C0  all 14 counties measured               294    N clusters = 14    13            13
    C7  7 counties measured, 7 all-missing y   147    N clusters = 7     13             6
    C12 2 counties measured                     42    N clusters = 2     13             1 (not run)

In C7 the run prints "N clusters = 7" in its own header and then "e(V) and e(df_r) carry the
standard-deviation family referred to the **t(13)** distribution"; the analytic run on the identical
data reports `e(df_r) = 6`.  t(13) = 2.160 against t(6) = 2.447, so every `bootdisplay(wald)`
interval, and everything that reads e(V)+e(df_r) without asking (the footer names esttab, coefplot
and honestdid), is about 12% too narrow.  In C12 only TWO clusters carry any data, 55 of 99
replications failed, and the run still reports 13 degrees of freedom and never prints the
few-cluster caution, because 14 >= 10.
This is the same root cause as FCW-02 (design sample vs estimation sample) on a different surface,
and it does not need a straddling cluster or a cell-weighted estimand to fire.

## SKEP-N2 — the unclustered leg T2 opened is mis-calibrated on panel data, by +38%.

**Own repro:** `n02_noclus_atti.do`.  Severity **medium**, realistic **yes** (direction is
conservative, which is why it is not higher).

`bootstrap(standard)` with a cell-weighted estimand and NO `cluster()` was refused r198 at old
ado:6353 and now runs, resampling observations.  Design: 12 hospital systems x 5 hospitals, T = 6,
onset at the system level, `i` = system; dependence is a hospital random effect (sd .8) plus iid
noise (sd 1); `tau_s` FIXED across the Monte Carlo, so the target is "the cohorts in the data".
Oracle = MC SD of att_i over 60 redraws of the noise.

    MC SD of att_i (60 draws)                  0.21252
    bootstrap cluster(hosp)  (keep-i branch)   .21398  .20873  .27769   mean .23346  ratio 1.099
    bootstrap NO cluster()   (observation)     .29362  .29409  .28238   mean .29003  ratio 1.365

Three datasets, reps(99).  The keep-`i` branch that T2 actually built scatters around the oracle
(ratios 1.01, 0.98, 1.31; the MC SD itself carries ~9% standard error at 60 draws), so it is
calibrated to within what this probe can resolve - a genuine negative result, since nothing
outside the project suite had measured that branch.  The unclustered leg the same ticket opened
sits at 1.33-1.38 on every dataset: about 36% too wide, consistently.
The mechanism is the one the DELETED refusal named at old ado:6353 - "observation-level resampling
cannot re-identify duplicated treatment cells": observation draws fragment units and cells rather
than redrawing them.  T2's own ticket text says of that refusal "revisit whether that refusal is
still necessary or can be relaxed to observation-level resampling"; it was relaxed, and this is
what the revisit measures.  `sthlp:241` states only the mechanism ("With no cluster() at all the
resampling is observation-level") and nothing about its calibration.

## SKEP-N3 — the retracted analogy survives in the file that records the measurement.

**Repro:** read-only.  Severity **low**, realistic **no** (internal record, not user-facing).

`test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md` (F3/F4 amendment 3) retracts the sthlp:255
analogy for the cluster()==i over-coverage, explicitly "so it does not propagate": "That passage
documents the ANALYTIC cluster-robust covariance at or above `i` and its direction is
anti-conservative ... Do not cite that passage for this finding."  The header of
`test/suite/mc/mc_fineclus.do`, lines 36-38, still says the opposite: "it is the same effect the
help already describes at 'What level to cluster at' - at or above the i level the covariance stops
seeing the treated cells' own variation. Here it lands on the conservative side."  That stale copy
is what FCW-03 found and spent its claim refuting.  Two records of one measurement disagreeing is
exactly the propagation the amendment tried to stop.

## SKEP-N4 (observation, not filed as a defect) — cross-cutting cluster() now runs silently.

**Own repro:** `n04_crosscut.do`.  40 units, T = 9, staggered onsets.  `cluster(t)` cross-cuts `i`
completely and was refused r459 before 2026-09-16.  It now returns rc 0 with 199/199 replications
complete, se(att) .2635 / se(att_i) .2629 against analytic cluster(unit) .1562 / .1527 and
bootstrap cluster(unit) .1833 / .1763, on `e(df_r) = 8`.  Nothing is wrong with the arithmetic -
resampling 9 time clusters is what was asked for - but nothing anywhere says that resampling the
time dimension re-draws the event-time structure the estimand is defined on.  Recorded, not filed.

---

## FCW-03 attribution, measured (`r03b_whichleg.do`)

`wdb_reid` is set only when a cell-weighted estimand is asked for, or `boottstat` with `cluster()`
(ado:378).  A run asking for `att` ALONE at `cluster() == i` therefore never enters the re-id leg;
it takes the count/reweighting path.  On the dgpE-shaped HET design (40 states, 2/4/6 units,
tau .4/1.2, noise sd .02, cluster(state) == i, reps(199), same rseed):

    re-id leg   (att att_i) :  se(att) = .0709192322   se(att_i) = .0886986123
    count path  (att alone) :  se(att) = .0709192322

Bit-identical.  The cohort-composition component that produces the recorded over-coverage is a
property of resampling cohorts as clusters at `cluster() == i`, present on both legs, and
independent of the re-id unit rule - which is also why the project's own `att` control over-covered
149/150 on the same arm.  FCW-03's mechanism is right; "the nesting re-id leg injecting" it is not.

---

## Verification of FCW-04's "regression" claim (I did not take it on trust)

`git show -U25 813c414 -- wooldidinstall/wooldid.ado`: the `if `reps' > 0` block is deleted from a
site whose 22 preceding lines are the CS-20 rationale comment and the `bootstrap() requires reps(#)`
error - top level, no cluster gate - and re-inserted after the two-cluster check inside
`if "`clus'" != "" {`.  So the move is real.  What also matters is the CS-20 comment's own reason
for the caution existing: "The footer that would have said so is conditional, not routine:
wooldid_pfloor prints only when a displayed p is exactly 0, and exits under bootdisplay(wald)
without boottstat - so a low-rep run whose p-values all happen to be non-zero was told nothing at
all."  The author ruled the caution in for that reason on 2026-09-07 and RF-2 confined it to
clustered runs on 2026-09-13 without saying so.  That is why this is a regression and not a
cosmetic move - but the ROUTINE disclosure line at ado:11785 (a different line from
wooldid_pfloor) does still print the resolution on every non-wald leg, which is why I put the
severity at low rather than medium.

---

## NEGATIVE RESULT (reported as such) — the branch T2 actually built is calibrated.

**Own repro:** `n01_keepi_oracle.do` (and `n02_noclus_atti.do`, same design).
12 hospital systems x 5 hospitals, T = 6, onset at system level, `i` = system, `tau_s` fixed,
hospital random effect sd .8 + iid sd 1.  Oracle = MC SD of att_i over 60 redraws.

    MC SD of att_i                                   0.21252   (ratio 1.000 by construction)
    analytic cluster(hosp), i.e. cluster BELOW i     0.22913   ratio 1.078
    bootstrap cluster(hosp)  = keep-i branch         0.23346   ratio 1.099
    bootstrap cluster(sys) == i = nesting branch     0.22641   ratio 1.066   (range .171-.331)
    bootstrap NO cluster()   = observation leg       0.29003   ratio 1.365   (range .282-.294)

The keep-`i` branch agrees with the legitimate analytic reference and with the oracle.  Three
datasets is too few to judge the `cluster(sys) == i` arm - its spread is 2x the others, which is
what 12 resampling clusters buys - so the composition claim rests on `r03_compvar.do`, where the
noise is turned down to 2% and the term is isolated cleanly.  I attacked the keep-`i` branch from
three directions (oracle calibration here, estimand invariance in r01a/r02 where b is bit-identical,
and the cross-cutting case in n04) and found nothing wrong with the branch itself.  What is wrong
is the PREDICATE that decides when to use it.
