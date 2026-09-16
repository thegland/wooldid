# Preregistered signature — clustering finer than i (tickets T1/T2/T3), 2026-09-16

Written BEFORE any line of `wooldid.ado` is changed for these tickets. The tickets are
`review_process/feature_tickets_2026-09-15.md` (T1 permute onset rebuild; T2 cell-weighted
estimand weights; T3 boottstat cluster re-id). T4 is closed by ruling and is not in scope.

These are FEATURE tickets, so the usual old/new axis is unavailable twice over: the
configurations they enable are hard errors at HEAD (`error 459` / `error 198`), and they do not
exist in v1.3 at all. What replaces the oracle is stated below, per ticket, before the fact.

## F0 — what is preregistered, and what is honestly not

- **PREREGISTERED HERE** — the three readouts in F2/F3/F4, the two invariance pins in F1, and
  the decision rule in F5.
- **INHERITED, NOT NEW** — the coverage judgment convention (central ~95% binomial region,
  normal approximation rounded outward, Wilson interval printed beside; `analyze_mc.py`), the
  MINN=30 not-judged gate, and the MARGINAL corridor [.90,.98] for a .95 coverage claim. All
  three come from `PREREG_STAGEC_2026-09-13.md` C0 and are reused unchanged, not re-derived.
- **NOT PREREGISTERED** — the tolerance on the 2-1/G ratio in F4 (set at +/-15% of the
  predicted value while writing this file, with no simulation in hand to calibrate it); it is
  stated now so it cannot later be presented as a prediction.

## F1 — the invariance pins (these run FIRST, before any check is relaxed)

The standing doctrine is that the nesting check is the LAST thing to go, after the machinery
that makes it unnecessary is in and verified. Every intermediate state must be one in which no
wrong number can print. Both pins are therefore falsifiable before the feature is reachable.

**F1a — the re-id leg is pinned by the battery under nesting.** Three battery configs reach the
physical re-identification leg, all at `cluster(i)` (nesting trivially):

| line | config | route to the leg |
|---|---|---|
| `runbattery.do:327` | `att cluster(i) gteffects bootstrap(standard) boottstat` | boottstat + cluster() |
| `runbattery.do:424` | `att cluster(i) subgroups(sg) adversarial jointtests bootstrap(standard) boottstat` | boottstat + cluster() |
| `runbattery.do:447` | `[aw=w] att att_it cluster(i) bootstrap(standard) bootdf(normal)` | cell-weighted estimand |

**PREDICTION:** after the two-regime `__wdbnewi` branch is implemented but before `ado:6353`
and `ado:6362` are relaxed, `battery` is IDENTICAL vs the committed baseline. A nesting run must
take the `group(__wdbns, i)` branch, which is byte-for-byte today's code path. Any movement here
means the branch predicate is wrong, and it is diagnosed before the feature is exposed.

**F1b — permute is pinned by its 19 battery configs with the option absent.** `runbattery.do`
carries 10 `permute()`, 7 `permute(shuffle)` and 2 `permute(probabilities)` configs.
**PREDICTION:** with the new T1 option not specified, `battery` is IDENTICAL. The rebuild is
reachable only under the user's explicit assertion, so an absent option must change nothing.

## F2 — T1: the rebuilt i must preserve sharp-null exactness (`mc_ri_ext.do`)

**Why this readout.** T1 rebuilds `i` inside every replication from the drawn onsets. If that
rebuild is wrong, the replication's design is not the observed design and the permutation
distribution is no longer the null distribution of the statistic. Exactness under the sharp null
is the one property that detects this and nothing else does: a wrong rebuild that is merely
*different* still produces plausible p-values.

**The arm.** `mc_ri_ext.do` gains one arm: the existing sharp-null DGP, run with `cluster()`
strictly finer than `i` under the new option, against the same statistic. Every existing arm
is retained unchanged as the control.

**P1 (falsifiable).** Under the sharp null, on the new arm: the mean of `main_att`'s p is inside
the band `analyze_riext.py` already applies to the existing arms, judged by that analyzer's
existing rule, at its existing MINN gate. The band is NOT re-derived for this arm; if the new
arm cannot be judged by the existing rule, that is itself the finding and is reported as such
rather than met with a fresh band.

**P2 (falsifiable, integrity).** `rc == 0` on every row; `p` in [0,1] on every row;
`e(perm_n)` equal to the requested replication count on every accumulated column. P2 is the
column-loss tripwire: `ado:7556` compares rep-0 `colfullnames` against the observed run's and
errors on mismatch, and `ado:15683`'s `select(draws, draws :< .)` silently drops blanks from the
denominator. A rebuilt `i` whose codes do not reproduce the observed column set fails P2 loudly
instead of producing a quietly thinner permutation distribution.

**P3 (falsifiable, negative control).** With `cluster()` at or above `i`, the new option
specified, and the same seed: results must be IDENTICAL to the same run with the option absent.
Where the rebuild is a no-op it must BE a no-op.

## F3 — T2: cell-weighted estimand coverage under finer-than-i clustering

**The DGP.** Heterogeneous cohort effects (so that `att_i` and `att` genuinely differ - a DGP
where they coincide cannot detect a weighting error), clustering strictly finer than `i`.

**P4 (falsifiable).** Coverage of the `att_i` bootstrap interval at nominal .95, judged by the
inherited binomial-band convention of F0. Reported beside it, never judged: coverage of `att` on
the same runs, as the control that the DGP's heterogeneity is live.

**P5 (falsifiable).** On the same DGP with `cluster()` at or above `i`, coverage is judged by
the same rule. This is the regression leg: the feature must not move the nesting case.

## F4 — T3: the 2-1/G negative control (the only proof the re-id worked)

**Why it is a negative control and not a coverage run.** A twice-drawn cluster must count as TWO
clusters in the replication's cluster-robust SE. The old count path could not express this:
copies pooled into one score and the meat took `m^2` where it needed `m`, giving a measured
inflation of `E[m^2]/E[m] = 2 - 1/G` in the variance (about sqrt(2) in the SE), which broke
percentile-t pivotality. The re-id leg exists to remove exactly that. A run that merely produces
plausible SEs does not demonstrate the leg is working; a run that can reproduce the defect on
demand and then not reproduce it does.

**The switch.** `ado:9135` performs the substitution `cluster(`clus')` -> `cluster(__wdbns)` as a
literal string operation whose own comment records that it "can silently do nothing", that
"nothing checks it", and that a no-op returns the full inflation. A diagnostic global is added
so the same seed can run with the substitution suppressed.

**P6 (falsifiable, the control).** Same seed, same data, `boottstat` + `cluster()`: the ratio of
mean `se_r` with the substitution SUPPRESSED to mean `se_r` with it ACTIVE equals
`sqrt(2 - 1/G)` for that run's G, within +/-15% (a tolerance set here without calibration, per
F0). If the suppressed arm does NOT inflate, the control is uninformative and P7 cannot be
credited - the test would be passing for want of sensitivity, not for correctness.

**P7 (falsifiable, the assertion that makes the substitution checked for the first time).** On
every replication that posts `e(N_clust)`, `e(N_clust)` equals the count of distinct `__wdbns`
on the replication frame. This closes the gap `ado:9126` admits in prose.

**P8 (falsifiable).** P6 and P7 hold on a finer-than-`i` cluster variable, not only under
nesting - a twice-drawn sub-unit cluster must still be two clusters.

## F5 — what a failure means, decided in advance

- **F1a or F1b moves** -> the branch predicate is wrong. Not a bless candidate under any
  reading. Fix the predicate; do not re-baseline.
- **P2 fails** -> the rebuilt `i` does not reproduce the observed column set. T1 must reuse
  observed `i` codes where the onset matches, or the accumulator must be made aware. Shipping a
  version that drops columns silently is excluded in advance.
- **P1 fails while P2 passes** -> the rebuild is mechanically sound and inferentially wrong.
  That is a design defect in T1, reported to the author, not patched by widening the band.
- **P4 fails while P5 passes** -> the per-replication estimand weights are wrong specifically
  under finer clustering. T2's core claim is refuted.
- **P6's suppressed arm does not inflate** -> report the control as UNINFORMATIVE. Do not
  report P7 as evidence the leg works.
- Any of these may be overturned by the author; none may be overturned by re-running with a
  different seed until it passes.

## F2 AMENDMENT - P2 was crediting a tripwire that cannot fire (recorded 2026-09-16, mid-implementation)

**What P2 claimed.** "A rebuilt {it:i} whose codes do not reproduce the observed column set fails P2
loudly instead of producing a quietly thinner permutation distribution", resting on `ado:7556`'s
rep-0 `colfullnames` check against the observed run.

**What was measured.** The T1 implementation measured that **no accumulated `e(b)` column name
embeds an `i` code**: relabelling `i` leaves `colfullnames` identical. The consequence is that
`ado:7556` can detect a change in the column SET - a cohort vanishing from a replication - but
**cannot detect a wrong onset-to-cohort map**, which is the failure mode T1 actually introduces.
A rebuild that assigns the right number of cohorts to the wrong onsets passes the tripwire
silently and produces a plausible, wrong permutation distribution.

**Why this is recorded as an amendment rather than a correction in place.** P2 was written before
the measurement and credited the tripwire with more than it can do. Editing the original text
would hide that the protection was assumed rather than verified. The original P2 stands as
written and is narrowed here.

**P2a (the original, narrowed).** `rc == 0` on every row, `p` in [0,1], and `e(perm_n)` equal to
the requested replication count on every accumulated column. This detects column-SET loss - the
`select(draws, draws :< .)` denominator drop at `ado:15683` - and nothing else.

**P2b (new, and it is the one that covers the actual failure mode).** An explicit assertion in
`wooldid_permrun` that rep 0's rebuilt `i` equals the observed `i` on the design sample. Its
falsifiability is established by construction rather than by hope: it was measured to fire r459
on a deliberately reversed onset-to-cohort map. A negative control that fires on a known-wrong
input is what P2 lacked.

**Standing lesson for the rest of this work.** A tripwire named in a preregistration has to be
demonstrated to fire on the failure it is credited with, not merely to exist at a line number.
The same question is open for T3's `e(N_clust)` assertion (P7) and is answered there by P6's
suppressed arm; if P6's control comes back UNINFORMATIVE, P7 inherits exactly this defect.

## F2 AMENDMENT 2 - P1 cannot detect the failure class F5 implicitly assigned to it (2026-09-16, post-run)

**What F5 said.** "P1 fails while P2 passes -> the rebuild is mechanically sound and inferentially
wrong." That phrasing treats P1 as the inferential readout and P2 as the mechanical one, and by
omission implies that P1 passing means the rebuild is inferentially right.

**Why that is wrong, established by the T1 implementation.** Sharp-null exactness detects a rebuild
that is **not a function of the drawn assignment** - one that leaks observed onsets into the
replication, or attaches a row to the wrong cluster's draw. It **cannot** detect a rebuild that is
a correct function of the draw but assigns the wrong codes. A relabelling of the rebuilt cohorts is
still deterministic in the draw, so the permutation distribution is still the null distribution of
the statistic and the p-values are still uniform. P1 passes on a relabelled map.

**So the converse case has no readout in F2 as written**, and the gap is not hypothetical: it is
exactly the failure mode the whole ticket introduces, since the rebuild is precisely a map from
drawn onsets to cohort codes.

**What covers it** is P2b from the first amendment, the rep-0 assertion that the rebuilt `i` equals
the observed `i` on the design sample, measured to fire r459 on a deliberately reversed map. P2b was
added as a mechanical tripwire; it is in fact the **only** readout either prediction has for the
relabelling class. P1 is not sufficient on its own.

**A second, smaller qualification on P1's weight, volunteered by the implementation.** On the R1 DGP
the cohort-`i` and unit-`i` designs are near-degenerate for `att`: at full `%21x` precision the
`att` draws coincide on only 1 of 199 replications, yet the observed `att` agrees to ~1e-13 and the
null SD to 12 significant figures. E6's `att` verdict is therefore close to inherited from E2's
rather than independent of it. The four non-`att` columns carry the independent content.

**Standing lesson, and it is the same one as the first amendment.** A preregistered readout must be
stated with the failure class it detects, not only the property it measures. "Exactness under the
sharp null" sounds total and is not: it is blind to any error that is still a function of the draw.

## F3/F4 AMENDMENT 3 - P5's disposition, and P7 finally falsified (2026-09-16, post-run)

### P5 FAILED, on the high side, and the failure is not this ticket's

**Measured.** NEST arm, `cluster(state) == i`, `e(N_clust)` 40: `att_i` covered 149/150 = .9933
against the band (137, 148). F5 anticipated "P4 fails while P5 passes"; this is the reverse, and
F5 had no disposition written for it.

**It is not a regression.** The arm's first 40 draws were re-run against pristine `bd41e01` on the
same seeds. Comparing field by field, **18 of the 19 columns are identical** - `rc`, `nclust`,
`ndone`, `nfail`, both point estimates and all eight interval endpoints. The only column that
differs is `secs`, wall clock. So the nesting leg behaves exactly as it did before the two-regime
rule landed, and T2's core claim is not refuted by this.

**It is not an `att_i` weighting effect either, which is what the prediction was really about.**
The `att` control, reported but never judged, covered **149/150 on the same arm** - identical to
`att_i`. An error in the per-replication cohort weights would move the equal-cohort-weighted
estimand and leave the observation-weighted one alone. Both moving together, by the same count,
says the interval is wide for reasons that have nothing to do with the weighting. For contrast the
FINE arm at `e(N_clust)` 160 sits at .9467 for `att_i` and .9800 for `att`.

**What it actually is, stated as a hypothesis and not as a measurement**: conservatism of the
percentile cluster bootstrap where every cohort is its own resampling cluster. This run is the
first time the nesting case's coverage has been measured at all.

**A wrong analogy, retracted here so it does not propagate.** The implementation's commit message
likened this to `sthlp:255` ("What level to cluster at"). That passage documents the ANALYTIC
cluster-robust covariance at or above `i` and its direction is **anti-conservative** - standard
errors measured at 0.66 to 0.75 of reference, joint tests oversized. This measurement is the
BOOTSTRAP and its direction is conservative. Same location in the design, different estimator,
opposite sign. Do not cite that passage for this finding.

**Disposition: recorded, not fixed.** It is pre-existing, it is on the leg this ticket did not
change, and it is conservative rather than anti-conservative. The author declined a G sweep on
2026-09-15 and nothing here reopens that.

### P7 was credited without ever being falsified - now it has been

**The defect in the prediction.** The implementation's table reads "P7: never fired - PASS". An
assertion that never fires is evidence of nothing until it is shown to fire on the failure it
names. F5 credited P7 on the strength of P6's suppressed arm inflating, but that proves the
SUBSTITUTION does something, not that the ASSERTION would catch it failing to.

**And the diagnostic switch structurally cannot test it.** `ado:9358-9361` sets `wdb_sub` to 1
only when the substitution is performed, and `WOOLDID_BOOTNOSUB` suppresses both together, so in
the single state where `e(N_clust)` and the distinct `__wdbns` count genuinely disagree, the
assertion is skipped by construction.

**So the P2b method was applied instead.** On a scratch copy the `subinstr` search string was
pointed at a token that cannot match, reproducing precisely the fragility `ado:9345-9356` warns
about in prose - the substitution silently does nothing - while leaving `wdb_sub` at 1 so the
assertion still runs. Measured: **`error 459`**, reporting 80 distinct resampled cluster
identities against the count the replication posted.

**P7 is now demonstrated to fire on the failure it is credited with**, and the `ado:9345` comment's
"nothing checks it" is no longer true.

### Standing lesson, third instance

P2 credited a tripwire that could not fire. P1 was credited with a failure class it is blind to.
P7 was credited on the strength of a different prediction's control. Each time the fix was the
same: construct the failure, and require the check to fire on it. A preregistration should name,
for every check, the input on which it is expected to fail.

## F4 POSTSCRIPT - the P7 assertion's cost, measured (2026-09-16)

The implementation flagged P7's per-replication cost as unmeasured: the assertion runs two gtools
calls on every `boottstat` + `cluster()` replication. It was also reported that a NEST-arm re-run
took 130.2 s against the original 69.2 s, which read as a possible regression.

Measured, paired, on the battery's `runbattery.do:327` config (`boottstat`, `cluster(i)`,
`gteffects`, reps 20, same seed), pre-merge build against merged build, one Stata session each so
neither pays the other's warm-up:

| build | rc | secs |
|---|---|---|
| pre-merge (`281ca06`) | 0 | 45.062 |
| merged | 0 | 45.800 |

Ratio 1.016. The assertion costs about 37 ms per replication and does not show above run-to-run
noise. The 130.2 s figure was machine load during a concurrent campaign, not the code. Recorded so
that a slow shard going red in a later full-suite run is not attributed here without evidence.
