# Feature tickets: decoupling `i` from the resampling unit

**Status 2026-09-16: T1, T2 and T3 are COMPLETE, merged and blessed. T4 CLOSED by ruling.**
All three are implemented, verified against a preregistration, merged to `bugfixes`, documented in
`wooldid.sthlp`, and covered in the suite from both directions - the nesting branch pinned so it
cannot move, the new keep-`i` branch pinned so it cannot regress, and each of the three retired
refusals given a config measured at r459/r459/r198 before the change and rc 0 after. The full
assessment suite is deferred to 2026-09-17 by the author; it is the only outstanding item, besides
the `fcol_sg` exposure recorded below.

**Superseded status line, kept for the record: T1 COMPLETE (worktree 84da86f, unmerged). T2+T3 IN PROGRESS. T4 CLOSED by ruling (documentation only).**
Implementation opened 2026-09-16; the survey below ("What the code actually assumes") was written
before any edit and changed the shape of the work: **T2 and T3 are one edit, not two**, and T2 is
substantially smaller than the ticket assumed. Preregistered readouts:
`test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md`, written before any line of `wooldid.ado` moved.

Opened 2026-09-15 at the author's instruction, after a help-file edit turned up the underlying
issue. These are FEATURE requests, not defects, so they do not belong in `open_defects.md` -
with the partial exception of T4, which the help already documents as a correctness limitation.

## The principle these all share

`i` is an **estimand-definition** choice: it says which cohorts the effects are averaged over.
`cluster()` is a claim about **dependence and assignment**: it says which observations are
correlated, and for the replication engines, what gets resampled or re-randomized. The two are
orthogonal, and the help now says so at the `i` entry (`8ca8ba0`): choose `i` by treatment
heterogeneity, not by the level at which you want fixed effects.

But three code paths use `i` as a **proxy** for the resampling unit and refuse anything finer:

| site | rule | gate |
|---|---|---|
| `wooldid.ado:5990` (`wooldid_permctl`) | `cluster()` must nest `i` | `permute()`, always |
| `wooldid.ado:6364` | same | `bootstrap(standard)` with a cell-weighted estimand |
| `wooldid.ado:6365` | same | `bootstrap(standard)` with `boottstat` and `cluster()` |

All three test the same thing - is the cluster id constant within each `i`? - so they refuse
BOTH cross-cutting clusters AND clusters strictly nested inside `i`. None can be bypassed with
`safety(off)`.

**Why that matters now.** The analytic path wants the opposite: cluster BELOW `i`, so several
clusters span each treated cohort-by-period cell (the help's "intended practice", and the
difference between 4.5 percent and 72.5 percent rejection at a nominal 5 percent in the
adversarial size study). At or above `i` the analytic covariance is "the wrong object rather
than a mis-scaled one" and there is no plug-in correction, so the help's stated remedies are a
finer `cluster()` **or** the replication engines. These tickets are what make the second remedy
available at the coarse `i` the help now recommends. Today a user who follows both pieces of
advice - coarse `i`, finer cluster - cannot use `permute()` at all.

**The refusals are not wrong today.** Each guards a real mechanical failure. The tickets are
about removing the failure, not the guard.

---

## What the code actually assumes (survey, 2026-09-16, before any edit)

The tickets were written from the three places that *check* nesting. Implementation began by
asking where nesting is *assumed* - a larger set, and the difference reshaped T2 and T3.

| site | what assumes cluster nests `i` | ticket |
|---|---|---|
| `ado:5990` (`wooldid_permctl`) | the permute check itself | T1 |
| `ado:6353` (`wooldid_bootctl`) | `error 198`: cell-weighted estimand with no `cluster()` | T2 |
| `ado:6362-6365` (`wooldid_bootctl`) | the bootstrap check itself, two distinct messages | T2/T3 |
| `ado:9107` (`wooldid_bootrun`) | **`__wdbnewi = group(__wdbns, i)`** - the actual mechanism | T2/T3 |
| `ado:9135` (`wooldid_bootrun`) | `subinstr` of the rep VCE onto `cluster(__wdbns)`; unchecked | T3 |
| `ado:7352` (`wooldid_bootctl`) | disclosure text asserting "clustering at or above i" | T3 |
| `ado:11968` (F13) | `subgroups()` constant within the `(i,t)` cell | T1 |

### Three findings that changed the plan

**1. The re-id leg already does most of T3.** `bsample, cluster(clusv) idcluster(__wdbns)`
already mints a fresh cluster id per draw, and `ado:9135` already substitutes the replication's
own VCE onto it. That IS the `2 - 1/G` fix, shipped 2026-09-04. What breaks under finer
clustering is a single line, `ado:9107`.

**2. T2 and T3 are therefore the same edit.** `__wdbnewi = group(__wdbns, i)` is correct under
nesting (each `i` sits in one cluster, so this mints one new unit per drawn copy) and wrong under
finer clustering, where a unit spanning two clusters is **split into two units in every
replication** whether or not anything was duplicated - unit fixed effects silently become
per-(cluster, unit). The rule that fixes both tickets is one two-regime branch:

- cluster nests `i` -> `group(__wdbns, i)`, exactly as today;
- cluster does not nest `i` (finer or cross-cutting) -> keep `i` unchanged, and let a twice-drawn
  sub-unit cluster appear as duplicate rows inside its unit's cell.

The nesting branch must stay `group()`. Keeping `i` there would make a twice-drawn cluster's
cells identical duplicates that the equal-cell weighting counts ONCE, collapsing the composition
variance of `att_i` toward zero. Detection is the `cmin`/`cmax` computation already standing at
`ado:6360`; it becomes a flag passed down rather than an error.

**3. T2's headline - "recompute the weights per replication" - is already true on this path.**
The re-id leg's replication is a full nested `wooldid y __wdbnewi t ttrevar, ...` call
(`ado:9137`), so `att_i`'s weights are rebuilt from whatever data the replication holds. T2 does
not need new weighting machinery. Its real content is the gate relaxation, the `__wdbnewi` rule
above, and the `error 198` relaxation for the unclustered case.

### Ruled: no display line for the rebuild (author, 2026-09-16)

The implementation proposed a one-line `part(perm)` footer saying that each replication
regrouped `i` by the onset it drew. **Ruled against.** Rebuilding `i` from the drawn
assignment is the natural behaviour of permutation inference in this setting rather than a
caveat about it, and the table does not need furniture saying so. `e(perm_irebuilt)` remains
posted for anyone who wants to check which of the two paths a run took, and the
`permonsetcohorts` help entry explains the rebuild in full. Not to be reopened.

### Five refinements from the T1 implementation (2026-09-16, measured)

**(a) The rebuild opens no column-loss channel at all.** The survey said column loss "IS new under
T1, and only under `permute(probabilities)`", with `ado:7556` as the net. Measured: no accumulated
`e(b)` column name embeds an `i` code (relabelling the cohort variable leaves `colfullnames`
byte-identical), and the realized relative-year set depends on the drawn onsets rather than on how
units are grouped. `e(gteffects)` is the only `i`-keyed surface and `permute()` blanks it outright.
Under `probabilities` what shrinks is the level set of `i`, which moves the per-replication
estimand, not the column set.

**(b) Under the assertion, "cluster nests i" collapses to "cluster == i".** `wooldid_permctl`
already requires every cluster to sit inside one onset group, and the assertion makes the onset
grouping BE `i`. So cluster is always a subset of `i`, and only two cases are constructible: equal,
where the rebuild is skipped, and strictly finer, where it is live. The ticket's "cluster at or
above `i`" has exactly one reachable case here, which is what makes the P3 no-op control exact
rather than approximate.

**(c) T1 needs nothing from T2.** The open question "estimands averaging over `i` must then be
recomputed per replication - see T2, which this shares machinery with" is already answered: a
permutation replication is a full nested `wooldid` call, so `att_i`'s weights are rebuilt from
whatever data the draw holds. The two tickets share a principle, not machinery.

**(d) Sharp-null exactness cannot distinguish a correct rebuild from a correctly relabelled one.**
A relabelling is still a deterministic function of the draw, so the p-values stay uniform. The
prereg's F5 named the case "P1 fails while P2 passes" but had no readout for its converse. The
rep-0 assertion is the only thing that covers it. Recorded as the prereg's second F2 amendment.

**(e) The rebuild's reserved never-treated code is unreachable.** A panel with no never-treated
units is a design `wooldid` rejects in data setup, which runs before the replication loop. Pinned
by cert Q1e, which asserts the option's screens do not pre-empt that rejection.

### One design question closed

T3 asked whether the estimator can carry two copies of the same `(i,t)` cell distinguished only
by cluster id. It already does: when `i` is a cohort every cell holds many rows, and the collapse
at `ado:1860` replaces the cell with its estimand-weighted mean. Duplicates are more rows
carrying more weight. **Resolved, not open.**

### T1's option name, as implemented

`permonsetcohorts` (declared `PERMONSETcohorts`, minimum abbreviation `permonset`), by author
ruling 2026-09-16. The implementation first used `permiisonsetgrouping`, encoding the assertion
"perm: {it:i} is onset grouping"; the name was put to the author because it is permanent
user-facing API, and the ruling preferred naming what {it:i} must BE over asserting a sentence
about it, the run-together "iis" reading as a typo. Registered
in the worker syntax line, the validate syntax line, the c_local blanket and the `permstray`
screen; `check_validate_sync.py` passes all three of its checks.

### One risk that is new only to T1

The accumulator matches columns by NAME against the observed `e(b)`, and a missing column becomes
a blank that `ado:15683`'s `select(draws, draws :< .)` drops from the denominator. Under the
bootstrap's keep-`i` regime this is not new - a unit vanishes only when all its clusters go
undrawn, against today's vanishing when its single cluster goes undrawn, so the accounting is
unchanged. Under T1 it IS new, and only under `permute(probabilities)`: `shuffle` permutes every
onset value across clusters and preserves the level set. The tripwire already exists at
`ado:7556` (rep-0 `colfullnames` vs the observed run).

**MEASURED 2026-09-16, and it narrows the paragraph above.** No accumulated `e(b)` column name
embeds an `i` code - relabelling `i` leaves `colfullnames` identical - so `ado:7556` detects a
changed column SET but cannot detect a wrong onset-to-cohort map. T1 therefore carries an explicit
assertion in `wooldid_permrun` that rep 0's rebuilt `i` equals the observed `i` on the design
sample, measured to fire r459 on a deliberately reversed map. See the F2 amendment in the prereg.


## T1 - `permute()`: recompute `i` after assignment when `i` is an onset grouping

**Current.** Each replication draws one onset per cluster and writes it to every row of that
cluster (`wooldid.ado:5955-5958`). With clusters nested inside `i`, a unit receives two onsets in
one replication, which violates the invariant wooldid enforces at parse time - `ttre` single-
valued within `i` - so the unit would be silently split across cohorts and each replication would
compute a different estimand from the observed statistic it is compared against.

**Proposed.** A user-declared option asserting that `i` was constructed as a treatment-onset
grouping. Under it, rebuild `i` inside each replication from the drawn onsets: after
reassignment, "units sharing the drawn onset" IS the onset grouping, so `ttre` is constant within
the rebuilt `i` by construction and the invariant holds. Then allow `cluster()` finer than `i`.

**Why an option rather than inference.** wooldid cannot tell from the data whether `i` is an
onset grouping or a substantive unit definition that happens to coincide with one. Getting this
wrong silently redefines the estimand, so it must be asserted.

**Design questions.**
- The number of levels of `i` can change between replications when no cluster draws a given
  onset. Estimands averaging over `i` must then be recomputed per replication - see T2, which
  this shares machinery with.
- Interaction with `permstrata()`: strata already constrain which onsets a cluster can draw.
- Does the rebuilt `i` need to be reported anywhere, or is it purely internal?

**Side benefit worth having on the record.** At `cluster() = i` the help cautions that
`permtstat` is unprotected, because the studentized statistic only guards against atypical
treated clusters when the cluster-robust standard error can see their variance, which at one
observation per `i`-period it cannot. Allowing a finer cluster makes `permtstat` actually
protective instead of nominally available.

---

## T2 - recompute cell-weighted estimand weights per replication

**Current.** `bootstrap(standard)` with `att_it`, `att_i` or `att_itime` requires `cluster()` to
nest `i` (`wooldid.ado:6364`), because the re-id leg duplicates drawn clusters as whole units and
a fragment of a unit cannot be re-identified as one.

**Proposed.** Recompute the estimand weights inside each replication rather than carrying the
original sample's. `att_i` is "the simple average of cohort-specific treatment effects"
(`wooldid.sthlp:226`), which is a functional of whatever data the replication holds - so
recomputing it over the cohorts present in that resample is the ordinary bootstrap principle, not
a compromise. With weights recomputed, finer-than-`i` clustering no longer breaks the estimand.

**Also covers.** `bootstrap(standard)` with a cell-weighted estimand and NO `cluster()` is
currently refused outright (`wooldid.ado:6353`, error 198: "observation-level resampling cannot
re-identify duplicated treatment cells"). Once weights are recomputed per replication, revisit
whether that refusal is still necessary or can be relaxed to observation-level resampling.

**Design questions.**
- Whether `permute()` needs the same treatment, or whether T1's rebuilt `i` already suffices
  there.
- Whether the recomputed estimand should be reported when it differs materially from the
  original sample's, so a user can see that the cohort set moved.

---

## T3 - `boottstat`: cluster re-identification for finer-than-`i` clustering

**Current.** Refused under `cluster()` (`wooldid.ado:6365`) because the physical expansion the
studentized bootstrap needs is undefined when a unit straddles clusters.

**The constraint any implementation MUST respect, and it is measured.** The re-id leg exists
because a twice-drawn cluster has to count as TWO clusters in that replication's own
cluster-robust standard error. Frequency weights cannot express "this cluster is here twice as
two clusters": under the old count leg the copies pooled into one score and the cluster-robust
meat took `m^2` where it needed `m`, giving a measured inflation of `E[m^2]/E[m] = 2 - 1/G`,
about `sqrt(2)` in the standard error, which broke percentile-t pivotality
(`wooldid.ado:367-375`, author ruling 2026-09-04). **Any scheme that allows finer-than-`i`
clustering must keep duplicate draws distinguishable at the cluster level or it reintroduces
exactly this.**

**Proposed.** Re-identify on the CLUSTER rather than on the unit: a duplicated sub-unit cluster
keeps its unit identity and receives a fresh cluster id, so the replication's SE clusters on
resampled cluster identity as required while units stay whole.

**Design questions, and this is the one with real unknowns.**
- Can the estimator carry two copies of the same `(i,t)` cell distinguished only by cluster id?
  That is the crux; if it cannot, the alternative is to relabel units too and accept partial
  units, which changes the `(i,t)` cell structure the estimator is built on.
- How unit fixed effects behave when a unit's rows are duplicated unevenly across its
  sub-clusters.
- Verification should reproduce the `2 - 1/G` result as a negative control: a deliberately
  wrong implementation should show the `sqrt(2)` inflation, and the right one should not.

---

## T4 - `impvalues`: refuse the incoherent combination rather than reweighting it

**REVISED 2026-09-15 after the author's assessment, which was right.** This ticket originally
proposed recomputing estimand weights under IM-style inference, by analogy with T2. That is the
wrong fix and the reason is worth keeping.

**Why reweighting cannot work here.** IM computes cluster-specific treatment effects, presents
their AVERAGE as the estimate, and t-tests the cluster-specific estimates against 0 on `q-1`
degrees of freedom (`wooldid.sthlp`, `impvalues` entry). The estimand is therefore *structurally*
the equally-weighted mean over clusters. That equal weighting is not incidental: Ibragimov and
Muller's validity result requires the group estimates to be independent, approximately normal
and to share a COMMON MEAN, with possibly heterogeneous variances, and it is the t-statistic of
their EQUALLY weighted mean that gets the conservative `t_{q-1}` reference at small `q`.
Reweighting the cluster estimates to reproduce `att`'s observation or cell weights would recover
the point estimand but forfeit that justification - it would need a Satterthwaite-style
approximation and would give up precisely the small-`q`, heterogeneous-variance robustness that
is the reason to use IM instead of plain cluster-robust standard errors. **Reweighting buys the
estimand by spending the inference.**

So IM coheres with a wooldid estimand only when `cluster()` = `i` and the estimand is `att_i`,
which is exactly what the help observes at the cohort-level case ("IM-style inference will cause
`att` to yield the same effect estimates as `att_i`"). That is a necessary preference, not an
accident to be engineered away.

**What IS actionable, and it is a defect rather than a feature.** The help says non-`att`
estimands under IM are "incorrect ... if clusters cross-cut cohorts or (worse) cohort-time period
cells", and wooldid computes and prints them anyway. Shipping numbers the documentation calls
incorrect, with no error and the caveat living only in the help, is the real problem here.

**RULED AND CLOSED 2026-09-15: documentation only. No refusal, and no runtime caution.**

Both code options were put to the author and both declined. The non-`att` estimands keep
printing under `impvalues` with cross-cutting clusters, and the help is the whole mitigation.
Recorded because a reader will otherwise assume this was an oversight: it was considered, the
options were costed, and the author chose the documented behaviour.

What was declined, so nobody re-proposes it:
- **Refusal.** There was precedent - `impvalues` is already refused with `subgroups()`,
  `jointtests`, `makecfxplots`, `adversarial`, `customvce()`, `fastnose`, `pretrendpower()`,
  `permute()` and `bootstrap()`, nine screens - and the author still declined a tenth. This is
  the standing position on not guarding specifications, applied to a case that sat near its edge.
- **A runtime caution at the table.** Declined even though it was free: all 15 suite configs
  using `impvalues` request `att` only, so nothing would have triggered it and no caution counts
  are pinned in `degeneracy_expected.csv`. Cost was not the reason.

**What WAS done instead, and it is the part that mattered** (`2d41a22`). The help entry said IM
"implicitly changes wooldid's target estimand" and then offered cohort-level clustering as an
example, which reads as a special case a user might avoid. It is the general rule. The entry now
leads with it - the procedure averages cluster-specific estimates with equal weight, so the
cluster-equally-weighted mean is what it reports in every run - names `cluster()` at `i` with
`att_i` as the pairing coherent by construction, and records why the equal weighting cannot be
reweighted away: it is what earns the `q-1` degrees of freedom, so trading it for the requested
estimand's weights would cost the inference that motivates the method.

So the transparency item is delivered in documentation rather than in `e()` or at the table,
which is where the author wanted it.

## Considered and NOT proposed, with reasons

**Multi-way `cluster()` under the replication engines.** (Author 2026-09-15: confirmed, not pursuing.) Refused by both today
(`wooldid.sthlp:1220`: "Both engines resample or relabel ONE unit, so there must be exactly
one"). This is a genuinely different problem - you cannot resample two crossed dimensions with a
single-dimension scheme - and the methods for it are their own literature, not an application of
the principle above. Worth its own investigation if wanted, but it does not belong in this group.

**`subgroups()` requiring the subgroup constant within `i`.** (Author 2026-09-15: confirmed, keeping the requirement as it stands.) The remedy the help gives is
defining `i` finer (`wooldid.sthlp:869`, `egen newi = group(i sg)`), which does pull against the
new coarse-`i` recommendation. But this is a real estimand constraint rather than a conflation:
if the subgroup varies inside an `i`-by-`t` cell, that cell's effect is not subgroup-specific and
no reweighting recovers it. Left alone deliberately.


## The `fcol_sg` exposure, and a correction to how I measured it (2026-09-16)

Retiring `fastcollapse` turned `fcol_sg` from a row pinned r198/r198 into an ordinary dual-build
row, so its `e()` reached content comparison for the first time since it was written. The suite's
comparer reported it DIFFERS. The author ruled to re-point the config off its straddling-subgroup
design rather than ledger the divergence.

**The re-point is sound and is not the problem.** `region` is constant within `fips` and straddles
`cohid`, `adm_sgsplit` already pins the straddle as its own r459 refusal, and the re-pointed command
runs without `safety(off)` at all, which the original needed only to bypass that refusal.

**What went wrong was my verification, and it is worth writing down.** I tried to confirm the
re-point by comparing `e(b)` and `e(V)` bitwise at `%21x` between the staged v1.3 build and HEAD.
It kept reporting a divergence, through two fixtures:

- `i=fips` at a thinned subsample: 120 of 3660 cells. That subsample turned out to give exactly
  2.00 counties per state, a thin-cluster degeneracy I had introduced myself while trying to cut
  runtime - the very thing the author ruled against testing by accident.
- `i=fips` at the deck's own `adm_sub` prep, 6.00 counties per state: 120 cells again.

The third measurement is the one that matters. Running the **exact `adm_base` command** through the
same harness also differs, 100 of 1225 cells - and `adm_base` carries no `b`/`V` ledger rows and is
clean in the suite. So the harness, not the config, was the outlier: a bitwise `%21x` comparison is
strictly stronger than the suite's comparer, which applies bounded storage-noise tiers and named
ledger classes. Every "still differs" I reported was measuring a difference the suite would have
tolerated.

**Consequence.** Nothing about `fcol_sg` should be decided on those three probes. The authoritative
reading is the `oldnew` component with the real comparer, which is deferred with the rest of the
full suite to 2026-09-17. `fcol_sg` is committed unchanged apart from the flag removal, so the
config is in a valid state meanwhile.

**The lesson, and it is the session's fourth instance of the same shape.** Three preregistered
checks were credited with catching failures they could not catch. This is the mirror image: a check
strong enough to report failures that are not there. Either way the fix is the same - establish what
the instrument does on a case whose answer you already know. One run of `adm_base` through my
harness at the start would have saved three rounds.

---

## `fcol_sg`: the authoritative reading, and why the ruled remedy does not address the cause (2026-09-16)

The full assessment suite ran on 2026-09-16 (569 configs, pool 8) and gives the authoritative
`oldnew` reading the section above deferred. **`fcol_sg` is the only DIFFERS in the deck**:

```
VERDICTS DIFFERS=1 EXPLAINED=414 MATCH=26 RC_EXPLAINED=128
COMPARED 569 of 569        SUMMARY unexplained=10 ... stale_ledger=0 unscoped=0 prov_fail=0
[fcol_sg] DIFFERS (+89 new-only objects)
    UNEXPLAINED scalar joint_att_it_allsgssame      0.0006310002949791 -> 0.0006310003791699
    UNEXPLAINED scalar joint_pre_att_it_allsgssame  0.6189793806821372 -> 0.6189793788099903
    UNEXPLAINED matrix b   30 cells, worst reldiff 5.04e-08
    UNEXPLAINED matrix V   30 cells, worst reldiff 1.40e-08
    UNEXPLAINED matrix pwsgs 72 cells, worst reldiff 1.60e-07
    UNEXPLAINED matrix wdidmainresults / wdidsgresults_1..4, worst reldiff 2.70e-07
```

So the divergence is real by the real comparer, not an artifact of the `%21x` harness that
misjudged it before. The comparer applies a **global** `--rtol 1e-10` with no absolute floor plus
per-config ledger rows; `fcol_sg` carries no ledger rows, and neither does the passing `adm_base`.

**But the cause is not the straddle, and so the ruled re-point cannot fix it.** Three measurements:

1. **The degeneracy guard never fires, on either build.** Running the config's exact command against
   the staged `build_old` and `build_new`: both return rc 0, both `e(wdid_fcguard)` **empty**, both
   `e(N) = 4162`. The premise behind re-pointing - that the straddle drives the two builds down
   different collapse branches - is not what is happening. Both builds collapse, identically.

2. **`adm_base` is the natural control and it passes.** Same dataset, same `adm_sub` prep, same
   `subgroups(region)`, same `cluster(state)`, same `safety(off)`, same `i = cohid`, same straddle.
   It differs from `fcol_sg` essentially in one respect: `fcol_sg` adds **`att_it`**. And every
   differing object above is an `att_it` surface - `joint_att_it_*`, `pwsgs`, `wdidsgresults_*`.

3. **The discriminating run settles it.** `att att_it` on `adm_sub` with **no `subgroups()` at all**,
   so no straddle is in play anywhere:

   | arm | old | new | verdict |
   |---|---|---|---|
   | A: `att att_it`, no subgroups | `+1.4b320b8e68b80X-001` | `+1.4b320b8828b20X-001` | **differs, ~1e-8** |
   | B: `att` only, with subgroups (the `adm_base` shape) | `-1.2916b1e875280X-004` | `-1.2916b1e875200X-004` | agrees to the last hex digit |

   A cell-weighted estimand diverges between the two builds with the straddle entirely removed;
   the straddle with no cell-weighted estimand does not diverge.

**Reading.** The cause is the numerical path for **cell-weighted estimands under the collapse** -
consistent with `1eea605`, "fastcollapse: evaluate the margins calls on one row per cell", which is
in HEAD and not in the published v1.3 that stages as `build_old`. That is a deliberate, already
shipped optimization, and a change of this shape moves cell-weighted results in the last few digits
by construction. `fcol_sg` is simply the only config in the deck that reaches content comparison
with `att_it` on data where the collapse is active.

**Two facts the author did not have when he ruled**, both of which bear on the remedy:

- Re-pointing `i: cohid -> fips` removes the straddle but **keeps `att_it`**, so by measurement (1)
  and (3) it should not remove the divergence. Tested directly rather than inferred; result recorded
  below.
- The re-point is not free. `i = fips` is **288 levels** on `adm_sub`, which trips wooldid's own
  "> 100 distinct values" runtime caution, and at Stata's default `maxvar` of 5000 the re-pointed
  command **fails outright with r(907)**. It only runs at all because `mkrun.py:137` sets
  `maxvar 20000`. It is also far slower - the runtime cliff that same caution describes - and this
  config runs twice per dual build. **Measured, paired, same session, same data:**

  | command | seconds |
  |---|---|
  | `fcol_sg` as it stands today (`i = cohid`) | **2.663** |
  | the ruled re-point (`i = fips`, 288 levels) | **522.189** |

  A factor of **196**. `oldnew` runs every config on both builds, so this one row would go from
  about five seconds of suite time to roughly seventeen minutes, on a component that already
  takes fifty. That is a large price for a change which, by measurements (1) and (3) above,
  does not address the cause.

**This wants a fresh ruling**, because the original one ("re-point rather than ledger the
divergence") was made on the understanding that the straddle was the cause. It is not.

### The remedy tested directly — it works, and that is the problem

I predicted above that the re-point would not remove the divergence, on the reasoning that it
keeps `att_it`. **That prediction was wrong, and testing it rather than trusting it is what
found the real objection.** Running the ruled re-point against both staged builds:

```
REMEDY old rc=0 ncols=60   b11=+1.5a450f9448eb0X-001  b12=+1.5a450f9448eb0X-001  b13=-1.44ff367212400X-004
REMEDY new rc=0 ncols=60   b11=+1.5a450f9448eb0X-001  b12=+1.5a450f9448eb0X-001  b13=-1.44ff367212400X-004
```

Bit-identical. The re-point does remove the divergence. **But it removes it by making the collapse
a no-op**, which is the one thing this config exists to exercise. Cell occupancy on `adm_sub`:

| `i` | rows per `(i,t)` cell | what the collapse does |
|---|---|---|
| `cohid` (as the config stands) | min 33, **mean 78.3**, max 109 | a real reduction - so `1eea605`'s "one row per cell" changes the numerical path, and the builds diverge |
| `fips` (the ruled re-point) | min 0, **max 1** | **identity**. There is nothing to collapse |

The data is one row per county-year, so an `i` at county level puts a single row in every
cell. `fcol_sg`'s own comment states its purpose - "the collapse WITH subgroups() (author-ruled
2026-09-01)" - and `test_fastcollapse.do` F8-F10 pin the screen, not the arithmetic. Re-pointing
to `i = fips` would leave a config that still *looks* like it covers the collapse, costs 196x the
runtime, and exercises an identity operation. That is worse than either alternative: it is silent
lost coverage, and lost coverage that reads as green.

### What this leaves

The divergence is a genuine, intended, already-shipped numerical change - `1eea605` evaluating the
margins calls on one row per cell - showing up on the only deck config that reaches content
comparison with a cell-weighted estimand on data where the collapse actually reduces anything.
Worst reldiff 2.70e-07.

**There is direct precedent for ledgering exactly this class.** `ledger.csv` already carries
`QR-sgjt-2,cont_sg,scalar,joint_pre_att_allsgssame,1e-9,"subgroup cell collapse class
(2026-09-01 default; blessed in test/baseline 2026-09-05 ...)"`. A bounded row is the house
convention (`A9-svmat-*` at 1e-6, `A17-part-b` at 1e-4) and bounds rather than blankets: a numeric
`rtol` only explains a finite numeric difference within it, so a real regression still fails.

Concretely, what a ledger entry would look like - **not applied, pending ruling**:

| id | config | kind | pattern | rtol |
|---|---|---|---|---|
| `FCOL-attit-b` | `fcol_sg` | `mat` | `b` | `1e-6` |
| `FCOL-attit-V` | `fcol_sg` | `mat` | `V` | `1e-6` |
| `FCOL-attit-pw` | `fcol_sg` | `mat` | `pwsgs` | `1e-6` |
| `FCOL-attit-res` | `fcol_sg` | `mat` | `wdidmainresults` | `1e-6` |
| `FCOL-attit-sg` | `fcol_sg` | `mat` | `wdidsgresults_*` | `1e-6` |
| `FCOL-attit-j` | `fcol_sg` | `scalar` | `joint_*att_it_allsgssame` | `1e-6` |

with 1e-6 bounding the observed worst of 2.70e-07 with roughly a factor of four of headroom,
matching the `A9-svmat-*` rows.
