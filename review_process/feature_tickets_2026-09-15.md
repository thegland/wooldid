# Feature tickets: decoupling `i` from the resampling unit

**Status 2026-09-16: T1 and T2+T3 IN PROGRESS. T4 CLOSED by ruling (documentation only).**
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

### One design question closed

T3 asked whether the estimator can carry two copies of the same `(i,t)` cell distinguished only
by cluster id. It already does: when `i` is a cohort every cell holds many rows, and the collapse
at `ado:1860` replaces the cell with its estimand-weighted mean. Duplicates are more rows
carrying more weight. **Resolved, not open.**

### T1's option name, as implemented

`permiisonsetgrouping` (declared `PERMIISOnsetgrouping`, minimum abbreviation `permiiso`) - the
spelling encodes the assertion "perm: {it:i} is onset grouping" rather than the effect. Registered
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
