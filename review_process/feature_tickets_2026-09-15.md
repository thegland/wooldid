# Feature tickets: decoupling `i` from the resampling unit

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

**Proposed.** Refuse non-`att` estimands under `impvalues` when the clusters cross-cut cohorts
(or cohort-by-period cells), or failing that print a caution at the table rather than only in the
help. **There is strong precedent for refusal**: `impvalues` is already refused in combination
with `subgroups()`, `jointtests`, `makecfxplots`, `adversarial`, `customvce()`, `fastnose`,
`pretrendpower()`, `permute()` and `bootstrap()` - nine screens. One more for a combination the
help already calls incorrect is consistent with how this option is governed, not a new policy.

**A second, cheaper item worth considering alongside it.** The help notes that IM "implicitly
changes wooldid's target estimand" even for `att`. wooldid could report which estimand IM
actually delivered rather than printing it under the requested name, so a user reading `att` in
the table knows it is the cluster-equal-weighted average. Pure transparency, no numbers move.

**Author ruling needed on one point.** The standing instruction is not to build guards against
arbitrary bad user specifications. This is arguably a different category - it is not a user
writing something silly, it is wooldid emitting a result its own documentation calls incorrect -
but the distinction is the author's to draw, not mine to assume.

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
