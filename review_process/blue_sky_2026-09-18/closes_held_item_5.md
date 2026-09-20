# The blue-sky review answers held-for-ruling item 5 (2026-09-13), which asked for exactly one diagnosis run

**This file makes no ruling and changes nothing.** It records that a question the author parked as
`undetermined` now has a measured answer, so the item can be classed. The ruling record itself
(`review_process/held_for_ruling_2026-09-13.md`) is untouched.

## What was asked

> **5 — `lattice(all 3)` returns an all-zero VCV on realistic multi-unit cells, under `robust` too.**
> **Kind: undetermined.** It is either the documented all-zero-VCE outcome or a lattice-construction
> defect, and the evidence does not separate them.

The recommended option was **(b)**: *"spend one diagnosis run to establish whether the 254
collinearity drops are the documented containment or a lattice-construction defect — the class
cannot be assigned without it."*

## The answer: neither. It is **exact fit of saturated lattice pieces.**

A lattice gives each treated (cell x bin) piece its own intercept **and** its own dose slope — two
free parameters at `contreatpoly(1)`. A piece holding **`1 + contreatpoly` observations or fewer is
therefore fitted exactly**, and its residuals are not small, they are **identically zero**.

Measured on the broken fixture (`regress` with explicit dummies, `vce(robust)`, model rmse 1.398):

| treated rows | n | mean abs residual | max abs residual |
|---|---|---|---|
| in pieces with <= 2 distinct dose values | **157** | **6.57e-11** | 4.96e-09 |
| in pieces with >= 3 | 690 | 0.741 | 4.18 |
| untreated / reference rows | 1353 | 1.110 | 4.55 |

157 is not an approximation. The piece-size distribution has 45 pieces of one observation and 56 of
two: `45 x 1 + 56 x 2 = 157` exactly. **18.5% of treated observations contribute a residual of zero.**

## Why that produces the all-zero VCV, and why two engines disagree about it

Every residual-based variance estimator — HC-robust, cluster-robust, and the bootstrap's own
analytic leg — assigns **zero** variance to precisely the directions in which the iid bread says the
variance is as large as **3.18e7**. A matrix holding both is what Stata calls "nonsymmetric or highly
singular". `reghdfe` zeroes it wholesale; `regress` posts a small rank-deficient one (`e(rank)` 381
of 537). That is the whole of the observed behaviour, including the "one path that reports a standard
error reports one 30x too small".

So the caution is firing for a real reason, but **the reason it names is wrong** — it blames
clustering on runs that clustered nothing.

## Why this is an extension of an accepted ruling rather than a new judgment call

It is the **R1 mechanism, one geometry finer.** `ado:12840-12847` already states it for cell
containment: *"the plug-in CRVE cannot see within-cell dose variation … the ame standard errors it
would deliver collapse toward the variance of the absorbed baseline dose control."* R1 fires when a
cell lies wholly inside one *cluster*. Here the cell is split by the *lattice* into pieces wholly
determined by their own parameters. Same disease; the package already ships the cure (the R1 scrub).

## The predicate, with no threshold to argue about

```
count(treated observations in a (cell, bin) piece) <= 1 + contreatpoly
```

On the broken fixture: **101 of 273 pieces**. On a healthy one: **0 of 234** (minimum 7 distinct dose
values per piece, median 13). The separation is clean, so no tolerance has to be calibrated.

## Where the detail and the proposed remedy live

`fixprop/proposals_LAT.md` — §0.6 for the diagnosis above, and proposals L1 (make the cautions state
the cause that actually applied; free, moves no number) and L3b (extend the R1 scrub to lattice
pieces). Both are still subject to the adversarial round; read them as proposals, not conclusions.

---

## CORRECTION, 2026-09-19 — the predicate above is wrong; the diagnosis is not

This file was written before the adversarial round on `proposals_LAT.md` returned. Two reviewers
then falsified the predicate it advertises. The correction is recorded here rather than by rewriting,
so the record shows what was claimed and when.

**What still stands — the diagnosis.** Saturated (cell x bin) lattice pieces are fitted exactly and
their residuals are identically zero; every residual-based variance estimator is therefore blind in
precisely the directions where the iid bread says the variance is large, which is what produces the
all-zero VCV and the engine-dependent disagreement. The measured residual table (157 rows at
6.57e-11 against 690 at 0.741) stands, and the `45 x 1 + 56 x 2 = 157` identity stands. The
false-clustering-blame defect was **independently confirmed**: the shipped all-zero-covariance caution
fires with text blaming clustering on runs that clustered nothing, on **14 of 24** configurations the
reviewer built. That is a real, reproducible defect and it is what proposal L1 fixes.

**What is wrong — the "threshold-free predicate".** The claim above was
`count(treated obs in a (cell,bin) piece) <= 1 + contreatpoly`. It is not right in general:

- A piece with **n = 2 and only ONE distinct dose value** has its slope dropped for collinearity, so
  it carries 1 free parameter for 2 observations and its residuals are **not** zero. The correct
  saturation condition is `n <= 1 + min(poly, nd - 1)`, where `nd` is the number of distinct dose
  values in the piece. That needs a `nunique()` pass, **not** the cheap `count()` pass this file
  priced as free.
- **The saturated-piece count does not order the outcomes.** Measured across 11 (fixture x lattice)
  combinations: a 20-unit panel at `lattice(all 3)` dies at **1.3%** saturated pieces and a 40-unit
  panel at `lattice(all 6)` dies at **0.9%**, while a discrete-dose 40-unit panel at `lattice(all 5)`
  **survives at 3.8%** with a healthy `e(V)` and an `ame` se of .0542 (t 6.13). So "101 of 273 versus
  0 of 234" is a clean separation *on those two fixtures*, not a general rule, and no threshold on
  this count is safe to ship.

**And the remedy this file pointed at is now DONT_SHIP.** Proposal `L3b` was rejected by both rounds:
the one-line implementation the author would have been ruling on **demonstrably undoes L1**. L1
itself (make the cautions state the cause that actually applied; free, moves no number) survives with
named changes — see `fix_proposal_status.md` §4.

**Net effect on held-for-ruling item 5.** The item can still be classed: the behaviour is neither the
documented all-zero-VCE outcome nor a lattice-construction defect, but exact fit of saturated pieces,
and the caution that reports it names the wrong cause. What this review does **not** deliver is a
cheap, safe screen for it. Option (b) is answered; a screen remains open work.
