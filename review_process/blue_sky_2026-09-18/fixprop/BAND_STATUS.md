# `simultaneousbands` / studentized inference — built, verified, and ruling-blocked

**Dated 2026-09-19.** Records the state of the BAND line after the author rejected BAND-1, plus a
finding that reframes what any fix here can achieve.

## The author's directive, which governs everything below

> "Absolutely not! No modifying SEs in this fashion, no! We want the actual SEs that get sampled,
> even if there is something wrong with them. If a given run is unidentified or cannot be estimated,
> then we should report all missing in the bootstrap case and let users know about non-viable draws
> that occurred, but absolutely never ever heuristically edit a SE or point estimate so that users
> get an edited one that looks real."

**BAND-1 is dropped.** It substituted a floor value. A package-wide audit found it would have been
the *only* such site: across 17,021 lines no standard error or variance is substituted anywhere —
no clamp, no floor, no winsorisation, no `max(se, eps)`, zero `makesymmetric`. Every degeneracy
guard blanks, refuses or cautions.

## What was built instead (BAND-2b)

One viability verdict per replication. Non-viable draws are withheld from the inference kernels by
masking a **copy**; the sampled values themselves are untouched. Build and patch at
`probes/rev6-band2fix/`.

**Verified:**
- `e(bootdist_se)` and `e(bootdist_pwse)` **bit-identical to pinned** on all eight arms, including
  ones where the earlier build had lost 373, 1486 and 627 cells.
- Printed output **identical to pinned** on every arm. (The earlier build added a false sentence to
  every healthy `boottstat` run — a bare `if` inserted between an `if`/`else if` pair re-bound the
  `else`. Fixed by moving the clause below the pair with its own test.)
- Clean runs: 30+ `e()` matrices per arm identical; the only change anywhere is `e(boot_novar)`
  absent -> 0.
- Its own first attempt was caught and corrected: v1 restored `e(bootdist_se)` but missed
  `e(bootdist_pwse)`, the pairwise-contrast draws.

**Not preserved:** `jointW_*`/`advW_*` columns in `bootsave()`, `e(boot_jointW)`, `e(boot_advW)`,
and the cell `_gtse`/`_cfx<k>se` accumulators — computed inside the replication and not maskable
afterwards. So "bootsave is untouched" would be wrong: the `se_*` columns are preserved exactly, the
joint-statistic columns are not.

## THE FINDING: this does not make the band usable

The 1e14 blowup is a **symptom, not the disease**.

- Design H (eight never-treated units): critical value **29.07 against 3.4-4.5 healthy, with ZERO
  non-viable draws.** Dropping draws could never have fixed it.
- The SC design has **no 1e14 anywhere** and its bands still fall 78.2/38.2/245.7 -> 12.5/17.0/14.6.
- Design D: 3.85e14 -> **155.5**. The absurdity is gone; 155.5 is still ~30x healthy, and 11
  retained draws carry two thirds of it.

### The mechanism, measured

A hypothesis that this was the documented analytic pairwise-SE undercoverage appearing inside the
studentization was **tested and refuted**. Arithmetic validated first by rebuilding the package's own
critical values from saved draws (29.074 vs 29.07; 155.460 vs 155.46).

`mean(se_r)/sd(b_r)` — the replications' own analytic SEs against realized dispersion — is **0.63 to
1.23 across ten arms, never the 0.1-0.25 undercoverage would require.** On the decisive arm it is
**1.226**, and the observed analytic SE is **1.54x** the bootstrap dispersion — conservative, the
opposite of undercoverage. The ratio does not track the critical value at all.

What it actually is: **a tail in the denominator, one column and one draw at a time.** Median
`max|t|` on design H is 0.96 — perfectly ordinary — while q95 is 29.1. Decomposing the draw that
*sets* the critical value:

| arm | max\|t\| | numerator | denominator |
|---|---|---|---|
| H, full window, cluster at `i` | 29.07 | 1.24 | **3.05e-02** (33x too small) |
| D (the 1e14 case) | 155.46 | 2.24 | **8.34e-03** (120x too small) |
| healthy control | 3.63 | 1.94 | 5.17e-01 |

Numerators are ordinary everywhere. Trimming draws whose worst column has `se_r/se_obs < 1e-1`
takes design H from **29.07 to 1.2** and design D from 155.46 to 32.3, while the healthy controls do
not move at all — there is nothing to trim.

So the residue is the **same per-column near-degeneracy that produces the exact-zero collapse,
operating continuously instead of at the 1e-14 limit.** An all-or-nothing *per-replication* verdict
structurally cannot reach it, because the affected draws are viable in every other column.

**Honest label for BAND-2b: it removes the numerical absurdity and it does not make the band usable
on a near-degenerate design.**

### Cluster geometry helps, in the documented direction, and does not explain it

Same data, full window: **below `i` 6.34 < at `i` 11.00 < above `i` 17.12.** The existing guidance
works. But even below `i` is 6.3 against 3.6-4.5 on healthy designs, and on healthy designs the three
geometries give 3.63/4.13/4.51 — geometry barely matters when the design is not near-degenerate. The
**full event window** is the larger lever: same data and geometry, windowed 17.7 vs full 29.1.

### Stated limit

This measures the symptom, not its origin. Whether the residue is a genuinely new defect or a new
face of the known analytic-SE pathology cannot be settled from outside the margins/pwcompare pass.

## The two rulings this needs

**1. Which half of the directive wins for a non-viable draw?**

| | kernels | sampled values | consequences |
|---|---|---|---|
| **(A)** as built | exclude non-viable draws | **every sampled value survives** in `e(bootdist_se)`, `e(bootdist_pwse)`, `bootsave()` | `jointW`/`advW` cells still drop; no per-row flag saying which draws were excluded; `e(bootdist_se)` no longer reconciles row-for-row with `e(boot_tstud)` |
| **(B)** | exclude | **missing**, sampled values included | none of the above arise |

Reverting to (B) is a one-token switch.

**2. The tolerance is 1e-8, and no tolerance-free predicate exists.** Five candidates were refuted by
measurement, including `e(r2) >= 1` (catches design D's 25 non-viable draws and **none** of design
E's 46). Bracketing: ~3.5 orders of headroom above 1e-8 before wooldid's own all-zero-VCV refusal
takes over, ~6 orders below it to the collapsed draws. Either accept 1e-8, or accept that no screen
ships.

**A third question the finding raises.** Because the damage is per-column and continuous, the shape
that would actually repair the band is a **per-column** exclusion, not a per-replication one — which
is the structure BAND-2 replaced. Note that excluding a column's contribution is not flooring, so it
is not what was rejected. Nothing has been built along those lines and it should not be until the
author says whether an unusable band on near-degenerate designs is worth further work at all.

## Untested
A collapse confined to one subgroup passing a max-over-the-diagonal verdict; partial collapse inside
one replication (some columns at 1e-14, others ordinary) — never observed by any agent, so the max
rule has never been exercised against it; `fweight` and `if` under collapse; poisson.
