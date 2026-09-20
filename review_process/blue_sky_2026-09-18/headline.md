# The headline: the documented remedy fixes the standard error and breaks the point estimate — and the complete recipe is one word longer

**Status line, dated 2026-09-18.** Every number below was measured by the coordinator directly, not
taken from an agent. Repros: `coord/remedy/rem.do`, `coord/femc/femc.do` in the review scratchpad.

## The chain

1. **Analytic SEs undercover at or above `i`.** Already known: the project's own 500-draw MC
   (`review_process/se_coverage_2026-09-14/`) measured 95% coverage 0.794-0.876 at `i` = unit with
   `cluster(i)`, against csdid/did_imputation at 0.926-0.939 on identical points.

2. **The documented remedy is to cluster finer.** sthlp:255 names "a finer `cluster()`" as the
   remedy; sthlp:67 states "The intended practice is to cluster BELOW the level of `i`". The whole
   finer-than-`i` feature (T1/T2/T3) shipped 2026-09-16 to make that reachable.

3. **The remedy genuinely works on the standard error.** ~1,500 draws over seven DGPs:
   remedy arm (`i`=cohort, `cluster(unit)`) median SE/SD **0.970**, coverage **0.935**; at-`i`
   control 0.589/0.739; cross-cutting cluster 0.999/0.950. Point estimates are **exactly**
   invariant to `cluster()` — mreldif 0 on every one of 774 draws. The harness reproduced the
   project's own positive control (0.852/0.796 against their 0.838/0.794).

4. **But reaching that geometry on a normal panel silently moves the point estimate.** On a
   one-row-per-unit-per-year panel, the only way to get `cluster()` below `i` is to coarsen `i` to
   the onset wave — which replaces the unit fixed effects with a single wave intercept, so sub-`i`
   permanent heterogeneity stops being absorbed. Zero-noise DGP, true att **exactly 1**:

   | panel | `i`=unit | `i`=wave (the remedy) | `i`=wave + `fe(unit)` |
   |---|---|---|---|
   | balanced | 1.0000000000 | 1.0000000000 | 1.0000000000 |
   | ragged (8% MCAR) | 1.0000000000 | **0.9838810126** | 1.0000000000 |

5. **`fe(unit)` repairs it exactly — and keeps the coverage gain.** This is the measurement nobody
   had made, and it decides what the author does next. 200 draws, true att = 1, ragged panel:

   | arm | mean bias | sd(b) | mean SE | SE/SD | coverage |
   |---|---|---|---|---|---|
   | `i`=wave, `cluster(unit)` (remedy, no `fe`) | -0.0055 | 0.2226 | 0.2204 | 0.990 | **0.970** |
   | `i`=wave, `cluster(unit)`, `fe(unit)` (full recipe) | -0.0026 | 0.1801 | 0.1888 | 1.048 | **0.970** |
   | `i`=unit, `cluster(i)` (known-bad control) | -0.0042 | 0.1824 | 0.1400 | 0.767 | **0.895** |

   `fe(unit)` costs nothing and buys a **24% smaller sampling SD**. The control reproduces the
   documented defect, which validates the harness.

## So the complete recipe is

    wooldid y cohort t ttre, att cluster(unit) fe(unit)

`i` = onset cohort for the estimand, `cluster()` = unit for the dependence claim, `fe()` = unit to
keep the fixed effects the coarse `i` gave up. **sthlp:255's remedy sentence names only the middle
one.** sthlp:64-70 does document the `i`/`fe()` pairing, and the SI-8 runtime caution at ado:1371
mentions `fe()` — but both present it as a modelling *preference*, and the Standard-Errors passage
that actually pushes users into the geometry omits it entirely.

## On the magnitude — report the range, not the headline number

The per-draw point error from omitting `fe()` is **not a systematic bias under MCAR attrition**: over
200 draws the mean bias is -0.6%, while the sampling SD inflates by 24%. The large single-draw
numbers an agent reported (+21%, -57%) came from attrition landing on a *treated* unit, where it is a
genuine bias rather than noise. The honest statement is: **omitting `fe()` costs efficiency always,
and costs unbiasedness when the raggedness is correlated with treatment.** A single zero-noise
demonstration overstates the typical case; this is the project's own recorded lesson about checking
order of magnitude against mechanism.

## What this implies for the author

The finer-than-`i` feature shipped one day before this review and its assessment suite was deferred.
On the **analytic** path it does what it was built to do. Three things remain open:

- The help's remedy sentence should name `fe()`. That is a one-sentence fix with a measured
  justification, and it is the single highest-value edit this review found.
- **`bootstrap()` and `permute()` at a finer-than-`i` cluster were never coverage-tested** — only the
  analytic path was. That gap is named in `unadjudicated.md`, and the `SMP` mechanism in
  `mechanisms.md` gives specific reasons to expect trouble there.
- `impvalues` **breaks** in exactly this geometry (see `critical.md`), so the recipe above must not be
  combined with it.
