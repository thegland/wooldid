# FIX-2 — the `contreatcontrols` reference-pool degeneracy screen. FINAL specification.

**Status: STABLE, pending two author rulings.** The author's open question — *would the `byTime`
screen fire on ordinary staggered panels at a rate I would accept?* — has a measured answer, and
it reframes the question. Closes CRIT-5.

## The answer to the open question

**Fire rate is not the right question; the false-positive rate is, and it is zero.**

- On the panel families where the screen fires, it fires on **100% of draws (20/20)** — the fire
  is *structural*, a deterministic function of panel geometry, not a tail event.
- Across 48 (dose-shape x type) cells at `e(r2) == 1`, plus 16 boundary designs:
  **0 run-level false positives, 0 false negatives.** Every run whose screen fired was an
  actually-wrong run; no run whose screen stayed silent was wrong anywhere.
- The predicate's count of dangerous groups **equals wooldid's own count of dropped
  treatment-interaction coefficients, exactly, in every case checked** (1/1, 7/7, 1/1, 4/4, 2/2,
  2/2, and 0 with no caution printed). Two independently written pieces of code agreeing on the
  count is what makes the fire table evidence rather than a model.
- `basic`, `byTreated`, `byTtre`: **0 fires in 17 of 20 panel families, 400+ ordinary draws.**
  Their false-positive exposure on ordinary panels is nil.
- `byTime` is the exposed one, and the exposure has exactly one cause: does the comparison pool
  *present in that calendar period* carry dose variation? It fails in every period at or after the
  last onset whenever the dose is zero or constant for never-treated units (a very common encoding
  of treatment intensity); in any period where the reference pool is a single row; and earlier as
  `esprelength()` grows, because ES pre-horizons are treated cells, not reference rows. On
  dose-zero-for-never-treated panels the count of dangerous periods is exactly
  `T - max(onset) + 1`, verified on four families.

**Decisive precedent:** the *shipped* `byCohort` screen already exits `r(459)` unconditionally on a
perfectly ordinary panel with dose defined everywhere, once `i` is unit-level and
`esprelength(4)` is requested — verified against the real build, not a replication. `ado:1383` has
**no `safety()` gate at all**. An unconditional `r(459)` that fires on ordinary data is therefore
already shipped precedent, not a new imposition.

**Ship FIX-2 as a hard `r(459)` refusal for all four types, `byTime` included** — with the four
changes below. R1 is a ship gate; without it a silently sign-flipped answer stays in place.

---

## R1 (CRITICAL) — the screens MUST fire under `contreatwithin`. Reverse the proposal.

The proposal recommends the new screens *not* fire under `contreatwithin`. That is backwards.
Measured, on a design with no never-treated units at all and dose zero pre-onset (truth
1.411111111111):

| spec | reported `ame` | cautions |
|---|---|---|
| `basic contreatwithin` | **-0.238888889** | **NONE. Zero cautions.** `se` 1.6e-15, `t` -1.5e14, `p` 0.000 |
| `byTreated contreatwithin` | -0.238888889 | none |
| `byTtre contreatwithin` | -0.544444444 | none |
| `byCohort contreatwithin` | **r(459)** | the shipped screen **already fires** under `contreatwithin` |
| `untreatedZero contreatwithin` | 1.411111111 ✓ | — |

A sign-flipped point estimate printed as overwhelmingly significant with no caution of any kind is
the worst user-visible state found anywhere in this study. `contreatwithin` does **not** rescue the
misspecification. And the carve-out would make the new screens inconsistent with their own
template, since `byCohort`'s screen is unconditional and fires here today.

The `ado:11918` steer the proposal cites addresses a *different* failure (a dropped treated
**dummy** when the dose is 0 in the untreated period). Correct reconciliation: fire under
`contreatwithin`, with a `contreatwithin`-aware remedy line pointing at
`contreatcontrols(untreatedZero)` — which returns the exact truth on precisely this design.

## R2 (should-fix) — `sd == 0` misses a band; the near-tier's TEXT is wrong regardless

Jitter sweep on the never-treated reference dose, truth 1.330555555556, `byTime`:

| ref-pool `sd` (relative) | hard predicate | near-tier | reported `ame` |
|---|---|---|---|
| 0 | FIRES | – | 0.752777777778 (wrong) |
| 3e-13 … 3e-7 | silent | fires | **0.752777777778 — the exactly-collinear wrong answer** |
| 2.9e-6 / 1.1e-5 | silent | fires | 0.752777777862 / 0.752777777997 |
| **3.2e-5** | silent | fires | 1.330555555509 (correct to 5e-11) |
| 3e-3 and above | silent | silent | exact ✓ |

Two consequences:

- **Text fix — do it regardless, it is free.** The `1e-3` near-tier message must stop calling the
  estimate *"computationally fragile: consider increasing `regtol()`"*. In the miss band the
  answer is **deterministic, not noisy, and unmoved by `regtol()`** — the flip sits three to four
  orders of magnitude above `regtol(9)`'s 1e-9, so tying the hard leg to `regtol()` would not
  close the gap either. The message must say the estimate **may be exactly the collinear answer**.
- **Threshold, if the author wants one.** The hard leg would need `sd <= tau * sd(dose|sample)`
  with `tau` around 1e-5. **Bracketed on only one design — bracket it on more before fixing the
  constant.** This gap is inherited from the shipped `bycoh` screen, so it is pre-existing; FIX-2
  is the moment to close it, and the `bycoh` call site should move with it.

## R3 — screen every requested type, including compound spellings

`contreatcontrols(basic byTime)` returns the *same* wrong number as plain `byTime` (0.752778 vs
truth 1.330556) while plain `basic` on the same data is exact. A screen that only looked at "the"
type would pass it through.

## R4 — the `byTime` message must not say "your dose has no variation"

`byTime` fires on panels where `basic` is exactly right, so the user's first reaction will be
"my data are fine" — and on those fires the dose usually *does* vary, just not inside the period's
own reference pool. The message must name the offending calendar periods and the two remedies
measured to work: (a) keep in-sample the units treated *after* the window — adding one
out-of-window-treated cohort took `byTime` from wrong (0.752778) to exact (1.330556) and silenced
the predicate 20/20; (b) drop `byTime` from the spec.

---

## The finding that changes the record (and one CRIT-5 correction)

**The fire is a property of the data; the damage is a property of the labels.** On family-D shape
with `byTime`, `grp2 = 8 - grp` is a pure relabel — identical partition, rows and onsets:

| `i` passed | `esprelength` | `main:ame` | `pretreatment:ame` (truth 0) |
|---|---|---|---|
| `grp` | 0 | 0.752777777778 | — |
| `grp` | 4 | 0.752777777778 | -0.000000000000 ✓ |
| `grp2` | 0 | **0.508333333333** | — |
| `grp2` | 4 | **-0.155555555556** | **-0.575000000000** |

1. **Label dependence REPRODUCES.** `critical.md` records the probe's label-dependence claim as
   "NOT reproduced by the coordinator … should not be repeated without its own repro". This is
   that repro. The coordinator's null was presumably on a design where the relabel did not change
   which cell sorted last; it is not a refutation. **Correct CRIT-5.**
2. **"Every pre-treatment row is bit-identical between the correct and the broken run" is NOT
   universal** — CRIT-5 fact 1 and the blast report's headline #4 both state it. Under `grp2` +
   `esprelength(4)` the pre-treatment `ame` is -0.575 against a truth of 0. A pre-trend check
   *can* catch this, on some labellings — which is worse than never, because the pre-trend
   evidence is itself label-dependent.
3. **Therefore a screen cannot be narrowed to "only the deficiencies that do damage."** Whether a
   dangerous period damages `main:ame`, the pre rows, or nothing is a function of Stata's factor
   sort order, not of the data. "Per-period false positive" is not a stable category; the
   run-level statement is the only defensible one.

`esprelength()` amplifies the fire count and the extra deficiencies are real — wooldid's own
dropped-coefficient counter agrees at every step (2,3,4,5,6 for `epre` 0..4). On ordinary dose the
count stays 0 at every window width, so `esprelength()` creates no new fires on data that was fine.

## Author rulings
1. Accept a hard `r(459)` for `byTime` given the above, or gate `byTime` alone behind `safety()`.
   Recommendation: **accept** — `byCohort`'s unconditional screen is the shipped template, and
   every measured fire was a wrong run.
2. R2's `tau`, if a threshold is wanted, after bracketing on further designs.
