<!-- Provenance: synthesized 2026-08-25 from three independent locked-prediction
theory reviews (score-structure, finite-sample, design-based) plus a grid-completeness
audit, all conducted WITHOUT running any estimation. The machine-readable grid and the
locked consensus bands live in test/suite/mc/advb_grid.json, which analyze_advb.py reads.
Committed before any AB rep fired. -->

# PREREG — Adversarial line-test SIZE BOUNDARY campaign (AB), 2026-08-25

**Status: LOCKED. Committed before any AB rep fires.** Three theorists produced LOCKED predictions and a methodologist audited the candidate grid; this file is the final synthesis. Predictions here are the scoring standard. Grid file: `test/suite/mc/advb_grid.json`. Driver: `test/suite/mc/mc_advb.do`. Runner: `test/suite/mc/run_advb.sh`. Analyzer: `test/suite/mc/analyze_advb.py`.

**Author's question (verbatim):** *"let's dig into whether there are any cases where adversarial works really well. Let's do some tests using coarse cohorts and units-within-cohorts level clustering, and using robust SEs, to see if we can find cases where adversarial line tests are correctly sized. If we can't find many cases where this works, we may wish to just abandon this for the honestdid style functionality."*

---

## 0. Facts established by reading before locking (corrections to the brief and to the audit)

These are verified against `wooldidinstall/wooldid.ado` and the committed record, not assumed. Each one moves a number in this file.

**F-1. The brief's "A5 → 48.4%" is wrong.** `test/suite/mc/MC_REPORT_2026-08-24.md:37` and `wooldid_old_bug_report.md:696` both give **A5 (G=10, singleton, `cluster(i)`) = 12/504 = 2.38% at .05, 4.17% at .10, chi2 = 77.5\*** — *under*sized and grossly non-uniform. 48.4% is **D2** (G=40, 8 units/cohort, `cluster(i)`); 65.7% is A5's *jointtests* number. All three theorists and the auditor caught this independently. Consequence: A5 is the campaign's proof that a structurally invalid V can be dragged *through* a PASS band by a lenient `F(q, G−1)` reference, which is why the uniformity chi2 is **co-primary** below and why AB5 is report-only.

**F-2. The event-study vector is 8 coefficients, q = 6 — not 9.** I overrule the audit's 9-coefficient claim. Trace: `wooldid.ado:4287-4288` bumps `esprelength` 4 → 5 (since `abs(esrelativeto)=1 <= 4`); with onset 6 and T=12, `min_t=5`, `max_t=7`, so `espre = min(5,5) = 5`, `espl = min(7,4) = 4`; `horizonstarter = 1 + min_t − espre = 1` (`:666`), `horizonend = 1 + 5 + 4 − 1 = 9` (`:1367`); `_t = relyr + min_t + 1` and the reference slot `_t == refperiod` is remapped to 0 (`:643`), so it is **not** a reported horizon. Reported set = relyr **{−5, −4, −3, −2, 0, 1, 2, 3}**, 8 coefficients, level+slope projected out → **q = 6**. Independently confirmed: `test/suite/mc/mc_amese_att.do:25` uses exactly `ES_att:pre5 pre4 pre3 pre2 post0 post1 post2 post3` on this identical geometry. Consequences: the b/SE dump is **16** columns (not 18); the shard-1 `V_es` lower triangle is **36** numbers (not 45); the post window ends at **+3**, not +4; and `altmethod_adversarial` has `nq = 4` post horizons, so its **q = 4** (T1 said 5, T2 and T3 said 4 — T2/T3 are right).

**F-3. `AB12` as originally specified dies on every rep.** The SI-A1 screen (`wooldid.ado:1123-1155`) raises `error 198` when `jointtests` meets a cluster variable at or above the level of `i`, unless `safetyinput == off`. Every AB config passes `jointtests`.

**F-4. `safety(off)` is inert on this DGP except for F-3.** I overrule T1's objection to using it. Every other `safetyinput`-gated screen is non-binding here: weights (`:4379` — none supplied), subgroups (`:741`, `:765` — none), lattice (`:4423`), `cfxlattice` (`:4448`), ttre-commensurability (`:4473` — integer `t` and integer onsets), F13 within-`i` ttre variation (`:341` — `ttre` is constant within `i` by construction). The two `safety`-gated *estimability* checks are a caution only (`:550`) and a pre/post-sufficiency assert that passes (`:630`); the actual support and reference-period drops at `:545` and `:583` run **unconditionally**. `safety(off)` therefore changes nothing but the SI-A1 error → caution, and it **preserves the `joint_horiz` secondary**, which is what makes AB12 comparable to D2's 36.0%. Dropping `jointtests` instead would forfeit that anchor and require a driver change. Guard: the per-rep `e(N)` and b/SE columns detect any sample change.

**F-5. `wooldid_amescreen` does not fire on any AB config.** `:5975-5996` fires only where a treated `icell × tcell` has `nunique(clusvar) == 1`. Under `cluster(unitid)` at U≥2, and under `cluster(cl2)` at U=6, every treated cell spans ≥2 clusters → `fire = 0` → `e(V_es_ame)` posts (`:3647`) and AB9/AB10 are measurable.

**F-6. `wooldid_advx` (`:6784-6866`) refers `W/q` to `Ftail(q, df_r, W/q)`,** with `df_r = e(df_r)` — cluster count −1 under `cluster()`, ≈ N−K under `robust`. Size is therefore jointly determined by whether V is right *and* whether the reference is right, and the two can cancel (F-1). This is the campaign's central threat to validity.

**F-7. `altmethod_adversarial` survives here.** The B16 guard (`:690`) disables it only when `espre <= 1`; `espre = 5`. It is a *modifier* of `adversarial` (`:4583`) and replaces the primary p, so AB11 costs a full config.

---

## 1. FINAL GRID

19 configs. `T = 12` and window `esprelength(4) espostlength(4) esfixedbaseperiod` everywhere; `att` everywhere; `adversarial jointtests` everywhere. `unitid = (i−1)·U + unit` (globally unique). `cl2 = (i−1)·2 + cond(unit <= ceil(U/2), 1, 2)` (two sub-clusters per cohort). *m* = number of clusters spanning each treated cohort-by-period cell.

| config | G | U | N | VCE | *m* | clusters | errors | mode | extra | M | role |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **AB1** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | iid | binary | — | 1000 | **CERT.** Central question: does clustering finer than the cell kill containment? |
| **AB2** | 30 | 6 | 2160 | `robust` | 6 | — | iid | binary | — | 1000 | **CERT.** Pure-binary isolation of A4's 5.75% — the only measured PASS |
| **AB2a** | 30 | 2 | 720 | `robust` | 2 | — | iid | binary | — | 1000 | **CERT.** Robust ladder low end; deconfounds leverage from inversion noise |
| **AB3** | 30 | 2 | 720 | `cluster(unitid)` | 2 | 60 | iid | binary | — | 1000 | **CERT.** Cell-size floor under fine clustering |
| **AB4** | 30 | 20 | 7200 | `cluster(unitid)` | 20 | 600 | iid | binary | — | 1000 | **CERT. FLAGSHIP.** Most favorable regime constructible |
| **AB5** | 5 | 4 | 240 | `cluster(unitid)` | 4 | 20 | iid | binary | — | 504 | report-only. Few clusters; the grid's largest theorist disagreement |
| **AB6** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | AR(1) ρ=.7 | binary | — | 1000 | **CERT.** Dependence the VCE *handles*, at the cross-horizon level the statistic lives on |
| **AB7a** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | cohortshock σ=**0.2** | binary | — | 504 | report-only. VCE violated, small dose |
| **AB7b** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | cohortshock σ=**0.5** | binary | — | 504 | report-only. VCE violated, larger dose |
| **AB8** | 30 | 6 | 2160 | `robust` | 6 | — | AR(1) ρ=.7 | binary | — | 504 | report-only. Sign test: is robust conservative or anti-conservative here? |
| **AB9** | 30 | 20 | 7200 | `robust` | 20 | — | iid | contreat | — | 504 | report-only. **The ame ruling**: do big cells rescue ame adversarial? |
| **AB10** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | iid | contreat | — | 504 | report-only. Does fine clustering help ame at all? |
| **AB11** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | iid | binary | `altmethod_adversarial` | 1000 | **CERT.** Does the shipped alternative behave in a good regime? |
| **AB12** | 30 | 6 | 2160 | `cluster(i)` | **1** | 30 | iid | binary | **`safety(off)`** | 504 | report-only. **NEGATIVE CONTROL** — calibrates the new driver against D2/D4 |
| **AB13** | 30 | 6 | 2160 | `cluster(cl2)` | **2** | 60 | iid | binary | — | 1000 | **CERT.** *m*-ladder middle rung: cell size 6 at *m*=2 |
| **AB14a** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | iid | binary | `pwr(0.05)` | 504 | report-only. **POWER** curve point |
| **AB14b** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | iid | binary | `pwr(0.20)` | 504 | report-only. **POWER** — decision trigger |
| **AB15** | 30 | 6 | 2160 | `cluster(unitid)` | 6 | 180 | iid | binary | `stag(1)` (onsets {6,7}) | 1000 | **CERT.** Staggered adoption — external validity |
| **AB16** | 10 | 6 | 720 | `cluster(unitid)` | 6 | 60 | iid | binary | — | 1000 | **CERT.** Few cohorts, several units each — the applied-user geometry |

**Certification arms (10):** AB1, AB2, AB2a, AB3, AB4, AB6, AB11, AB13, AB15, AB16. Only these can be CERTIFIED and only these feed the gate arithmetic in §4. Everything else is report-only by design.

### The three orthogonal ladders this grid buys

1. ***m*-ladder at fixed everything else** — AB12 (*m*=1) → AB13 (*m*=2) → AB1 (*m*=6), all G=30, U=6, N=2160, K identical, **and run on bit-identical datasets via a shared base seed**. This is the sharpest experiment in the campaign and it is exactly the predicate a runtime gate would be written on.
2. **Cell size at fixed *m* and fixed cluster count** — AB3 (*m*=2, cell 2, 60 clusters) vs AB13 (*m*=2, cell 6, 60 clusters). Separates "cell size" from "clusters per cell" with cluster count held constant.
3. **VCE channel at fixed geometry** — AB3 vs AB2a (U=2) and AB1 vs AB2 (U=6), each pair on identical datasets. Under `robust` the cluster-count/inversion-noise channel is structurally absent, so a U-gradient that survives there is leverage, and one that vanishes was inversion noise.

### Overrules (one line each)

- **vs the audit — 9 ES coefficients:** wrong; it is 8, q=6 (F-2). Changes the dump width, the V-triangle size, the post window, and altmethod's q.
- **vs the audit — add AB2c (robust, U=20) and AB4b (G=10, U=18):** dropped. The *m*-ladder plus AB3↔AB13 achieves the same deconfounding at a third of the cost, AB9's att side already covers robust+big cells, and no theory places the boundary at the big-cell end. AB4b's G=10 confound (K and treated-cohort count move too) is avoided entirely by AB16.
- **vs T1 — run AB12 without `jointtests` rather than `safety(off)`:** overruled. F-4 verifies `safety(off)` is inert here, it needs no driver change, and it preserves the `joint_horiz` secondary that makes AB12 comparable to D2.
- **vs the audit — cut AB10 and AB11:** kept. AB10 resolves a named 12-point disagreement (T2 predicts 44%, T1/T3 predict 32%) and is the direct answer to the author's "units-within-cohorts clustering" for the ame; AB11 decides whether the shipped alternative lives or dies alongside the primary.
- **vs the audit — run AB7 at σ ∈ {0.2, 0.5}:** adopted, and the σ=1.0 dose is dropped as a foregone conclusion. All three theorists' headline AB7 points assumed σ=1 and are recorded as out-of-grid; their own stated sensitivity values at σ≈0.25 and σ=0.5 are what gets scored.
- **vs the task brief — "contreat at M=500 (8×63)":** 8×63 = **504**, and `analyze_advb.py:67` requires `n == m` exactly. M is locked at **504**, matching the A4/A5 house precedent.
- **vs T2 — add G=10/U=6:** adopted as AB16. It is the only LOCKED theorist prediction sitting in the gap between AB13 and AB1, and it is the geometry applied users actually have.
- **vs all three theorists — no power arm in their grids:** T1 and T3 both named one as decision-changing; adopted as AB14a/AB14b. A correctly sized test with no power is also an abandon.

---

## 2. LOCKED PREDICTION TABLE

**LOCKED** = the theorist predicted this exact config. **INFERRED** = I derived the cell from that theorist's stated framework because they did not predict this config; inferred cells are marked `~` and **do not count as theorist hits or misses in the scorecard**. Only the CONSENSUS band goes into `advb_grid.json` (it is what `score_prediction` reads).

### Config-name mapping (three different proposals were all called "AB13")

| theorist's label | their spec | this grid |
|---|---|---|
| T1 "AB13" | G=30, U=6, cluster on `ceil(unit/3)` within cohort, *m*=2 | **AB13** |
| T2 "AB13" | G=10, U=6, `cluster(unitid)` | **AB16** |
| T2 "AB5b" | G=5, U=20, `cluster(unitid)` | not run — see FOLLOWUP in §4 |
| T3 "AB13" | cohort-period shock under `cluster(i)` | not run — see FOLLOWUP in §4 |
| T3 "AB9-att" | free att secondary of the AB9 run | **AB9:adversarial** (secondary) |
| T3 "AB14" | power arm on AB4 geometry | **AB14a/AB14b** (AB1 geometry — cheaper, and `cluster(unitid)` is the author's named regime) |
| T3 "AB15" | staggered onsets 5/7/9 | **AB15** (onsets {6,7} — preserves the identical 8-horizon reported set, so it is a single-factor twin of AB1) |
| audit "AB13" | staggered onsets {6,7} | **AB15** |

### Predictions (rejection % at nominal .05)

| config | T1 | T2 | T3 | **CONSENSUS** | band | verdict |
|---|---|---|---|---|---|---|
| **AB1** | 5.8 | 5.4 | 5.5 | **5.6** | [4.0, 8.0] | correctly sized |
| **AB2** | 5.6 | 5.5 | 6.0 | **5.7** | [4.0, 8.0] | correctly sized |
| **AB2a** | ~6.0 | ~5.3 | ~5.5 | **5.8** | [4.0, 8.5] | correctly sized (all INFERRED) |
| **AB3** | 7.0 | 7.5 | **10.5** | **8.0** | [4.5, 14.0] | mildly oversized — **D-6** |
| **AB4** | 5.4 | 5.1 | 5.2 | **5.2** | [3.8, 7.0] | correctly sized (unanimous) |
| **AB5** | **9.0** | **19.5** | **22.0** | **17** | [3.0, 35.0] | uncertain — **D-1**, the grid's widest split |
| **AB6** | 7.0 | **5.8** | 7.0 | **6.6** | [4.2, 10.5] | mildly oversized — **D-5** |
| **AB7a** (σ=.2) | ~11 | ~11 | ~13 | **12** | [7, 22] | oversized (interpolated from their σ-curves) |
| **AB7b** (σ=.5) | 51 | 46 | ~47 | **48** | [30, 70] | severely oversized |
| **AB8** | 1.0 | 0.5 | 0.4 | **0.6** | [0.0, 3.0] | **undersized** (unanimous, contrarian) |
| **AB9** ame | 8.5 | 8.0 | 9.0 | **8.5** | [5.5, 13.0] | mildly oversized (unanimous) |
| **AB9** att (sec) | ~5.3 | ~5.2 | 5.2 | **5.3** | [3.8, 7.0] | correctly sized |
| **AB10** ame | 32 | **44** | 32 | **34** | [22, 50] | severely oversized — **D-4** |
| **AB10** att (sec) | ~5.8 | ~5.6 | ~5.8 | **5.8** | [4.0, 9.0] | correctly sized (all INFERRED) |
| **AB11** | 6.0 | 5.2 | 6.0 | **5.7** | [3.8, 8.5] | ≈ AB1; alt buys nothing in a good regime |
| **AB12** | 46 | **58** | 45 | **49** | [35, 65] | severely oversized — **D-3** |
| **AB13** (*m*=2, cell 6) | **14.5** | ~7.5 | ~10.5 | **11** | [5.5, 20.0] | severely oversized — **D-2, the decisive falsifier** |
| **AB14a** (pwr .05) | ~ | ~ | ~ | **75** | [45, 95] | power, all INFERRED |
| **AB14b** (pwr .20) | ~ | ~ | ~ | **99.5** | [95, 100] | power, all INFERRED |
| **AB15** (staggered) | ~5.6 | ~5.4 | **6.2** | **6.0** | [4.0, 9.5] | ≈ AB1 or mild degradation — **D-7** |
| **AB16** (G=10, U=6) | ~7.0 | **6.8** | ~7.5 | **6.8** | [4.5, 10.0] | mildly oversized |

### Material disagreements — the campaign's highest-value output

These are not papered over. Each names the exact quantity that resolves it.

- **D-1 (AB5, 9 vs 19.5 vs 22).** T1 leans low, anchoring on A5's *measured* 2.38% and arguing that at 20 clusters the `F(6,19)` reference plus the 1e-12 whitening tolerance dropping near-null eigenvalues over-corrects the inversion inflation. T2/T3 predict a blowout from ν ≈ 13–15 against q=6. **Turns on:** whether the reference-distribution leniency and rank collapse can flip a structurally noisy V to *conservative*. **Resolved by:** AB5's rate *plus* its chi2 *plus* the shard-1 `V_es` rank. A low rate with a catastrophic chi2 confirms T1's mechanism and is scored as a failure, not a pass.
- **D-2 (AB13 vs AB3 — the decisive one).** T1 predicts AB13 (14.5%) is *worse* than AB3 (7.0%) despite a 3× larger cell, because CR1's `c = (N−1)/(N−K)` is 1.42 at U=2 (K/N = 0.29) and over-corrects, but only ≈1.11 at U=6 (K/N = 0.08) and cannot offset `ρ_cell = 1 − 1/m = 0.5`. T2 and T3's frameworks both give AB13 ≈ AB3 (identical ν ≈ 44, identical 60 clusters). **Turns on:** whether size is governed by *clusters per cell* (T1) or by *cluster count and cell size* (T2/T3). **This is the single result that determines what predicate a runtime gate is written on.** T1 named it his own named falsifier: "If AB13 comes back at AB1's level, my score-structure account is wrong."
- **D-3 (AB12, 58 vs 46/45).** T2 argues A1's singleton design had K/N = 0.56, so CR1's c = 1.98 *accidentally* restored the halved scale, and that no such rescue exists at K/N = 0.08 — so AB12 should be much worse than A1's 20.5%. T1/T3 anchor directly on D2/D4 (48.4/44.8) and argue ρ is U-invariant under containment. **Turns on:** whether the CR1 correction interacts with K/N in the containment regime. Either way AB12 fails; the *magnitude* adjudicates.
- **D-4 (AB10, 44 vs 32/32).** T2 predicts fine clustering makes the *ame* **worse** than A4's robust 31.55%: under `robust` the ame cross-horizon off-diagonals are structurally zero, which happens to match the truth, whereas `cluster(unitid)` makes each unit's time-invariant dose weight generate nonzero, PR-limited-df off-diagonals that are pure noise. T1/T3 predict a near-replication of A4 (CR1 1.21–1.26 ≈ HC1 1.21–1.25, one obs per unit per cell). **Turns on:** whether unit-clustered off-diagonals in `V_es_ame` are noise or signal. Directly answers the author's "units-within-cohorts clustering" for the continuous case.
- **D-5 (AB6, 5.8 vs 7.0/7.0).** T2 says AR(1) inside the cluster is exactly scale/shape invariant, so AB6 ≡ AB1 and the measured A1→A6 gap (20.5 → 23.2) was MC noise. T1/T3 say demeaning removes more error variance per observation under AR(1) than under iid (E[e²] ≈ .935σ² vs .975σ²), which CR1's uniform correction cannot know, costing ~1.5pp. **Resolved by:** the paired AB1−AB6 difference (different CRN groups, so unpaired; use the binomial MC SE of ≈0.75pp each).
- **D-6 (AB3, 10.5 vs 7.0/7.5).** Magnitude of the U=2 penalty. Combined with D-8 it splits leverage from inversion noise.
- **D-7 (AB15, 6.2 vs ≈5.5).** T3 expects modest degradation from staggering (thinner per-horizon cell counts, time-varying control pool) and flags that *no* result in the record generalizes without this arm; T1's and T2's frameworks are onset-blind (*m* and ν do not change). **Turns on:** whether pooling cells from two calendar periods into one horizon changes the score structure or only the estimand mix.
- **D-8 (the AB3 − AB2a gap: T1 ≈1pp, T2 ≈2.2pp, T3 ≈5pp).** Under `robust` there is no cluster-count channel. A gap near zero says the U-gradient is inversion noise (T2); a large gap says it is cell-size leverage (T3). AB2a exists solely for this.

---

## 3. PASS CRITERIA (house conventions, reused verbatim)

Central ~95% binomial region, normal approximation rounded outward — exactly `analyze_advb.py:58-61`, exactly the M2 convention that produced A1–A6 and A4.

**M = 1000** (8 shards × 125):
- rej@.05 PASS band = **[36, 64]** rejections = [3.6%, 6.4%]
- rej@.10 PASS band = **[81, 119]** rejections = [8.1%, 11.9%]

**M = 504** (8 shards × 63):
- rej@.05 PASS band = **[15, 35]** rejections = [2.98%, 6.94%]
- rej@.10 PASS band = **[37, 64]** rejections = [7.34%, 12.70%]
- (Sanity check against the record: A4 att measured 29 and 61 → PASS. ✓)

**Uniformity (CO-PRIMARY, not secondary).** 10 equal-probability bins on (0,1]; chi2 with df = 9; critical value **21.666** at p = .01. This is co-primary because of F-1: A5 sat inside a rate band while its p-value distribution was destroyed (chi2 = 77.5, median p = .658, 15.7% of p-values above 0.9). `analyze_advb.py:147-149` already implements this correctly — a config is certified only if `verdict == "PASS" and not chi2_flag`. Endorsed as written.

**MARGINAL corridor** (house): rate@.05 ∈ [.025, .085] and rate@.10 ∈ [.070, .135]. **MARGINAL is not CERTIFIED.**

**CERTIFIED(C)** ≡ C is one of the 10 certification arms **AND** verdict = PASS (both bands) **AND** chi2 ≤ 21.666 **AND** integrity clean (`n == M`, zero `rc != 0`, zero missing p).

**Report table additions (mandatory).** For every config the report must print the p-value CDF at α ∈ {.01, .025, .05, .10, .20}, not just the two band counts. Pass-by-cancellation shows up as a *shape*, not a point.

---

## 4. THE DECISION RULE (mechanical; applied to the results table without further judgment)

Let `rate(C)` = 100 × rej@.05 / n for C's primary adversarial quantity. Apply in order; **the first trigger that fires decides.**

**V0 — CAMPAIGN VALIDITY.** If any certification arm has integrity FAIL (`n != M`, any `rc != 0`, or any missing p) → **NO RULING**; fix and re-run that arm before applying anything below.

**V1 — NEGATIVE-CONTROL CALIBRATION.** If `rate(AB12) < 25` → **NO RULING**. The driver is not reproducing the known containment failure (D2 = 48.4%, D4 = 44.8%), so no PASS anywhere in the grid is trustworthy.

**A1 — FLAGSHIP.** If **NOT CERTIFIED(AB4)** → **ABANDON.** AB4 is the most favorable regime constructible (20 clusters per cell, 600 clusters, CR1 nearly inert, VCE exactly matching the DGP). All three theorists staked the feature on it. If it cannot be sized there, no runtime gate can rescue it.

**A2 — POWER.** If AB14b rejection < **90%** → **ABANDON.** A correctly sized but powerless test is an abandon; honestdid-style bounds, which report a magnitude rather than a reject/accept, are then strictly better.

**A3 — ORDINARY GEOMETRY.** If `rate(AB1) > 8.0` → **ABANDON.** Six units per treated cell with unit-level clustering is the author's named regime and the median applied panel. A gate that refuses it refuses most real datasets, which is a worse user experience than not offering the test (T1's abandon-on-practicality, adopted as a hard trigger).

**A4 — REALISTIC ERRORS.** If `rate(AB6) > 8.0` → **ABANDON.** AB6 clusters at exactly the level of the dependence. If the test cannot survive ordinary within-unit serial correlation *with the correct VCE*, "use fine clustering" is not deployable advice for panel data.

**A5 — STAGGERING.** If `rate(AB15) > 8.0` → **ABANDON.** wooldid is a staggered-adoption estimator. A test valid only under a single common onset is valid only outside the designs the package exists to serve.

**K — Otherwise: KEEP, GATED.** Gate parameters set mechanically from the results:

- `m* = min{ m(C) : C a certification arm with a cluster() VCE and CERTIFIED(C) }` — candidates AB13/AB3 (*m*=2), AB16/AB1/AB6/AB11/AB15 (*m*=6), AB4 (*m*=20).
- `n* = min{ clusters(C) : C a certification arm with a cluster() VCE and CERTIFIED(C) }`.
- `cell* = min{ U(C) : C a certification arm with robust VCE and CERTIFIED(C) }` — candidates AB2a (2), AB2 (6).
- **Override:** if `rate(AB1) ∈ (6.4, 8.0]` (below the abandon trigger but outside the band), set `m* = 6` regardless of what the formula returns, and the release note must state that ~6-unit cells are mildly oversized.

**Ship form.** Mirror the SI-A1 screen (`wooldid.ado:1123-1155`), reusing the `wooldid_amescreen mode(detect)` predicate (`:5975-5996`), which already computes clusters-per-cell:
1. **Hard error** when any treated cohort-by-period cell lies inside a single cluster (*m* = 1) — for `adversarial`, `altmethod_adversarial`, **and the honestdid `e(V_es_*)` handoff**. The shared object at risk is `V_es`, not the test built on it; gating one and not the other relocates the failure rather than fixing it.
2. **Caution** when min-over-cells *m* < `m*`, or cluster count < `n*`, or (robust) min cell size < `cell*`.
3. **Silent** above.
4. **Documented caveats** per DOC1/DOC2 below.

**AME — INDEPENDENT SUB-RULING** (does not affect the binary att ruling; deciding on the ame would be reading the wrong evidence). Keep the ame-side adversarial iff `AB9:ame_adversarial` has rej@.05 ∈ **[15, 35]** and chi2 ≤ 21.666. Otherwise **WITHDRAW** ame-side `adversarial` and `joint_Choriz` and route ame joint inference to RI / cluster bootstrap / IM p-values. On the consensus prediction (8.5%) this fires as WITHDRAW.

**ALT — ALTMETHOD SUB-RULING.** Keep `altmethod_adversarial` iff CERTIFIED(AB11). If CERTIFIED(AB11) and NOT CERTIFIED(AB1), **promote** altmethod to the documented default for the line test. If CERTIFIED(AB1) and NOT CERTIFIED(AB11), **withdraw** altmethod. Note when reading: AB11 tests q=4, AB1 tests q=6 (F-2) — they are not the same statistic and a small gap is expected.

**DOC1 — COHORT-SHOCK RIDER.** If `rate(AB7b) > 20`, the sthlp design-conditions section must state that the line test is invalid under cohort-by-period common shocks *even when clustering finer than the cell*, that this is undetectable at runtime, and that the containment escape route trades one failure for another. If additionally `rate(AB7a) > 15` (a shock contributing only ~3.8% of error variance), the ruling is downgraded from KEEP-GATED to **KEEP-GATED-WITH-PROMINENT-WARNING**: the escape from containment then requires assuming away exactly the dependence that motivates cohort-level clustering in the first place.

**DOC2 — AR(1)/ROBUST RIDER.** If `rate(AB8) < 3.6`, the documentation must state that robust SEs make the line test **conservative, not safe**, under within-unit serial correlation, and that AB14's power result does not transport to that case. Undersizing is not a pass and must never be reported as one.

**FOLLOWUP — BOUNDARY GAP.** If NOT CERTIFIED(AB13) and NOT CERTIFIED(AB3) and CERTIFIED(AB1), the *m*/cluster boundary lies in (2, 6] and a second, cheap campaign must map it (*m* ∈ {3, 4} at U = 6 and 12) before `m*` is documented as 6; until then the gate ships at `m* = 6`. Two further follow-ups are pre-registered as *not* run here and are triggered only by their conditions: T2's G=5/U=20 arm (run iff AB5 fails and AB16 certifies, to separate few-cohorts from few-clusters), and T3's cohort-shock-under-`cluster(i)` arm (run iff DOC1 fires, to establish that no analytic VCE recovers an absorbed shock).

---

## 5. DGP SPECIFICATION

Exact construction, per rep, matching `mc_advb.do`. `G`, `U`, `T`, `VCE`, `ERR`, `MODE`, `SCS`, `STAG`, `PWR`, `EXTRA` are the config's arguments.

```stata
local ONSET = ceil(`T'/2)                        // = 6 at T=12
set seed `BASE' + `SHARD'*10000 + `rep'
clear
qui set obs `G'
gen i = _n
gen ttre = cond(mod(_n,2)==0, `ONSET' + `STAG'*mod(_n/2,2), .)   // STAG=0 -> all onsets 6
                                                                  // STAG=1 -> onsets alternate 7,6,7,6,...
gen a    = rnormal()                             // cohort intercept: absorbed by the cohort FE (inert)
qui expand `U'
bysort i: gen unit = _n
gen unitid = (i-1)*`U' + unit                    // GLOBALLY UNIQUE
gen cl2    = (i-1)*2 + cond(unit <= ceil(`U'/2), 1, 2)   // 2 sub-clusters per cohort
gen ure  = cond("`ERR'"=="unitre", rnormal(), 0) // 0 in every AB arm (see below)
gen dose = exp(0.8*runiform() + 0.2*mod(i,5))    // unit-level, varies WITHIN each cell
qui expand `T'
bysort i unitid: gen t = _n

if "`ERR'" == "ar1" {
  bysort unitid (t): gen e = rnormal()
  bysort unitid (t): replace e = .7*e[_n-1] + sqrt(1-.49)*rnormal() if _n > 1
}
else gen e = rnormal()

if "`ERR'" == "cohortshock" {                    // correlated ACROSS units within a cohort-period
  tempvar cs
  bysort i t: gen `cs' = rnormal() if _n == 1
  bysort i t (`cs'): replace `cs' = `cs'[1]
  qui replace e = e + `SCS'*`cs'
}

gen y = a + ure + .05*t + `PWR'*cond(ttre < ., (t - ttre)^2, 0) + e

wooldid y i t ttre, att `VCE' `CONTOPT' esprelength(4) espostlength(4) ///
        esfixedbaseperiod adversarial jointtests `EXTRA'
// CONTOPT (mode==contreat only) = contreat(dose) contreatcontrols(basic)   [as in A4]
```

**Treated/never-treated split.** Cohorts with even `i` are treated; odd `i` are never-treated and form the comparison pool in every period. G=30 → 15 treated / 15 control; G=10 → 5/5; G=5 → **2 treated / 3 control** (AB5 confounds few clusters with few treated cohorts; this is why it is report-only and why AB16 exists).

**Onset timing.** `STAG=0`: all treated cohorts adopt at t=6, so relyr = t−6 ∈ [−5, 6]. `STAG=1`: onsets alternate 7 and 6 (8 cohorts at 7, 7 at 6 for G=30), so relyr ∈ [−6, 6] pooled. **Verified:** both give `min_t=5,6`, `espre=5`, `espl=4`, and the identical reported horizon set relyr {−5,−4,−3,−2,0,1,2,3} (8 coefficients, q=6). Every horizon is supported by both onsets (relyr −5 at t=1/t=2; relyr +3 at t=9/t=10), the reference relyr −1 is observed for both (t=5/t=6), and never-treated units supply control support in all 12 periods, so the `controlcheck` at `:545` is non-binding. **AB15 is therefore a strict single-factor twin of AB1.**

**Why the null is exactly true, arm by arm.**
- All arms except AB14a/AB14b: `y` contains **no treatment effect and no dose effect at all** (`PWR = 0`). The true event-study profile is flat zero, which lies on a line — so the composite null "all relative-time coefficients on some line" holds exactly, for both the att and the ame quantities.
- `a` is cohort-level → absorbed by the cohort FE, inert.
- `.05*t` is period-level → absorbed by the period FE, inert.
- `ure` is identically 0 in every AB arm. An unabsorbed unit random intercept **differences out** of contrasts under `esfixedbaseperiod` while inflating residuals, making robust *conservative* for this statistic rather than invalid — so it is not a usable violation arm and is deliberately unused. (The `cond()` still burns one draw per observation in every arm, so AB streams are not comparable to A-series streams. Harmless for size.)
- `dose` is generated after `expand U` so it varies across units within every cohort-period cell → the ame is identified and the SI-A4 void screen cannot fire.
- **AB14a/AB14b are the deliberate exception: the null is FALSE by construction.** `PWR ≠ 0` plants τ(h) = PWR·h² on treated units at every relative time h. Over the 8 reported horizons the deviations from the best-fit line are PWR × {8.467, 1.475, −3.517, −6.508, −6.492, −3.483, 1.525, 8.533}, with Σr² = 258.0·PWR². Design targets, using SE(b_h) ≈ 0.21 at the AB1 geometry: PWR = 0.05 → λ ≈ 15 → power ≈ 0.78; PWR = 0.20 → λ ≈ 235 → power ≈ 1.00. The SE figure is approximate, which is why the **abandon trigger A2 uses only PWR = 0.20**: even a 3× SE mis-estimate leaves λ ≈ 26 and power > 0.97. The realized λ̂ is recomputed at analysis time from the b/SE dump and the shard-1 `V_es`, so a miss is diagnosable rather than arguable. **Power is interpretable only if AB1 clears A3**; the rule's ordering already makes an ABANDON verdict insensitive to this.

**Error structures.** `iid`: e ~ N(0,1) i.i.d. `ar1`: within-unit AR(1) with ρ = 0.7 (the A6 precedent), stationary-scaled. `cohortshock`: e + SCS·v_{i,t} with v ~ N(0,1) i.i.d. across cohort-periods — SCS = 0.2 (3.8% of error variance) and SCS = 0.5 (20%). Note the shock inside a *treated* cell is collinear with that cell's dummy (`wooldid.ado:1186` fits `i._i#i._t#i.targetflag`), so it passes wholly into the estimate and leaves zero trace in the residuals: invisible to every analytic VCE. That absorption is the point of the arm, not a defect.

### Driver, runner and grid edits required before launch

All in `test/suite/mc/` unless noted. **No `wooldid.ado` change.**

1. **`cl2`** — one `gen` line after `unit` exists (above). Consumes no randomness, so it does not perturb any seed stream.
2. **`STAG`, `SCS`, `PWR`** — three new positional args (defaults 0, 1, 0), threaded through `advb_grid.json`, the `plan()` emitter in `run_advb.sh:35-38`, and **both** `read` loops (`run_advb.sh:50` and `:64`).
3. **Per-rep b/SE dump** — 16 columns `b1..b8, s1..s8` from `e(wdidmainresults)` rows `ES_att:pre5 pre4 pre3 pre2 post0 post1 post2 post3`, using the proven `rownumb()`-guarded idiom at `mc_amese_att.do:45-58`. Use `e(wdidmainresults)`, **not** `e(b_es_att)`/`e(V_es_att)` — those post under `cap` (`:3643-3649`) and can silently fail to post, which would look like missing data rather than an error. Contreat arms additionally dump the 8 `ES_ame:` rows into `ab1..ab8, as1..as8`; all arms write `.` where a row is absent.
4. **`n` and `dfr` columns** — `cap`-guarded `e(N)` and `e(df_r)`. `e(N)` must equal `G*U*T` in every arm; any deviation is a silent-drop detector and an integrity FAIL.
5. **Shard-1 `V_es` dump** — on shard 1 only, write the 36-number lower triangle of `e(V_es_att)` (and `e(V_es_ame)` for contreat arms) to `out/<cfg>_V_s1.csv`. This is what turns "AB3 failed" into "AB3's V is understated 1.4× uniformly, consistent with leverage 1/U, and the adversarial subspace is not differentially hit."
6. **`cl2` sanity assertion** — before the first `wooldid` call in each shard, assert that `nunique(cl2)` within each treated `i × relyr` cell equals the intended *m*; abort the shard otherwise. A silently-wrong cluster variable is exactly the failure mode the driver header at `:20-22` was written to prevent.
7. **`advb_grid.json` naming (load-bearing for `analyze_advb.py:147-149`).** The primary att quantity is named exactly **`adversarial`** in the 10 certification arms and **`adversarial_ro`** in the 9 report-only arms (AB5, AB7a, AB7b, AB8, AB9, AB10, AB12, AB14a, AB14b). The ame primary is named **`ame_adversarial`**. This makes the analyzer's `k.endswith(":adversarial")` filter select exactly the certification set with **no analyzer code change**. `decision_rule.cluster_type_keys` must list the literal strings `AB1:adversarial, AB3:adversarial, AB4:adversarial, AB6:adversarial, AB11:adversarial, AB13:adversarial, AB15:adversarial, AB16:adversarial` (AB2/AB2a are robust and are excluded from that line by design).
8. **`AB12` extra = `safety(off)`**, via the existing `extra` slot (`mc_advb.do:18` already anticipates it).

### Seeds (published, not left to the implementer)

`seed = base_seed + shard*10000 + rep`, shards 1..8, reps 1..125 (M=1000) or 1..63 (M=504). Bases are spaced by 100000 > max offset 80125, so no two configs can collide.

| base_seed | configs | CRN meaning |
|---|---|---|
| **8200001** | AB1, AB2, AB10, AB11, AB12, AB13, AB14a, AB14b, AB15 | **bit-identical datasets** (G=30, U=6, iid; `cl2`/`ttre`/`PWR` consume no draws). The *m*-ladder, the cluster-vs-robust contrast, the binary-vs-contreat att contrast, the staggering contrast and the power arms are all **paired** on the same data. |
| **8300001** | AB6, AB8 | identical AR(1) datasets — cluster vs robust, paired |
| **8400001** | AB7a, AB7b | identical shock draws, scaled — the σ-gradient is paired |
| **8500001** | AB2a, AB3 | identical U=2 datasets — robust vs cluster, paired |
| **8600001** | AB4, AB9 | identical U=20 datasets |
| **8700001** | AB5 | — |
| **8800001** | AB16 | — |

**Reading rule:** differences *within* a CRN group are paired comparisons and must be reported with the paired MC standard error (much smaller than the unpaired one). Differences *across* groups are unpaired; use the binomial MC SE (≈0.75pp at 5% and M=1000, ≈1.0pp at M=504). An undocumented shared stream is what invalidates a "difference" claim after the fact — hence the table.

### Cost control (mandatory, pre-committed)

Run `pilot_median.py` at 10 reps each on the three cost archetypes — **AB1** (N=2160 binary), **AB4** (N=7200 binary), **AB9** (N=7200 contreat) — before launching. If projected wall-clock at `STATA_JOBS=8` exceeds **6 hours**, apply this demotion ladder in order until the projection clears, amending `advb_grid.json` **before any real rep fires** and recording the amendment here:

**D1** AB2a → 504 · **D2** AB11 → 504 · **D3** AB16 → 504 · **D4** AB2 → 504 · **D5** AB4 → 504 · **D6** AB15 → 504.

Never demoted, never dropped: **AB1, AB3, AB13** (the *m*-ladder and cell-size contrast), **AB6**, **AB9**. Nothing goes below M=504. Every demoted config is judged against the M=504 bands in §3, which are already stated, so a demotion is a preregistered move and not a post-hoc one. `shards × reps_per_shard == M` must hold exactly for every config (`analyze_advb.py:67` treats a mismatch as an integrity failure).

---

## 6. RISKS AND GUARDS

| # | What could make this campaign uninformative | Guard |
|---|---|---|
| **R-1** | **Pass-by-cancellation.** A structurally blind V dragged through the PASS band by a lenient `F(q, df_r)` reference — A5 is the proof this happens (F-1, F-6). A "PASS" would then certify a broken test. | Uniformity chi2 is **co-primary**, not secondary (`analyze_advb.py:147-149` endorsed as written). Plus the full p-value CDF at 5 alphas in every row, plus the b/SE dump: a config in the band whose reported SEs are 15% too small is not a pass, and only the dump can see that. |
| **R-2** | **Silent rank collapse.** The 1e-12 whitening tolerance (`:6812`) drops near-null eigenvalues, q shrinks, and a broken config looks artificially *good* — exactly A5's trap. | Shard-1 `V_es` dump gives the realized eigenvalue spectrum and rank; the b/SE dump shows a vanished horizon as a `.` immediately. AB5, the arm most at risk, is report-only and cannot enter the certified set (naming rule, edit 7). |
| **R-3** | **A "FAIL" that cannot be decomposed.** One p-value per rep can only ever say pass/fail, never why — and per R-1 both are ambiguous. | The 16-column b/SE dump plus shard-1 `V_es` yields empirical cov(b_es) across reps versus mean reported V, **by direction**, including the small-eigenvalue cross-horizon contrasts where the damage lives. Every failure becomes a variance ratio, not an adjective. |
| **R-4** | **AB12 errors on every rep** (SI-A1, F-3), destroying the negative control and with it the campaign's calibration. | `safety(off)` pinned in the grid for AB12 only (F-4 verifies it is otherwise inert); rule V1 makes a missing or non-firing negative control a NO RULING rather than a silent pass elsewhere. |
| **R-5** | **AB7 unfalsifiable** — all three theorists' points swing 70pp on an unstated σ. | σ pinned at two levels (0.2, 0.5) and both predicted from the theorists' own sensitivity curves. σ is a fixed constant, so no per-rep variance decomposition is needed. |
| **R-6** | **A driver edit silently changes the sample or the cluster variable.** Six edits are being made before launch, including a new cluster construction. | `e(N)` per rep must equal `G*U*T`; the `cl2` assertion (edit 6) aborts a shard whose clusters-per-cell is not the intended *m*; a 2-rep smoke run of **every** config with the dump inspected precedes the full launch. This is also the guard against `safety(off)` doing anything unexpected in AB12. |
| **R-7** | **Prediction scoring becomes unauditable** — three configs were all called "AB13", and many cells are my synthesis rather than a theorist's commitment. | Config-name mapping table in §2; every synthesized cell marked `~` INFERRED; **only LOCKED cells count as theorist hits or misses**; only consensus bands enter `advb_grid.json` (which is what `score_prediction` reads). |
| **R-8** | **The power arm's realized noncentrality lands far from design**, making A2 fire (or not fire) for the wrong reason. | Two magnitudes trace the curve; the abandon trigger uses only the large one, where a 3× SE mis-estimate still leaves power > 0.97; λ̂ is recomputed ex post from the dump. |
| **R-9** | **CRN seed sharing mistaken for independence**, inflating confidence in a cross-config difference. | Published CRN table; within-group differences reported paired, across-group unpaired, with the MC SE stated. |
| **R-10** | **Cost overrun** starves the flagship arms or forces unpreregistered cuts. | Mandated pilot on three cost archetypes + the D1–D6 demotion ladder above, with both band sets already stated so a demotion is not a judgment call. |
| **R-11** | **A certification arm certifies for a reason outside its intended contrast** (e.g. AB5 or AB9's att side sneaking into the certified set and setting `n*` = 20). | The certified set is restricted to the 10 designated arms by the `adversarial` / `adversarial_ro` naming rule; **MARGINAL is explicitly not CERTIFIED**; gate arithmetic (`m*`, `n*`, `cell*`) draws only from that list. |
| **R-12** | **The ruling is read off the wrong evidence** — abandoning the binary att test because the ame arm fails, or keeping it because a violation arm happened to look conservative. | The ame ruling is a separate, independently-triggered sub-ruling; DOC2 forbids scoring undersizing as a pass; AB7/AB8 are report-only and enter the decision only through DOC1/DOC2 documentation riders, never through the KEEP/ABANDON trunk. |
| **R-13** | **Everything certifies, and the result still does not generalize** — all 12 existing measurements and 18 of 19 AB configs use a single common onset, in a package that exists for staggered adoption. | AB15 is a strict single-factor twin of AB1 (verified horizon-set identity, F-2/§5) and its failure is a **hard ABANDON** at A5, not a footnote. |
| **R-14** | **KEEP is ruled, but the escape from containment is illusory** — the fix for *m*=1 is to cluster finer than the cell, and clustering finer than the cell is invalid precisely when there is cell-level dependence. | DOC1 makes AB7a/AB7b decide the *strength* of the shipped warning, and a failure at σ=0.2 (3.8% of error variance) downgrades the ruling to KEEP-GATED-WITH-PROMINENT-WARNING. If the honest release note has to say "valid only when you may assume no cohort-by-period shocks", the author should see that stated before the results arrive, not after. |

**One framing note that belongs in the record before results exist:** honestdid-style functionality is not an escape hatch from this problem. It consumes the same `e(V_es_*)`. If `V_es` is directionally understated — the containment family, or any design with cohort-by-period shocks — honestdid's sensitivity bounds are wrong in the same directions by the same factors. Whatever design condition this campaign establishes must be imposed on the honestdid handoff too, which is why the ship form in §4 gates both.