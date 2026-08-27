# Adversarial line-test size boundary — REPORT OF RECORD (AB campaign, 2026-08-26)

Preregistration: `test/suite/PREREG_ADVBOUND_2026-08-25.md` (grid, locked per-config
predictions and the mechanical keep-vs-abandon rule committed at `263c37c`, before any rep
fired; author-requested amendment 1 at `293b9dc`). Harness: `mc_advb.do`, `run_advb.sh`,
`analyze_advb.py`, `advb_grid.json`, all committed at `e478711` before results existed.
Machine-readable output: `out/advb_report.md` (out/ is gitignored; this file is the record).

**25 configs, 17,560 reps, zero nonzero return codes, zero missing p-values, every config at
its exact preregistered rep count.** Locked consensus bands contained the measured rate in
**25 of 26** scored quantities.

---

## 1. The rulings

| ruling | outcome | basis |
|---|---|---|
| **Binary `adversarial`** | **KEEP, GATED** | no abandon trigger fired; 5 configs certified |
| **`altmethod_adversarial`** | **KEEP — and it strictly dominates** | better than standard at *every* boundary; certifies once where standard fails (at the band edge) |
| **Continuous `adversarial` (ame)** | formal: **WITHDRAW**; recommended: **GATE ON CELL SIZE** | 31.55% at 6 obs/cell → 9.92% at 20 → **5.95% PASS at 50** (AB22, amendment 2) |
| Gate parameters | m\* = 6 clusters per cell, cell\* = 6 obs (robust), n\* ≥ 180 clusters (not localized) | minima over the certified set |

Decision rule applied mechanically: V0 integrity clean → V1 negative control 56.35% ≥ 25 (campaign
valid) → A1 flagship certified → A2 power 100% ≥ 90 → A3 4.80% ≤ 8.0 → A4 5.60% ≤ 8.0 →
A5 7.20% ≤ 8.0 → **KEEP, GATED**.

---

## 2. Both methods, side by side (the author's question)

Nominal 5%. PASS band [3.6, 6.4]% at M=1000 and [2.98, 6.94]% at M=504; uniformity chi2 must
also be ≤ 21.666 (co-primary — see §5). Altmethod arms are CRN-paired to their standard twins:
identical datasets, so the differences below are paired.

| regime | standard | altmethod | who wins |
|---|---|---|---|
| 6-unit cells, unit clustering (AB1/AB11) | **4.80% ✓certified** | **4.40% ✓certified** | both work |
| 2-unit cells, unit clustering (AB3/AB17) | 9.10% ✗ | **6.94% ✓certified (band edge)** | **altmethod only** |
| few clusters, 20 clusters (AB5/AB20) | 17.46% ✗ | 11.31% ✗ | altmethod better, both fail |
| cohort-period shock σ=0.5 (AB7b/AB19) | 44.84% ✗ | 36.71% ✗ | altmethod better, both fail |
| containment: cluster(i) (AB12/AB18) | 56.35% ✗ | 37.70% ✗ | altmethod better, both fail |
| robust + AR(1) (AB8/AB21) | 0.00% ✗ | 0.60% ✗ | both inert (see DOC2) |

**Altmethod is uniformly better than standard — at every single boundary measured, without
exception.** The margin is large where the failure is structural (56.35 → 37.70; 44.84 → 36.71;
17.46 → 11.31) and it is decisive exactly once: at 2-unit cells it crosses from failure into
certification, which is what the preregistered ALT-INDEPENDENT rule required to keep it as an
independently valuable method rather than a footnote.

**The knife-edge, stated plainly.** AB17's certification is the weakest possible pass: 35
rejections out of 504, and the band's upper limit is *exactly* 35. One additional rejection
would have failed it. The exact paired comparison on the 504 identical datasets is standard 47
vs altmethod 35 (discordant 31 vs 19), a **+2.38pp** paired improvement with a paired MC standard
error of **1.40pp** — about 1.7 SE, suggestive but not decisive on its own. **Recommendation:
treat "2-unit cells are supported under altmethod" as provisional and confirm at M ≥ 2000 before
documenting it.** Everything else in the altmethod column is a magnitude reduction that does not
change a verdict.

Why altmethod does better is mechanical and worth recording: it tests q = 4 dimensions (post
horizons) where the standard test's projection leaves q = 6. Fewer dimensions inverted through
a directionally-starved covariance means less amplification. That predicts exactly the observed
pattern — same answer where V is honest, uniformly smaller damage where V is not.

---

## 3. The boundary map (standard `adversarial`, binary att)

**Works (certified):**

| config | geometry | VCE | rate | chi2 |
|---|---|---|---|---|
| AB1 | 30 cohorts × 6 units | cluster(unitid) | 4.80% | 3.4 |
| AB2 | 30 × 6 | robust | 4.60% | 7.2 |
| AB4 | 30 × 20 | cluster(unitid) | 5.00% | 6.8 |
| AB6 | 30 × 6, AR(1) ρ=.7 | cluster(unitid) | 5.60% | 8.7 |
| AB11 | 30 × 6, altmethod | cluster(unitid) | 4.40% | 9.4 |

**Fails:**

| config | what breaks it | rate | chi2 |
|---|---|---|---|
| AB13 | **only 2 clusters per cell** (60 clusters, 6-unit cells) | 14.50% | 238.4\* |
| AB3 | 2-unit cells | 9.10% | 25.5\* |
| AB2a | 2-unit cells, robust | 7.40% | 14.4 (marginal) |
| AB5 | 20 clusters | 17.46% | 231.1\* |
| AB7a | cohort-period shock σ=0.2 (**3.8% of error variance**) | 12.30% | 72.7\* |
| AB7b | cohort-period shock σ=0.5 | 44.84% | 1213.1\* |
| AB12 | cohort-level clustering (containment) | 56.35% | 1674.3\* |
| AB15 | staggered onsets {6,7} | 7.20% | 24.2\* (marginal) |
| AB16 | 10 cohorts (60 clusters) | 6.70% | 15.3 (marginal) |
| AB8 | robust under AR(1) — **inert, not safe** | 0.00% | 1072.4\* |

**Power is not a concern** where size holds: planted quadratic deviations rejected at 98.81%
(PWR=.05) and 100.00% (PWR=.20).

**Free secondary — the binary joint pre-trend test (`jointtests`) tracks the adversarial result
closely** in every regime (AB1 3.80%, AB4 4.90%, AB6 5.10% where adversarial certifies; 38.29%
at AB12, 33.53% at AB7b where it fails). This closes a gap flagged during the SI-A1 discussion:
binary joint tests under fine clustering were "expected correct but unmeasured." They are now
measured, and correct — which independently validates the steer text shipped in the SI-A1 screen.

---

## 4. The continuous (ame) side — WITHDRAW, and the advice inverts

| config | geometry | VCE | ame rate | att twin |
|---|---|---|---|---|
| (A4, prior campaign) | 30 × 6 | robust | 31.55% | 5.75% ✓ |
| AB9 | 30 × 20 | robust | **9.92%** (chi2 35.6\*) | 4.96% ✓ |
| **AB22** | **30 × 50** | **robust** | **5.95% ✓ PASS** (chi2 8.2) | 5.16% ✓ |
| AB10 | 30 × 6 | cluster(unitid) | **48.41%** (chi2 1491.8\*) | 5.16% ✓ |

Two findings:

1. **Large cells do NOT merely help — they cure it.** The completed ladder under robust is
   **31.55% (6 units) → 9.92% (20 units) → 5.95% (50 units, chi2 8.2, PASS)**. Excess over
   nominal falls 26.55pp → 4.92pp → 0.95pp, i.e. by ~5.3× for each ~3× increase in cell size,
   a strikingly regular decay (roughly excess ∝ cell size^−1.4) that confirms the leverage
   mechanism quantitatively rather than merely directionally. **The ame adversarial test is
   therefore not irredeemable: it is correctly sized under robust SEs once cohort-period cells
   contain on the order of 50 observations.**

   *Status of the ruling.* The preregistered sub-ruling keyed on AB9 (20-unit cells) and
   mechanically returns WITHDRAW; AB22 is a disclosed post-hoc amendment and by its own terms
   cannot move a trigger. So the formal ruling stands as written. But the honest scientific
   reading is that the withdrawal is **conditional, not absolute**, and the recommendation to
   the author is a GATED KEEP for the ame side rather than removal: refuse/caution below a cell
   size threshold and permit it above one. Caveats: this is a single config at M=504, the pass
   is comfortable (30 of 504, band limit 35) but the threshold between 20 and 50 units is
   unmapped, and it holds for **robust only** — see finding 2.
2. **Clustering finer makes the ame WORSE, not better** — 48.41% (unit-clustered) vs 31.55%
   (robust) on the same 6-unit geometry. **The advice therefore inverts between the two sides:**
   "cluster finer than the cell" rescues the binary att test and actively harms the continuous
   ame test. One theorist predicted this against the other two; the mechanism he gave is that
   unit clustering manufactures noisy, few-df off-diagonal terms in the ame covariance where
   robust structurally has zeros that happen to match the truth. This belongs in the help file.

The att twins in both continuous configs stayed correctly sized (4.96%, 5.16%), so the binary
result is robust to having `contreat()` in the model.

---

## 5. Two methodological findings that changed the campaign

**F-1 — a band pass is not proof of correctness.** The theory review caught an error in the
task briefing (A5 was cited at 48.4%; the record says **2.38%**, undersized, with uniformity
chi2 = 77.5). 48.4% belongs to D2 and 65.7% to A5's jointtests figure. The lesson: a structurally
blind covariance matrix can be dragged *through* a plausible rejection rate by a lenient
F(q, df_r) reference. Uniformity was therefore made **co-primary** and every row reports the
full p-value CDF. This campaign confirms the risk is real: AB8 sits at 0.00% rejection — a rate
a naive reader would call "conservative, therefore safe" — with chi2 = 1072.

**Corollary for the earlier SI-A1 work: coarse cohorts do not rescue cohort-level clustering.**
D2 (8-unit cohorts, cluster(i)) measured 48.4% and AB12 (6-unit cells) measures 56.35% — both
*worse* than the singleton-cell A1 at 20.5%. Making cells bigger while clustering at the cohort
level makes containment worse, not better.

**F-2 — the tested geometry is 8 coefficients, q = 6.** Derived from the code before locking and
confirmed empirically by the 36-entry V_es lower triangle in the smoke run. This fixed the
altmethod comparison (q = 4 vs 6, §2) and the b/SE dump width.

---

## 6. Theory scorecard — locked predictions vs measurement

25 of 26 scored quantities inside the locked consensus band (the one miss: AB14a power,
predicted 75% [45,95], measured 98.81% — power was underestimated, which cannot affect any
ruling). Named disagreements, all resolved by measurement:

| # | question | resolved | detail |
|---|---|---|---|
| **D-2** | is size governed by **clusters per cell** or by cluster count/cell size? | **clusters per cell** | T1 predicted AB13 = 14.5 vs the others' ~7.5–10.5; measured **14.50**. He named it his own falsifier and it went his way. **This is the finding a runtime gate must be written on.** |
| D-1 | AB5 few clusters: conservative flip or blowout? | blowout | 17.46%, chi2 231 — T2/T3 right, T1 wrong |
| D-3 | AB12 magnitude | T2 right | predicted 58, measured 56.35 |
| D-4 | does fine clustering help the ame? | **no, it hurts** | T2 predicted 44 (worse than robust); measured 48.41 |
| D-5 | does AR(1) degrade size under matching clustering? | no | AB6 5.60% ≈ AB1 4.80% — T2 right |
| D-6 | AB3 magnitude | between | measured 9.10 (T3 10.5 closest) |
| D-7 | does staggering degrade? | **yes, mildly** | AB15 7.20% + chi2 24.2\* vs AB1 4.80% — T3 right; the others' frameworks were onset-blind |
| D-8 | is the U-gradient leverage or inversion noise? | mostly inversion noise | AB3 − AB2a = 1.70pp (T1/T2 ~1–2.2pp right; T3's 5pp wrong) |

---

## 7. What should ship (recommendations — no code written)

1. **Gate on clusters-per-cell, not cluster count or cell size** (D-2). Hard error when any
   treated cohort-by-period cell lies inside a single cluster (m = 1) — the containment case,
   56.35%. The `wooldid_amescreen mode(detect)` predicate already computes this.
2. **Caution below the certified floor:** m < 6 clusters per cell, or fewer than ~180 clusters,
   or (robust) cells smaller than 6 observations.
3. **The gate must cover the `e(V_es_*)` honestdid handoff, not just the test.** The object at
   risk is the covariance matrix; gating the test alone relocates the failure.
4. **Ame-side `adversarial`: GATE, do not remove** (§4, revised by AB22). It is correctly
   sized under robust at ~50 observations per cell (5.95%) and severely oversized at 6 (31.55%).
   Suggested shape: permit above a cell-size threshold, caution/refuse below it, and route the
   refused cases to RI / bootstrap / IM. The exact threshold in (20, 50] needs one cheap
   follow-up ladder before it is documented.
5. **Document the inverted advice** for att vs ame clustering (§4, finding 2).
6. **DOC1 (fires): cohort-period shocks.** Even a shock worth 3.8% of error variance gives
   12.30%. Escaping containment by clustering finer assumes away precisely the dependence that
   motivates cohort clustering, and this is undetectable at runtime. The rule's
   KEEP-GATED-WITH-PROMINENT-WARNING downgrade does **not** fire (AB7a 12.30% < 15), but the
   caveat is mandatory.
7. **DOC2 (fires): robust under serial correlation is inert, not safe** — 0.00% with chi2 1072.
   Never report undersizing as a pass.
8. **`altmethod_adversarial` earns its keep** (§2) and is the better default wherever cells are
   small — subject to confirming AB17 at higher M first.

**Open, not answered by this grid:** nothing certified below 180 clusters while 60 clusters came
in marginal (AB16 6.70%), so the cluster-count boundary lies somewhere in (60, 180] and n\* ships
conservatively until a cheap follow-up localizes it. The preregistered FOLLOWUP clause also calls
for mapping m ∈ {3, 4} before documenting m\* = 6 as tight.


---

## 8. AMENDMENT 3 (2026-08-26): binary results at 50 observations per cell

Six CRN-paired configs at 30 cohorts × 50 units × 12 periods, varying only the VCE and the error
structure. Report-only, predictions locked before running (5 of 6 inside band).

| config | clustering | clusters/cell | errors | rate | chi2 | verdict |
|---|---|---|---|---|---|---|
| AB24 | robust | — | iid | 4.37% | 12.6 | passes |
| AB23 | unit level | 50 | iid | 4.37% | 14.0 | passes |
| AB27 | unit level | 50 | AR(1) | 4.96% | 6.5 | passes |
| AB25 | half-cohort | **2** | iid | **17.86%** | 235.0\* | fails |
| AB26 | cohort level | **1** | iid | **59.33%** | 2009.5\* | fails |
| AB28 | unit level | 50 | cohort shock σ=.5 | **98.02%** | 4476.2\* | fails |

### Cell size read across every regime

| regime | 6/cell | 20/cell | 50/cell | effect of larger cells |
|---|---|---|---|---|
| binary, clustering finer than the cell | 4.80% | 5.00% | 4.37% | none — already correct, stays correct |
| **continuous, robust** | 31.55% | 9.92% | **5.95%** | **cures it** — the one place size is the remedy |
| binary, 2 clusters per cell | 14.50% | — | 17.86% | worse |
| binary, cells inside one cluster | 56.35% | — | 59.33% | worse (20.5% at singleton cells) |
| binary, cohort-period common shock | 44.84% | — | **98.02%** | **catastrophically worse** |

**Finding 1 — the clusters-per-cell account is now cleanly separated from cell size.** AB13
(m=2, 6-unit cells) could not distinguish "two clusters per cell is fatal" from "small cells are
fatal", because both moved together. AB25 holds cells at 50 and drops to m=2: **17.86%, worse
than AB13's 14.50%**. Clusters per cell governs; cell size is not protective and is mildly
aggravating. **A guard written on cell size would pass AB25 and be wrong.**

**Finding 2 — larger cells make the common-shock failure catastrophic: 44.84% → 98.02%.** This is
the campaign's one badly missed prediction (locked at 65%, measured 98.02%); the predicted
*direction* was right and the magnitude was not. Mechanism: the cohort-period shock is absorbed
into each cell's own dummy and leaves no residual trace, so no analytic VCE can see it, while the
idiosyncratic component that the VCE *does* see shrinks like 1/U. The ratio of true to estimated
variance therefore grows with cell size. **Practical consequence: "use large cells and cluster
finely" is sound only absent cohort-period common shocks, and is close to worst-case advice when
they are present.** This strengthens the DOC1 rider from mandatory-caveat to
prominent-warning status for any guidance that recommends enlarging cells.

**Finding 3 — the binary test is indifferent to cell size once clustering is fine enough**
(4.80% → 5.00% → 4.37%), the exact mirror of the continuous side, where cell size is everything.
Both facts belong in the design-conditions documentation together.

Campaign totals after amendment 3: **31 configs, 20,584 reps, 0 bad rcs, 32 of 34 locked
predictions inside band.**
