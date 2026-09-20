# hv-cwm — FINDINGS

---
## hv-cwm-01  (SEEDED, half 1: CONFIRMED)
### customweightmultiplier() silently reweights the treatment cells, so att_it / att_i / att_itime stop being the estimands their help entries define

**severity** medium · **confidence** certain · **realistic** true
**surface** `customweightmultiplier()` x `att_it` / `att_i` / `att_itime`
**repro** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-cwm/repro_cwm_estimand.do`
(`stata-mp -e do repro_cwm_estimand.do`; read the VERDICT block in the `.log`)

**claim.** With `customweightmultiplier(v)` the att_it estimate is not "the simple, unweighted mean
of all of the cohort-period-specific estimates" (sthlp:220-221) but a mean of those estimates
weighted by `sum_cell(w*v)/sum_cell(w)` — i.e. by `v` itself when `v` is constant within the cell.
Same for att_i (sthlp:227, "the simple, unweighted mean of the cohort-specific treatment effects") and
att_itime. Nothing in the output names the substitution.

**mechanism (verified in code and against a raw-data oracle).** `customweighttext` = `* v`
(ado:1581) is appended to the NUMERATOR of all four estimand weights (ado:1630, 1637, 1644, 1650)
while the normalizers are built from the sample weight alone — `W_it = total(esamp*w) by(i,t)`
(ado:1619) and `W_i = total(w*esamp) by(posttreatment,i)` (ado:1620). So each cell's total
aggregation weight is `sum(w*v)/sum(w)` instead of 1. Measured: the `w_att_it` column of
`e(gteffects)` matches that raw-data quantity to <= 2.9e-16 in every arm of `p1_cwm_estimand.do`,
including a `v` that varies WITHIN the cell (cell weight = the w-weighted mean of `v`, arm D).

**observed** (8 counties, 4 staggered adopters, 10 years, 5 establishments per county-year,
`pop` = 1 for six counties and 50 for two of the treated ones — an ordinary 50:1 population spread):

| | no cwm | `customweightmultiplier(pop)` |
|---|---|---|
| att       | 2.7147886 | 3.7412488 |
| **att_it**| **2.74890010** | **3.74381423** |
| att_i     | 2.9520313 | 3.8071510 |
| att_itime | 2.9518714 | 3.8078950 |
| simple mean of the gteffects cells | 2.74890010 | 2.74890010 (unchanged — the cells are regression coefficients) |
| `w_att_it` per cell | 1 … 1 | 1 … 50 |
| per-cohort `w_att_i` (pre and post blocks) | 1,1,1,1 | 1,1,50,50 |

att_it is **36.2% away** from the estimand its help entry defines, and the reported
pop-weighted mean of the cells reproduces it exactly (3.74381423, bit-for-bit).

**expected.** The att_it / att_i / att_itime entries should cross-reference
`customweightmultiplier()` and say that it reweights the cells (respectively the cohorts), and
the runtime note should name the estimands whose definition it changed. As written, sthlp:218
tells the `att` reader the multiplier applies, sthlp:1107 says the estimand weights are
"multiplied by this variable", and sthlp:220-221 / 226-227 / 232-234 tell the cell-weighted reader the
estimand is a simple unweighted mean — with no cross-reference between them, so whichever entry
the user reads is the one they believe.

NOT proposed as the fix: renormalizing the cell (`w*v/sum_cell(w*v)`). That would make `v` an
exact no-op on att_it in a binary design, and it would destroy a real and otherwise unreachable
estimand — cells weighted by `v` with the sample weights acting only WITHIN the cell, which is
what distinguishes arm B's att_it (3.7438, cells weighted by pop) from arm B's att (3.7412,
cells weighted by pop x the cell's sample-weight mass). The problem is the silence, not the
arithmetic.

**consequence surface.**
1. **The one comparison the four estimands exist to support is silently voided.** att vs att_it is
   the package's own weighted-vs-equal-cell heterogeneity check. Under cwm(pop) they read
   3.7412 vs 3.7438 — a 0.07% gap — so a user concludes "the result is insensitive to cell
   weighting". The true equal-cell number is 2.7489, **27% lower**. That is a wrong substantive
   conclusion produced from correct-looking output.
2. **In a plain binary design the ONLY thing the multiplier does to att_it is the cross-cell
   reweighting.** There is exactly one coefficient per (i,t) cell, so within-cell weights cannot
   move a cell effect. The user trades the documented equal-cell property for a different
   estimand (cells weighted by `v`) without being told the trade happened.
3. **No escape inside one run.** All four estimands take the multiplier, so a user who wants a
   population-weighted `att` cannot get an equal-cell `att_it` from the same call; they must run
   wooldid twice and know to.
4. **A 0/1 cwm used as a subpopulation selector is worse, not better** (p2 arm G): with a mixed
   cell the cell weight becomes the within-cell subgroup SHARE (measured range 0.20-0.86), so
   cells are weighted by how female/rural/… they are, and the cell effects being averaged are
   still the whole-cell effects.
5. The `e(gteffects)` `w_att*` columns are honest — they publish the real weights, and the
   sthlp:1136 identity (weighted mean of `estimate` reproduces the aggregate) holds at ~1e-16 in
   every arm. A user who requests `gteffects` and reads the weight column can see this. A user who
   does not, cannot.
6. The only runtime signal is `Note: all Estimates Shown Apply pop Weights Post-Estimation`
   (ado:11369/11529), which reads as confirmation of what was asked for, not as "att_it is no
   longer an equal-cell mean".
7. Secondary, **measured in the repro's ARM D**: `W_i` is `total(w*esamp) by(posttreatment, i)`
   (ado:1620), so with a `v` that varies over time within a cohort (`pop*1.15^year`) the
   PRETREATMENT att_i row and the main att_i row weight the cohorts differently — the placebo row
   printed directly above the estimate is not the same weighted functional as the estimate.
   Measured per-cohort share of its block's total `w_att_i`:

   | cohort | main row | pretreatment row | difference |
   |---|---|---|---|
   | 1 | 0.008400 | 0.007422 | -11.6% |
   | 2 | 0.009057 | 0.008119 | -10.4% |
   | 3 | 0.479530 | 0.463144 |  -3.4% |
   | 4 | 0.503014 | 0.521315 |  +3.6% |

   With no cwm both blocks are exactly 1 per cohort, so the two rows are the same functional.

---
## hv-cwm-02  (SEEDED, half 2: **REFUTED**)
### The bootctl cell-weighted gate's reasoning is correct; "either cwm is a no-op or the gate is wrong" is a false dichotomy

**severity** low · **confidence** certain · **realistic** n/a · **surface** `bootstrap(bayes)` gate, ado:6430-6442
**repro** same file, ARM C (and `p1_cwm_estimand.do` arm E)

**claim under test.** "bootctl's cell-weighted gate comment reasons about cwm in a way that is
inconsistent with this, so either cwm is a silent no-op on cell-weighted estimands or that gate's
reasoning is wrong."

**observed.** Neither. The gate (ado:6437-6442) says "att_it/att_i/att_itime aggregate with equal
weight per treatment cell, so a within-cluster weight multiplier cancels out of them". It is
talking about the Dirichlet draw, which `wooldid_bootrun` applies as the **sample weight** —
`qui gen double repw = basew*drawv` then `wooldid … [aweight=repw]` (ado:9448-9451) — so it lands
in `W_it` as well as in the numerator and cancels **exactly**. Measured: arm C
(`[aw=sw*pop]`, no cwm) gives `w_att_it` min=max=1, and p1 arm E ([aw=w*m] with `m` varying WITHIN
the cell, i.e. a cluster finer than `i`) also gives min=max=1 — cancellation is exact in both
regimes, since `sum_cell(w*m)/sum_cell(w*m) = 1` identically. `customweightmultiplier()` is a
different path: it is deliberately NOT folded into `w` (ado:624-626 says so explicitly, because
folding it in would break the fweight demotion), so it appears in the numerator only and does not
cancel. Both statements are true simultaneously; the gate needs no change.

**small residual (not filed as a defect).** With a cluster finer than the cell the bayes draw does
perturb within-cell composition even though the cell TOTAL stays 1, so the refusal message's
rationale ("cluster reweighting cannot resample treatment cells") is narrower than the refusal.
No wrong number: the combination is refused under default safety anyway.

---
## Probed and CLEAN
- **Default cell collapse vs `nofastcollapse` under cwm**: binary and `contreat(dose)
  contreatcontrols(basic)`, with `v` varying within the cell and with a `v` that is exactly **0 on
  the cell's first row** — the case ado:5636 names as having once produced a wrong number.
  Max relative difference on `e(b)`: 1.3e-13 (contreat), 3.2e-16 (binary). The invariant asserted
  at ado:5627-5644 ("every estimand weight is w x cwm over a cell-constant denominator") holds.
  (`p3_cwm_collapse.do`)
- **Cross-path weight consistency**: `e(escomp)`'s `wshare_att_it` column reproduces the horizon
  shares rebuilt from `e(gteffects)`'s `w_att_it` column exactly (102 at each of post0-3, 101 at
  pre3 where only 3 cohorts are estimable). The F8 block's independent reconstruction of the
  weights (ado:3832-3834) agrees with the main construction. (`p3_cwm_collapse.do` arm L)
- **The replication engines under cwm** (`p4_cwm_engines.do`): `bootstrap(standard) reps(25)
  cluster(county)` and `permute(shuffle) permreps(25)` both leave the reported point estimates
  exactly where the plain run put them (max relative difference 0 and 1.19e-16), so the
  sthlp:532 invariant holds with a multiplier in play; `bootstrap(bayes)` with `att_it` + cwm is
  refused r(198) by the cell-weighted gate as designed; `bootstrap(bayes)` with `att` alone runs
  and its point estimate is unchanged (0).
- **The `e(gteffects)` identity under cwm**: the weighted mean of `estimate` under each `w_*`
  column over `relyr>=0` reproduces the corresponding reported aggregate at <= 4.8e-16 in all
  five p1 arms. The published table does not lie about the weights it used.

---
## Read-only observations (no probe run; recorded so nobody re-derives them)
- `wooldid_ptp` builds the same asymmetry deliberately and says so: `ptpwa = wvar*cwm` but
  `ptpwpit = total(wvar)` by (i,t), with the comment "cell totals use the raw weight, without any
  custom weight multiplier - the documented asymmetry" (ado:10392, 10401-10403). So
  `pretrendpower()`'s hbar scalars are consistent with the estimands as implemented, and would
  have to move with them under any fix.
- A **missing** `customweightmultiplier()` value zeroes `ifvar` (ado:704), i.e. the row leaves the
  ESTIMATION SAMPLE and therefore changes the stage-1 coefficients — it is not merely excluded
  from the aggregation. sthlp:1107's "missing values drop the observation" covers this, but
  sthlp:218's "that multiplier enters the margins weights only" reads against it.
