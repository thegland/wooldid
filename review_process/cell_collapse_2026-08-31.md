# Collapse-to-cells for the margins calls: feasibility, boundary, and measured payoff (2026-08-31)

Investigation of the Tier 2 "collapse to cells before running margins" proposal from
`efficiency_review_2026-08-31.md`, run as a 10-agent workflow: source analysis, a shared
harness, a 12-cell correctness matrix across three VCE regimes x continuous treatment x
joint tests, an adversarial refutation pass, and a dedicated sequential timing pass.

**Verdict in one line: the idea works and is worth real money at scale, but its stated
premise is false on more configurations than anyone expected, and the safe gate is much
narrower than "no continuous treatment".**

---

## 1. What was being tested

Every wooldid estimand is a weighted average of cell-level treatment effects: the
regression carries one dummy per `(_i, _t, targetflag)` cell, and the reported effect is a
weighted average of those cell coefficients. The proposal's premise was that *because
targetflag does not interact with controls, its derivative is constant within a cell* -
so margins can be handed a dataset collapsed to one row per cell with per-cell summed
weights, and return the same answer far faster.

The correctness question is not whether the algebra is right. It is **on which
configurations the premise actually holds**.

---

## 2. The premise is false in three distinct ways

Two of these were found by adversarial attack and were in nobody's plan.

### 2a. `interactivecontrols()` - the serious one. NEW vector, no `contreat()` involved.

`wooldid.ado:1058` builds `i._i#i._t#i.targetflag#c.(interlist)`. That interacts
targetflag with cohort-demeaned continuous controls, so the discrete 0-to-1 effect inside
a cell is `b_cell + g_cell * ic_obs` and **varies observation by observation**. The
collapse key carries no `ic`, so the cell row gets one arbitrary observation's value
standing in for the whole cell's distribution.

Measured on `dta_s_admin.dta`, non-degenerate at 556:1, and verified independently by the
orchestrator from the stored dumps:

| quantity | full data | collapsed | error |
|---|---:|---:|---|
| headline ATT | 0.646391 | 0.667649 | 3.3% |
| ATT se | 0.028539 | 0.034027 | 16% |
| ATT variance | 8.145e-4 | 1.158e-3 | 1.42x |
| horizon pre1 | **-0.010943** | **+0.011223** | **sign flip** |
| horizon post0 | 0.453528 | 0.349616 | 23% |
| horizon post1 | 0.513273 | 0.701422 | 37% |
| V[8,8] | 1.700e-3 | 8.998e-3 | 5.29x |

This is mainline usage - the source comment at line 1040 calls interactive controls
"closer to the object that wooldridge covers in his paper for covariates". A fast path
gated only on `contreat()` being empty would silently corrupt these runs.

The A/B partner isolates the cause exactly: plain `controls()` on identical data, estimand
and sites is **ULPCLASS** (b 3.9e-15). Only the interaction term breaks it.

### 2b. `contreat()` - breaks the *binary* site, not the continuous one

With `contreat()` the regression carries a targetflag-by-contreat interaction, so the same
mechanism applies. Measured (cluster-multiple-per-i, `continuous_control_compact.dta`):
b 7.0e-4, V 1.9e-2, se 9.8e-3 at the `dydx(targetflag)` site - 10-11 orders above the
ULPCLASS ceiling.

The instructive part: **both contreat-operator sites are BITEXACT**, including the
contreat joint test's p-value. What breaks is the binary targetflag margins, purely
because `contreat()` is present. Because the damage lands on the point estimate, it is
predicted to be VCE-regime-independent.

**A repair was measured, not speculated:** carrying the contreat variable into the cell
row as its weight-weighted cell mean `sum(w*contreat)/sum(w)` - rather than an arbitrary
value - restores ULPCLASS at all four sites. This follows from the algebra: the AME
gradient entry for the cell's dose-slope coefficient is exactly `W_cell * meandose_cell / W`.
Extending the collapse *key* with the contreat variable is a dead end: dose is unit-level,
so the key extension yields `ncell == nfull`, i.e. no collapse at all.

### 2c. The operator matters, independently of the model - `semielasticity` (eydx)

`controls()` are not absorbed (line 1285 puts them in reghdfe's main regressor list), so
they are visible to margins and `yhat` varies within a cell. For `dydx`, the controls'
contribution appears in both `xb(tf=1)` and `xb(tf=0)` and **cancels**, leaving the cell
coefficient - constant within cell, collapse exact. For `eydx`, the division by `yhat`
retains them, so the averaged quantity varies within the cell. Measured: b 7.5e-4,
V 1.3e-3, 15 difference groups at rtol 1e-13.

So the premise depends on the **operator**, not only on the model.

### 2d. `unconditionalse` - structural, and must fail closed

`vce(unconditional)` builds its variance from observation-level influence functions and
**verifies that the visible row count equals `e(N)`**. Collapse changes that count by
construction, so the check can never pass: `rc 0` on full data, **`rc 459` collapsed**.
Isolated in a self-contained 800-row probe with no wooldid involved: under
`vce(unconditional)` the collapsed call errors; under the delta VCE the *same collapse* is
bit-exact in both b and V. The VCE regime is the entire difference. This is not a
tolerance question and is not fixable by marking `e(sample)` on collapsed rows.

---

## 3. The 12-cell matrix: VCE regime is irrelevant, `contreat()` is decisive

Three VCE regimes x continuous treatment x joint tests. Cluster structures were verified
empirically in-band, in the run that produced the dumps.

| VCE regime | contreat | jointtests | worst rel. b | adjudicated |
|---|---|---|---:|---|
| robust | no | no | 1.47e-13 | **ULPCLASS** |
| robust | no | yes | 1.47e-13 | **ULPCLASS** |
| robust | **yes** | no | **4.54e-01** | **BROKEN** |
| robust | **yes** | yes | **4.54e-01** | **BROKEN** |
| cluster, multiple per i | no | no | 1.47e-13 | **ULPCLASS** |
| cluster, multiple per i | no | yes | 1.47e-13 | **ULPCLASS** |
| cluster, multiple per i | **yes** | no | **7.04e-04** | **BROKEN** |
| cluster, multiple per i | **yes** | yes | **4.54e-01** | **BROKEN** |
| cluster, one per i | no | no | 1.47e-13 | **ULPCLASS** |
| cluster, one per i | no | yes | 1.47e-13 | **ULPCLASS** |
| cluster, one per i | **yes** | no | **7.04e-04** | **BROKEN** |
| cluster, one per i | **yes** | yes | **4.54e-01** | **BROKEN** |

Realized structures: multiple-per-i used `dta_s_admin` with `i=cohid`, `cluster(fips)` -
distinct clusters per cohort min 194 / mean 500 / max 814. One-per-i was verified two
independent ways at exactly 1 (min = mean = max). Robust cells confirm `e(clustvar)`
empty and `e(N_clust)` missing in every dump.

**Two conclusions, both clean:**

1. **The VCE regime does not affect correctness at all.** Every no-contreat cell lands on
   the identical 1.47e-13 / 1.776e-15 quantum regardless of regime; every contreat cell
   lands on the identical 7.04e-04 and 4.54e-01. Robust, many-clusters-per-cohort and
   one-cluster-per-cohort are indistinguishable. This is expected in hindsight - the
   clustering lives entirely in `e(V)`, which the collapse never touches - but it was the
   open question that motivated the matrix, and it is now measured rather than assumed.
2. **`contreat()` is decisive, and the damage is in the point estimate.** The same numbers
   appear in all three regimes because the defect is in `b`, not `V`. This confirms, across
   three independent regimes, the prediction made from the mechanism.

Joint tests add no failure mode of their own: they inherit whatever the underlying site
does. Where the site is clean the joint p moves 1.4e-14; where the site is broken the
joint test is broken with it (in one contreat cell the joint p moved by a factor of
11,200 and F moved 70.1 to 39.4 - a genuine test-statistic change, not tail amplification).

### Adjudication of the disputed verdicts

Sibling agents split on the no-contreat cells: some returned ULPCLASS, others BROKEN, on
*the same measured number*. The pre-registered gate said "max relative difference <= 1e-13
on b, V and se", and 1.47e-13 breaches it by 1.47x. **I adjudicate these ULPCLASS**, and
the gate wording - mine - was the thing that was wrong. Four independent pieces of
evidence:

1. The absolute error is **exactly 2^-49**, one quantum, *identical* on two cells with
   different values and different subpops. Algebraic error does not produce a fixed
   power-of-two discrepancy. **CORRECTION (2026-09-01):** 2^-49 is not a general constant -
   it is `ULP(|xb|)` on this dataset, where |xb| sits in binade 2^3. A later amplification
   ladder tracked the predicted binade exactly across three decades of |xb|. The scale-free
   statement is that the collapse perturbs the answer no more than permuting the input rows
   does; see `cell_collapse_followup_2026-09-01.md` section 4.
2. Every large coefficient is **bit-exact**. The same quantum on 0.4575 would read 3.9e-15;
   it appears only where cancellation makes the result near-zero (-0.0199, -0.0121).
3. **The partial-collapse ladder.** Refining the collapse key with an arbitrary bucket
   sweeps K from 60 to 27,913 against N=29,067. Every rung is algebraically the same
   estimand, so all movement is pure floating point. Measured: the gap is a
   0-or-one-quantum flip, **discrete rather than a monotone decay** (G=16 lands bit-exact,
   G=64 is off by one quantum again), bounded by one quantum at *every* rung including
   maximal collapse, and **identically zero at every rung with K >= 12,316**.
4. **Calibration against the estimator's own floor:** pristine, uninstrumented wooldid with
   only the *row order* permuted moves the same coefficient by **2.8e-12** - 19x more than
   the collapse does, and 28x the gate itself. For these near-zero quantities the 1e-13
   gate sits below wooldid's own reproducibility floor, so no implementation could pass it.

`V` and `se` pass the gate outright everywhere (worst 4.3e-16 and 1.9e-16).

---

## 4. Where the premise holds - attacks that failed to break it

Reported explicitly, because negative results are the evidence that the gate is drawable.
All non-degenerate, all with asserted stressors:

- **Weights x all four estimands.** `att`, `att_it`, `att_i`, `att_itime` under `aweight`,
  with within-cell weight variation in 60 of 60 cells and `pop` spanning 300 to 3,161,661.
  All ULPCLASS (worst b 2.7e-14). All four have different normalisers and all four
  collapse, because margins forms `sum(w*f)/sum(w)` with `f` constant in the cell.
- **`pweight` + `customweightmultiplier()`** - ULPCLASS (b 2.7e-15). `fweight` cannot
  reach wooldid at all (the syntax line accepts only `aw pw`).
- **Ragged/unbalanced panels** - a whole relative period missing for one cohort, a
  singleton cell, a cell thinned to 3 rows - ULPCLASS (b exactly 0).
- **Subpop realisation**, the edge predicted to be sharpest: an `if` removing 68% of rows
  and cutting 60 of 60 cells partially - ULPCLASS (b 2.7e-15). **Important correction to
  the natural explanation:** having the if-variable in the collapse key is a *vacuous*
  defence, since it is constant on every surviving row. The real defence is `e(sample)`
  containment - the user's `if` is folded into `e(sample)` before either arm runs.
- **Plain `controls()`** - ULPCLASS, isolating 2a's interaction term as the entire cause.

---

## 5. A methodology finding that invalidates naive benchmarks

**The overall collapse ratio is the wrong non-vacuity metric.** `nfull/ncell` measures the
collapse over all of `e(sample)`, but each margins site averages only over its **own
subpop**. Where subpop cells are singletons, the collapse is a no-op *exactly where it
matters* and bit-equality is free and meaningless.

This invalidated three passes inside the investigation itself, and it retroactively
qualifies any result taken on `test/data/dta_binary.dta`, whose main sites are singleton
(subpop rows 275/100/200 equal their weight-sums). `test/data/dta_path_relyrgap.dta` - the
canonical "ragged" fixture - is degenerate everywhere at one row per unit-period and
**cannot test collapse at all**.

**Any future collapse benchmark must publish its per-site subpop ratio.**

---

## 6. Performance - measured, and the headline number does not survive

The prior lens's "~12x" is not wrong so much as not a constant. Timing was taken by a
dedicated agent, one Stata at a time, arms interleaved within a session, medians of 5,
with every run screened for machine contention.

**The collapsed arm's cost is flat in N; the full arm's scales with it.** Admin panel,
K=60, p=32, robust + jointtests:

| site | N=29,067 | N=116,268 | N=465,072 |
|---|---|---|---|
| main att | 0.114 / 0.027 = **4.2x** | 0.331 / 0.027 = **12.3x** | 1.116 / 0.027 = **41.3x** |
| event study | 0.197 / 0.064 = 3.1x | 0.574 / 0.064 = 9.0x | 2.518 / 0.063 = 40.0x |
| joint test | 0.197 / 0.059 = 3.3x | 0.506 / 0.059 = 8.6x | 2.283 / 0.059 = 38.7x |

Measured cost law over 49 clean site-rows (N from 576 to 465,072, p from 32 to 683):
**saving per call ~= 1e-07 x (N-K) x p seconds.**

**Two results that contradict expectations, including mine:**

1. **VCE regime and joint tests barely matter.** Robust and clustered agree within
   rep spread at every size. Only N/K and p move the ratio. This directly answers the
   question that motivated the matrix.
2. **The collapse bites LEAST at high K, not hardest.** `p` grows with K and *both* arms
   pay the p-work. The 150k/high-K panel has the largest absolute saving per call (up to
   6.5s) and one of the smallest ratios (2.46x).

**End-to-end gain, the number a user feels** (all inputs measured; net charges one
collapse against the saving):

| config | N | margins share | net gain |
|---|---:|---:|---:|
| 960-obs binary | 960 | 76.6% | **3.3%** |
| admin, robust | 29,667 | 63.2% | **36.3%** |
| admin, clustered | 29,667 | 66.8% | **38.1%** |
| admin x4 | 118,668 | 60.8% | **40.0%** |
| admin x16 | 474,672 | 62.7% | **50.5%** |
| 150k high-K | 150,000 | **21.8%** | **12.0%** |

The collapse itself is essentially free (0.013-0.122s) and **one collapse serves the whole
run** - the by-key is identical at every site, only the weight column differs. Break-even
is inside the first margins call at admin scale and above.

### A correction to the efficiency review's headline

The review states that essentially all wall time is inside margins, with a ~3% ceiling on
everything else. **That is config-specific and fails at high K.** On the 150k/40-cohort
panel margins is only **21.8%** of a 147-second run; the other 115 seconds is the
683-coefficient regression and wooldid's setup. Collapsing every margins call to zero
would cap the gain there at 21.8%. The 95% figure came from an admin-shaped config, where
the measured share is 63-67%.

---

## 7. The corrected gate

Not "`contreat()` is empty". The conjunction of:

1. the derivative operator is **`dydx`** (not `eydx`/`eyex`/`dyex`), **and**
2. **no model term interacts `targetflag` with anything varying within a cell** - which
   excludes `contreat()` and `interactivecontrols()`, unless the measured weighted-cell-mean
   aggregator is used for the interacting variable, **and**
3. the VCE is the delta method - **`unconditionalse` must fail closed** to the full path, and
4. every `over()` variable at that site is in the collapse key - which the current key does
   **not** satisfy for the `over(cluster)` sites (reachable under `impvalues`) or for
   `over(cfxlatticeplot)`.

---

## 8. Open items

- **Poisson is unmeasured.** `poisson` + `dydx` reaches `predict(xb)`, the linear index, so
  it is reasoned-safe by the same argument that clears plain `controls()` - but it was not
  measured, and it deserves its own check before being cleared.
- **Clustered VCE under adversarial stressors.** The stressor runs used robust or
  unweighted VCE. Cluster collapse-safety rests on the matrix cells, which is good evidence
  but is not the same as having been attacked.
- **The near-zero-coefficient classification question.** Event-study pre-trend cells show
  large *relative* differences from tiny *absolute* ones (1.47e-13 relative from 1.78e-15
  absolute on a coefficient of -0.012, while every large coefficient is bit-exact). A
  calibration run found that **pristine wooldid with only the row order permuted moves the
  same coefficient by 2.8e-12** - 19x more than the collapse does, and 28x the 1e-13 gate.
  For these quantities the gate sits *below the estimator's own reproducibility floor*, so
  the pre-registered threshold was too strict, not the result. Adjudicated as ULPCLASS.

---

## 9. Recommendation

**Do not adopt this as a general fast path. Adopt it, if at all, as a narrowly gated one -
and only after deciding the baseline question.**

The case for: at admin scale and above the net end-to-end gain is 36-50%, the collapse
costs essentially nothing, one collapse serves the whole run, and correctness is
indifferent to the VCE regime. That is a real win on exactly the panels where wooldid is
slow enough to hurt.

The case against, and it is substantial:

- The gate is a **four-way conjunction** (section 7), not a one-line check, and two of its
  four clauses were discovered only by adversarial attack. A gate that subtle is a
  maintenance liability in a package where an option interaction that looked harmless has
  already produced one shifted-coefficient bug this month.
- **`interactivecontrols()` is mainline usage** and fails catastrophically - sign flips and
  37% coefficient errors, silently, at rc 0. Any gate that misses it corrupts results
  rather than merely being slow.
- The gain is **smallest exactly where users hurt most**: 3.3% on a small panel, and 12% on
  the high-K panel where runtime was 147 seconds. High K is the axis the review identified
  as the real scaling problem, and collapse is weakest there because `p` grows with K and
  both arms pay the p-work.
- Even where it works it is **ULPCLASS, not bit-exact**, so it needs a versioned change and
  re-blessed baselines.

**If it is adopted, the safest shape** is not a general fast path but a targeted one: apply
the collapse only where the gate is provably satisfied, verified per-site at runtime rather
than inferred from options at parse time, defaulting to the full-data path whenever any
clause is unproven. `unconditionalse` must fail closed. And the contreat and
interactivecontrols cases should be left on the full path rather than repaired with the
weighted-cell-mean aggregator - the repair is measured and works, but it adds a second
correctness argument to maintain for the configurations that are hardest to reason about.

**The question to settle first is not technical.** Both this and the direct-Mata option
move low-order digits on every reported number. If the project is willing to re-baseline
the suite for a performance change, collapse is clearly the better of the two - it keeps
margins as the engine of record and so inherits all of its bookkeeping instead of
re-implementing it. If the project is not willing, neither option is available, and the
remaining efficiency work is the conventional kind already completed.
