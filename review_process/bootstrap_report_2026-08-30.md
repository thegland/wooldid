# Bootstrap standard errors — chunk report (2026-08-30)

Work-file chunk: *"implement a facility for producing bootstrap clustered SEs … standard
issue cluster bootstrapping … also the Bayesian bootstrap … borrow reps(), ties, level(),
dots(), and rseed() … everything does get bootstrapped, even the jointtest p-values …
ideal would be working with the built in stata commands … DIY is fine if we can prove we
match."*

Plan of record: `review_process/bootstrap_plan.md` (sections 1-10 design, 11 author
rulings, 12-13 review amendments). This report is the closeout: what shipped, what it was
measured against, and what remains open.

---

## 1. What shipped

| Commit | Contents |
|---|---|
| `4b08962` | **A — prefix hardening.** Worker-entry `set coeftabresults on`; shell `properties(svyr)`. |
| `8cfd915` | **B — option surface + `wooldid_bootctl` gates** (loud not-yet-enabled tail). |
| `1ad2f65` | **C — Bayesian engine + replace-mode `e()` rewiring + deferred display.** |
| `d3d9e99` | **D — standard mode**: count path + physical re-id leg. |
| `8d8316f` | **E — joint tests, subgroup contrasts, CI trio, `ties`.** |
| `58e4467` | **bayes default alpha flipped to `equal`** (author-ruled on MC evidence). |
| `6693e2f` | **SI-17** — `permute()` refuses `cleanmatrices` (adjacent pre-existing defect). |
| `133b11b` | **F — docs, suite configs, calibration MC, this report.** |
| `4906bd3` | **`bootdf()`** — t critical values as a user choice (author-requested from the MC finding in section 3). |

### The user-facing surface

`bootstrap(standard|bayes)` with `reps(#)` (mandatory, no floor — author ruling),
`ties`, `level(#)`, `dots(#)`, `rseed(#)`, and `bootalpha(equal|size)`.

**Replace-mode inference** (author ruling Q1): the reported standard errors, z
statistics, p-values, confidence intervals, joint tests and subgroup contrasts *are* the
bootstrap quantities — `vce(bootstrap)` style — rather than an appended block. Point
estimates never change (asserted bit-exactly in the cert suite). The analytic covariance
is preserved in hidden `e(V_analytic)` / `e(V_es_analytic_*)`.

**Two engines, one rep loop.** Both modes inject a per-replication weight rather than
physically rebuilding data, except where the estimand demands otherwise:

- `bootstrap(standard)` resamples clusters with replacement. For the observation-weighted
  estimands (`att`, `ame`) it converts the draw to per-cluster multiplicities — provably
  the same statistic as physical duplication — and the draws come from `bsample,
  cluster() idcluster()` on a map frame, i.e. the official prefix's own draw form.
- When a **cell-weighted** estimand (`att_it`, `att_i`, `att_itime`) is requested, the
  replication instead physically resamples and re-identifies duplicated clusters, so a
  twice-drawn cluster contributes its cells twice (author ruling Q2). This leg requires
  `cluster()` to nest `i`; otherwise it is refused.
- `bootstrap(bayes)` draws cluster-level Dirichlet weights and multiplies them into any
  existing weights. No cluster ever leaves the sample, so replication estimability is
  stable. With a cell-weighted estimand it is refused under default safety and
  cautioned under `safety(off)` (author ruling Q2), because reweighting cannot resample
  cells.

**"Everything gets bootstrapped."** Replications capture the *complete named* coefficient
vector — pooled, pre-treatment, every event-study horizon, every subgroup block — stacked
by column name (widths vary across replications when a resample extinguishes a horizon).
From that matrix: standard errors per column; a listwise complete-case covariance for
`e(V)` and for the event-study `e(V_es_*)` handoff (so HonestDiD consumes bootstrap
Sigma); chi-squared Wald joint tests under the same `e(joint_*)` names the analytic path
uses, including the per-subgroup ones and the subgroup-equality chains; pairwise subgroup
contrasts; and normal, percentile and bias-corrected confidence intervals.

---

## 2. "Prove we match the built-ins"

Everything below is measured, reproducible from the cert suites, and re-derived from a
build extracted at the commit under test.

| Claim | Oracle | Result |
|---|---|---|
| Standard-mode replications match Stata's `bootstrap` prefix | official `bootstrap _b, cluster() idcluster() saving(,double)`, same seed | max reldif **2.6e-15** over 15 reps × 4 coefficients |
| Bayes draws match `bayesboot`'s generator | manual rep loop driven by `rwgen bayes` weight variables, same seed | **bitwise identical** (maxreldif 0), both alpha conventions |
| Count path == physical re-id path for `att` | same seed, both legs | **8.0e-16** |
| Weight-injection == physical resampling | reghdfe / wooldid / ppmlhdfe, cluster duplicated + dropped | ≤ **9e-16** across the full estimand vector |
| Joint test == Wald on the bootstrap covariance | hand computation from `e(bootdist)` | **reldif 0** |
| Percentile and BC confidence intervals | `bstat` fed our own replication matrix | **exact** (0) |
| BC median-bias `z0`, plain and `ties` | `bstat` | **exact** (0) |

Two convention subtleties were caught by testing at *odd* replication counts, and are
worth recording because an even-count-only test masks both: Stata's percentile rule is
the classic average-at-integer-*np* rule rather than interpolation, and its BC endpoints
are order statistics, not interpolated.

### The built-in prefix route

Working with Stata's own commands turned out to be *partly* possible, and the blocker was
worth fixing on its own merits. Every Stata replication prefix runs its loop under
`set coeftabresults off`, which silently strips `margins`' `r(table)`; wooldid harvests
that at 16+ sites and died r(503) on every replicate. Two one-line changes (commit A) fix
it, and `bootstrap _b, cluster() idcluster(): wooldid …, fastnose` now runs end to end on
unweighted data and prints the official bootstrap table. It remains a secondary route —
the prefixes reject weighted commands outright, `bayesboot` cannot drive wooldid at all
(it injects `iweight`; wooldid takes `aw`/`pw`), and the prefix silently *discards*
replications whose coefficient vector shrank, which biases event-study horizon SEs toward
full-support draws. The integrated option avoids all three.

---

## 3. Monte Carlo validation

Script: `test/suite/mc/mc_bootse.do` (feature-scoped; not part of `run_full.sh`). Design
A is the primary one — 24 clusters strictly **above** the unit level, onsets assigned at
cluster level, unit effects plus within-unit AR(1) plus a cluster-by-period shock (the
dependence that survives unit fixed effects). Design B is the cluster-at-*i* contrast.
Both are unbiased for the known true att = 2 over 1000 draws. `ratio` is mean reported SE
over the true sampling SD of the estimate, so 1.0 is calibrated; the SD denominator comes
from a 1000-draw plain pass on the same seeds, contributing ~2.2% error.

**Design A (primary, clusters above i):**

| estimator | n | ratio | coverage |
|---|---|---|---|
| `bootstrap(standard)` | 47 | **0.967** | 89.4% |
| `bootstrap(bayes)`, default `equal` | 47 | 0.877 | 85.1% |
| `bootstrap(bayes)`, `size` | 10 | 0.275 | 40.0% |
| analytic clustered VCE | 1000 | 0.869 | 89.8% |
| Stata `bootstrap` prefix | 10 | 0.953 | — |

**Design B (cluster-at-i contrast):** standard 1.004, bayes-equal 0.941, bayes-size
0.389, analytic 0.897.

Four things worth carrying forward:

1. **Standard mode is calibrated, and what it fixes is the analytic SE.** On the same
   data the bootstrap SE runs 18.5% above the analytic one (paired, se 0.035), and the
   analytic one is 0.869 of truth.
2. **The `size`-alpha shortfall reproduces across three geometries and tracks
   1/√(mean cluster size)** — measured 0.312 / 0.415 / 0.372 against predicted 0.289 /
   0.408 / 0.354, consistently just above the asymptotic prediction, which is what a
   finite-cluster correction looks like. This is what the default flip was ruled on.
3. **Bayes with the default `equal` alpha runs 6–10% below standard mode** (paired
   0.909 on A, 0.936 on B) and is mildly anticonservative. It is a usable default, and
   resampling remains the better-calibrated choice where it is available.
4. **Four exact identities held on every single draw** (max absolute difference exactly
   zero, not merely small): point estimates unchanged by bootstrap inference; the
   official prefix's observed estimate equals wooldid's; hidden `e(V_analytic)` equals
   a separate plain run's SE; and `e(boot_N_V)` was 199 of 199 on every arm, so no
   replication was ever silently dropped.

**One flag the MC raised, since settled.** Comparing *independent* seed streams it
measured the internal SE ~6% above the official prefix's and correctly noted that a
Monte Carlo cannot settle an exact per-replication claim. Tested directly at a common
seed, the internal and prefix **reported standard errors agree to 1.3e-09** with
identical point estimates — that residual is the prefix's own float replicate storage,
which is documented behavior unless `saving(, double)` is used.

**One question this raised, since answered by the author** (see open item 5). Coverage sits below
nominal even where the SE is calibrated — design B standard has ratio 1.004 with 81%
coverage. The reported interval is the estimate plus or minus 1.96 bootstrap standard
errors, and at 24 clusters the bootstrap SE is itself noisy across draws, so a
mean-correct SE still undercovers. A t(G−1) critical value (2.069 at G=24, ~5.6% wider)
would recover a few points. That is a critical-value convention choice, not a bootstrap
defect — Stata's own bootstrap prefix is z-based too — and the author ruled it in as the
optional `bootdf()`, leaving the default alone.

---

## 4. Author decisions recorded

| # | Decision | Basis |
|---|---|---|
| Q1 | Replace reported inference rather than append a block | author ruling |
| Q2 | Split engines; bayes + cell-weighted refused under default safety, cautioned under `safety(off)` | author ruling, with the safety-gating modification the author specified |
| Q3 | *(revisited)* `bootalpha(equal)` is the default | author ruling on MC evidence, see below |
| Q4 | `reps()` mandatory, **no** floor | author ruling — low-rep runs are legitimate for timing tests |

**The Q3 revisit.** The initial ruling took `bayesboot`'s own convention (α_c = cluster
size) as the default for compatibility. The MC then measured that convention at **0.35**
of the true sampling standard deviation with **55%** coverage, against **0.94 / 90%** for
the Rubin equal-α convention. The mechanism is exact rather than noise: with α_c = n_c the
Dirichlet's total concentration is the *observation* count, so cluster-weight dispersion
shrinks by roughly the square root of mean cluster size (measured 0.353 against a
predicted 0.363). `bayesboot`'s `id()` semantics are a row-duplication-multiplicity
device — correct when duplicate rows are independent evidence, wrong for panel clusters.
The author flipped the default to `equal`; `size` remains available for exact bayesboot
comparison and is documented with this caveat.

---

## 5. Review and verification record

Two adversarial design reviews of the plan (the second blind to the first), then two
verifiers on the C+D engine layers, then one on E.

- **R1** (plan) — found one blocker: the re-id leg was undefined when `cluster()` does not
  nest `i`, and the observation-level composition was reachable. Gate added. Nine further
  amendments.
- **R2** (plan, blind) — verdict amend-then-implement, six substantive findings, all
  adopted. The load-bearing one: a pairwise-complete covariance over
  availability-varying columns can be **full-rank yet indefinite**, at which point
  `invsym` silently drops a constraint and returns a confident wrong p-value. Everything
  inverted now uses a listwise complete-case covariance with a PSD guard, and a failed
  guard posts missing plus a caution. R2 also found SI-17.
- **v1** (C+D, adversarial) — PASS-with-observations. Rewiring coherence exact, keep/strip
  fidelity across the option surface, sentinel semantics safe, state clean, HonestDiD
  works on the rewired handoff. All observations addressed in commit E.
- **v2** (C+D, do-no-harm) — PASS-with-observations. Off-state **bit-identical on 41/41
  configs at rtol 0**, **9/9 console goldens byte-identical**, zero on-state residue
  across 37 sentinels, five cert suites green.
- **E verifier** — *pending* (see section 7).

---

## 6. Test assets added

| Asset | Purpose |
|---|---|
| `test/test_prefixharden.do` | the commit-A hardening, incl. the naked official prefix end to end |
| `test/test_bootstrap_gates.do` | the option surface and every gate — grading **messages**, not just return codes |
| `test/test_bootstrap_bayes.do` | the bayes engine, incl. the bitwise `rwgen` oracles |
| `test/test_bootstrap_std.do` | standard mode, incl. the official-prefix matching proof |
| `test/test_bootstrap_joint.do` | joint tests, contrasts, CI trio, `ties` — all oracle-checked |
| `test/test_permute2.do` (extended) | SI-17 |
| `test/suite/mc/mc_bootse.do` | the SE-calibration Monte Carlo of section 3 |
| 8 suite configs (`boot` component + `patho_boot`) | dual-build coverage: both engines, the re-id leg, weighted, `ties`, subgroups, joint tests, observation-level, and one pathological-data probe |

---

## 7. Open items

1. **The commit-E verifier pass** has not run. It was launched twice and killed by API
   spend limits both times, and the author elected to finish the remaining work first.
   E's own oracle suite passes, but E has not had the independent adversarial pass that
   C and D received.
2. **Deliberate v1 scope calls, documented rather than fixed:** plots, saved plot data,
   `gteffects` cell standard errors, histogram outputs and adversarial line-test p-values
   are produced from the analytic covariance during estimation, before the rewiring. A
   caution says so at runtime. Re-rendering plots from rewired `e()` is filed as an
   enhancement; `cfx` plots can never be re-rendered from `e()` alone.
3. **`pretrendpower()` is refused with `bootstrap()`** (plan A13): the Roth diagnostics
   compute and post from the analytic event-study covariance inside the worker, before
   any rewiring, so the combination could only have produced a mislabeled number. The
   refusal message points at feeding the rewired `e(b_es_*)`/`e(V_es_*)` to the
   `pretrends` package directly.
4. **`e(df_r)` keeps the analytic degrees of freedom.** User-run `test`/`lincom` after a
   bootstrap run are therefore t/F-metric on the bootstrap covariance and will not
   exactly match the reported z-based p-values. Documented rather than changed —
   `ereturn repost` cannot unset it.
5. **Critical-value convention — RESOLVED** (`4906bd3`). Ruled in by the author as a
   user choice: `bootdf(normal)` (default, unchanged, z-based like the built-in prefix),
   `bootdf(clusters)` for t with G−1 degrees of freedom, or an explicit integer. It moves
   only the reference distribution — p-values, normal-theory intervals, and the joint
   tests, which become F(k, df) — never the standard errors, and it aligns `e(df_r)` so a
   hand-run `test`/`lincom` reproduces the reported p-values.
6. **Coverage of the percentile and bias-corrected intervals was not measured.** Both are
   already stored from the same replications, so measuring them costs no additional
   refits — noted at the foot of the MC script.
7. **Untested corners** carried from the verifiers: `Break` mid-loop, `r907` maxvar
   re-id replications, several caution texts that no probe managed to trigger, weighted
   on-state residue, and the all-columns-extinct listwise corner.
