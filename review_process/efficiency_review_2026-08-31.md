# wooldid efficiency review (2026-08-31)

A large-scale search for efficiency improvements, run as five independent measuring
lenses: phase profiling across a size ladder, an instrumented hunt for discarded and
repeated work, a scaling study along five axes, a Stata-idiom pass over the hot paths,
and a study of the bootstrap/permute replication engines.

**Nothing here is implemented.** This is a report for you to rule on. Every claim below
is measured unless labelled otherwise, and every proposal carries an explicit read on the
risk to reported numbers, since that is the project's cardinal constraint.

---

## The headline

**One line of code is worth 2–10× on every poisson run.** Line 1767 sets
`expression(predict(xb))`, which forces `margins` onto its *numerical* derivative path —
it then issues one interpreted pass over the estimation sample **per coefficient**.
Plain `predict(xb)` reaches the analytic path instead. Two lenses found this
independently and measured it in isolation and end to end:

| config | pristine | one-line change | speedup |
|---|---:|---:|---|
| poisson basic (dta_binary, 960 obs) | 23.3 s | 3.0 s | **7.8×** |
| poisson + event study | 16.2 s | 3.4 s | 4.7× |
| poisson + att_it + subgroups | 138.9 s | 13.5 s | **10.3×** |
| poisson + jointtests | 25.1 s | 3.9 s | 6.4× |
| poisson + unconditionalse | 20.0 s | 2.2 s | 9.2× |
| the isolated margins call | 23.75 s | 0.43 s | 55× |

Point estimates agree to 16+ significant digits (`mreldif` on b: 0 to 3.8e-16) and
variances to 13+ (0 to 7.5e-13).

**It cannot be a blanket swap.** With `semielasticity` (i.e. `eydx`), the pristine build
runs clean but `predict(xb)` **fails outright** — margins aborts on discontinuous
numerical derivatives, surfacing as r(459). The fix is therefore a per-site split: use
`predict(xb)` where the operator is plain `dydx`, keep `expression(predict(xb))` wherever
it is `eydx`/`eyex`/`dyex`. The `contreatelasticitytype()` forms under poisson were not
tested, and the measured eydx failure makes all of them suspect until they are.

One caveat to weigh: those variance differences, while invisible at display precision,
are ~100× larger than the package's own measured run-to-run wobble (1.4e-15 — the
pristine build's `e(V)` is not bitwise reproducible on at least one configuration). So
this is a real change to low-order digits and needs judging against the suite's
comparison tolerance, not waved through as noise.

---

## Where the time actually goes

Three lenses agree on the shape, which is worth internalising before optimising anything:

- **Essentially all wall time is inside `margins`, not in wooldid's own code.** On a rich
  admin-panel config, 44 margins calls account for ~95% of runtime; wooldid's macro and
  matrix loops, the F8 diagnostics and the display together are ~3%. *Optimisation effort
  anywhere outside the margins calls has a ceiling of about 3%.*
- **Of each margins call, 61–84% is the delta-method variance calculation** (99.7% on the
  poisson path). Every reported estimand is a fixed weighted average of cell means, so
  this is the delta method re-deriving, numerically, a quantity that is structurally a
  linear functional of the coefficient vector.
- **Cost is driven by the cell-dummy count K = levels(i) × relative periods, not by N.**
  Measured decisively: 29,067 regression rows at K=29 cost 1.4 s, while 4,152 rows at
  K=2257 cost 86.8 s — 7× fewer rows, 63× more time. Growth is ~K^1.3 to K^1.7; N,
  subgroup count and event-study width are all *linear*. The worst axis is `lattice()`, at
  ~K^1.9, where reghdfe rather than margins becomes dominant.

That last point also explains something you flagged earlier: halving the observations
barely helped the bootstrap loop (57.1 s → 47.2 s) because halving rows leaves the cell
count, and therefore K, unchanged.

---

## Tier 1 — worth doing, low risk

**1. The poisson `predict(xb)` split.** Above. Biggest single win in the report.
*Risk:* point estimates unchanged to 16 digits; variances move in the 13th. Requires the
dydx-only gate and testing of the elasticity forms.

**2. Joint pre-trend tests — REFUSED on investigation (2026-08-31), and it surfaced a bug
ticket.** The structural-identity premise below is false under `interactivecontrols()` and
under `lattice()` + `contreatcontrols(basic)`: 20 of 193 probes diverge, at rc 0 with no
caution, with p moving 0.146 vs 0.615 and 3.6e-3 vs 7.3e-44. Divergence happens exactly
when the reference-cell margin is not a structural zero. The author declined the change
and filed the residue as **SI-20**: the reference period can carry a real estimated effect
in those configurations, which matters not only for the joint tests' null but for the
`e(b_es_*)`/`e(V_es_*)` handoff to HonestDiD and pretrends, which assume a normalised
reference. The original proposal is retained below for the record.

*Original proposal.* Each
`margins … contrast(overjoint)` joint test is a full margins re-computation of a Wald
test whose inputs — the event-study b and V — the run has *already captured* moments
earlier. Measured four independent ways at **23–25% of any `jointtests` run** (15.0 s of
59.5 s on a 45k-obs panel; 40.7 s of 165.5 s at 150k; A/B with and without the option:
94.2 s vs 72.3 s). The replacement costs 0.0078 s per test. Equivalence was verified two
ways, and the `post`/`test`/`restore` pattern it needs is already established in the file
by the E1 fix. *Risk:* low but not zero — the estimands are structurally identical rather
than coincidentally equal, so it wants the usual verification wave.

**3. Subgroup level margins are computed twice.** The `pwcompare(effects)` pass already
returns the level margins in `r(table)`, **bitwise identical** to the ones the code then
recomputes with a separate `over(subgroups)` call (measured `mreldif` exactly 0). Worth
~9% of a subgroups run. *This is the one saving in the report that cannot perturb a
reported digit* — it is reading a number already computed instead of computing it again.
Two ordering constraints to handle.

**4. Stop computing F8 composition diagnostics inside replications.** Bootstrap and
permute replicates consume only `e(b)` (or `_b`/`_se`), yet each one computes the full
event-study composition block and discards it: **8–16% of rep-loop time** on the admin
panel. Found by three lenses independently. *Implementation note:* gate it on an internal
rep-only flag, **not** on `fastnose` — a fastnose gate would silently remove
`e(escomp_*)` for users who call fastnose directly, which is a surface change.

**5. Permute replications compute a whole event study that is thrown away.** 16% of
permute wall time without `permtstat`, **43% with it**. *Measured constraint:* the event
study window must stay in the *regression* — dropping it moves `main_att` by 1–10% — so
only the aggregation may be skipped, which cannot move the reported statistic. The naive
one-line gate does not work; the lens says so explicitly.

**6. `summarystats` uses an e-class `total`**, which destroys `e()` and forces an
`est restore` of the full model — an O(K²) operation. Small today (0.01–0.11 s) but grows
with K. Replaceable with `gstats sum`.

---

## Tier 2 — real, but structural

**Collapse to cells before running margins.** Every margins call evaluates over every row,
but every quantity it computes is a cell-level weighted average. Benchmarked inside the
live model, collapsing to one row per (i, t, targetflag) cell made the same calls **~12×
faster** (event study 1.504 s → 0.119 s; the four-call subgroup loop 5.639 s → 0.461 s).
The algebra holds because `targetflag` does not interact with controls, so its derivative
is constant within a cell.

**Or replace margins at the ~16 event-study sites with direct Mata.** The margin is a
linear functional of b, so the estimate and its variance are two matrix products, versus a
delta-method pass that is superlinear in K.

Both are core rewrites that change summation order and would therefore move low-order
digits on **every** reported number. They are viable only as a versioned rewrite with a
re-baselined suite. I record them because they are what the ~90% is being spent on, not
because I recommend doing them now.

**Documentation, zero risk.** The granular-`i` caution is worded around Stata's variable
limit, but the binding cost is superlinear in levels(i) × periods and bites far below that
limit. Users' single biggest lever is coarsening `i`; the K-curve belongs in the help
text. Same for `lattice()` (~K^1.9) and for the re-id bootstrap leg, which roughly doubles
per-replication cost for `att_it`/`att_i`/`att_itime` — inherent to making those estimands'
bootstrap correct, but worth telling people.

---

## Verified negatives — recorded so they are not hunted again

- **Consolidating the per-subgroup `over(_t)` loop into one `over(subgroups _t)` call does
  not pay.** Two lenses measured it independently: one found ~5%, the other found it
  *slower*. See the adjudication note below.
- The discarded half of every margins result (the `0b.targetflag` columns) is **free**
  under `set emptycells drop`.
- `subpop()` expression evaluation is not a cost, and `if` is worse.
- **The gtools migration is already complete** — there are no Stata-idiom stragglers left
  to find in the hot paths.
- Option re-parsing and validation are **not** the per-replication overhead.
- `bootstrap(bayes)` costing more per replication than `standard` is a sample-size
  mechanism, not an inefficiency.
- `gteffects` is effectively free.
- Pathological small panels show no hot-loop blow-up.
- `set emptycells drop` is worth 2.24× — which also means `emptycellsoverride` is an
  expensive user option, worth a documentation note.

## One disagreement, adjudicated

The profiling lens proposed consolidating the per-subgroup event-study margins into a
single `over(subgroups _t)` call and estimated a saving of (S−1)/S — but marked it
explicitly as *not verified*. Two other lenses **built and measured it**: one found ~5%,
the other found it slower than the status quo. Measurement beats inference, so the
proposal is **refuted** and the subgroup block's cost stands as inherent. I am recording
this rather than quietly dropping it, because the inferred version is superficially
compelling and someone will propose it again.

---

## Method and honest caveats

Five agents, each measuring on real workloads with `timer` instrumentation on scratch
copies of the ado (the repo was never modified), across a ladder from `dta_tiny` through
`dta_binary` to the realistic ~30k-row admin panel, plus purpose-built synthetic panels
for the scaling sweeps. Instrumented builds were validated bit-identical to the repo build
where claims depended on it.

- The agents ran concurrently, so the machine was under load. Within-run phase splits and
  interleaved A/B comparisons are clean; cross-config wall-clock multipliers are
  directional. Where a claim rests on a small difference, medians of ≥3 runs were used.
- Scaling exponents come from single runs per point in some sweeps and are indicative
  rather than precise.
- The poisson elasticity forms (`eyex`/`eydx`/`dyex` under `contreatelasticitytype()`)
  were **not** tested and must be before the tier-1 change ships.
- No fix in this report has been implemented or verified against the assessment suite.
