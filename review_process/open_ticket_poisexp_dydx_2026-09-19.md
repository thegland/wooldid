# TICKET — `poisson poisexpresults` returns a truncated standard error on the plain ATT

**Opened 2026-09-19 on the author's instruction:** *"The problem in the dydx is a subject of immense
concern. Let's make a separate ticket for this and make measuring the scale of this effect across
scenarios a big priority."*

**Status: CLOSED 2026-09-20 - MOOTED BY REMOVAL, not by repair.**
**Headline as measured: on the then-shipped build, a realistic count panel reported a standard
error 23% too small, with nothing in the output saying so, and the distortion did not shrink
with sample size.**

> **Close-out.** The author ruled on 2026-09-20 that `poisexpresults` be removed outright rather
> than repaired, together with `semielasticity` and the `eyex`/`eydx` forms of
> `contreatelasticitytype()`: *"Efforts to bail out these features I think would be too much,
> especially when some of the issues are sort of upstream-ish."* The truncated standard error
> this ticket tracks is reachable **only** through `poisexpresults`, so it is unreachable on the
> current build - every margins call wooldid now makes under poisson is a plain `dydx` or `dyex`
> of the linear index, where the Poisson score equations pin `mean(mu-hat) = mean(y)` and the
> Jacobian row on the absorbed part is zero.
>
> **This file is kept, not deleted.** Its campaign results are the evidence base for the removal:
> the distortion's failure to shrink in N or T, the refutation of the analytic-SE-undercoverage
> hypothesis, and the measurement that exonerates ppmlhdfe's own `margins` claim for its modal
> use. The execution record is
> `review_process/blue_sky_2026-09-18/execution_2026-09-20.md`.

## What it is

Under `poisson poisexpresults` **with no elasticity option typed**, the reported standard error on
the ordinary ATT does not account for the sampling variation of the absorbed fixed effects.

This is NOT an elasticity defect. It reaches users who never asked for an elasticity, and no ruling
about `semielasticity` touches it.

## Why it happens

`ppmlhdfe_p.ado` computes `mu = exp(xb + d)`, reading `e(d)` - the stored variable `_ppmlhdfe_d`.
When `margins` perturbs the coefficient vector to build its Jacobian, `xb` moves and **`d` does
not**. So the Jacobian carries no rows for the absorbed effects.

Standalone mechanism proof: `ppmlhdfe y i.cellg#i.D, absorb(i.coh i.t) vce(robust) d` followed by
`margins, dydx(D) subpop(if D==1)` returns an `r(Jacobian)` with **53 columns - exactly ppmlhdfe's
`e(b)`, and no FE columns**. SE 0.363297 against the explicit-dummy twin's 0.296566.

The **point estimate is correct and always was** - ppmlhdfe's predict already includes `d`. Only the
variance is affected. It is therefore **not** repairable by the wrapper that fixes the linear-scale
elasticity path.

## Measured so far - one design family only

8 cohorts x 12 units x 10 periods, `mu = exp(1.0 + .1*coh + .05*t +/- .40*D)`, `y ~ Poisson(mu)`,
`mu` in 2.5-11, row `main_att`, `nofastcollapse` on both engines:

| spec | SE ratio, positive effect | SE ratio, negative effect |
|---|---:|---:|
| `poisson` | 1.00000 | 1.00000 |
| `poisson semielasticity` | 1.00000 | 1.00000 |
| `poisson poisexpresults semielasticity` | 1.00000 | 1.00000 |
| **`poisson poisexpresults`** | **1.22501** | **0.80824** |

The first three are clean for a reason: without `poisexpresults` wooldid remaps `eydx -> dydx` and
`eyex -> dyex` and hands margins `predict(xb)`, so every reported quantity is a difference on the
linear index and `d` cancels. With `poisexpresults` **and** `semielasticity`,
`ln(exp(xb + d)) = xb + d` and `d` cancels again.

## The mechanism, stated generally

Any margin whose Jacobian has a **nonzero row on the absorbed effects** has a truncated standard
error when those effects live outside `e(b)`. On the linear scale a difference of predictions
cancels `d`, so `d*` operators are clean and `ey*` are not. **On the exponentiated scale the roles
invert**, because `exp(xb1+d) - exp(xb0+d) = e^d * (...)` does not cancel while `ln(exp(xb+d))`
does.

What is omitted is `Var(F-hat) + 2*Cov(theta-hat, F-hat)`. The variance term is always positive; the
covariance term carries the sign of the effect and dominates at moderate effect sizes - which is why
the ratio sits above 1 for a negative effect and below 1 for a positive one.

## What the campaign must establish

1. **Scale across scenarios** - panel shape, N, T, effect size and sign, FE dispersion, cluster
   geometry, weights, subgroups, event-study horizons, `att_i`/`att_it`/`att_itime`.
2. **Asymptotics.** Does the distortion vanish in expectation on realistic data? The decisive
   question is whether it shrinks with **T** (each fixed effect better estimated) while staying put
   as **N** grows with T fixed - the incidental-parameters shape. If it does not vanish in T, it is
   not a small-sample artefact.
3. **Does it reach `e(V)` off-diagonals**, and therefore the joint tests and simultaneous bands?
4. **Which reported rows** carry it, and whether null rows are inert as the linear-scale analysis
   found (there, every pre-period row measured 0.9985 because the omitted Jacobian row vanishes
   where there is no effect).

## Repair options, pending the measurement

The only repair is to put the fixed effects in `e(b)` - the explicit-dummy refit that
`customvce()`/`unconditionalse` already run, where the delta method is complete by construction.
`bootstrap()` and `permute()` never touch the analytic Jacobian, so they are unaffected.

## Provenance
Found while diagnosing the linear-scale elasticity standard error
(`review_process/blue_sky_2026-09-18/`). Evidence and do-files in
`scratchpad/elasse/` (`R1.do`-`R4.do`, `NOTES.md`).


---

# CAMPAIGN RESULTS, 2026-09-20

Six agents, 0 errors. Adversarially verified: **SOUND_WITH_CORRECTIONS** (corrections folded in
below). Evidence in `scratchpad/elasNT/`, `scratchpad/elasse/`, `scratchpad/econexp/`.

## 1. It does not vanish with sample size, in any dimension

The predicted incidental-parameters shape - shrink in T, not in N - is **wrong**. It shrinks in
neither, because `Var(F_hat)` and `Var(theta_hat)` fall at the same rate.

| rows | se_trunc | se_full | ratio |
|---|---|---|---|
| 1,536 | 0.01650 | 0.01964 | 0.8508 |
| 6,144 | 0.00861 | 0.01026 | 0.8379 |
| 24,576 | 0.00435 | 0.00516 | 0.8429 |

Both SEs fall by 3.79x for 16x the rows; the ratio is flat (3-seed range at the largest cell
0.84229-0.84329). Sixteen times more incidental parameters, each estimated no better: ratio flat
(0.828 -> 0.839 over G = 4..64).

**Lengthening the panel the ordinary way makes it WORSE**: linear 0.831 -> 0.797 and exponentiated
1.199 -> 1.346 as T goes 8 -> 50 with the pre-window fixed. Only an unusual long-pre/short-post
shape improves, and it **plateaus at ~0.87** - 400 periods, 396 of them pre-treatment, 19,200 rows,
still 13% wrong.

**Why it plateaus.** The margin loads on *post-period* time effects, which are identified only from
the control cohorts' rows at that period - 32 rows here - and get no more precise however long the
panel runs. Decomposed: absorbing the cohort FE alone gives 0.981, the time FE alone 0.895, both
0.8755 (~= the product).

`se_trunc/se_full` is therefore a **shape constant**, set by the effect size relative to the
baseline level and by the treated cohorts' post-treatment share. Measured rule of thumb,
**corrected upward by the verifier**: `|1 - ratio| ~= w * rho` with **w in 0.45-1.04** across
designs (the campaign's own 0.42 understates it by about 2x).

## 2. The realistic case, on the SHIPPED build

`wooldid nevent hid year ttre, att cluster(hid) poisson poisexpresults` - a 120-hospital count panel
with a 28.6% reduction in adverse events. No elasticity option typed.

| | value |
|---|---|
| reported se | 0.34315 |
| complete delta method | 0.44629 |
| ratio | **0.7689** (Jacobian-clean normalizer exactly 1.0000, so all of it is truncation) |
| published 95% CI | [-4.279, -2.934] |
| correct 95% CI | [-4.481, -2.732] |
| actual coverage of that "95%" | **~87%** |

**The dangerous direction on the count scale is a NEGATIVE effect** - counts falling - which is the
inverse of the linear-scale rule: 0.84 at -30%, 0.56 at -70%, 0.45 at -70% with a single treated
cohort (correct SE 2.2x reported, coverage 62%; worst constructed 0.36, coverage ~51%). A positive
effect is merely conservative and **saturates** near 1.55-1.6. Outcome scale is irrelevant - mean
counts from 0.60 to 167 move the ratio only 0.800 -> 0.777.

## 3. Worst surface: the across-subgroup equality test

Contrast SE ratio **0.665-0.732** (t inflated ~1.4x); the reported p crossed 0.01/0.05 in **3 of 7**
linear draws - 0.00092 reported against 0.0817 from the complete engine. Off-diagonals are
distorted *more* than diagonals (covariance ratios 0.548-2.477 against variance ratios 0.597-1.426).

Note the plumbing: `e(V)` is diagonal except for the within-family event-study blocks spliced in by
`wooldid_vsplice`, so no off-diagonal touches a `main_*` or `pretreatment_*` row.

## 4. CORRECTION - null rows are NOT inert

The earlier finding that pre-treatment rows measure 0.9985 and are therefore inert **does not
generalise**. It was the `Var(F_hat)` floor on a draw whose realized pre effects happened to be ~0.
On a design with a genuine pre-trend the pre rows carry the truncation exactly like post rows, in
proportion to their own realized effect: ES_att_pre1 0.9581, pre2 0.9746, pre3 0.9741,
pretreatment_att 0.9603.

## 5. This is a general Stata exposure, not a wooldid invention

`margins` assumes the prediction is a function of `e(b)` and the data. Any `predict` program that
reads an *estimated* quantity out of a dataset variable silently truncates the Jacobian, and
`margins` cannot detect it.

Measured with **no wooldid in the session**: `ppmlhdfe y i.x z, absorb(i.id i.t) vce(robust) d` +
`margins r.x` against the explicit-dummy `poisson` twin - points bit-equal (rel 2e-10), SEs 0.5435
vs 0.5758, **ratio 0.9439**; negative twin 0.9201. That is 5.6% and 8.0% on a design nobody tuned.

**Verifier's scope correction:** this requires a covariate besides the contrasted one. The bare
two-way-FE Poisson DiD - `ppmlhdfe y i.treat, absorb(i.unit i.t)` then `margins r.treat`, no
controls - does **not** show it.

## 6. Corrections the verifier made to the campaign's own claims

- The rule-of-thumb coefficient understates by ~2x (above).
- "The correct SE is invariant to the effect size" is **design-specific**, not a property of
  `se_full`; reproduced exactly on one design and broken on two others.
- "Fine `i` is not catastrophically worse" was an artefact of comparing different designs. On
  identical data: **0.9716** (coarse `i`, no county FE) -> **0.9130** (coarse `i` + `fe(i.cty)`) ->
  **0.8741** (`i` = county). Absorption depth carries the first step; `i` granularity a further ~1.45x.
- Cross-engine point agreement was quoted as 1e-8; the worst measured is **8.2e-6**.
- The worst-realistic 0.45 is an 8-seed mean over range [0.382, 0.570], not a design constant.
- **Contradicting the brief's own framing:** the zero-null test is correctly sized (0.0463
  two-sided). Only a correct p in [0.0500, 0.0528] is pushed below .05, because the omitted
  covariance is proportional to the *realized* effect. "Too-small SE produces false significance" is
  wrong at the null; the damage is to intervals and to tests of non-zero nulls.

---

# WHY ppmlhdfe's DOCUMENTATION IS NOT WRONG (2026-09-20)

ppmlhdfe advertises `margins` compatibility and posts `e(marginsok) = "default"`, blessing the
exponentiated prediction. That is **correct for the case its users mostly occupy**, and the
distinction is exact rather than approximate.

**The protected case: any margin the Poisson first-order conditions pin.**

With a full set of fixed effects the score equations force `sum(y - mu_hat) = 0` within every FE
group, hence `mean(mu_hat) = mean(y)` **exactly**. So for a continuous regressor

    AME = mean(d mu / d z) = mean(beta_z * mu_hat) = beta_z * mean(y)

and `mean(y)` is a sample mean with no sampling variation of its own. The true Jacobian row on `d`
is **zero at the optimum**, so nothing is omitted and the truncated variance is exactly right.

**The exposed case: any margin requiring a counterfactual prediction no FOC pins.** Setting `D = 0`
for treated observations produces `exp(xb_0 + d)`, which is not the fitted value of any observation.
Nothing forces it to match anything in the data, and `d`'s estimation error enters in full.

**Measured on one fit, one dataset** (`scratchpad/foc/foc.do`), 60 units x 12 periods, continuous
`z` and dummy `D`, `ppmlhdfe y z i.D, absorb(i.id i.t) vce(robust) d` against the explicit-dummy
`poisson` twin:

| quantity | point agreement | SE ratio vs twin |
|---|---|---|
| `mean(mu_hat)` vs `mean(y)` | rel -1.9e-14 | - |
| `margins, dydx(z)` vs `beta_z * mean(y)` | rel 2.6e-9 | - |
| **`margins, dydx(z)`** (continuous) | exact | **1.000000** |
| **`margins r.D, subpop(if D==1)`** (contrast) | rel 1.3e-10 | **1.2608** |

Same fit, same data: one margin exact, the other 26% wrong. The only difference is whether a first
order condition pins the quantity.

**Consequence for wooldid.** A DiD ATT is a counterfactual contrast *by construction*, so wooldid is
in the exposed case on every `poisexpresults` run. This is not a defect in ppmlhdfe and any write-up
should say so plainly: their claim and this finding are statements about different margins.

**Revise the general-exposure claim accordingly.** The earlier line - "every applied paper running
`ppmlhdfe ..., absorb()` then `margins, dydx()` is reporting a truncated SE" - is **wrong** and must
not ship. The continuous-regressor AME, which is the modal postestimation use after a gravity or
count regression, is exactly correct.

## Separate, unresolved, and NOT this defect

A **bare** two-way-FE Poisson DiD with no other covariate - `ppmlhdfe y i.treat, absorb(i.unit i.t)
vce(robust) d` then `margins r.treat` - returns a silently **wrong point**: measured 1.594 vs 3.613
(-56%), 4.678 vs 16.877 (-72%), 1.714 vs 4.435 (-61%) across three panels. It matches neither the
correct value nor the no-`d` nor the mean-`d` value. Adding any second covariate makes it exact.
`d(name)`, `noestimcheck`, `predict(mu)` and `[pw=1]` do not fix it. The only visible tell is the
table header reading "Contrasts of adjusted predictions" rather than "Contrasts of predictive
margins". Cause unresolved; it does not affect the wooldid numbers in this ticket, which all carry
additional covariates.

---

# RULING, 2026-09-20: REMOVE

The author has ruled that `poisexpresults` and the failing elasticity cases be **removed** rather
than repaired: *"Efforts to bail out these features I think would be too much, especially when some
of the issues are sort of upstream-ish."* Explicitly: do **not** criticise ppmlhdfe in the help
files.

## What goes

| surface | why |
|---|---|
| **`poisexpresults`** (whole option) | its `dydx` cell is the truncated-Jacobian defect; the option is the only route to the exponentiated scale |
| **`semielasticity` on the NON-poisson path** | eydx on the linear index - measured 0.79 to 1.53 |
| **`contreatelasticitytype(eydx)` / `(eyex)` on the NON-poisson path** | same mechanism - 1.003 to 1.052 and 0.995 to 1.028 |

## What stays, because it measured clean

| surface | evidence |
|---|---|
| plain binary `dydx` (the default) | SE bit-identical cross-engine at every effect size |
| `contreatelasticitytype(dydx)` / `(dyex)` | 1.00000 exact on both designs |
| the elasticity options **under `poisson`** | BS-6 (author ruling 2026-09-07) already remaps them to xb-scale `dydx`/`dyex`; they relabel the header and do not change the estimate |

## The end state is the reason this is clean rather than merely safe

BS-6's own comment records that under `poisson` without `poisexpresults` **"no ey* operator reaches
margins any more."** Removing `poisexpresults` removes the one path where ey* still did, and
removing the linear-path ey* cases removes the other. **After both removals no `ey*` operator ever
reaches `margins` anywhere in wooldid**, and every remaining margins call is a plain `dydx`/`dyex`
on the linear index - where the difference cancels `d` identically and the Jacobian row is exactly
zero rather than merely absent.

That eliminates the defect class at the root rather than screening for it.

## Scoping consequence to confirm when implementing

`semielasticity` survives as a **poisson-only** option. On the poisson path it is clean, delivers
the quantity its header promises, and is pure relabelling of the xb-scale operator. Its removal
target is the non-poisson path only. Same for `contreatelasticitytype(eydx|eyex)`.

## Implementation notes

- `poisexpresults` is currently implemented by ABSENCE: the gate at `ado:1894` populates
  `poismarginstext`/`poismarginstextC` with `predict(xb)` / `expression(predict(xb))` only when
  `poisexpresults` is empty. Removing the option means that gate's condition is always true.
- The BS-6 remap block (`ado:5081-5090`) and its caution text reference `poisexpresults` by name and
  will need rewriting, not just deleting - its caution currently ends "Specify poisexpresults to
  obtain the same quantities from the exponentiated prediction instead."
- `sthlp:246-248`, `301`, `372`, `382` and the option entry at `sthlp:715-716` all reference
  `poisexpresults`.
- `e(gteffects)` carries a note saying its cells stay on the xb scale "whether or not
  poisexpresults is specified" - that note becomes vacuous.
- Decide what a user who types the removed options should get: an error naming the removal, or
  silent acceptance. An error is consistent with the package's other retired-option handling.
