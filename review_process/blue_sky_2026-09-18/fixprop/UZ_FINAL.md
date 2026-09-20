# `contreatcontrols(untreatedZero)` — FINAL fix specification

**Status: STABLE. Do not re-synthesize.** Three proposal probes converged, three adversarial
reviewers returned SHIP_WITH_CHANGES, and the changes they named are folded in below. The only
open items are two author rulings, marked as such. Further iteration would re-litigate a
question the 368-run adversarial grid already settled.

Closes: `fp-contreat-E` (critical), `sk-untreatedzero-regresspath` (critical),
`FP-CTZERO-ENGINE` (high — with the refuter's inversion: the *explicit* arm is the wrong one),
`SK-UZ-POISSON` (high). These are CRIT-2 and CRIT-3 in `critical.md`.

---

## The defect, in one paragraph

The gate at `ado:1234` builds a pooled baseline dose control (`ccringerhead`) whenever the run is
on an explicit-regressor engine (`unconditionalse`, `customvce()`, `poisson`). That manipulation
presupposes a pooled `c.contreat` column exists in the margins-visible block to partial out.
Under `untreatedZero` there is none — `ado:1314` spells the injection with a single `#` — so the
ringer is not a numerical twin but a real column, exactly collinear with the cell#dose set
whenever the dose is 0 off treated cells, which is *always* when `untreatedZero`'s own
precondition holds. Rank does not rise; the engine drops a cell slope to pay for the ringer;
`margins` never toggles the ringer, so `att` and every post-treatment `ES_att` row silently lose
the pooled slope and `ame` becomes `mean(beta) - beta(dropped cell)`. On the absorbing ppmlhdfe
path the ringer sits inside `absorb()`, the pooled slope is partialled out first, and under a
homogeneous effect the visible cell slopes — hence `ame` — go to **exactly zero**.

So `untreatedZero`, the option whose entire purpose is to control for nothing, was acquiring the
one control it exists to exclude.

---

## What ships

**The primary AND Form C, together. Form C is required, not optional.**

### Primary — prefer **Form B**; Form A is exactly equivalent

Form B: clear the three inserter locals immediately after `ado:1314`, inside the `untreatedZero`
branch. Form A: add the exclusion to the `ado:1234` gate —

```stata
    if strpos("`contreatcontrols'",  "basic") != 0 | ( ( "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" ) & "`contreatcontrols'" != "untreatedZero" ) {
```

**Form A and Form B are byte-identical over 368 config x build runs** in four independent grids —
`0` differing keys, including every configuration the proposers never ran (all six `lattice()`
types plus a premade lattice variable, `aw`/`pw`/`fw`, `if`, `timetrends`,
`interactivecontrols()`, `ccc_absorb()`, `controls()`, `contreatwithin`, `subgroups()`, `fe()`,
`customweightmultiplier()`, `semielasticity`, `contreatelasticitytype()`, `poisexpresults`,
`fastnose`, `nofastcollapse`, `impvalues`, `adversarial`, `bootstrap(standard)`,
`contreatpoly(2/3/4)`). Prefer **B**: it is the one a future edit to the 1234 gate cannot reopen.
If A is chosen, keep its comment pointing at `ado:1313` — the two `== "untreatedZero"` predicates
must never drift apart.

Ship the companion comment block on the DF-1 paragraph either way; it states the rule that makes
the exclusion non-obvious.

### Form C — required (`ado:1314`, the `contreatpoly` loop)

Each polynomial degree becomes its own `#`-joined, cell-interacted term, so no pooled main effect
is ever created; `c.x#c.x` stays legible to `margins` for `dydx`, preserving `ado:1100`'s stated
requirement. At `contreatpoly(1)` the loop body never runs and the macro reduces to **textually
the same string as today**, so poly-1 cannot move — confirmed by string construction *and* by
measurement.

**Why C is required and not optional** — the adversarial evidence is stronger than the proposers':

| arm | A alone | A+C | truth |
|---|---|---|---|
| `contreatpoly(2) customvce(robust)` (d1) | 0.380010603028 | exact | 1.188350230862 |
| `contreatpoly(3)` / `contreatpoly(4)` (d6) | 0.401894023222 / 0.401894022058 (= `basic`'s answer) | 1.314259241731 | 1.314259241731 |
| **`contreatpoly(2) poisson` — the DEFAULT ppmlhdfe engine** (d7) | **4.878692729346** | 0.579867572686 | 0.581271768823 (poly-1) |

The proposal explicitly denied the last row ("the poisson arm's `ame` is repaired by form A
alone"). It is false: **A alone leaves the default poisson path wrong by a factor of eight.**
Form C's derivative is correct at poly 2, 3 and 4 on real curvature (`ame` and `att` exact to
3e-15 on every engine).

---

## Ledger text — four corrections that MUST be made before this is written up

The synthesis's supporting claims were falsified in four places. The fix is sound; these
sentences are not.

1. **Strike "0 return-code changes."** There are **five**, four from the primary and one from
   Form C. All are error -> success, so the direction is benign, but each means a user who saw an
   error yesterday sees a number tomorrow:
   `untreatedZero lattice(all 3) contreatpoly(2) poisson` r(3201) -> rc 0; `untreatedZero
   customvce(ols)` / `customvce(robust)` / `unconditionalse robust` on a cohort-by-period dose
   r(322) -> rc 0; and `contreatpoly(2) customvce(robust)` r(322) -> rc 0 under Form C only.
2. **Strike "Form C moves ONE default-path number."** It moves **six** materially, by O(1), plus
   six more at 1e-15 — all on `contreatpoly(2/3/4)` with precondition-violating or near-violating
   data. On precondition-*satisfying* data C is inert on the default path at every
   `contreatpoly()` tried, confirming the proposers.
3. **Strike "the number it moves is already the as-if-`basic` answer."** It does not reproduce:
   at `contreatpoly(2)` with real dose curvature the pre-C number belongs to **neither** model
   (-0.076915721698, against `basic`'s 0.237498926680). Replace with: *"it is solver-dependent
   and belongs to no model — it changes with `contreatpoly()` on the same data."* This
   **strengthens** the case for C.
4. **Qualify "0 fit invariants move."** `e(ll)`, `e(rss)` and `e(rmse)` do move — 1 ulp on the
   poisson arms under the primary, materially on the six default-path poly arms under C
   (`e(ll)` 14059.16 -> 13616.99). `e(N)`, `e(df_r)` and `e(r2)` do not.

## Pre-existing residual to record in the same ledger entry (DEF-3)

Not caused by this fix, but it is the reason the fix does not achieve full cross-engine identity,
and it must not be discovered later as a regression. When the dose is **constant inside a treated
cell** (e.g. every unit in the 2010 cohort got the same first-year subsidy),
`customvce()`/`unconditionalse` return `att` exactly right and `ame` short by
`beta(dropped cell) x weight(cell)` — the same *truth-minus-one-dropped-cell-slope* signature as
CRIT-2, arriving by a route the fix does not close — while the default engine refuses with
`r(322)`. It is **not silent**: every such arm prints wooldid's existing generic multicollinearity
caution. It is nonetheless easy to dismiss, because the `att` beside it is exactly right and
`e(r2) = 1`.

Optional follow-on, outside this fix's scope: extend the drop check so `r(322)` fires when a
treated cell's `contreat` *slope* is dropped, not only when the cell *dummy* is.

## Verification design — one leg is unsafe as proposed (DEF-6)

**Do not adopt the proposers' cross-engine-identity suite leg unguarded.** "Under `untreatedZero`,
every engine's point-estimate vector must agree with the default `reghdfe` path" is **flaky**: on
thin cells (8 obs/cell) the post-fix `ame` still differs between engines by 3-9% for
`lattice(byCohort 2)`, `lattice(baseTreatedPost 3)` and a premade lattice, because the
per-(cell x lattice) slope is not identified there. At 40 obs/cell every lattice arm agrees to the
last bit and equals the truth. Gate the leg on a minimum obs per (cell x lattice) — 40 is
measured-adequate — or restrict it to lattice-free runs, or use an exact-fit oracle instead.

## Two author rulings

1. **DEF-2 (high).** Post-fix, on an ordinary cohort-by-period intensity dose (minimum wage,
   subsidy rate, capacity cap), the **default** engine refuses with the informative `r(322)` while
   `customvce()`/`unconditionalse` now return `att` = truth with `ame` missing and a generic
   caution. Before the fix all seven engines refused. The new behaviour is defensible — the `att`
   is correct, the missing `ame` is honest — but it is new, it is engine-dependent in the
   *opposite* direction from CRIT-3, and it is outside the proposal's stated envelope. A user
   refused by the default engine who adds `customvce(robust)` on a referee's advice now gets a
   table where they used to get an error.
2. **The six moved default-path numbers** under Form C (correction 2 above) need a re-bless.

## Untested
`makeplots`/`saveplots` output files; `permute()` (refused by design under `untreatedZero`);
`boottstat`/`bootdisplay`; `honestdid`/`pretrends` consumers; `contreatpoly() > 4`; panels large
enough for the cell collapse to matter; `oldsyntax`; multi-way `cluster()`; designs with more
than six cohorts.
