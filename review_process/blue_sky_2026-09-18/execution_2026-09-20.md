# Execution record — the 2026-09-20 rulings

Companion to `fix_proposal_status.md`. That file records what was *ruled*; this one records
what was *built*, with the measurements that justify each step. Everything below is applied to
the working tree (uncommitted) unless a line says otherwise.

## 1. The removals (author ruling 2026-09-20)

`poisexpresults` and `semielasticity` removed outright; `contreatelasticitytype()` narrowed to
`dydx` and `dyex`, with `eyex`/`eydx` refused by name.

**Why the collapse of the poisson margins gate is numerically inert.** The gate now reads

    if "`poisson'" == "poisson" { local poismarginstext predict(xb) ; local poismarginstextC predict(xb) }

where it used to open with `expression(predict(xb))` and fall back to `predict(xb)` for
`dydx`/`dyex`. The `expression()` branch was **unreachable on the surviving path**:
`semielasticity` was effectively boolean (`eydx` when given, `dydx` otherwise) and BS-6 remapped
`eydx` -> `dydx` and `eyex` -> `dyex` under poisson-without-`poisexpresults` *before* the worker
saw them. So the branch was live only under `poisexpresults`, which set the whole gate empty.

Checked empirically as well as by construction: pinned `c1d8979` against the patched build, in
**separate Stata processes** (per `agent_brief.md`, since `adopath ++` does not reload a loaded
program), over four poisson paths - `att`, `att`+contreat, `att`+contreat+`dyex`, and a poisson
event study - dumping `b` and diagonal `V` at `%21x` and comparing as strings: **bit-identical on
every cell**.

**Build verification.** 38-config battery, one differing line across the whole comparison
(`M wdid_semielas = dydx`, an undocumented macro that posted a constant); `check_validate_sync.py`
43/98 rc 0; SMCL braces balanced; `help wooldid` loads. Working-tree `.ado` is byte-identical to
the battery-verified artifact.

**Author ruling on stored estimates:** replaying a pre-removal stored estimate prints a plain
"Effects" header over numbers computed under the removed operator. Ruled **accept** - stored
estimates from a removed feature are unsupported. No guard added.

## 2. The impvalues IM-geometry screen: caution -> r(459)

Predicate unchanged from the 2026-09-19 correction (counts over `targetflag == 1`, not over all
of `ifstatement`).

**DISCLOSURE - a design choice that was not ruled.** The ruling was "escalate to r(459)
refusal"; what is built is a refusal *with a `safety(off)` escape hatch*. The reason is that the
predicate is **sufficient, not necessary**: a panel where a few units change cluster mid-panel
trips it while its inference is barely affected (measured se ratio 0.968), and that false-positive
class is ordinary data. Most 459 screens are gated this way (`ado:784`, `:918`, `:959`); the
byCohort screen at `:1386` is the unconditional alternative, and it detects an actual
identification failure rather than a sufficient condition for one. **Say the word and it becomes
absolute.**

**A defect found in the screen before it shipped.** The first version keyed on `` `safety' ``, the
EFFECTIVE value. CS-98 (`ado:4856`) rewrites that to `off` whenever `contreatwithin` accompanies
`contreat()`, and nothing refuses `impvalues` alongside `contreatwithin` - so adding
`contreatwithin` to the command silently skipped the refusal. Measured on a crossing geometry:
r(459) without it, **rc 0 with it**. This is exactly the failure the `ado:4944` note predicts. The
screen now reads `safetyinput`, the user's own value, as the neighbouring screens at `:716`,
`:770` and `:1353` already do.

Verified on nine geometries, all matching expectation:

| case | geometry | rc |
|---|---|---|
| A | `cluster(unit)` finer than `i` | **459** |
| B | `cluster(cohid)` == `i` | 0 |
| C | never-treated pooled into one `i`, treated cohorts nested | **0** (the false positive the corrected predicate must not fire on) |
| D | case A under `safety(off)` | 0 |
| E | finer cluster, no `impvalues` | 0 |
| W1 | crossing + `contreat()` | **459** |
| W2 | crossing + `contreat()` + `contreatwithin` | **459** (was 0 before the `safetyinput` fix) |
| W3 | `cluster == i` + `contreat()` + `contreatwithin` | 0 |
| W4 | W2 under `safety(off)` | 0 |

## 3. Suite and docs

- `README.md`: dropped the `semielasticity` and `poisexpresults` bullets; narrowed the
  `contreatelasticitytype` bullet to `dydx`/`dyex`.
- `wooldid.sthlp`: the `impvalues` entry now states the refusal and names `safety(off)`; the
  "What level to cluster at" paragraph names `impvalues` as the one enforced exception.
- `test/suite/configs_ext.csv`: dropped `ext_bin_ext_semielasticity`,
  `ext_cont_ext_semielasticity`, `ext_contract10_poisexp_without_poisson` (their only purpose was
  a retired option NAME, and `coverage_check.py` makes an unresolvable token FATAL - author
  ruling 2026-09-11 - so the token cannot stay and stripping it leaves a duplicate of a base row);
  stripped the flag from `ext_bin_ext_poisson` and `ext_cont_ext_poisson`, keeping them for their
  poisson coverage per the `option_consolidation_2026-09-16.md` precedent (verified: creates no
  duplicate signature); flipped `ext_cont_ext_elasticity_eyex` and
  `ext_cov10_cont_elasticity_eydx` to `rc_new 198`, which they CAN host because the option NAME
  survives even though the VALUES do not.
- `test/suite/configs.csv`: `adm_impv` re-pointed to `rc_new 459`. It was measured to have **3**
  treated `cohid` levels crossed with **48** states (144 `i`-by-cluster groups) - the exact defect
  the screen exists to catch. There is no valid IM geometry on that fixture at a meaningful
  cluster count: `i=fips` is nested in state and onset-constant but demands 17,159 variables
  (r(907) preflight), and onset is not constant within `state` or `county3`. Clean IM coverage
  remains on `real_impv`, `impv`, `trig_imdegen`, `trig_imdeg_part` and two ext rows.
- `test/test_si89.do`: T1/T2 rewritten as the retirement pins. This file is the cert **host** for
  them, because `configs_ext.csv` structurally cannot carry a retired option name. **SI89_ALL_PASS.**
- `test/test_poismargins.do`: rewritten for the unconditional gate; P2/P4/P7 (all
  `poisexpresults`/`semielasticity`/BS-6) removed, P5 re-pointed at `dydx`/`dyex`.
  **POISMARGINS_ALL_PASS** (P1 = .12400474, the plain-poisson value BS-6 documented as the target).
- `test/test_fastcollapse.do`: F4 and F5 reduced from six gates to five, with `lattice()` taking
  the `semielasticity` slot; F7's `poisexpresults` half removed; F11b/F11c reduced from six legs
  to three, which also removes the BS-6 special case that drove all their conditional logic.

`coverage_check.py` returns **COVERAGE_OK** after all of the above, with `impvalues` still at
15 configs / 6 clean.

## 4. Open / carried forward

- **THE BATTERY RE-BLESS IS NOT DONE, and the removal is not landed until it is.** Every config's
  macro dump loses one line (`wdid_semielas`), so the next `test/runbattery.do` produces a diff on
  essentially every config. Per `battery-diff-adjudication.md` and `wooldid-ruling-batches.md`:
  name a new cell CLASS in the comparer - "macro `wdid_semielas` absent" = expected - rather than
  blanketing the matrix in the ledger, run the baseline pristine-twice, then re-bless. Without the
  named class the first post-removal battery run is unadjudicable.
- ~~A clean cert re-run is owed.~~ **DONE.** All three touched cert files re-run against the
  final build (removals + IM refusal + the `safetyinput` fix): **SI89_ALL_PASS**,
  **POISMARGINS_ALL_PASS**, **FASTCOLLAPSE_ALL_PASS**, no rendered failures in any of them.
  F15's crossed-cluster leg still reproduces `nofastcollapse` at 7.25e-15 over 2,878 cells under
  `safety(off)`, so standing the screen down there preserved exactly the property it tests.

- **F15 in `test_fastcollapse.do` still needs re-specification.** Its leg (a) is
  `wooldid y cohid year ttre, att cluster(state) impvalues` on the admin panel - the geometry the
  new screen refuses - and it is a bare `qui wooldid`, so it will hard-abort the file. Intended
  fix: add `safety(off)` to the crossed-cluster calls, since the collapse-exactness property it
  tests is orthogonal to whether IM inference on that geometry is sound. Not yet applied because
  the file was being executed (never edit a running suite do-file).
- The `test_fastcollapse.do` run launched before the refusal was inserted, so it exercised the
  pre-refusal build. It needs a re-run after the F15 edit.
- Pre-existing, not introduced here: `fcol_pois` in `configs.csv` has an unquoted comma in its
  trailing `note`, so it parses ragged. Harmless (note is the last column) but worth tidying.
- Cosmetic: F14's `di` line labels a `customvce()` check as `with_impvalues`.
- `review_process/open_ticket_poisexp_dydx_2026-09-19.md` is **mooted** by the removal - the
  23%-truncated SE it tracks is reachable only through `poisexpresults`. Close with a pointer to
  this record rather than deleting it; its campaign results are the evidence for why the feature
  went.

---

# The strongest result in this slate: the reghdfe-vs-reg paranoia, answered with a number

R6 (`ELAS`) was ruled "ship it" and then **mooted by R2**, which is a quiet-looking ledger cell for
what is actually the most reassuring measurement in the whole review, so it is written out here
rather than left in a table.

The worry the review was built to test is that wooldid's default absorbing engine and its
explicit-dummy engine might not agree — that `reghdfe` with fixed effects in `absorb()` could
return something different from plain `reg` with the same effects written out as dummies. ELAS
looked like evidence that they do not agree: it carried a **13.12% relative gap on `timetrends`**.

The defect turned out to be narrower than that, and specific. `reghdfe`'s `predict` omits the
absorbed part of the fit. Any margins operator that **divides by the prediction** therefore divides
by the wrong thing — and every row in ELAS P1's table is exactly such an arm: a `semielasticity`,
`contreatelasticitytype(eydx)` or `contreatelasticitytype(eyex)` run. The control rows in the same
table, `dydx` and `dyex`, showed error of exactly zero, because those operators never divide by the
prediction at all. When the author ruled to remove the dividing operators outright rather than
repair them, the defect left with them.

Measured on the build that ships, default `reghdfe` against `unconditionalse robust`, worst
relative gap over matched row NAMES rather than row indices (comparing by index is the documented
verification trap here, since the two engines can drop different columns):

| surface | worst relative gap |
|---|---|
| binary, plain | 1.6e-15 |
| `timetrends` | **4.1e-14** (was 13.12% on the pinned build) |
| event study | 2.1e-13 |
| `subgroups()` | 1.6e-15 |
| `contreat()` `dydx` | 6.7e-14 |
| `contreat()` `dyex` | 6.7e-14 |
| `dyex` + `timetrends` | 5.5e-14 |

Seven surviving surfaces, and the worst disagreement anywhere is 2.1e-13 — double-precision
accumulation noise, not a difference in what is being estimated. **The two engines agree.** The one
place they did not was an operator family that has now been removed, and no code shipped for R6
because there is nothing left to repair.
