# uz-blast — full blast radius of the `ado:1234` ringer insertion, and fix proposals

Agent angle: MAP THE FULL BLAST RADIUS. Everything below is measured on exact-fit designs
(R2 = 1, no noise) where the truth is known in closed form, on the pinned build (sha c1d8979).
Workdir `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uz-blast`; every probe do-file and log is there;
`NOTES.md` carries the running record.

---

## 0. HEADLINE — the defect is larger than CRIT-2/CRIT-3 state, in three ways

1. **`contreatcontrols(basic)` is silently wrong on the DEFAULT reghdfe path** whenever the dose has
   no variation outside the treated cohort-periods — the exact data shape `untreatedZero` is
   documented for. Same wrong number, same size, no engine opt-in required.
2. **`byTime` / `byTreated` / `byTtre` are wrong on every engine** on such data, by a *different*
   mechanism (the `ccringertail`, which no gate guards), and are therefore **not fixed by any change
   to the `ado:1234` gate**. `byCohort` is the only member of that family already protected
   (`wooldid_amescreen mode(bycoh)` refuses r(459)); the other three have no screen.
3. **The reported number is a function of the arbitrary FACTOR ORDERING, not of the data.** Renaming
   the cohort identifier — same partition, same rows, `i2 = 9 - i` — moves the reported `ame` from
   **0.373750000** to **-0.106250000** (truth 0.973750000, a sign flip). Both shifts equal minus the
   slope of whichever treated cell happens to sort last.

And two smaller ones:

4. **Every pre-treatment and event-study pre row is EXACTLY identical between the correct and the
   broken run** (all zero on these designs). A pre-trend / placebo check cannot detect this.
5. **`interactivecontrols()` shows the same engine-dependent divergence** on a control that is zero
   off treated cells: `att` 1.754034491 (default, = truth) vs 0.890310304 (`customvce(robust)`).

---

6. A **null ringer column is harmless**: a group whose dose is identically zero everywhere (a
   never-treated-only group) contributes an all-zero column that drops without displacing anything.
   That is the difference between a screen that works and a screen that refuses good runs — see
   FIX-2.

## 1. THE UNIFIED MECHANISM (one sentence, then the invariant)

The ringer set is placed **before** the margins-visible treatment terms precisely so Stata's sweep
drops a margins-visible column; that is correct only when the column it displaces is the
**redundant pooled base-dose carrier**. It misfires whenever a ringer column lies in the span of
the margins-visible **cell#dose** columns, because then the sweep drops a cell slope — or the
`c.contreat` carrier the AME needs — and every reported quantity comes back shifted.

**THE INVARIANT A FIX MUST RESTORE (and the cheapest regression test):**
> for a given `(contreatcontrols, lattice, contreatpoly)` specification and a given dataset, the
> point estimates must not depend on which ESTIMATION ENGINE runs — reghdfe, `reg` under
> `unconditionalse`/`customvce()`, or ppmlhdfe under `poisson`.

The pinned build violates engine-invariance on exactly the `untreatedZero` rows. It *satisfies*
engine-invariance and is still wrong on the `basic`/`byX` rows, which is why engine-invariance is
necessary but not sufficient — the closed-form truth is the second leg.

Two triggers:

- **(i) head-induced, `untreatedZero`.** `contreatinj` at `ado:1314` uses a single `#`, so the
  margins-visible set `i._i#i._t#i.targetflag#c.contreat` is **saturated over the whole sample** —
  the `(_i=0,_t=0,tf=0)` cell carries its own dose slope. The pooled ringer `c.(ringer)` is in that
  span **by construction, on any data**, not only when the dose is zero off treatment. Design N
  (dose nonzero everywhere) proves it: 1.343750000 -> 0.373750000.
- **(ii) span-induced, `basic` / `byX`.** These use `##`, so the pooled column is redundant only
  when the reference pool (`_t == 0` rows, *within the byX group*) carries dose VARIATION. With no
  variation there, the ringer again lies in the cell span, and the sweep takes out the
  `c.contreat` carrier (head) or a cell slope (tail). Fires on **all engines**.

**`basic` is broken by the HEAD, not by a rank deficiency of its own spec.** Under `##` the AME is
`b_c.contreat + b_cell`, which is identified whichever column drops — *unless* the ringer forces
`b_c.contreat` to drop, which is exactly what the head does. Proven by patchB (head never built):
`dz`/`basic` returns the exact truth 0.973750000 on all seven engines. **But patchB is not a fix**:
on ordinary data it moves `dn`/`basic` from 0.973750000 to 1.343750000 on every engine — the
baseline slope g0 = 0.37 leaks straight into the AME. The head IS load-bearing (this confirms the
DF-1 comment at `ado:1258-1299`).

**`byX` is broken by the TAIL, which no gate touches.** patchB leaves `dz`/`byTime` at 0.200000000
on the default path (the tail `i.t#c.(ringer)` is still there and still collinear).

---

## 2. THE BLAST-RADIUS TABLE — the maintainer's regression checklist

Design Z = dose EXACTLY 0 outside treated cohort-periods (what `untreatedZero` documents).
Design N = dose defined everywhere, pooled baseline slope g0 = 0.37 (ordinary continuous treatment).
Both exact-fit, 8 cohorts x 12 units x 10 periods = 960 obs, 16 treated cells,
**true att-weighted ame = 0.973750000000**, true att (dz) = 1.754034490530.
Engines: E1 `robust` (reghdfe), E2 `cluster(unit)`, E3 `customvce(robust)`, E4
`customvce(cluster unit)`, E5 `unconditionalse robust`, E6 `ccnoabsorb robust`,
E7 `ccnoabsorb customvce(robust)`; plus P1 `poisson robust`, P2 `poisson customvce(robust)` on the
log-linear twins (true semi-elasticity 0.347500000000).
Raw: `gridA_pinned.csv`, `gridA_patchA.csv`, `gridA_patchB.csv`, `pois_pinned.csv`, `pois_A.csv`.

| # | contreatcontrols | data | E1 E2 E6 (reghdfe) | E3 E4 E5 E7 (reg) | P1 P2 (ppmlhdfe) | verdict | fixed by |
|---|---|---|---|---|---|---|---|
| 1 | `untreatedZero` | Z | **0.973750000 OK** | **0.373750000 (-0.600)** | **0.127500000 (-0.220)** | CONSEQUENTIAL, head | FIX-1 |
| 2 | `untreatedZero` | N | 1.343750000 (its own spec) | **0.373750000 (-0.970)** | **0.127500000 (-0.370)** | CONSEQUENTIAL, head, **on ordinary data too** | FIX-1 |
| 3 | `basic` | Z | **0.373750000 (-0.600)** | **0.373750000 (-0.600)** | **0.127500000 (-0.220)** | CONSEQUENTIAL, head, **default path** | FIX-2 (or FIX-2alt) |
| 4 | `basic` | N | 0.973750000 OK | 0.973750000 OK | 0.347500000 OK | REDUNDANT — this is the case the head exists for | — (must not move) |
| 5 | `basic byTime`, `byTime` | Z | **0.200000000 (-0.774)** | **0.200000000 (-0.774)** | **0.087500000 (-0.260)** | CONSEQUENTIAL, tail | FIX-2 only |
| 6 | `basic byTreated`, `byTreated` | Z | **0.373750000 (-0.600)** | **0.373750000 (-0.600)** | **0.127500000 (-0.220)** | CONSEQUENTIAL, tail | FIX-2 only |
| 7 | `basic byTtre`, `byTtre` | Z | **-0.046250000 (-1.020)** | **-0.046250000 (-1.020)** | **-0.085000000 (-0.433, SIGN FLIP)** | CONSEQUENTIAL, tail | FIX-2 only |
| 8 | `byCohort`, `basic byCohort`, `basic byTime byCohort` | Z | r(459) | r(459) | r(459) | **already protected** (`amescreen mode bycoh`) | — |
| 9 | every `basic`/`byX` spelling | N | 0.973750000 OK | 0.973750000 OK | 0.347500000 OK | REDUNDANT | — (must not move) |

Modifier cross (design Z, `gridB_pinned.csv`):

| modifier | effect on the defect |
|---|---|
| `lattice(all 2)` | **makes it worse.** `untreatedZero` + lattice + `customvce`: `ame` still 0.373750000 AND `att` 0.763210 (vs 1.754034) — a second drop. `basic` + `lattice(all 2)` on Z errors **r(322)** (protected by accident). Only `all 2` tested. |
| `contreatpoly(2)` | changes the geometry — see FIX-6, held as UNCERTAIN |
| `subgroups(sg)` | no effect on the defect (identical numbers) |
| `contreatwithin` | **no effect** — does not suppress the ringer, even though it declares there is no untreated comparison pool |
| `timetrends` | no effect |
| `[aw=w]` | **no effect on the defect** (`gridW_pinned.csv`): dz `untreatedZero` 0.973750000 default / 0.373750000 customvce, dz `basic` 0.373750000 both, dn all correct — identical to the unweighted grid |
| `bootstrap()` | cannot inherit the `reg` door: `ado:6410/6414` refuse `bootstrap()` + `unconditionalse`/`customvce()`. It CAN inherit the **poisson** door, so `bootstrap(...) poisson contreat() contreatcontrols(untreatedZero)` bootstraps a distribution centred on the wrong point. |

**Per-estimand damage** (design Z, `untreatedZero`, default vs `customvce(robust)`, probe03):

| family | shift |
|---|---|
| `main:ame`, `ame_it`, `ame_i`, `ame_itime` | **exactly -0.600000000 = -beta(dropped cell)**, all four |
| every `ES_ame` post horizon | **exactly -0.600000000** (post0..3: 0.750->0.150, 0.930->0.330, 1.020->0.420, 1.120->0.520) |
| `main:att` and `att_it` | 1.754034491 -> 1.083210355 (**-38%**) |
| `att_i`, `att_itime` | 1.526760312 -> 0.867407530 |
| `ES_att` post0..3 | -0.661824951 / -0.629990583 / -0.697763924 / -0.702665094 |
| **every pretreatment row and every ES pre row** | **0.000000000 in BOTH runs — bit identical** |

**The `att`-side identity is verified, not inferred** (probe08): the shift on every `att` row equals
`-beta(dropped cell) x mean(dose)` over the relevant observations, to 9 decimals —
0.6 x 1.118040225087 = **0.670824135052** against a measured main-`att` shift of 0.670824135, and
0.6 x the per-horizon mean dose gives **0.661824951443 / 0.629990582880 / 0.697763924497 /
0.702665093663** against measured 0.661824951 / 0.629990583 / 0.697763924 / 0.702665094. A
maintainer can use this identity as an assertion.

**Label dependence** (probe05, `RELAB`): dz with the cohort variable renamed `i2 = 9 - i` (identical
partition, identical rows). `untreatedZero`+`robust` 0.973750000 (unchanged, correct);
`untreatedZero`+`customvce(robust)` **-0.106250000** (was 0.373750000); `basic`+`robust`
**-0.106250000** (was 0.373750000). Both equal truth minus 1.08, the slope of the cell that now
sorts last. The error is a pure artifact of factor ordering.

**Predicate validation across three data shapes** (probe05 `MIX`, truth 0.973750000 everywhere):

| design | `basic` | `byTime` | `byTreated` |
|---|---|---|---|
| Z: dose 0 off all treated cells | **0.373750** wrong | **0.200000** wrong | **0.373750** wrong |
| P: dose 0 off treatment, except one never-treated cohort carries a dose | 0.973750 OK | 0.973750 OK | **0.373750** wrong |
| Q: dose nonzero for every treated unit incl. its pre-periods, 0 for never-treated | 0.973750 OK | **0.823750** wrong | 0.973750 OK |

Nine cells, and the outcome is predicted **exactly** by "does the dose vary over the `_t == 0`
reference rows **within the grouping that control type uses**": whole sample for `basic`, calendar
period for `byTime`, ever-treated indicator for `byTreated`. Every zero-variation group gives a
wrong number; every positive-variation group gives the exact truth. That correspondence is what
FIX-2 is built on.

---

## 3. FIX PROPOSALS

### FIX-1 — stop inserting the ringer under `untreatedZero` (CLASS i). **Ship this.**

Two equivalent sites; I recommend **1b**, with **1a** as the minimal diff.

**FIX-1b (recommended) — enforce the invariant where the semantics are declared.**
At `wooldidinstall/wooldid.ado:1313-1315`, the `untreatedZero` branch becomes:
```stata
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
      * untreatedZero declares g() = 0, so no ringer of ANY spelling may survive. The pooled
      * i.`null'## set exists to displace the margins-visible pooled base-dose column; under
      * untreatedZero the margins-visible dose set is saturated over the whole sample (single `#',
      * so the (_i=0,_t=0) cell carries its own slope), the pooled ringer is inside its span by
      * construction, and the sweep displaces a TREATED CELL SLOPE instead - every ame comes back
      * as truth - beta(dropped cell). `wooldid_validate' (ado:4797) already guarantees
      * untreatedZero is never paired, so `ccringertail' is empty here and these three locals hold
      * nothing but the head.
      local contreatcontrolinserterFE
      local contreatcontrolinserterABS
      local contreatcontrolinserter
    }
```
Why 1b over 1a: it is immune to a future disjunct being added to the `ado:1234` gate, and it sits
next to the comment that already says "the only thing to do is not control for anything!".

**FIX-1a (minimal diff)** — at `ado:1234`, wrap the existing disjunction:
```stata
    if ("`contreatcontrols'" != "untreatedZero") & (strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson") {
```

**Evidence (patchA = FIX-1a, full 154-cell grid, `gridA_patchA.csv`):** exactly **8 cells move**, all
of them `untreatedZero` x {`customvce(robust)`, `customvce(cluster unit)`, `unconditionalse robust`,
`ccnoabsorb customvce(robust)`} on both designs; **146 cells bit-stable at rtol 1e-11; 0 return-code
changes**. dz 0.373750000 -> 0.973750000 (= truth); dn 0.373750000 -> 1.343750000 (= the default
engine's own answer for that spec). On the poisson twins (`pois_A.csv`) `pz`/`untreatedZero` goes
0.127500000 -> **0.347500000 = truth** on both poisson engines, and `pn` becomes engine-invariant at
0.497500000. **FIX-1b (patchC) was run on 32 verification cells — the 8 moved ones plus 24 controls spanning `basic`, `basic byTime`, `byTtre` on both designs — and returns numbers IDENTICAL to FIX-1a everywhere** (`gridC_patchC.csv` vs `gridA_patchA.csv`). The two sites are interchangeable, as the reading predicts: under `untreatedZero`, `wooldid_validate` forbids pairing, so `ccringertail` is empty and the three inserter locals hold nothing but the head.

**Blast radius / what could regress.** Only `untreatedZero` runs on a non-reghdfe engine change at
all, and they change from "silently the `contreatcontrols(basic)` answer" to "the answer the default
engine already gives". Three second-order consequences a maintainer must expect and not flag as a
regression: (a) `se`, `t`, `p`, the interval, `e(df_r)`, `rmse`, `r2_a` all move on those 8 cells,
because one explicit column is gone and the residual df rises by one; (b) the generic
multicollinearity caution stops printing there; (c) **`poisson` + `untreatedZero` no longer trips
ppmlhdfe's rc-3200 fallback** (measured: `FBBOX` in probe04_pinned.log vs probe04_A.log), so the
three-line "initial estimation of ppmlhdfe failed" caution disappears and the run gets faster.
`e(cmdline)`, `e(N)`, the sample, and every binary-treatment path are untouched.

**Risk: low.**

**Verifies how.** (1) On design Z, `untreatedZero` must return 0.973750000000 on all seven OLS
engines and 0.347500000000 on both poisson engines, to rtol 1e-10. (2) Engine-invariance sweep: for
every `contreatcontrols` value and both designs, all engines must agree to rtol 1e-10 — the only
rows that may differ across engines after the fix are none. (3) Re-run `gridA_pinned.csv` on the
patched build and confirm the diff is exactly the 8 cells and nothing else.

---

### FIX-2 — a reference-pool identification screen for the other control types (CLASS ii)

`wooldid_amescreen mode(bycoh)` (`ado:12885-12938`) already implements exactly the right test for
`byCohort`: within each treated cohort, over its `_t == 0` reference rows, `sd(contreat)` must be
`> 0` with at least 2 observations, else hard-refuse r(459). **Generalize its grouping variable and
call it once per requested control type.**

Proposed edit — replace the single `bycoh` call site at `ado:1378-1390` with a loop, and give
`wooldid_amescreen` a `gvar()` argument in place of the hard-wired `iorig()`:

| requested type | grouping variable passed as `gvar()` | reference pool |
|---|---|---|
| `basic` | a constant (one group) | all `_t == 0` rows |
| `byTreated` | `evertreated` | that group's `_t == 0` rows |
| `byTtre` | `modified_ttre` | that group's `_t == 0` rows |
| `byCohort` | `iorig` (unchanged, today's behaviour) | that cohort's `_t == 0` rows |
| `byTime` | `t` | that period's `_t == 0` rows |
| `untreatedZero` | — no screen: there is no g() to identify | — |

Message, per type: name the type, name the offending group, and state the remedy. For the common
case — dose identically zero outside treated cohort-periods — the remedy is one line:
`contreatcontrols(untreatedZero)`. The existing `bycoh` text is the template; keep its
`< 1e-3 * sd(sample)` near-degeneracy caution tier too.

**Cost.** One `gegen sd(), by()` plus one `gstats sum` per requested type on the estimation sample —
the same shape and the same order of cost as the screen that already runs for `byCohort`, i.e. O(N)
with gtools. Well inside "very computationally cheap even with large datasets".

**The predicate needs ONE refinement that the measurements forced, and it matters.** The raw
`bycoh` predicate — `min over groups of sd(contreat | _t == 0)` — has a **false positive**, and I
measured it: design Q, `byTreated`. Its never-treated group has `sd = 0` over its reference rows
(the dose is identically zero for never-treated units), yet `byTreated` returns the exact truth
0.973750000 there. The reason is that a group whose ringer column is **identically zero** is dropped
as an all-zero column and displaces nothing. The damage needs a group whose ringer column is LIVE
(nonzero somewhere) but has no independent variation on its reference rows — then the column lies in
the span of that group's treated cell#dose columns and the sweep takes a cell slope.

So the screen must be `bycoh`'s predicate **restricted to groups that carry a live treated dose**,
which is exactly what `mode(bycoh)` already gets for free by looping only over TREATED cohorts
(`bctr = (tvar < .)`). Written out:

> group `g` is DANGEROUS iff `g` contains at least one treated cohort-period row with a nonzero
> `contreat`, AND `sd(contreat)` over `g`'s `_t == 0` reference rows is 0 (or fewer than 2 such rows).

Note the `sd == 0` rather than "all zero": a nonzero CONSTANT reference dose is also dangerous,
because the constant times the reference-cell indicator is itself in the model's span. That is
already `bycoh`'s convention and it should be kept.

**Measured predicate values** (probe06 `PRED2`, `min` over groups of `sd(contreat | _t == 0)`):

| design | whole sample (`basic`) | by `t` (`byTime`) | by evertreated (`byTreated`) | by onset (`byTtre`) | by cohort (`byCohort`) |
|---|---|---|---|---|---|
| Z | 0.0000000000 | 0.0000000000 | 0.0000000000 | 0.0000000000 | 0.0000000000 |
| N | 0.5954769677 | 0.4963359987 | 0.5658205075 | 0.5658205075 | 0.4756605546 |
| P | 0.4245073339 | 0.3157090146 | 0.0000000000 | 0.0000000000 | 0.0000000000 |
| Q | 0.6692288946 | **0.0000000000** | 0.0000000000* | 0.0000000000* | 0.0000000000* |

(`modified_ttre` really does separate the onsets — on design N it takes 5 values, 4/6/8/10 and the
never-treated sentinel, 120/120/120/120/480 rows; its min coincides with the ever-treated min because
the minimising group is the never-treated pool, which both groupings contain.)

and the measured outcomes (section 2): Z wrong everywhere; N right everywhere; P right for
`basic`/`byTime`, wrong for `byTreated`; Q right for `basic`/`byTreated`, wrong for `byTime`.
**With the live-group restriction, the predicate reproduces all twelve outcomes exactly.** The
starred cells are the false positives the restriction removes (a never-treated-only group with a
null ringer column).

**Blast radius / what could regress.** This turns silently-wrong runs into refusals, so it can only
break a run that is producing a number today. Two real false-positive risks:
- **`byTime` is the exposed one.** Per-period reference pools are thin in ordinary staggered panels,
  and the last periods of a panel often contain no treated unit in its pre-period at all (design Q
  is exactly that case). A hard r(459) there would refuse specifications users run today. I would
  propose `byTime` follow the **SI-A1 pattern**: error under default `safety()`, caution and proceed
  under `safety(off)`, gated on `safetyinput` — while `basic`/`byTreated`/`byTtre` keep the
  `bycoh` precedent of an unconditional r(459). This is the one point in this proposal that is a
  judgement call for the author rather than a measurement.
- **`contreatwithin` must be reconciled with the screen, and the package's own text says so.**
  `contreatwithin` declares there is no untreated comparison pool, and `ado:11918` — the EXISTING
  multicollinearity error — steers users to it in exactly these words: *"Specify contreatwithin if
  you have a continuous treatment intensity variable and wish to partial out any effect of treatment
  when the treatment intensity is 0."* A `basic`-mode screen firing under `contreatwithin` would
  therefore refuse the run the package just recommended. **My recommendation: the new screens do NOT
  fire under `contreatwithin`** — but that is a ruling, not a measurement, and if the author prefers
  them to fire, `ado:11918`'s steer has to change with it. (`contreatwithin` is measured to have no
  effect on the defect itself: it does not suppress the ringer.)
- Because it is a *data* screen, the predicate is **necessary but not sufficient**: reference-pool
  dose that varies across but not within the absorbed groups would pass the screen and still be
  collinear. That is why FIX-2 and FIX-3 are a pair, not alternatives.

**Risk: medium** (it changes r(0) into r(459) for some live specifications).

**Verifies how.** (1) design Z: `basic`, `byTime`, `byTreated`, `byTtre` must all exit r(459) on all
seven engines, matching `byCohort`'s existing behaviour. (2) designs N, P, Q: every combination that
returns the exact truth today must still return it, bit-identical — the screen must fire nowhere
(`gridA_pinned.csv` dn rows and the `MIX` P/Q rows are the frozen expectation). (3) A thin-but-
nonzero design (a single reference observation per period) to confirm the `< 2 obs` leg behaves.

**FIX-2alt for `basic` only, offered because the author should see it.** patchB proves that
`basic`-without-the-head returns the exact truth on design Z. So a second remedy exists for `basic`
alone: **skip the head when the reference pool has no dose variation**, using the same predicate, and
carry on. It would produce a number rather than a refusal — but that number is `untreatedZero`'s
answer under a `basic` request, i.e. wooldid silently substituting a different model, which is the
pattern the `wooldid_validate` pairing refusal at `ado:4797` exists to prevent, and which the
project's "users debug by subtraction" doctrine argues against. **No such option exists for the byX
types — there the tail *is* the control, and removing it removes the requested specification.** I
recommend the refusal for uniformity; the choice is the author's.

---

### FIX-3 — escalate the partial-drop detector that already exists (defence in depth, zero cost)

`ado:11929-11945` already detects that a treatment-period interaction coefficient was dropped for
multicollinearity and counts it in `catastrophe`; `ado:11950-11980` (SI-A4/R6) already classifies the
**all-dropped** case, voids the ame block and explains it. The **partial-drop** case — the one this
whole report is about — falls through to a caution that says only *"This could imply your results
are not identified, depending on exactly what was dropped."*

Proposed edit: in the same `ado:11950` block, alongside the `amelive == 0 & amedrop > 0` branch, add
`amelive > 0 & amedrop > 0`, with text that states the consequence rather than the possibility:
every reported `ame` row is shifted by minus the slope of the dropped cell; every `att` and `ES_att`
row by that slope times the cell's mean dose; pre-treatment rows are unaffected, so a pre-trend check
will not reveal it; and if the dose is zero outside treated cohort-periods, the specification wanted
is `contreatcontrols(untreatedZero)`.

Keep it a **caution, not an error**, and say honestly what it cannot distinguish: the existing
`amelive`/`amedrop` scan sees *that* a dose-slope term was omitted, not *why*, so it cannot separate
"a ringer displaced a cell slope, every ame is shifted" from "this cell has no within-cell dose
variation, only its own contribution is affected". Under `lattice()` drops are routine and the
existing text already says so — so gate any stronger wording on `"`lattice'" == ""`.

**Cost:** zero. The loop and both counters already run.
**Blast radius:** display only; no number moves. **Risk: low.**
**Verifies how:** design Z `basic` on the default engine must print the new text (it prints the
generic caution today, with `catastrophe` = 1); design N `basic` must print nothing new; a
`lattice(all 3)` run with expected cell drops must keep today's "likely innocuous" wording.

---

### FIX-4 — `poisson` always pays a failed ppmlhdfe attempt when the head is the only absorbed ringer

Measured (probe04, `FBBOX`): ppmlhdfe's first attempt at `ado:1462` fails r(3200) **every time**
`ccringerheadabs` survives into `absorb()` — i.e. for plain `basic` (and, before FIX-1, for
`untreatedZero`) — because `i.`null'##c.(`ringer')` is a single-level absorbed fixed effect. With any
`byX` present, DF-1 blanks `ccringerheadabs` and the first attempt succeeds. So DF-1's "absorb the
finest partition" optimisation is **dead code on the poisson path for exactly the specification it
was written for**, and every such run prints a three-line caution about the user's data that is
really about wooldid's own `i.null`.

Proposed edit at `ado:1455-1470`: when `"`poisson'" == "poisson"` and `ccringerheadabs` is non-empty,
build the absorb spelling without the head and carry the head in `contreatcontrolinserterFE` from the
start — which is precisely what the rc-3200 fallback already does.

**Blast radius:** none numerically *by construction* — the fallback is already what produces every
reported number on those runs, so skipping the attempt that always fails cannot change an estimate.
Verified on `pois_pinned.csv` (28 cells) and the five `FBBOX` traces only; **the maintainer should
confirm bit-identity on a wider poisson sweep** before shipping. It removes one wasted full
estimation and one misleading caution per run.
**Risk: low.** **Verifies how:** every `pois_pinned.csv` row reproduced bit-identically while
"initial estimation of ppmlhdfe failed" disappears from the log for `basic`; `byX` rows unchanged.

---

### FIX-5 — documentation (small, and independent of the code)

- The `untreatedZero` help section should state the identification requirement of the *other* types,
  since that is the choice a user is making: `basic`/`byTime`/`byTreated`/`byTtre`/`byCohort` all
  need the dose to VARY over the untreated/reference observations of the relevant group; if the dose
  is zero (or constant) there, `untreatedZero` is the only estimable type.
- `sthlp:1076`'s "on some of these designs" family should gain the fact that the generic
  multicollinearity caution, when it names a `co.<contreat>` term, means the ame rows are shifted —
  not merely "could imply".

---

### FIX-6 — `contreatpoly(>1)` under `untreatedZero` emits a bare pooled `c.contreat`. **HELD — uncertain.**

Structural part, **confirmed**: `contreatinj` at `ado:1314` concatenates to
`i._i#i._t#i.targetflag#c.contreat##c.contreat` for `contreatpoly(2)`. Stata's factor-variable
expansion of `A#B##C` emits a **bare `C` main effect**, and a direct column dump of that exact term
(probe05, `POLYCOLS`) shows a bare `dose` column in the design matrix. Under `untreatedZero`, whose
documented definition is "by not including a g() function at all", that is a margins-visible pooled
baseline dose term the option says is not there.

Numeric part, **now measured**: I dumped the omitted columns of a hand replica of the
`customvce` varlist on design N with `contreatpoly(2)` (probe06 part C). Of 42 omitted columns,
**twelve are structurally empty cells** (cohort 4 exists only at `_tt == 1`, cohort 3 only to
`_tt == 3`, cohort 2 only to `_tt == 5`) and the other three are `o.dose`, `1o.nullv#co.ringer` and
`1o.nullv#co.ringer2`. **No live cell slope is dropped.** So on design N the bare pooled `dose` term
that the `A#B##C` expansion emits *absorbs the drop that would otherwise take a cell slope* — which
is why `ame` comes back at the truth 0.973750000 there, and it also explains why that run's `att`
equals `basic`'s. On design Z the same specification is still wrong (`ame` 0.373750000), because the
dose is zero off treatment and there are two independent dependencies to resolve rather than one.

So `contreatpoly(>1)` does not fix anything; it **changes which model `untreatedZero` estimates** —
it adds a pooled baseline dose term the option says is not there, and on ordinary data that term
silently converts `untreatedZero` into something closer to `basic`. **I am still not proposing a code
edit**, because the right remedy depends on a ruling: (a) respell the `untreatedZero` injection so
every polynomial power is cell-interacted and no bare term can appear — e.g.
`i._i#i._t#i.targetflag#(c.contreat c.contreat#c.contreat ...)` — which changes numbers on every
`untreatedZero` + `contreatpoly(>1)` run that exists today, or (b) document the term. Severity if
confirmed as a defect: `untreatedZero` + `contreatpoly(>1)` is not the model the help describes.

---

## 4. DOES THE MECHANISM REACH ANY SURFACE OTHER THAN `contreat`?

- **`interactivecontrols()` — YES, confirmed engine divergence.** `interinject` at `ado:1171` is
  `i._i#i._t#1.targetflag#c.(interlist)  i.t#c.(interactivecontrols)`; the second, pooled-by-period
  set is collinear with the first whenever the control is zero off treated cells. On design Z with
  `interactivecontrols(xz)`, `xz = dose` (binary-treatment run, true att 1.754034490530):
  default `robust` returns **1.754034490530** (exact), `customvce(robust)` returns
  **0.890310303808**. **The divergence is confirmed; the MECHANISM is NOT established** — and the
  contreat analysis does not transfer directly, because `interinject` interacts the cell terms with
  `1.targetflag` (treated cells only) where `contreatinj` uses `i.targetflag` (all cells), so the
  collinearity geometry is different, and on a reading of the varlist order the cell terms precede
  the `i.t#c.x` terms on BOTH engines, which is not what a simple sweep-order story predicts. I did
  not chase the dropped column. What is solid: on an exact-fit design where the truth is
  1.754034490530, one engine returns it exactly and another is 49% low, with no caution that
  distinguishes them. Reproducible from `probe05.do` (`SURF` lines). This surface has its own screen
  for the *all*-dropped case (`wooldid_icscreen`, `ado:12706`) and nothing for the partial case, so
  **FIX-3's escalation would at least surface it**.
- **`coarsecohortcontrols()` — not reached by this mechanism, and untested.** `ado:1138` builds
  `i._i#i._t#(coarsecohortcontrols)` with **no pooled counterpart**, so there is no ringer to
  displace anything. My smoke call returned r(198) (a syntax mismatch on my part, not a finding).
- **`ccc_absorb()`, `timetrends`, `controls()`** — no pooled-twin construction of this shape;
  `timetrends` measured to have no effect on the defect.

---

## 5. UNTESTED / OPEN (honest list)

- `lattice()` beyond `all 2`: the other eleven lattice types are untested. The mechanism does not
  depend on lattice type, but I do not claim coverage.
- `contreatpoly(3+)`, `semielasticity`, `contreatelasticitytype()`, `poisexpresults`.
- `permute()` propagation (the `reg` door is refused for `bootstrap()`; I did not check `permute()`'s
  own refusals beyond the `unconditionalse` one at `ado:5802`).
- `fe()` with a unit-level fixed effect alongside a coarse `i` — the configuration the help now
  recommends as the default — crossed with these control types.
- Whether the `byTime` screen of FIX-2 would fire on real staggered panels at a rate the author finds
  acceptable. This needs a realistic-panel sweep, not an exact-fit design.
- The exact dropped column on the `interactivecontrols()` divergence.
- Cluster-robust SE behaviour under the fix: I verified point estimates throughout and only noted
  that SEs move on the 8 fixed cells.
- **`i` at UNIT level with no `fe()`** — the configuration SI-8 steers users away from but which is
  the historical default — was never run. `untreatedZero`'s saturation argument (a single `#`, so the
  `(_i=0,_t=0)` cell carries its own dose slope) holds whatever the granularity of `i`, so FIX-1
  transfers; FIX-2's live-group restriction interacts with cell size and I have no measurement there.
