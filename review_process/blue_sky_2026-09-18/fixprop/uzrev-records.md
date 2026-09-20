# uzrev-records — the record correction, and whether the untreatedZero fix is complete

**Agent:** uzrev-records. **Workdir:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uzrev-records`
**Pinned build:** `.../bluesky/build/wooldid.ado` (c1d8979). **Patched copy under test:** `probes/uzrev-records/buildB/wooldid.ado` (FORM B, the three-line clear).
**Repo untouched.** Everything below is text for the author to paste; nothing was edited in `/Users/hegland/Documents/Research/wooldid_dev`.

---

## VERDICT IN ONE PARAGRAPH

**Ship the primary fix (either form — they are provably, not just measurably, identical).** It is
strictly larger than advertised: besides the 8 grid cells, it repairs `poisson` on its *default*
engine, `lattice()` under `customvce`, and `contreatpoly(2)`, and on every design I built it lands
the non-reghdfe engines exactly on the default-path answer. **But the fix set is NOT complete**, and
two of the three open items are now settled rather than open:

1. **CRIT-5 is untouched — confirmed by measurement, not inference.** A 60-cell pinned-vs-patched
   grid moves **exactly 2 rows**, both `untreatedZero x customvce(robust)`. Every `basic`, `byTime`,
   `byTreated`, `byTtre`, `byCohort` row is **bit-identical in hex**, including the wrong ones.
2. **The `interactivecontrols()` divergence mechanism is now ESTABLISHED** (it was "unestablished"),
   and it is the *same shape* of defect on a *different* varlist with **no gate at all** — so no
   form of the untreatedZero fix can reach it. It also **contradicts a shipped help-file paragraph**
   (`sthlp:1074-1078`) and its twin code comment (`ado:1160-1170`), which say SI-20 removed this
   order-dependence. It removed one instance of it; a second instance is live today.
3. **The factor-ordering claim is REAL, and the coordinator's non-reproduction is also real.** Both
   are explained by one exact identity I verified on six runs across two designs. The claim should
   go back into the record, with the identity, not be dropped.

**What to change before shipping (the short list):**
1. **Ship the primary fix, FORM B**, and put **R4** (the code comment directly above the branch) in
   the same commit — otherwise the next refactor recreates the defect.
2. **Hold FIX-UZ-C / Form C** out of that commit. It is the only piece that moves a default-path
   number, and its failure mode did not reproduce (§2E).
3. **Land R1, R2, R3, R5** as record corrections. R5 is a shipped help-file paragraph and a code
   comment that are now measurably wrong; it is not bookkeeping.
4. **Open two tickets**, both confirmed unmoved by this fix: **CRIT-5** (`basic`/`byTime`/
   `byTreated`/`byTtre`, every engine, `byTtre` sign-flipped) and an
   **`interactivecontrols()` zero-off-treatment screen**.
5. **Add the missing battery cells**: `untreatedZero` x {`customvce`, `unconditionalse`, `poisson`}
   and `interactivecontrols()` x `customvce`. The battery has never run either off the default engine
   (R6).

**FIX-UZ-C / Form C: I could not reproduce the failure it repairs.** On an exactly-saturated
quadratic design (`e(r2) = 1.00000000`), Form B *alone* took `contreatpoly(2) + customvce(robust)`
from 0.872496406950 to the design truth 1.703627665609 exactly. I am not refuting Form C — its
`##`-emits-a-pooled-main-effect reading of Stata's grammar is correct — but the author should know
that its necessity is not demonstrated on my designs, and that it is the only piece of the proposal
that moves a default-path number.

---

# PART 1 — THE RECORD CORRECTION

Six records. R1 and R2 are the two the task names; R3, R4 and R5 are new; R6 is the coverage gap
that let the battery stay green. For each: the **verbatim** false sentence, why it is false, and
**replacement text**.

---

## R1 — `review_process/cell_collapse_contreat_2026-09-02.md`, §1 "What the prototype does", lines 115-120

### The false sentence (verbatim, lines 118-120)

> `contreatcontrols(untreatedZero)` **alone** puts the
>   ringer nowhere - it only rewrites `contreatinj` (1152-1154) - so the prototype's own comment,
>   which names `untreatedZero` as a ringer path, is wrong on that point (section 7).

### Why it is false

`untreatedZero` puts the ringer **in the model on every engine except reghdfe**. The gate that
builds the ringer (`ado:1234` on the pinned build) is

```
if strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
```

— a **disjunction with no `untreatedZero` exclusion**. The `untreatedZero` branch at `ado:1313-1315`
rewrites `contreatinj` and never clears `contreatcontrolinserterFE` / `ABS` / the plain inserter, so
under `customvce()`, `unconditionalse` or `poisson` the ringer is carried into the regression exactly
as it is under `basic`.

**The document refutes itself, twice.** Two sentences earlier, at lines 115-117, it correctly
writes that the ringer reaches the regression "**and unconditionally under `unconditionalse`,
`customvce()` or `poisson`**". And row (vi) of its own proviso table, at **line 151**, states the
rule correctly a third time: *"its level enters `xb` wherever `contreatcontrolinserterFE` is
non-empty (1123-1141: `basic` without a lattice, `byTreated`/`byTtre`/`byCohort`/`byTime`, or any of
`unconditionalse`/`customvce()`/`poisson`)"*. The refuting evidence is adjacent to the false
conclusion on both sides; the error is that the sentence at 118 scoped "alone" to the wrong thing.
"alone" is true of the *default reghdfe path* and false of the option. **Line 151 needs no edit** —
it is already right, and it is the sentence a reader should be sent to.

**Measured** (`probes/uzrev-records/grid_pinned.csv`, exact-fit design, `e(r2) = 1.0000`, truth
ame 0.973750000000):

| engine | `untreatedZero` ame |
|---|---|
| default (reghdfe) | 0.973750000000 |
| `robust` | 0.973750000000 |
| `ccnoabsorb` | 0.973750000000 |
| **`customvce(robust)`** | **0.373750000000** — which is exactly `contreatcontrols(basic)`'s answer |

### Replacement text (paste over lines 118-120)

```
  the absorb into the explicit regressor list. `contreatcontrols(untreatedZero)` **does not** stop
  the ringer being built: the gate one clause above is a disjunction, so whenever `unconditionalse`,
  `customvce()` or `poisson` is also given, `contreatcontrolinserterFE`/`ABS` are populated exactly
  as under `basic` and the `untreatedZero` branch (1152-1154 in the prototype, 1313-1315 shipped)
  rewrites `contreatinj` without clearing them. The ringer therefore reaches the model on every
  engine except plain reghdfe, and because the dose is identically zero off treated cohort-periods
  the pooled column equals the sum of the cell#dose columns, so one cell slope is dropped and every
  reported quantity comes back as `truth - beta(dropped cell)`. `untreatedZero` puts the ringer
  nowhere ONLY on the default reghdfe path, which is the one path that never reads
  `contreatcontrolinserterFE`. The prototype's own comment naming `untreatedZero` as a ringer path
  is therefore RIGHT for the non-reghdfe engines and wrong only for the default one; it should be
  narrowed, not deleted (see section 7). CORRECTED 2026-09-19 after blue-sky CRIT-2/CRIT-3.
```

---

## R2 — `review_process/cell_collapse_contreat_2026-09-02.md`, §7 item 2, lines 802-806

This is the operative half of R1: it is a **build instruction** derived from the false claim, telling
the author to delete a comment that was correct.

### The false sentence (verbatim, lines 802-803 and the first half of 804)

> 2. **The prototype's own comments** need two corrections before they ship: the collapse-site
>    comment names `untreatedZero` as one of the paths that put the ringer in the model, which it is
>    not (section 1);

### Why it is false

Same mechanism as R1. The prototype comment was correct for `customvce()`/`unconditionalse`/`poisson`;
this item instructed its removal. The item's **second** correction (the worker comment's "measured:
att-type reads moved by up to 0.045 absolute", lines 804-806) is unaffected and must be preserved —
hence the replacement below spans the whole item, lines 802-806.

### Replacement text (paste over lines 802-806, the whole of item 2)

```
2. **The prototype's own comments** need two corrections before they ship: the collapse-site
   comment names `untreatedZero` as one of the paths that put the ringer in the model, which is
   correct for `unconditionalse`, `customvce()` and `poisson` and wrong only for the default reghdfe
   path, so it should be narrowed to name the engines rather than removed (section 1, corrected
   2026-09-19); and the worker comment's "measured: att-type reads moved by up to 0.045
   absolute" should carry the sharper V1b(f) measurement of section 2.2 (a pre-trend horizon moving
   p 0.167 -> 0.0022).
```

---

## R3 — `review_process/amese/verifier1_adversarial.md`, §1 "False-fire attacks", lines 21-22

### The false sentence (verbatim)

> **P4** dose zeroed outside treated cells + contreatcontrols(untreatedZero) + safety(off) +
> cluster(unitid): no fire, ame se .0628, V_es_ame posted. **OK**

### Why it needs correcting

The finding itself ("no false fire of the AMESE R1 predicate") is **still true** — that is not what
is wrong. What is wrong is the **verdict token `OK` applied to the configuration**. P4 is the only
place in the AMESE campaign where `untreatedZero` was attacked, and it ran on the **default reghdfe
path only** (`cluster(unitid)` keeps the run on reghdfe). Under `customvce()`, `unconditionalse` or
`poisson` the same configuration returns the `contreatcontrols(basic)` point estimates with the same
`e(N)`, the same `e(r2)` and the same `e(df_r)`. A scrub-predicate probe that never leaves the
default engine cannot see it, and the `OK` was read downstream as clearance for the option.

### Replacement text (paste over lines 21-22)

```
**P4** dose zeroed outside treated cells + contreatcontrols(untreatedZero) + safety(off) +
cluster(unitid): no fire, ame se .0628, V_es_ame posted. **OK for the R1 predicate, DEFAULT ENGINE
ONLY.** Scope correction 2026-09-19: this arm keeps the run on reghdfe, which is the one path that
never reads `contreatcontrolinserterFE`. `untreatedZero` under `customvce()`, `unconditionalse` or
`poisson` silently returns `contreatcontrols(basic)`'s point estimates (blue-sky CRIT-2/CRIT-3); this
`OK` is about false-firing of the scrub, and is NOT clearance for `contreatcontrols(untreatedZero)`.
```

---

## R4 — `wooldidinstall/wooldid.ado`, lines 1311-1312 (code comment)

### The false sentence (verbatim)

> ```
>     * if always zero when not treated, the only thing to do is not control for anything! this results in no contreat controls
>     * it is implemented via modifying the contreatinj statement and just not having contreatcontrols of any kind
> ```

### Why it is false

"this results in no contreat controls" is true only when the gate 79 lines above did not fire.
Modifying `contreatinj` is **not** sufficient to have "no contreatcontrols of any kind": the three
inserter locals are already populated by then on three of the four engines. This comment is the
in-code statement of the same wrong belief as R1, and it sits immediately above the branch the fix
edits — so it must change in the same commit or it will re-create the defect at the next refactor.

### Replacement text (paste over lines 1311-1312; assumes FORM B is taken)

```
    * if always zero when not treated, the only thing to do is not control for anything! this results in no contreat controls
    * it is implemented via modifying the contreatinj statement AND clearing the three ringer
    * inserters. Clearing them is load-bearing, not tidiness: the gate above is a disjunction
    * (`basic` OR unconditionalse OR customvce() OR poisson) with no untreatedZero exclusion, so on
    * every engine that is not plain reghdfe the pooled i.`null'##ringer set was already built and
    * would be carried into the regression -- and because the dose is identically 0 off treated
    * cohort-periods that pooled column EQUALS the sum of the cell#dose columns, so the engine drops
    * one cell slope and every reported quantity comes back as truth - beta(dropped cell), i.e.
    * exactly contreatcontrols(basic)'s answer, with no caution and the same e(N)/e(r2)/e(df_r).
    * `ccringertail' is always empty here (wooldid_validate forbids pairing untreatedZero with any
    * other control type), so clearing these three is equivalent to excluding untreatedZero from the
    * gate above; this spelling is preferred because it keeps the exclusion next to the branch that
    * needs it. Blue-sky CRIT-2/CRIT-3, 2026-09-19.
```

---

## R5 — `wooldidinstall/wooldid.sthlp`, lines 1074-1078 (and its code twin, `ado:1160-1170`)

**This is a CONTRADICTION of a shipped, ruled record, and it is the reason the `interactivecontrols()`
divergence has no owner.**

### The false sentence (verbatim, sthlp lines 1074-1078)

> The treatment-effect interactions are built for the treated cohort-by-period cells only: the pooled reference cell (untreated observations) carries no interaction slope, so its
> de-meaned-control variation is partialled out through the time-period interactions and the reference-period effect is a structural zero. Before 2026-09-05 the reference cell also received a
> slope column that was exactly collinear with the time-period interactions plus the treated-cell slopes, so the estimation engine's collinearity drop decided the reported horizons: the
> horizon estimates, their standard errors, the pre-treatment joint tests and {it:e(b_es_*)}/{it:e(V_es_*)} were solver- and regressor-order-dependent wherever the control varies within
> cohort over time, under weights, or on an unbalanced panel. Results under this option changed on that date.

The same claim is in the code at `ado:1160-1170` ("With 1.targetflag the reference-cell slope is
partialled out through i.t#c.x, the reference margin is a structural zero (0bno.) and the horizons
are estimable functions of the model.").

### Why it is false

SI-20 removed **one** instance of the collinearity — the reference cell's. A **second** instance is
live today, and it does not involve the reference cell at all. `interinject` (`ado:1171`) is

```
local interinject i.`_i'#i.`_t'#1.`targetflag'#c.(`interlist')  i.`t'#c.(`interactivecontrols')
```

Whenever the interactive control is **zero outside treated cohort-periods**, the pooled
`i.t#c.(control)` column at each calendar period `t` is *identically* the sum of the treated-cell
slope columns at that `t` (modulo the cell dummies, which are already in the model, because
`interlist` is the within-`i` demeaned copy). So there is one exact collinear relation **per
post-onset calendar period**, and the horizons are still "solver- and regressor-order-dependent".

**Measured** (`probes/uzrev-records/p2_interctl.log`, design Z, `xz = dose`, true att 1.754034490530):

- Span check, single period: regressing the raw control at `t = 7` on the cell dummies + cell slopes
  returns `r2 = 1.0000000000000000`, `rmse = 0.0000000000000e+00`. The relation is exact.
- **reghdfe spelling** drops **all ten** `t#c.xz` columns, keeps every live cell slope, returns
  **1.754034490530** — the truth.
- **plain-`reg` spelling** (the `customvce()`/`unconditionalse`/ppmlhdfe-fallback path) keeps
  `t#c.xz` for `t = 4..10` and instead drops **seven live treated-cell slopes** —
  `1#1, 1#2, 2#1, 2#2, 3#1, 3#2, 4#1`, i.e. the first treated cell at each post-onset period — and
  returns **0.890310303808**, 49% low.
- Control: the same `reg` **without** the pooled `i.t#c.xz` set reproduces reghdfe's cell
  coefficients to 10 dp. The pooled set is the whole story.
- **Scope of that last bullet, stated honestly.** The omitted-column lists come from a hand
  reconstruction of `ado:1476`/`ado:1485` (same span, same spelling, same order), not from inside
  wooldid. Its corroboration is that the reconstruction's seven extra drops match wooldid's own
  *"triggered 7 times total"* exactly, and that the pooled-set-removed control reproduces the correct
  engine's coefficients. I did **not** show the reconstruction reproduces 0.890310303808 itself, so
  the shipped help sentence above is written to what is proven — *which* columns the explicit paths
  drop is corroborated, not certified.
- What the user sees: `Caution: a treatment period coefficient interacted with something ... was
  dropped for multicollinearity` / `Note: the multicollinearity caution above was triggered 7 times
  total.` The count is exactly the number of collinear relations. `wooldid_icscreen` (`ado:12720`)
  covers the all-dropped case only and does not fire.

### Replacement text (paste over sthlp lines 1074-1078)

```
The treatment-effect interactions are built for the treated cohort-by-period cells only: the pooled reference cell (untreated observations) carries no interaction slope, so its
de-meaned-control variation is partialled out through the time-period interactions and the reference-period effect is a structural zero. Before 2026-09-05 the reference cell also received a
slope column that was exactly collinear with the time-period interactions plus the treated-cell slopes, so the estimation engine's collinearity drop decided the reported horizons: the
horizon estimates, their standard errors, the pre-treatment joint tests and {it:e(b_es_*)}/{it:e(V_es_*)} were solver- and regressor-order-dependent wherever the control varies within
cohort over time, under weights, or on an unbalanced panel. Results under this option changed on that date. That change removed the REFERENCE-CELL instance of the collinearity only. A second
instance remains and is not repaired by it: if an interactive control is zero (or constant) outside treated cohort-periods, the pooled time-period interaction at each post-onset period is
exactly the sum of that period's treated-cell slopes, so one column per post-onset period must be dropped, and which one is dropped still depends on the engine. On the default reghdfe path
the pooled set is dropped and the treated-cell slopes survive; on the {it:customvce()}, {it:unconditionalse} and poisson-fallback paths treated-cell slopes are dropped instead, which moves
the reported {it:att} and every event-study horizon. A multicollinearity caution is printed, once per dropped column, but it does not say which engine kept which set. Do not use an interactive control that is zero outside treated cohort-periods; specify it as a continuous treatment ({cmd:contreat()}) instead. {p_end}
```

(The code twin at `ado:1160-1170` needs the same qualification appended: SI-20's `1.targetflag`
removed the reference-cell column, **not** the per-period relation among the treated-cell slopes.)

---

## R6 — the coverage statements that are true but read as clearance

Not falsehoods; name them so the next reader does not take them for coverage.

- **`test/suite/configs.csv:15`** — the *only* `untreatedZero` battery config is
  `contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2)`.
  `cluster(i)` keeps it on reghdfe. **The battery has never run `untreatedZero` off the default
  engine.** Suggested ledger note: *"cont_uz is default-engine only; the untreatedZero x
  customvce/unconditionalse/poisson cells were added 2026-09-19 after CRIT-2/CRIT-3."*
- **`review_process/assessment_grid.md:25-29` (D2) vs `:35` (D4)** — D2 enumerates the
  `contreatcontrols()` values and D4 enumerates the estimator paths, but nothing requires D2 x D4 to
  be crossed. That un-crossed product is the hole. Suggested amendment to D4: *"D2 x D4 must be
  crossed for every `contreatcontrols()` value: the ringer gate at `ado:1234` is engine-conditioned,
  so a control type can be correct on reghdfe and wrong on `reg`."*
- **`review_process/research_2026-09-06/F032.md:292-296`** is the near-miss and should be cited, not
  corrected — it already recorded, honestly, that *"`contreatcontrols(basic)` and `untreatedZero`
  were not both swept"* and that *"`unconditionalse`, `customvce()`, `impvalues` and
  `poisexpresults` paths were not run"*. Both halves of the defect are named in one untested-list.

### Checked and found NOT false (negative results, so nobody re-opens them)

- `wooldid_standing_issues.md:320-326` (the AMESE untreatedZero SE answer) — "Under robust/fine
  clustering untreatedZero inherits only the standard small-cell condition". Plain `robust` stays on
  reghdfe; measured ame 0.973750000000, correct. **Stands.**
- `review_process/research_2026-09-08/engine_parity.md:56` (permute refuses `untreatedZero`,
  bootstrap forwards it) — the refusal is live at `ado:5871-5872`. **Stands.** Read in full: despite
  its title this document is about the *bootstrap-vs-permute* engines, not the *estimation* engines,
  and contains no reghdfe-vs-`reg` parity claim anywhere. The name is a trap for the next searcher
  and is worth one clarifying clause in its opening line.
- `wooldidinstall/wooldid.ado:11985` (R6 amescreen "Scoping on targetflag==1 exempts untreatedZero's
  expected zero-dose reference periods") — true; the screen covers the all-dropped case, and the
  ringer defect is a partial drop. **Stands.**
- `README.md:103`, `sthlp:98/986/993-995` — descriptions of the option, not of the code path. Stand.
- `review_process/blue_sky_2026-09-18/unadjudicated.md` — read in full; it lists neither the
  `interactivecontrols()` divergence nor the factor-ordering claim, so nothing there needs closing
  out against this report. Its six open items are unrelated.
- **Search method, stated so it can be audited.** Pass 1 was the token `untreatedZero` across
  `*.md`/`*.sthlp`/`*.ado`/`*.do`/`*.doh`/`*.csv` in the whole repo. Pass 2 was token-free, hunting
  engine-invariance claims scoped to continuous treatment:
  `grep -rn -iE "customvce|unconditionalse" review_process/*.md review_process/*/*.md | grep -iE "contreat|ringer|identical|agree|parity"`.
  Pass 2 surfaced no additional false record — only the two `cell_collapse_contreat` sites already
  listed, and correct statements elsewhere.
- **`sthlp:1033-1035`** ("the control type {it:untreatedZero} can be used to specify g() = 0") is
  false as written on three of four engines, but it becomes true the moment the primary fix ships.
  **Do not paste a correction; ship the fix.** If the fix is deferred, it needs a one-sentence
  engine caveat.

---

# PART 2 — IS THE FIX SET COMPLETE?

## 2A. CRIT-5 — CONFIRMED UNTOUCHED (measured, not inferred)

**Static argument.** Both forms live behind `if "\`contreatcontrols'" == "untreatedZero"` (Form B) or
inside the `ado:1234` gate qualified by `!= "untreatedZero"` (Form A). `basic`, `byTime`,
`byTreated`, `byTtre` never take that branch. Nothing in either form can reach them.

**Measurement.** `grid.do`, 60 cells = {design Z, design N} x {untreatedZero, basic, byTime,
byTreated, byTtre, byCohort} x {default, robust, customvce(robust), unconditionalse, ccnoabsorb},
pinned vs Form-B patched, compared in `%21x` hex:

```
$ diff grid_pinned.csv grid_patchB.csv
4c4    dz,untreatedZero,customvce(robust) ...
34c34  dn,untreatedZero,customvce(robust) ...
```

**Exactly 2 of 60 rows move. Both are `untreatedZero`.** `e(N)`, `e(df_r)` and `e(r2)` are unchanged
on both. Every `basic`/`byTime`/`byTreated`/`byTtre`/`byCohort` row is bit-identical.

**And here is what stays broken, on the pinned build, design Z, `e(r2) = 1.0000`, truth ame
0.973750000000 / att 1.754034490530 — identical before and after the fix:**

| `contreatcontrols()` | default | robust | customvce(robust) | ccnoabsorb |
|---|---|---|---|---|
| `basic` | **0.373750000000** | 0.373750000000 | 0.373750000000 | 0.373750000000 |
| `byTreated` | **0.373750000000** | 0.373750000000 | 0.373750000000 | 0.373750000000 |
| `byTime` | **0.200000000000** | 0.200000000000 | 0.200000000000 | 0.200000000000 |
| `byTtre` | **-0.046250000000** | -0.046250000000 | -0.046250000000 | -0.046250000000 |
| `byCohort` | r(459) | r(459) | r(459) | r(459) |

`byTtre` is **sign-flipped** against a truth of +0.97375 on the default path with no opt-in. This is
CRIT-5, independently reproduced on a third design, and **the primary fix does not move it by one
bit.** It needs FIX-2 (the generalised reference-pool screen) or its own ruling.

One caution about the grid for anyone reading it: `unconditionalse` returns r(198) on every row of
mine, because the option requires `robust` or a cluster spec alongside
(*"Error: the option unconditionalse requires specifying either robust or clustered..."*). That is
**not** a refutation of CRIT-3's `unconditionalse` arm — CRIT-3's run supplied the companion.

## 2B. `interactivecontrols(xz)` — MECHANISM ESTABLISHED, AND NOT COVERED BY ANY FORM OF THE FIX

Full mechanism in **R5** above. Three points for the fix decision:

1. It is **the same class of defect** — a pooled twin that is an exact linear combination of the
   per-cell terms, with an engine-dependent drop — on a **different varlist** (`interinject`,
   `ado:1171`).
2. There is **no gate to patch**. `contreatinj` has `ado:1234` and `ado:1313`; `interinject` is built
   unconditionally. Neither Form A, nor Form B, nor Form C touches `ado:1171`.
3. **The remedy is not "delete the pooled set."** `i.t#c.(interactivecontrols)` is load-bearing:
   `ado:1160-1170` records that it is what partials out the reference-cell slope, which is precisely
   what SI-20 relied on. The right shape is a **screen**: when an interactive control has no
   variation outside treated cohort-periods, either refuse (matching `wooldid_amescreen` mode(bycoh),
   which is the only member of the contreat family already protected) or steer the user to
   `contreat()`. That is FIX-3's escalation, and this measurement is the evidence for it.

Confirmed numbers, design Z, binary run, true att **1.754034490530**:
`robust` (reghdfe) **1.754034490530** exact; `customvce(robust)` (plain `reg`) **0.890310303808**.
Both `rc 0`. Both print the same generic 7-times multicollinearity caution.

## 2C. FACTOR ORDERING — SETTLED. The claim is REAL; the non-reproduction is also real; one identity explains both

**The identity.** On a dose-zero-off-treatment design, the reported ame is

> **reported_ame = truth − beta(the treated cell that sorts LAST in the factor ordering of the cohort id)**

Verified on six runs, two designs, three labellings each, all at `e(r2) = 1.000000`
(`p3_relab.log`, `p4_pinned.log`).

**Design A = the probe's design Z** (8 cohorts, onsets 4/6/8/10; `beta = 1.00 − 0.10i + 0.15r − 0.02r²`).
Last-cell slope by cohort, measured: **1.08, 1.08, 0.92, 0.60**.

| cohort labelling | last treated cohort in factor order | predicted = 0.97375 − beta_last | `basic`, default path (measured) |
|---|---|---|---|
| `i` (1..8) | 4 -> beta_last 0.60 | 0.373750000000 | **0.373750000000** |
| `9 - i` | 8 (= original 1) -> 1.08 | -0.106250000000 | **-0.106250000000** |
| a third permutation (1->3, 2->1, 3->4, 4->2) | 4 (= original 3) -> 0.92 | 0.053750000000 | **0.053750000000** |

Three labellings, three different answers, one sign flip, identical data. **The claim is real.**

**Design B = the coordinator's own `repro/c5b.do`** (4 cohorts, onsets 4/6/8;
`beta = 0.55 + 0.20·i + 0.10·r`). Measured `basic` = **-0.226666666667 under all three labellings**,
matching the coordinator exactly.

**Why it does not move there.** Last-cell slope by cohort, measured: **1.35, 1.35, 1.35**. It is an
arithmetic accident of the parameterisation: the design gives cohort `i` onset `2 + 2i`, so its last
relative period is `r_max = 8 − 2i`, and
`beta_last(i) = 0.55 + 0.20i + 0.10(8 − 2i) = 1.35` for **every** cohort — the `+0.20` per cohort is
exactly cancelled by the two periods of exposure each later cohort loses. Whichever cohort sorts
last, the dropped cell's slope is 1.35. And indeed `1.123333333333 − 1.35 = −0.226666666667`, the
measured value to 12 dp.

**Ruling this settles:** the coordinator's line in `critical.md` — *"So that claim is design-specific
at best; it should not be repeated without its own repro"* — was right that it is design-specific and
right to demand a repro. But "design-specific" understates it: **the claim is generically true and
the coordinator's design is the special case.** Cross-cell slope heterogeneity *at the last period of
each cohort* is what the effect needs, and that is the ordinary case. Suggested replacement for that
paragraph of `critical.md`:

```
**Reproduced and explained 2026-09-19** (uzrev-records, `p3_relab.do`): the probe's factor-ordering
claim is correct and generic. The reported ame equals `truth - beta(last treated cell in the factor
ordering of the cohort id)`, verified to 12 dp on three labellings of the probe's design (0.373750000
/ -0.106250000 / 0.053750000 against measured last-cell slopes 0.60 / 1.08 / 0.92). It does not move
on MY design because that design's last-cell slope is 1.35 for every cohort by arithmetic accident
(beta = .55 + .20i + .10r with r_max = 8 - 2i), so relabelling cannot change the dropped cell's
slope. My run was a true negative on a design that could not have shown the effect, not a refutation.
```

## 2D. What the primary fix covers that the proposal did not claim (all measured, all landing exactly on design truth)

| arm | pinned | Form-B patched | design truth |
|---|---|---|---|
| **`poisson`, DEFAULT engine**, exact-fit log-linear, true semi-elasticity 0.25 | ame **1.5e-16** (a structural zero) | **0.250000000000** (hex `+1.0X-002`) | 0.25 |
| `lattice(all 2)` + `customvce(robust)` | ame 0.373750000000 / att **0.763210355478** (a *second* drop) | 0.973750000000 / **1.754034490530** | 0.97375 / 1.754034490530 |
| `contreatpoly(2)` + `customvce(robust)`, exactly-saturated quadratic design, `e(r2) = 1.00000000` | **0.872496406950** | **1.703627665609** | 1.703627665609 |
| `contreatpoly(1)` + `customvce(robust)` on the same design | 0.874075896226 | 1.831450761426 = the default path, to the bit | — |

The poisson row is the strongest single argument for shipping: CRIT-2's claim that poisson is broken
on its *default* engine with no opt-in reproduces, and the fix takes it from a structural zero to the
exact truth.

**Two notes on that row.** (a) The pinned `1.5e-16` is **2C's identity at work**, not a separate
pathology: that design has a homogeneous true slope of 0.25 in *every* cell, so
`truth − beta(dropped cell) = 0.25 − 0.25 = 0`. The poisson breakage and the factor-ordering
behaviour are one mechanism seen from two sides — and a homogeneous-effect design is exactly where
this defect hides as a clean-looking zero rather than a visibly wrong number. (b) I did not establish
which scale wooldid's poisson `ame` reports (a level-scale `dydx` reading of this design would be
0.7496), but the contrast is scale-independent: a structural zero before, an exact hit on a design
constant after.

Post-fix, `untreatedZero` under `customvce(robust)` agrees with the default reghdfe path to
**~3e-15 relative** (`+1.f28f5c28f5c21X-001` vs `+1.f28f5c28f5c15X-001`) — ordinary cross-solver
noise, not a residual.

## 2E. FIX-UZ-C / Form C — necessity NOT reproduced (stated as a caution, not a refutation)

- **The grammar reading is correct.** `contreatinj` at `ado:1314` appends
  `#\`contreatpolypart'` to `i._i#i._t#i.targetflag`, and at `contreatpoly(2)`
  `contreatpolypart` is `c.dose##c.dose`, so Stata expands the top-level `##` into
  `(cells#dose) + (dose) + (cells#dose#dose)` — a **pooled `c.dose` main effect** really is emitted.
  At `contreatpoly(1)` there is no `##` and no pooled term, which is why Form C is inert at the
  default. That much I confirm.
- **The pooled term really does survive the primary fix — I measured it.** Hand-rebuilding
  `ado:1485`'s varlist under Form B at `contreatpoly(2)` (`p7_formc.log`), the bare pooled `dose`
  column is **not** dropped: it is estimated live at **0.600000000003**, with
  `e(r2) = 1.0000000000000000`. So Form C's premise is correct — the `##` spelling leaves a pooled
  dose main effect in the model that the primary fix does not remove.
- **But the failure it repairs did not occur.** On that same exactly-saturated quadratic design
  (`e(r2) = 1.00000000`, truth 1.703627665609), Form B *alone* returned the truth **exactly** at
  `contreatpoly(2)` on the default, `robust` and `customvce(robust)` arms. The surviving pooled
  column is margins-VISIBLE — it is a term in `contreat` itself, not a differently-named ringer copy —
  so `margins dydx()` differentiates through it and the reparameterisation cancels. That is the
  structural reason it did not bite here, and it is also why the two are not the same defect: the
  ringer is invisible to `margins`, the `##` main effect is not.
- **So the open question is narrower than either side has stated it.** Not "does the pooled term
  survive the primary fix" (it does — measured) but "is there a design on which the surviving,
  margins-visible pooled term biases a reported number". I could not build one.
- **Consequence for the decision.** Form C is the **only** piece of the proposal that moves a
  default-path number (att 2.444765194869 -> 4.298562435389), and that move is on
  *precondition-violating* data, where the pooled term is not collinear with anything — so there it
  is a model change, not a defect repair. Combined with the fact that on precondition-*satisfying*
  data I could not make the surviving term bias anything, the cheap, risk-free ordering is:
  **ship the primary fix now, hold Form C** until the cell that proposed it supplies a design where
  Form B alone still returns a wrong number at `contreatpoly(>1)`. If Form C is taken anyway, its
  default-path movement should be certified on its own, not folded into the primary fix's
  evidence.
- Form A vs Form B is **not** a judgement call: `wooldid_validate` (`ado:4827-4828`) forbids pairing
  `untreatedZero` with any other control type, so `ccringertail` is *always* empty in that branch and
  the two forms are **provably** identical, not merely measured so. I recommend **Form B**, because
  the exclusion then sits next to the branch that needs it, and because it survives a future edit to
  the `ado:1234` gate. Its only cost is a stray `tempvar null` constant (one variable against the
  `si8demand` preflight) that Form A would not create.

## 2F. Residual gaps — what neither the fix nor this audit covers

1. **CRIT-5** (`basic`/`byTime`/`byTreated`/`byTtre` on dose-zero-off-treatment data, every engine,
   `byTtre` sign-flipped). Needs FIX-2 or a ruling. **Unmoved by this fix.**
2. **`interactivecontrols()` on a control that is zero off treated cells.** Mechanism now
   established; no gate exists; needs a screen. **Unmoved by this fix.**
3. **`controls(xz)` with the same control** returned 0.650127060512 on the `robust` arm where the
   binary truth is 1.754034490530. I did **not** establish whether that is the defect or the correct
   answer for a different model (controlling for a post-treatment variable legitimately changes the
   estimand). Flagged, not diagnosed.
4. **`contreatwithin`, `subgroups()`, weights, `if`, and the bootstrap/permute replication frames**
   crossed with the fix: not run here. `permute()` already refuses `untreatedZero` (`ado:5871`);
   `bootstrap()` does not.
5. **Standard errors under the fix.** I verified point estimates and `e(N)`/`e(df_r)`/`e(r2)`
   throughout. I did not audit whether the post-fix `customvce`/`poisson` SEs are right, only that
   they are no longer computed off a rank-deficient design.
6. `lattice()` beyond `all 2`; `contreatpoly(3+)`; `semielasticity`; `contreatelasticitytype()`.

---

# EVIDENCE INDEX (all absolute, all self-contained)

| file | what it shows |
|---|---|
| `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uzrev-records/buildB/wooldid.ado` | the Form B patched ado (three added lines at 1315-1317) |
| `.../uzrev-records/grid.do`, `grid_pinned.do`, `grid_patchB.do` | the 60-cell pinned-vs-patched grid |
| `.../uzrev-records/grid_pinned.csv`, `grid_patchB.csv` | `%21x` hex results; `diff` shows exactly 2 moved rows |
| `.../uzrev-records/p2_interctl.do` / `.log` | the interactivecontrols span check (`r2 = 1`, `rmse = 0`) and the two engines' omitted-column sets |
| `.../uzrev-records/p3_relab.do` / `.log` | three labellings x two designs, the factor-ordering settlement |
| `.../uzrev-records/p4.do`, `p4_pinned.log`, `p4_patchB.log` | last-cell slopes (`ALAST`/`BLAST`), the icscreen check, `contreatpoly(2)` on design Z |
| `.../uzrev-records/p5.do`, `p5_pinned.log`, `p5_patchB.log` | quadratic exact-fit design, `lattice(all 2)`, poisson |
| `.../uzrev-records/p7_formc.do` / `.log` | Form C: the `##`-emitted pooled `c.dose` survives Form B, live at 0.600000000003 |
| `.../uzrev-records/NOTES.md` | running state |

**Reproduce:** `cd` to the workdir, then
`"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do grid_pinned.do` (then `grid_patchB.do`,
`p2_interctl.do`, `p3_relab.do`, `p4_pinned.do`, `p4_patchB.do`, `p5_pinned.do`, `p5_patchB.do`).
`grid.do` reads `dz.dta`/`dn.dta` from `.../bluesky/probes/uz-blast/`; `p3`/`p5` build their own data.
