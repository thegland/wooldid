# CRIT-2 / CRIT-3 — `contreatcontrols(untreatedZero)` on non-reghdfe engines
# FIX PROPOSALS: design and comparison of the alternatives

Agent: `uz-altfix`. Workdir `.../bluesky/probes/uz-altfix` (NOTES.md, probe_*.do, probe_*.log,
res_*.csv, build_{A,B,C2,D,E}/wooldid.ado). **Nothing in the repo was edited.** Every variant was
built as a private copy of the pinned `wooldid.ado` and loaded with `adopath ++`.

---

## 0. Two corrections to the brief, up front

**(a) The angle brief's premise about the detector is wrong, and the synthesis must not carry it
forward.** `wooldid_amescreen , mode(detect)` is **not** a colname scanner and does not look at
`co.`. It is the AMESE **R1 cell-containment** predicate: `gegen nunique(clustervar)` within
treated `(_i,_t)` cells, firing when a treated cell lies wholly inside one cluster. The `co.` /
`o.` colname scanner lives in **`wooldid_ccscreen`** (`ado:11816-11985`), called from `ado:1713`
and **gated on `safety != off`**. Its R6 branch is the `amelive == 0 & amedrop > 0` one.

**(b) The defect is engine-conditioned, not `contreatcontrols`-conditioned.** Under
`untreatedZero` on the **default** engine the gate at `ado:1234` never fires at all —
`strpos("untreatedZero","basic") == 0` — which is precisely why the default path is correct. The
gate fires under `untreatedZero` only through its three engine words (`unconditionalse`,
`customvce()`, `poisson`). This matters for the fix: the bug is the engine clause reaching a
control type that asks for no control set, not a mis-specified control set.

---

## 1. The mechanism, restated so the edit is obviously right

`contreatinj` is normally `` `latticepart2hash'##`contreatpolypart' ``, so the regression term
`` i.`_i'#i.`_t'#i.`targetflag'`contreatinj' `` expands (via the top-level `##`) to

    cell dummies  +  cell#dose slopes  +  A POOLED, MARGINS-VISIBLE `c.contreat` MAIN EFFECT.

That pooled main effect is what the ringer set `i.`null'##c.(ringer)` exists to knock out: the
ringer is a bit-exact double twin of the dose that appears **earlier** in the regressor list, so
the engine's collinearity drop takes the margins-visible column and the baseline dose slope is
partialled out invisibly. That is the manipulation the long comment at `ado:1225-1233` describes,
and `ado:1290-1299` (DF-1) is explicit that on the explicit-regressor paths the `i.`null'##`
spelling is load-bearing **for the drop order**.

`untreatedZero` rewrites `contreatinj` at `ado:1313-1315` to

    `latticepart1hash'#`contreatpolypart'  i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'

— a single `#` at the top level plus a bare cell-dummy term. **There is no pooled
margins-visible `c.contreat` main effect under `untreatedZero`.** Nothing is left to knock out.
This is exactly `sthlp:1033`'s "the control type `untreatedZero` can be used to specify g() = 0",
and `ado:1311`'s own comment, "the only thing to do is not control for anything!".

So under `untreatedZero` the ringer set is pure surplus rank. And because the dose is exactly 0
off treated cohort-periods (the option's premise), `c.ringer` **equals the sum of the cell#dose
columns**: exact collinearity, one treated-cell slope dropped, every reported quantity returned as
`truth - beta(dropped cell)`.

**This also settles the scope of the fix.** The engine clause is DESIGNED FOR every other
`contreatcontrols` value, including a `byX` with no `basic`: those all keep the `##` spelling and
therefore still have the pooled margins-visible dose column the ringer must take out.
`untreatedZero` is the **only** value that removes it. That is the argument against refusing the
combination (option C below) rather than fixing it.

---

## 2. What I measured (own DGP, own probes — the project suite was not used)

Design D1: 4 groups (onsets 4/5/6 plus one never-treated), 5 units each, t = 1..8, N = 160;
dose exactly 0 off treated cohort-periods, with within-cell dose variation; **exact fit**, so the
truth is analytic. `TRUTH_AME = 1.025000000000`, `TRUTH_ATT = 2.265500000000`,
`TRUTH_AME(poisson semi-elasticity) = 0.273333333333`.

### 2.1 The pinned build, reproducing CRIT-2/CRIT-3 on a fourth independent design

| run | `main:att` | `main:ame` |
|---|---|---|
| `untreatedZero`, default reghdfe | **2.265500000000** | **1.025000000000** (exact truth) |
| `untreatedZero customvce(ols)` | 0.954666666667 | **-0.275000000000** |
| `untreatedZero customvce(robust)` | 0.954666666667 | -0.275000000000 |
| `untreatedZero unconditionalse robust` | 0.954666666667 | -0.275000000000 |
| `contreatcontrols(basic)`, **default engine** | 0.954666666667 | **-0.275000000000** |
| `untreatedZero poisson` | 0.076333333333 | **-0.066666666667** (truth **+0.273333333333**) |

Implied dropped slope `1.025 - (-0.275) = 1.3000000000`, and this DGP's
`beta(g = 3, rel = 2) = 0.5 + 0.2*3 + 0.1*2 = 1.3` **exactly**. The broken number is
byte-for-byte `contreatcontrols(basic)`'s answer — the pairing `wooldid_validate` (`ado:4797`)
refuses when a user asks for it. The poisson arm is **sign-flipped** with no opt-in.
`untreatedZero + lattice(all 2) + customvce(ols)` is worse still on the att side (**-0.545333**).
`untreatedZero + ccnoabsorb` on the default engine is **correct** — the gate never fires there.

### 2.2 Why the project's suite could not have caught it (verified read-only)

`test/suite/configs.csv` + `configs_ext.csv` contain exactly four `untreatedZero` configs:
`cont_uz`, `cont_uz_ss` (default engine, `cluster(i)`, `safety(off)`),
`ext_cont_ext_untreatedzero` (default engine, `robust`, `safety(off)`) and
`ext_val_untreatedzero_combo` (a parse-refusal test).
`grep untreatedZero ... | grep -Ei 'customvce|unconditionalse|poisson'` returns **nothing**.
There is no `untreatedZero` x non-reghdfe-engine cell anywhere in the suite. The project's own
`test/suite/degeneracy_waste.py` already calls `untreatedZero` "THE ONE THIN SPOT ... a SINGLE
config with live inference". Note also that three of the four run `safety(off)`, which skips
`wooldid_ccscreen` outright — so even the generic caution would not have printed in-house.

---

## 3. The five shapes, designed concretely and (four of five) built and run

I built each as a private copy of `wooldid.ado` and ran an identical 24-call probe under each.

### Option A — exclude `untreatedZero` from the gate at `ado:1234`  **[RECOMMENDED]**

**Exact edit.** `ado:1234`, currently

```stata
    if strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
```

becomes

```stata
    * untreatedZero rewrites contreatinj at the foot of this block to a single-# spelling that
    * carries NO pooled margins-visible c.contreat main effect, so there is nothing for the
    * ringer set to take out and g() = 0 (sthlp:1033). The engine words above force the FE
    * spelling of a set that basic/byX requested; untreatedZero requests none, and on the
    * explicit-regressor paths the pooled ringer is then exactly the sum of the cell#dose
    * columns, so the engine drops a treated-cell slope and every reported quantity comes back
    * as truth - beta(dropped cell). [CRIT-2/CRIT-3, 2026-09-18]
    if (strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson") & strpos("`contreatcontrols'", "untreatedZero") == 0 {
```

and a one-line pointer goes into the `untreatedZero` branch at `ado:1313`:

```stata
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
      * BELT AND BRACES (they are already empty: the gate above is excluded on this value, and
      * wooldid_validate refuses untreatedZero paired with any byX, so ccringertail is empty too).
      * The invariant is that untreatedZero carries NO ringer set on ANY engine - CRIT-2/CRIT-3.
      local contreatcontrolinserterFE
      local contreatcontrolinserterABS
      local contreatcontrolinserter
    }
```

Optionally spell the predicate once, above the block, as
`local ccuz = strpos("`contreatcontrols'","untreatedZero") != 0` and use `& !`ccuz'`. That is
strictly more robust than the exact `==` comparison already used at `ado:1313`, which would fail
silently if `contreatcontrols` ever arrived with surrounding whitespace.

**Measured result.** All broken rows land on the analytic truth to the last digit:

| run | base `main:ame` | A `main:ame` | truth |
|---|---|---|---|
| `untreatedZero customvce(ols)` | -0.275000000000 | **1.025000000000** | 1.025000000000 |
| `untreatedZero customvce(robust)` | -0.275000000000 | **1.025000000000** | 1.025000000000 |
| `untreatedZero unconditionalse robust` | -0.275000000000 | **1.025000000000** | 1.025000000000 |
| `untreatedZero poisson` | -0.066666666667 | **0.273333333333** | 0.273333333333 |
| `untreatedZero poisson customvce(robust)` | -0.066666794455 | 0.273333358513 | 0.273333333333 (1e-7, ML tol) |
| `untreatedZero lattice(all 2) customvce(ols)` | -0.275000000000 | **1.025000000000** | 1.025000000000 |
| noisy, `customvce(cluster unit)` | 0.308808587696 | **1.081671067033** | = the default path's own answer, exactly |

**Blast radius, measured at hex precision.** Of 66 recorded quantities, **exactly 18 move, in 9
tags, and every one of them is `untreatedZero` x a non-reghdfe engine.** Everything else is
byte-identical to the pinned build: `contreatcontrols(basic)` on all four engines, `basic byTime`
on two engines, `basic + lattice(all 2)` (r322 both builds), every default-engine
`untreatedZero` run, `untreatedZero + ccnoabsorb`, `untreatedZero + lattice(all 3)`, and the
noisy default-engine run.

**`e(N)` and `e(df_r)` do not move anywhere** — `untreatedZero customvce(ols)` reports 160 / 125
before and after. The removed column was exactly collinear, so the rank, and with it every
df-dependent quantity, is unchanged; only the basis the engine reports changes. This is the
cleanest evidence that the edit is a reparameterization on premise-satisfying data, not a model
change.

**Two side benefits, both measured.** (i) The generic multicollinearity caution stops firing on
`untreatedZero` runs (13 occurrences in the base log, 7 in A's; the 7 survivors are the
legitimate `basic`/lattice arms). (ii) The **spurious `ppmlhdfe` rc-3200 fallback disappears for
`untreatedZero + poisson`** (2 fallback cautions in base, 1 in A): the single-level `i.`null'`
absorb term guaranteed ppmlhdfe's first attempt would fail, so every `untreatedZero + poisson`
run was paying for a failed fit plus a four-line caution it did not need.

**Risk: LOW.** Confined to one boolean; the disabled branch's only consumers are
`contreatcontrolinserterFE`/`ABS` (grep: `ccringerhead`/`ccringerheadabs`/`ccringertail` appear
nowhere else in the file).

### Option B — clear `contreatcontrolinserterFE`/`ABS` inside the `untreatedZero` branch at `ado:1313`

```stata
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
      local contreatcontrolinserterFE
      local contreatcontrolinserterABS
      local contreatcontrolinserter
    }
```

(The third line is belt-and-braces: `contreatcontrolinserter` is still empty at this point; the
`ccnoabsorb` block at `ado:1318-1325` runs afterwards and copies the now-empty FE macro.)

**Is it equivalent to A? PROVABLY YES TODAY, AND MEASURABLY SO.** `wooldid_validate` (`ado:4797`)
refuses `untreatedZero` paired with any other control word, so under `untreatedZero`
`ccringertail` is necessarily empty; hence `contreatcontrolinserterFE == ccringerhead` and
`contreatcontrolinserterABS == ccringerheadabs`, and clearing the composed pair is the same as
never building the head. **`res_A.csv` and `res_B.csv` are byte-identical on all 66 rows**,
including the hexadecimal standard errors. There is no path on which they differ.

**They differ only in where the invariant lives and how it fails later.**
- A states it where ringers are *born*, beside the comment that explains why the engine words are
  in the gate; a future control type added to the gate inherits the reasoning.
- B states it where `untreatedZero` is *special-cased*, and is the more defensive shape: it
  neutralises anything the block above may build, including a future `ccringertail` term someone
  adds without re-reading the gate. Its failure mode is the opposite one — it silently discards
  work that a future `untreatedZero`-compatible control set might legitimately want.
- B also leaves `tempvar null` created and `i.`null'` composed and thrown away, which is harmless
  but leaves a reader wondering why.

**Recommendation: ship A, and put B's clearing lines in as well if the author wants belt and
braces.** They cost nothing and cannot change a number (measured). Do **not** ship B alone: it
leaves the gate looking as though `untreatedZero` still wants a ringer set.

**Risk: LOW** (identical to A).
### Option C — REFUSE `untreatedZero` + a non-reghdfe engine  **[built, works, and I argue against it]**

**Exact edit** (prototype `build_C2`; site chosen per the existing precedent, immediately after the
`untreatedZero cannot be paired with other control types` refusal at `ado:4797-4800` inside
`wooldid_validate`, so it fires at parse time before any lattice is built):

```stata
    if strpos("`contreatcontrols'","untreatedZero") != 0 & ("`unconditionalse'" == "unconditionalse" | `"`customvce'"' != "" | "`poisson'" == "poisson") {
      dis as error "Error: contreatcontrols(untreatedZero) cannot be combined with unconditionalse, customvce() or poisson."
      error 198
    }
```

**Measured.** It refuses exactly the nine defect tags with r(198) and leaves every other arm
byte-identical to the pinned build. It works.

**Why I would not ship it.**
1. **It refuses `poisson`, which is `untreatedZero`'s most natural companion.** The option exists
   for a dose that "is unknowable or undefined outside of treated cohort-periods"
   (`sthlp:993-995`) — count outcomes with an undefined off-treatment intensity are the canonical
   case. `poisson` is the DEFAULT-engine poisson (ppmlhdfe); the user never opts into an explicit
   regressor path, so C refuses a combination the user cannot see is special.
   Measured on a realistic count panel (`ycnt = rpoisson(mu)`, `cluster(unit)`): the pinned build
   reports `ame = -0.399459850135` where the fix reports `+0.209670460377` — a sign flip C would
   answer with a refusal instead of a number.
2. **The underlying model is perfectly estimable on those engines.** Under A they return the
   default path's answer to the last bit. Refusing an estimable specification to avoid a bug in
   how it is spelled is the wrong trade when the spelling is a two-token fix.
3. **It is not analogous to the existing refusals.** `ado:4797` refuses `untreatedZero` + another
   control type because the two are substantively contradictory; `ado:5842` refuses `permute()`
   because reassigning onset while the dose path stays fixed *violates the premise*. Neither is
   true here: nothing about `customvce()` or `poisson` contradicts "dose is zero off treatment".
4. **It silently widens.** `bootstrap()` is allowed with `untreatedZero` today, and
   `bootstrap()` + `poisson` is allowed, so C also kills that combination; so does
   `impvalues` + `poisson` + `untreatedZero`.

**Where C IS the right tool:** as a *narrow* companion for the residual `contreatpoly(>1)` case
described in §5 — if the author does not want to re-spell the polynomial injection now.

**Risk: LOW to implement, MEDIUM in user cost.**

### Option D — detect the collinearity at runtime and refuse/caution  **[built, measured, NOT an enforcement layer]**

**What exists now (corrected).** `wooldid_ccscreen` (`ado:11816-11985`), called at `ado:1713`
**only when `safety != off`**, tokenizes `colnames e(b)` when `o.` is present. Two branches reach
the defect's neighbourhood: the `wordcountproblem > 3` branch (`ado:11929`) prints the GENERIC
caution the user actually sees today, and the R6 branch (`ado:11969`) fires only on
`amelive == 0 & amedrop > 0`. R6's own comment says it "exempts untreatedZero's expected
zero-dose reference periods by construction" — `untreatedZero` was considered here and waved
through.

**Prototype.** `build_D` passes `cc(`contreatcontrols')` into `wooldid_ccscreen`, records the
dropped dose tokens, and adds
`strpos("`cc'","untreatedZero") != 0 & `amedrop' > 0 & `amelive' > 0` as a second predicate.
`build_DA` is the same instrumentation on top of fix A. Measured counts (`amelive`/`amedrop`):

| run | pinned build | after fix A |
|---|---|---|
| `untreatedZero`, default engine | 12 / 0 (silent) | 12 / 0 |
| `untreatedZero customvce/unconditionalse/poisson` | **11 / 1 -> fires** | **12 / 0 -> nothing left to detect** |
| `contreatcontrols(basic)`, **default engine**, ordinary correct run | **11 / 1** | 11 / 1 |
| `contreatcontrols(basic byTime)`, correct run | **7 / 5** | 7 / 5 |
| `untreatedZero lattice(all 5)`, default engine | 24 / 0 | 24 / 0 |
| `untreatedZero lattice(all 5) customvce(ols)` | — | **23 / 1 -> still fires after A** |
| `untreatedZero customvce(ols)`, one treated cell with no within-cell dose variation | 10 / 2 | **11 / 1 -> still fires after A** |

**Verdict: widening the screen is not the fix, for five measured reasons.**
1. **An UNSCOPED widening false-positives on ordinary correct runs.** Plain
   `contreatcontrols(basic)` on the default engine already shows `amelive = 11, amedrop = 1`, and
   `basic byTime` shows `7 / 5`. Both are correct runs. So "a treated-cell dose slope was dropped
   while others live" is simply not a defect signature.
2. **A SCOPED widening (to `untreatedZero`) needs no scan at all.** The condition it would be
   testing for is "is a ringer set in the model under `untreatedZero`", which is known for free at
   `ado:1234`, before the regression. Detecting after the fit what you could have prevented before
   it is strictly worse.
3. **It fires after the expensive estimation**, so it cannot save the runtime, only the
   conclusion.
4. **It is skipped under `safety(off)` — and three of the project's four `untreatedZero` suite
   configs use `safety(off)`.** A user following the package's own fixtures would not be covered.
5. **After A there is nothing left to detect** for this defect: `amedrop` falls to 0 on every
   `untreatedZero` arm. The two cases where it still fires after A are a genuinely unidentified
   cell and an over-parameterised lattice — i.e. it degenerates into a generic
   "some treated-cell dose slope is not identified" caution, which is a separate (and defensible)
   feature request, not this fix.

**If the author wants something here anyway**, the cheap, honest version is to sharpen the
EXISTING generic caution at `ado:11929`: when `contreat()` is in force and the dropped token is a
dose slope (`co.<contreat>`), say so by name and name the cell, instead of "a treatment period
coefficient interacted with something (likely, a term implementing a continuous treatment
variable ...)". That is a message change with no new predicate and no false-positive cost.

**Risk: LOW as a message change; MEDIUM-HIGH as an enforcement layer (measured false positives).**

### Option E — keep the ringer but move it to the END of the explicit-regressor list  **[built, run, DISPROVED]**

The idea: if the collinearity drop takes the *later* column, putting the ringer last would make
the engine drop the ringer instead of a treated-cell slope. `build_E` moves
`` `contreatcontrolinserter' `contreatcontrolinserterFE' `` after the treatment-cell block on both
explicit spellings (`ado:1485` and the ppmlhdfe rc-3200 fallback at `ado:1468`).

**It fails on both halves.**
- **It does not fix `untreatedZero` on the `regress` path.** `customvce(ols)`, `customvce(robust)`
  and `unconditionalse` all still return `main:ame = -0.275000000000` (the broken value) — the
  drop does NOT move when the ringer is one column against a twelve-column cell#dose set. Only the
  ppmlhdfe fallback arm moved.
- **It breaks `contreatcontrols(basic)`**, which is the whole installed base of the explicit
  paths: `basic customvce(ols)` goes -0.275000000000 -> 1.025000000000, `basic unconditionalse`
  the same, `basic poisson` -0.066666666667 -> 0.273333333333, `basic byTime customvce(ols)`
  -0.083333333333 -> 1.216666666667. With the ringer last, the pairwise collinearity between it
  and the pooled margins-visible `c.contreat` resolves the other way: the base dose slope stops
  being partialled out and every ame-family row is off by it. **That is exactly the defect SI-20
  fixed on 2026-09-05** (`sthlp:933`, `sthlp:1075-1078`).

So E reintroduces a shipped-and-closed defect while fixing nothing. It does, however,
**independently confirm the DF-1 comment's claim** that the `i.`null'##` spelling is load-bearing
for the collinearity drop order on the explicit paths — that ruling is now verified in code AND in
measurement.

**Risk: HIGH. Do not ship.**

### Option G (companion, not an alternative) — assert `untreatedZero`'s premise, cheaply

`untreatedZero` is an unverified user assertion. `sthlp:996` only says "If this option is
specified when the always 0 condition is not met, wooldid will return poorly identified
estimates." **Measured how poorly**: on a panel where the dose is nonzero off treatment with a
true baseline dose effect of 0.6, the default `untreatedZero` path returns
`main:ame = 1.625000000000` where the within-cell truth is `1.025000000000` — off by exactly the
0.6 baseline slope, silently.

Two `count`s over the estimation sample, no sort, no `gegen`, O(N) and free even on large data:

```stata
    * untreatedZero premise screen: the dose must be zero outside treated cohort-periods
    qui count if `contreat' != 0 & `contreat' < . & `ttre' >= . & `ifstatement'
    local uzbadC = r(N)
    qui count if `contreat' != 0 & `contreat' < . & `ttre' < . & `t' < `ttre' & `ifstatement'
    local uzbadP = r(N)
    if `uzbadC' + `uzbadP' > 0 dis as text "Caution: contreatcontrols(untreatedZero) assumes `contreat' is zero outside treated cohort-periods, but `=`uzbadC'+`uzbadP'' estimation-sample observation(s) violate that (`uzbadC' never-treated, `uzbadP' pre-onset). The reported effects are then off by the baseline dose slope; use contreatcontrols(basic) or another control type instead."
```

Caution rather than refusal, and skipped under `safety(off)`, matching the package's convention.
Why it belongs with this fix rather than separately: the stray ringer was, on premise-violating
data, an *accidental* partial control. A removes it, which is right (it makes every engine agree
with the documented g() = 0 model — measured: after A, `untreatedZero customvce(ols)` on
premise-violating data returns the default path's `2.870500000000 / 1.625000000000` exactly), but
it also means the documented "poorly identified" warning is now the only thing standing between a
user and a silently wrong number.

**Risk: LOW.**
---

## 4. Comparison table and recommendation

| | fixes the defect | breaks anything | cost to the user | fires before the fit | risk | ship? |
|---|---|---|---|---|---|---|
| **A** gate exclusion at `ado:1234` | **yes, exactly** (18 of 66 quantities move, all of them `untreatedZero` x non-reghdfe) | no (65 of 66 non-defect quantities hex-identical; `e(N)`/`e(df_r)` unchanged everywhere) | none; also removes a spurious caution and a spurious ppmlhdfe refit | yes | **low** | **YES** |
| **B** clear inserters at `ado:1313` | yes, **byte-identical to A on all 66 rows** | no | none | yes | low | as a belt-and-braces addition to A, not alone |
| **C** refuse `untreatedZero` + engine | yes (by refusal) | no | **high** — kills `poisson`, `customvce()`, `unconditionalse`, and with them `bootstrap()`/`impvalues` under poisson | yes (parse time) | low code / medium user | no, except as a narrow companion for `contreatpoly(>1)` |
| **D** widen the colname screen | only as a post-hoc caution | measured false positives on ordinary `basic` runs (11/1) and `basic byTime` (7/5) | a caution the user cannot act on | **no** (after the fit; skipped under `safety(off)`) | medium-high | no; keep the *message* sharpening only |
| **E** reorder the ringer last | **no** (reg path unchanged at -0.275) | **yes — reintroduces the SI-20 defect** on every `basic` explicit path | severe | yes | high | **no** |
| **G** premise caution (companion) | n/a (different failure) | no | none | yes | low | yes, alongside A |

**Ship A, with B's three clearing lines added inside the `untreatedZero` branch as belt and
braces, plus G.** A is the smallest edit that removes the cause; B alone would leave the gate
reading as though `untreatedZero` still wanted a ringer set; C pays a large user cost to avoid a
two-token fix; D cannot be the enforcement layer and has nothing left to detect once A ships; E is
a disproof.

**I do not disagree with the narrow fix.** What I would add to it is: (i) the invariant should be
written at BOTH sites, because the two halves of this bug are 79 lines apart and each looks
locally innocent; (ii) the fix is *necessary but not sufficient* — see §5.

---

## 5. What A does NOT fix (found while testing it; report as separate tickets)

### 5.1 `contreatpoly(>1)` + `untreatedZero` reintroduces a pooled dose main effect — `g()` is not 0

`contreatpolypart` is `c.<dose>##c.<dose>` at `contreatpoly(2)`. Concatenated into the
`untreatedZero` injection, the model line reads
`i.`_i'#i.`_t'#i.`targetflag'#c.dose##c.dose`, and `##` binds at the LOWEST precedence, so
Stata expands it as `(cells#c.dose) ## (c.dose)`. **Measured with `fvexpand`: the expansion
contains exactly one bare pooled `c.dose` term.** `untreatedZero` therefore does NOT deliver
`g() = 0` (`sthlp:1033`) whenever `contreatpoly() > 1`.

Measured on an exact-fit quadratic design (analytic truth `att 2.726766666667`,
`ame 1.831666666667`):

| run | `main:att` | `main:ame` |
|---|---|---|
| default engine | **2.726766666667** (exact) | **1.831666666667** (exact) |
| `customvce(ols)`, pinned build | 0.954666666668 | -0.275000000000 |
| `customvce(ols)`, **after A** | 1.415933333337 | **1.831666666667** (exact) |

**Why the ame is exact and the att is not.** With the pooled `c.dose` column surviving on the
reg path and a cell slope dropped instead, `margins` SEES the pooled column, so its `dydx` read
sums `beta_pooled + (beta_cell - beta_pooled)` and the ame is invariant to the reparameterization
- hence the exact 1.831666666667. The att counterfactual instead carries `beta_pooled * dose` on
one side only, so it does not cancel and the att is understated by roughly
`beta_pooled * mean(dose | treated)`: `2.726766666667 - 1.415933333337 = 1.310833`, against
D1's dropped slope 1.3 times a mean treated dose of about 1.008 - the same signature as CRIT-2's
att damage.
Note also that the DEFAULT path's exactness at `contreatpoly(>1)` rests on reghdfe happening to
drop the pooled column rather than a cell slope. That is precisely the solver- and order-dependent
fragility SI-20 removed for `interactivecontrols()` (sthlp:1075-1078), so the `#`-only re-spelling
below fixes a fragility on ALL FOUR paths, not only the explicit ones.

So A fixes the ame exactly and leaves the att path-dependent. **Proposed follow-up:** under
`untreatedZero`, build the injection without any `##` — e.g. emit the polynomial as explicit
`#`-joined terms, `i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'#c.x`,
`...#c.x#c.x`, ... , plus the bare cell term — so that no pooled main effect is created at any
polynomial order. Until that is done, the narrow companion refusal (option C, restricted to
`contreatpoly(>1) & untreatedZero & non-reghdfe engine`) is defensible.

### 5.2 `untreatedZero` + a fine `lattice()` still diverges between engines after A

`lattice(all 2)` agrees exactly across engines after A (`2.265500000000 / 1.025000000000`
everywhere). `lattice(all 3)` and `lattice(all 5)` do not: after A, `lattice(all 5)` on a noisy
panel gives `main:ame` 1.0863082 on the default engine and 1.0400148 under `customvce(ols)`, with
the explicit path still dropping one dose slope (`amelive = 23, amedrop = 1`). This is the
over-parameterised-lattice family already open in the project's §7 record ("`lattice(all 3)`
returning an all-zero VCV ... classified undetermined"), not something A causes or should fix;
A strictly reduces the divergence (`lattice(all 3)` noisy: base -0.2078 vs default 1.5224;
after A 0.9667 vs 1.5224).

### 5.3 One treated cell with no within-cell dose variation: default refuses, explicit path proceeds

On a design where one treated cohort-period cell has a constant dose, the **default** path dies
`r(322)` ("a treated cohort-period dummy coefficient was dropped for multicollinearity") while
`customvce(ols)` returns a number in **both** the pinned build and under A (ame -0.183333333333
before, 0.916666666667 after). Pre-existing, orthogonal to this fix, worth its own line: the same
data is fatal on one engine and silent on another.

---

## 6. How to prove A did not regress anything

1. **The evidence I already have.** `probe_base.do` vs `probe_A.do` (and `probe_B.do`) on 24
   specifications: 48 of 66 recorded quantities byte-identical in hex, and the 18 that move are
   all `untreatedZero` x non-reghdfe. `e(N)` and `e(df_r)` identical on every arm. Re-runnable:
   `stata-mp -e do probe_A.do` in the workdir; `res_A.csv` is the artefact.
2. **The project battery.** Expect **zero** baseline movement, because `test/suite/configs.csv`
   and `configs_ext.csv` contain no `untreatedZero` x `customvce`/`unconditionalse`/`poisson`
   cell (verified by grep). If a battery config does move, the fix is wrong and should be
   stopped — that is the sharpest available falsifier.
3. **New cert legs to add with the fix** (they are the reason the defect survived three years):
   - `cont_uz_cv` : `contreat(x) contreatcontrols(untreatedZero) customvce(ols)`
   - `cont_uz_unc`: `... untreatedZero unconditionalse robust`
   - `cont_uz_pois`: `... untreatedZero poisson`
   and, as the load-bearing assertion, a **four-engine agreement test** on a fixture whose dose is
   exactly zero off treated cohort-periods: default, `customvce(ols)`, `unconditionalse robust`
   and (on a count fixture) `poisson` must return the same `main:att`/`main:ame` to rtol 1e-10.
   That assertion is what no existing config makes, and it fails loudly on the pinned build.
4. **A ringer-presence assertion**, cheap and direct: under `verbose`, the echoed regression line
   for an `untreatedZero` run must contain no `i.<null>##c.(` term on any of the four engines.

---

## 7. Documentation that must ride with the fix

- **`sthlp:993-996`** (the `untreatedZero` synopt): add a dated results-changed note in the exact
  style of `sthlp:933` and `sthlp:1075-1078`: before this date, `untreatedZero` combined with
  `unconditionalse`, `customvce()` or `poisson` silently estimated the
  `contreatcontrols(basic)` model instead — every att, ame and event-study row was off by the
  dropped cell's dose slope, and under `poisson` it happened on the default engine with no opt-in.
  Also state plainly that `untreatedZero` sets `g() = 0` **on every estimation path**.
- **`sthlp:1033-1035`** (the g() = 0 paragraph): note the `contreatpoly()` caveat of §5.1 once it
  is settled.
- **`review_process/cell_collapse_contreat_2026-09-02.md`** recorded the opposite conclusion —
  that `untreatedZero` "puts the ringer nowhere". It needs a correction line, because it is the
  document a future reviewer will find first.
- The AMESE adversarial verifier's **P4** probed `untreatedZero` on the default engine only; its
  record should say so, or the same gap regrows.

---

## 8. Honest limits of this work

- All estimation runs are on ONE panel shape (N = 160, 4 cohorts, 8 periods, 5 units per cohort)
  plus three variants of it. The defect's mechanism is exact algebra and does not depend on size,
  but I did not check a large or an unbalanced panel.
- I did not run `bootstrap()`, `impvalues` or `subgroups()` with `untreatedZero`. Their point
  estimates go through the same worker, so A moves them the same way, but I have not measured it.
- I did not test `untreatedZero` under weights (`aweight`/`pweight`/`fweight`).
- Option D's false-positive measurements are from the token counts my instrumented build printed,
  not from a shipped screen; a real screen would need its own message and its own gate.
- §5.1's proposed re-spelling of the polynomial injection is designed but **not built or tested**.
