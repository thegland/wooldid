# uz-verify — the readout that would prove an `untreatedZero` fix correct

**Angle:** design the test, not more diagnosis. Nothing here should ship without this readout
agreed in advance, and the harness is written so the author can run it himself against any
candidate patch.

| deliverable | absolute path |
|---|---|
| the harness (runnable, build-parameterised) | `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/fixprop/uz_oracle_harness.do` |
| the reader for its CSVs (the three verdicts) | `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/fixprop/uz_compare.py` |
| probe workdir: logs, CSVs, NOTES.md | `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uz-verify/` |
| patched build A — gate exclusion at `ado:1234` | `.../probes/uz-verify/build_A/wooldid.ado` |
| patched build B — clear in the branch at `ado:1313` | `.../probes/uz-verify/build_B/wooldid.ado` |
| spelling + ppmlhdfe-fallback side probe | `.../probes/uz-verify/spell.do` |

Run it:

```stata
global WDBUILD "/path/to/dir/with/wooldid.ado"
global WDTAG   "fixB"
global WDOUT   "/path/uz_fixB.csv"
global WDDESIGNS "1 2 3 4 5 6 7 8 9 10"
do uz_oracle_harness.do
```
then `python3 uz_compare.py uz_unfixed.csv uz_fixB.csv`.

---

## 1. The three verdicts the harness produces

**(O) ORACLE.** Zero-noise designs where `y` is fitted EXACTLY by `i`-FE + `t`-FE (+ unit FE)
+ cell intercept + `beta_cell * dose`, with the dose exactly 0 outside treated cohort-periods.
wooldid's own model spans that, so the truth is a *weighted mean of the structural parameters* —
no regression, no external package:

```
truth(ame, weighting e) = weighted mean of beta over treated rows
truth(att, weighting e) = weighted mean of (kappa + beta*dose) over treated rows
every pre-treatment row  = EXACTLY 0
```

with the four estimand weights read off `ado:1594-1650` (`att`: w; `att_it`: w/cell-w;
`att_i`: w/post-w-in-i; `att_itime`: w/(cell-w x #post cells of i)), computed on wooldid's own
realized `e(sample)`.

**(X) CROSS-ENGINE IDENTITY.** Every engine of a family must return the SAME point estimate on
the SAME specification. This needs no oracle at all and is the sharpest test available: it is a
property of the *specification*, not of the data, so it holds on any design, and it is what the
defect violates. Threshold `rtol 1e-10`, the project's own indifference level — do not tighten it,
because plain `reg`/`poisson` clustered paths are not bit-reproducible run to run.

**(N) BIT IDENTITY / NEGATIVE CONTROLS.** Point estimates are dumped in `%21x`, so "did this move"
is string equality. `uz_compare.py` splits every moved cell into *allowed* (`untreatedZero` on an
engine the `ado:1234` gate reaches) and *must not move* (everything else: the default reghdfe path,
`robust`, `cluster()`, `ccnoabsorb`, every `contreatcontrols(basic)` cell, the whole binary path).
It also checks `e(N)`, `e(df_r)` and `e(r2)`/`e(ll)` did not move ANYWHERE, including on the
cells whose estimates are supposed to change — because the fitted span is unchanged by the fix
(the pooled column is identically the sum of the cell columns), only the labelling moves.

## 2. The designs, and why each earns its place

| # | design | what only it can catch |
|---|---|---|
| 1 | balanced, `i` = cohort group of 5 units, ES rows requested | the headline case; the ES and pre-treatment rows |
| 2 | design 1 + unit-level intercepts, run with `fe(i.unit)` | the author's own recommended `i`=group + `fe()` usage |
| 3 | ragged panel (rows dropped, >=5 kept per unit) | where cell weights stop being proportional |
| 4 | `[aw=w]` | weights enter the estimand weights and the truth |
| 5 | `[fw=w]` | the other weight algebra |
| 6 | design 1 run with `lattice(all 2)` | `ado:1240` injects a SECOND ringer twin `i.lattice##c.(ringer)`; a fix touching only the pooled head would leave it |
| 7 | poisson twin, HOMOGENEOUS semi-elasticity 0.25 | a convention-free target for the one family with no unbroken wooldid path to validate against |
| 8 | poisson twin, cell-varying gamma | the same mechanism with heterogeneity, so a dropped cell is visible |
| 9 | dose varies EVERYWHERE, `contreatcontrols(basic)` | proves a fix does not damage `basic` where `basic` is the RIGHT estimator — an oracle leg, not just a bit-control |
| 10 | binary mode, no `contreat()` | proves the edited block is unreachable from the binary path |

Engines (trimmed to one per behaviour class because eight Stata instances share the machine;
the four dropped — `cluster(i)`, `ccnoabsorb robust`, `customvce(cluster unit)`,
`unconditionalse cluster(unit)` — were measured on design 1 and behave exactly like their
class-mates, and the harness header says to restore them for a release gate):

* linear: `(default)`, `robust`, `cluster(unit)`, `ccnoabsorb`, `customvce(ols)`,
  `customvce(robust)`, `unconditionalse robust safety(off)`
* poisson: `poisson`, `poisson robust`, `poisson ccnoabsorb`, `poisson customvce(robust)`

`unconditionalse` needs `safety(off)` alongside non-`att` estimands (`ado:4617`) and is an outright
error without `robust`/`cluster()`/`customvce()`/pweights (`ado:4676`); `customvce()` is refused
with `robust`/`cluster()` (`ado:4667`). Those three constraints fix the spellings above.

`safety(off)` there cannot be confounding the result, and the data says so rather than the
argument: on the unfixed build `unconditionalse robust safety(off)` is hex-identical to
`customvce(ols)` and `customvce(robust)` on every broken row of design 1, and after the fix it
lands on the same truth to the same 5.33e-15. It gates screens, not estimation.

---

## 3. BASELINE ON THE CURRENT (UNFIXED) BUILD — THE FAILING MATRIX

Build `sha c1d8979`, one Stata invocation per design, CSVs in
`probes/uz-verify/out/uz_base_d<k>.csv`, logs `probes/uz-verify/base_d<k>.log`.
`pass`/`fail` count ROWS against the design truth.

| design | engine | pass | fail | worst rel. error |
|---|---|---:|---:|---|
| **1** base-balanced (20 rows: 8 main, 2 pre, 6 ES post, 4 ES pre) | `(default)`, `robust`, `cluster(unit)`, `ccnoabsorb` | 20 | 0 | **2.41e-15** |
| | `customvce(ols)`, `customvce(robust)`, `unconditionalse robust` | 6 | **14** | **8.67e-01** |
| **2** unit FE via `fe(i.unit)` (8 main rows) | the four correct engines | 8 | 0 | 5.75e-16 |
| | the three gate engines | 0 | **8** | **7.91e-01** |
| **3** ragged panel | the four correct engines | 8 | 0 | 1.10e-15 |
| | the three gate engines | 0 | **8** | **7.99e-01** |
| **4** `[aw=w]` | the four correct engines | 8 | 0 | 4.60e-16 |
| | the three gate engines | 0 | **8** | **7.91e-01** |
| **5** `[pw=w]` | the four correct engines | 8 | 0 | — |
| | `customvce(robust)`, `unconditionalse robust` | 0 | **8** | **7.91e-01** |
| | `customvce(ols)` | — | — | r(198): `reg [pw=w], vce(ols)` is not legal |
| **6** `lattice(all 2)` | the four correct engines | 8 | 0 | 1.59e-15 |
| | the three gate engines | 0 | **8** | **8.48e-01** |
| **7** poisson, homogeneous gamma = 0.25 | `poisson`, `poisson robust`, `poisson ccnoabsorb`, `poisson customvce(robust)` | **0** | **8** | **1.00e+00** |
| **8** poisson, cell-varying gamma | all four poisson engines | **0** | **8** | **7.91e-01** |
| **9** dose everywhere, `contreatcontrols(basic)` | ALL SEVEN engines | 8 | 0 | 7.32e-15 |
| **10** binary, no `contreat()` | ALL SEVEN engines | 4 | 0 | 5.67e-15 |

### What the matrix says

**1. The predicted partition is exactly right, and it is a partition of ENGINES, not of data.**
On every linear design the same four engines are exact to ~1e-15 and the same three are wrong by
80-87%. The three are exactly the arms of the `ado:1234` disjunction that an `untreatedZero` run
can reach (`customvce()` in any spelling, `unconditionalse`, and `poisson`); `ccnoabsorb` — which
also moves the FE set out of `absorb()` and into the regressor list, but is not in the gate — is
CORRECT. That is the cleanest possible evidence that the gate, not the engine, is the cause.

**2. The mechanism signature is an exact identity, and it is weighting-free.** On design 1,
(correct − broken) is the SAME CONSTANT on all seven ame-family rows:

```
main:ame        0.720000000000181      ES_ame:post0  0.720000000000200
main:ame_i      0.720000000000174      ES_ame:post1  0.720000000000200
main:ame_it     0.720000000000180      ES_ame:post2  0.720000000000122
main:ame_itime  0.720000000000174
```

and the design's `beta` at cell (i = 3, relyr = 2) is exactly **0.72** — agreement 2.5e-13.
The same constant 0.72 appears on designs 2 (`fe()`), 3 (ragged), 4 (`[aw=w]`) and 5 (`[pw=w]`) —
weights and raggedness do not touch it, because the shift is a regression-coefficient
relabelling, not an averaging artefact. On design 8 (poisson, cell-varying) the constant is
exactly **0.18**, which is that design's `gamma` at the same cell (`beta/4 = 0.72/4`). The
`att` rows shift by `-0.72 x (that row's weighted mean dose)` instead, so they all differ
(`main:att` 0.9234367868, `main:att_i` 0.9538434559, `ES_att:post0` 1.0154537005). A fix must
therefore move every ame row by exactly one constant and every att row by a row-specific amount;
the harness predicts both, so "did the fix do the right thing or merely a different thing" is
decidable.

**3. Every pre-treatment row is EXACTLY zero on the broken engines too.** Six of design 1's twenty
rows pass on the broken engines, and they are the six placebo rows. The part of the table a
careful reader checks first is precisely the part the defect does not touch.

**4. THE POISSON CASE IS WORSE THAN THE WRITE-UP STATES, AND ONLY THE ORACLE CATCHES IT.**
On design 7 the true semi-elasticity is a homogeneous **+0.25** everywhere. Every poisson engine,
default included, returns

```
main:ame  = -0.0000000000     (truth 0.25)   — EXACTLY ZERO
main:att  =  0.3000000000     (truth 0.6206) — EXACTLY the cell intercept kappa
```

and `contreatcontrols(basic) poisson` returns the same numbers to the bit. The reason it is exactly
zero is the mechanism: with a homogeneous slope the dropped cell's slope IS 0.25, so every surviving
cell slope becomes `0.25 − 0.25 = 0`. A constant dose effect is reported as a precise null.

The second half of that finding is a fact about TESTING, not about wooldid: **cross-engine identity
finds nothing on design 7**, because all four poisson engines are broken identically and therefore
agree with each other. `uz_compare.py` prints "no split anywhere" for design 7. Only the exact-fit
oracle catches the poisson family. Any fix-verification plan that leans on cross-engine agreement
alone would pass this build.

(It also settles a convention: the broken poisson `att` is exactly `kappa` = 0.30, which it could
not be if wooldid reported on the exp scale. So wooldid's poisson `att` is the xb/log-scale
quantity and the harness's truth formula for it is right. The patched build confirms it directly.)

**5. The lattice arm needs the fix too.** Design 6 breaks on the same three engines, and its `att`
rows break HARDER than without a lattice (`main:att` 0.4193 against a truth of 2.0527, relspread
0.796 vs 0.450 on design 1), consistent with the second ringer twin `i.lattice##c.(ringer)` at
`ado:1240` costing a further dropped column per lattice level. A fix that clears only the pooled
head would leave design 6 broken. Both candidate fixes clear the whole block, so both cover it —
but this is the leg that would catch a narrower fix.

**6. The negative controls already pass, so a fix has something to break.** Design 9
(`contreatcontrols(basic)` on data where the dose varies everywhere — the world `basic` is FOR)
is exact to 7.3e-15 on ALL SEVEN engines including the three gate engines. Design 10 (binary, no
`contreat()`) is exact on all seven. Any candidate fix must leave both bit-identical.

**7. Design 2 is the author's own recommended usage.** `i` = cohort group with unit fixed effects
supplied through `fe()` is what `sthlp` guidance and the SI-8 caution steer users toward, and it is
broken on the three gate engines exactly like the rest.


---

## 4. Fix proposals

All line numbers are the pinned build (`sha c1d8979`), `wooldidinstall/wooldid.ado`.

### FIX-B — RECOMMENDED. Clear both inserter spellings inside the `untreatedZero` branch

**Site:** `ado:1313-1315`. **Becomes:**

```stata
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
      * "not having contreatcontrols of any kind" is what the comment above promises, so
      * clear BOTH spellings of the ringer control set here, in the same if that rewrites
      * contreatinj. The gate at :1234 builds them whenever the engine leaves reghdfe
      * (unconditionalse / customvce() / poisson), and under untreatedZero the dose is
      * exactly 0 off treated cohort-periods, so the pooled i.null##c.(ringer) column
      * EQUALS the sum of the cell#dose columns: an explicit-regressor path drops one cell
      * slope and every reported quantity comes back as truth - (slope of the dropped cell).
      * Nothing is given up - under untreatedZero contreatinj carries no pooled dose main
      * effect (it is spelled #, not ##) for the ringer head to partial out.
      local contreatcontrolinserterFE
      local contreatcontrolinserterABS
      local contreatcontrolinserter
    }
```

**Why this site.** It is the same `if` that already rewrites `contreatinj`, it sits AFTER the two
assignments at `:1300`/`:1302` and BEFORE the `ccnoabsorb` copy at `:1318-1325`, so all three
macros the four engine calls read (`:1461/1462/1468`, `:1475/1476`, `:1484/1485`) are cleared
whatever the gate above did. It cannot drift apart from the `contreatinj` rewrite under a future
edit to the gate — which is exactly how this defect arose: `fc17051` (2023-08-15) shipped the
disjunction and nothing ever re-checked it against the branch three dozen lines below. It also
covers the lattice twin at `:1240` with no extra code.

**Blast radius.** `contreatcontrolinserterFE/ABS/inserter` are read at `:1461, :1462, :1468, :1475,
:1476, :1484, :1485` and nowhere else (verified by grep over the whole ado). The `ringer` variable
itself is still created at `:1104-1105`, so the cell collapse's `contreat(`contreat' `ringer')` at
`:1900` is untouched, and so is the DF-1 absorb-thinning logic for every other control type.
Nothing outside `untreatedZero` is reachable from this edit, because the enclosing test is equality
with the single string `untreatedZero` and `wooldid_validate` (`:4796-4800`) already refuses to
pair it with anything.

**What could regress.** (i) If any real design needed the pooled head under `untreatedZero` to
partial out a base dose level, removing it would break that design. It cannot: under
`untreatedZero`, `contreatinj` is spelled with `#` and carries NO pooled `c.dose` main effect at
all, so there is nothing to partial out — and the DEFAULT reghdfe path has never carried the head
and is the path the project has always called correct. (ii) The residual-df bookkeeping: under
`untreatedZero` the pooled column is exactly collinear, so the engine was already dropping one
column, and `e(df_r)`/`e(N)`/`e(r2)` should NOT move. That is checked, not assumed — see the
`meta_identity` leg of `uz_compare.py`. (iii) Poisson: clearing `contreatcontrolinserterABS`
removes the single-level `i.null` from `absorb()`, so ppmlhdfe's first attempt at `:1462` should
now SUCCEED and the rc-3200 fallback at `:1464-1468` should stop firing — a behaviour change the
user sees (four caution lines disappear). Intended, and it must be confirmed rather than assumed.

### FIX-A — equivalent, keep as the named alternative. Exclude `untreatedZero` from the gate

**Site:** `ado:1234`. **Becomes:**

```stata
    if "`contreatcontrols'" != "untreatedZero" & (strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson") {
```

**For A:** the `null` tempvar and the head are never built at all. `` `null' `` is referenced only
at `:1236-1238` (verified by grep), so nothing downstream can miss it.
**Against A:** the gate's four arms are about ENGINES; hanging a control-type exclusion off it
states in a second place a fact `:1313` already knows, which is the shape that produced the defect.

**A and B are NOT trivially identical and the harness is what settles it.** Under B the gate still
fires, so `tempvar null` and `qui gen `null' = 1` still execute: a variable is created that no
model uses, and — more to the point — every tempvar allocated after it gets a different internal
name than under A. The project's own record (NFFAST-02) says `margins`-inherited row stripes embed
wooldid tempvar names, which is why stable keys were introduced. If any surviving output is
name-sensitive, A and B will differ. **Test:** run the harness under both and require every `hex`
cell to be string-equal. If they are equal, ship B (it is the more robust site) and note A as
equivalent; if they are not, the difference is itself a finding.

### FIX-C — MANDATORY whatever code fix ships: correct the records that caused this

1. `ado:1225-1233`. The offending clause is the one that justifies the engine arms of the gate:
   *"and also initializes it this way if there is any chance that this FE object might get kicked
   outside of an HDFE absorb, in which case this section being included is necessary in order to
   correctly manipulate the behavior of reg and poisson's colinearity drop rules"* — together with
   the sentence below it about "ensur[ing] that the base level of contreat gets partialled out".
   Neither qualifies by control type, and under `untreatedZero` there is no base level of
   `contreat` in the model to partial out: `contreatinj` is spelled with `#`, not `##`, so it
   carries no pooled `c.dose` main effect. The comment must say that the manipulation applies only
   to control types that actually put a pooled dose column in the model.
2. `review_process/cell_collapse_contreat_2026-09-02.md:112-120` contradicts itself in adjacent
   sentences, and the second one is the record that matters. It first gets the mechanism right —
   the FE inserter is built "under `contreatcontrols(basic)` without a lattice, under
   `byTreated`/`byTtre`/`byCohort`/`byTime`, **and unconditionally under `unconditionalse`,
   `customvce()` or `poisson`**" — and then, at **line 118-120**, writes:

   > "`contreatcontrols(untreatedZero)` **alone** puts the ringer nowhere - it only rewrites
   > `contreatinj` (1152-1154) - so the prototype's own comment, which names `untreatedZero` as a
   > ringer path, is wrong on that point (section 7)."

   The word doing the damage is **"alone"**: `untreatedZero` puts the ringer nowhere *by itself*,
   but the three engine arms listed one sentence earlier put it there regardless of the control
   type, and that combination is never considered. The record also over-rules the prototype
   comment that had it RIGHT. Both the sentence and the "(section 7)" ruling it supports need
   correcting, or this is a standing instruction to re-introduce the defect.
   Note the line numbers there are the 2026-09-02 prototype's (1123-1141 / 1152-1154), which are
   `ado:1234-1302` / `ado:1313-1315` on the pinned build.
3. `sthlp:1033-1035` describes `untreatedZero`'s model without saying it is engine-invariant.
   One sentence: whichever engine the other options select, `untreatedZero` inserts no dose
   control of any kind. `sthlp:1086` already names the refusal to pair it with `basic`; that
   refusal and the silent acquisition of `basic`'s own answer are the same fact and should be
   stated together.

### FIX-D — REJECTED alternative, named so the adversarial review can dismiss it explicitly

Reorder: move `contreatcontrolinserterFE` to AFTER `contreatinj` at `:1485` (and at `:1468`) so
Stata's "collinear terms get dropped in the order in which they appear" rule drops the redundant
pooled column rather than a cell slope.
**Rejected:** `ado:1454` hedges that rule with "(mostly)" in the author's own words; it would leave
an omitted column in `e(b)` that every downstream name-walk (`:1671`'s regexm over `bvars`, the cfx
key builders) has to keep ignoring; and it fixes nothing about the `i.null`-in-`absorb()` poisson
fallback. It trades a wrong number for a fragile one.

### FIX-E — REJECTED alternative: refuse the combination

Error out on `untreatedZero` + (`customvce()` | `unconditionalse` | `poisson`).
**Rejected:** `untreatedZero` + `poisson` is a natural and documented pairing (a count outcome
whose dose is undefined outside treated cohort-periods), and the help already advertises
`untreatedZero` in an example with plotting options (`sthlp:1551`). Refusing removes a capability
instead of restoring one. Hold it only as the fallback if the author will not move any numbers.

### FIX-F — ship a cross-engine certification leg alongside the code fix

The cheapest permanent guard is not a runtime check but a test: on ONE small exact-fit design,
assert every engine returns the same point estimate under `untreatedZero`. That is seven wooldid
calls and it would have caught this in 2023. `uz_oracle_harness.do` with
`global WDDESIGNS "1"` is that leg, already written; the (X) verdict needs no oracle, no baseline
file and no external package, so it never goes stale. The (O) oracle leg is the stronger
statement and costs nothing more on the same design.

### ADJACENT — a separate finding, not part of this fix

Under `contreatcontrols(basic)` + `poisson`, `contreatcontrolinserterABS` is
`i.`null'##c.(ringer)` with a single-level `i.null` inside `absorb()` at `:1462`, so ppmlhdfe's
first attempt fails, the rc-3200 branch at `:1464-1468` prints four caution lines the user cannot
act on, and the whole estimation is run a second time. Under `basic` the fallback's answer is the
one `basic` wants, so this is a cost and an alarming message, not a wrong number — but it fires on
every such run, with no opt-out. Fixing `untreatedZero` removes it there and leaves it under
`basic`. Measured by `probes/uz-verify/spell.do` part 2.

---

## 5. THE TWO CANDIDATE FIXES, RUN THROUGH THE HARNESS

Both patches were applied to **private copies** of the ado (the repo was not touched):
`probes/uz-verify/build_A/wooldid.ado` and `probes/uz-verify/build_B/wooldid.ado`.
Each build was run through the same ten designs and read with the same `uz_compare.py`.

### FIX-B — full pass on every leg

```
(O) ORACLE   : no engine fails any row on any design.  Worst relative error anywhere
               5.33e-15 (linear), 7.36e-11 (poisson customvce, MLE convergence).
(X) CROSS-ENGINE : "no split anywhere: every engine agrees within rtol on every row."
(N) BIT IDENTITY unfixed -> fixB, over the merged 10-design matrix:
       rows compared            : 732
       moved, ALLOWED to move   : 224
       moved, MUST NOT MOVE     :   0
       present in base, missing :   0
       of the allowed moves, 224/224 now hit the design truth.
       fit-invariants (e(N), e(df_r), e(r2)/e(ll)) that moved : 0
```

Read that last line again, because it is the claim the fix rests on and it was measured, not
assumed: **not one of `e(N)`, `e(df_r)`, `e(r2)`, `e(ll)` moves anywhere**, including on the cells
whose point estimates change by 80%. That is exactly what the mechanism predicts — the pooled
column is identically the sum of the cell columns, so the fitted span, the residuals and the rank
were always the same; only which column carried the slope changed. A fix that moved the residual
count would be doing something else, and this readout would say so.

Per-design highlights:

| design | before | after |
|---|---|---|
| 1 base-balanced | 3 engines x 14 rows wrong (rel. err 0.867) | all 7 engines 20/20, worst 5.33e-15 |
| 6 `lattice(all 2)` | 3 engines x 8 rows wrong (0.848) | all 7 engines 8/8, worst 2.44e-15 |
| 7 poisson homogeneous | **all 4** poisson engines 0/8, `ame` exactly 0 against a truth of 0.25 | all 4 engines 8/8, worst 1.22e-15 |
| 9 `basic`, dose everywhere | already exact on all 7 | **bit-identical** |
| 10 binary | already exact on all 7 | **bit-identical** |

### FIX-A — run on designs 1, 7, 6, 9 for the A-vs-B equivalence question

FIX-A was run on designs 1, 7, 6 and 9 and compared to FIX-B cell by cell:

```
rows in A: 456   in B: 456   common: 456
hex differences A vs B            : 0
META (e(N), e(df_r), e(r2)/e(ll)) : 0
```

**A and B are bit-identical across 456 rows on four designs, four engine families and both
control types.** The concern that B leaves `tempvar null` allocated (shifting every later
tempvar's internal name) does not reach any output. So the choice between them is a maintenance
judgement, not a numerical one, and the recommendation stands: **ship B**, because it lives in the
same `if` as the `contreatinj` rewrite and cannot drift from it; keep A as the documented
equivalent.

FIX-A also passes the oracle on its own: design 7 (`poisson`, homogeneous gamma) returns
`att = 0.6206377732`, `ame = 0.250000` on all four poisson engines against a truth of
`0.6206377732 / 0.25`.

### The spelling question the fix should close

`wooldid_validate` accepts `contreatcontrols()` by looping `foreach word in ...` and counting
words, so surrounding whitespace survives validation; the branch at `ado:1313` tests
`"`contreatcontrols'" == "untreatedZero"` exactly. `probes/uz-verify/spell.do` runs the default
engine under `untreatedZero`, ` untreatedZero`, `untreatedZero `, `  untreatedZero  `,
`untreatedzero` and `basic` and prints the `ame` each returns.

**Result: not a defect — every whitespace variant returns the identical hex.**

```
SPELL [untreatedZero]        rc=0  ame=0.966153846154  hex +1.eeabb78845512X-001
SPELL [ untreatedZero]       rc=0  ame=0.966153846154  hex +1.eeabb78845512X-001
SPELL [untreatedZero ]       rc=0  ame=0.966153846154  hex +1.eeabb78845512X-001
SPELL [  untreatedZero  ]    rc=0  ame=0.966153846154  hex +1.eeabb78845512X-001
SPELL [untreatedzero]        rc=198                      (correctly refused: case-sensitive)
SPELL [basic]                rc=0  ame=0.246153846154  hex +1.f81f81f81f84dX-003
```

Stata's `syntax` trims the option text before `ado:1313` sees it, so the exact-equality test is
safe. The leg is kept in the harness as a standing guard rather than as a finding. The same run
under FIX-B returns the identical hexes, so the fix does not disturb the default path.

**One correction to the CRIT-2/CRIT-3 wording, which an adversarial reviewer would otherwise
find.** Those write-ups say the broken number is *exactly* `contreatcontrols(basic)`'s answer. On
design 1 it is the same answer to **1.8e-12 relative at worst**, but only 3 of 20 rows are
bit-identical (the three that are exactly zero). It is the same estimator reached through a
different column arrangement, so the agreement is floating-point, not bitwise. The substantive
claim survives intact; the word "exactly" does not, and the report should say "to 1e-12".

**The ppmlhdfe rc-3200 fallback: confirmed, and the fix removes it under `untreatedZero`.**
`spell.do` part 2 runs `poisson` non-quietly and greps for the caution at `ado:1464`:

| build | `untreatedZero poisson` | `basic poisson` |
|---|---|---|
| unfixed | fallback FIRES; `ame` = -0.000000000000, `att` = 0.300000000000 | fallback FIRES; `ame` = 0.000000000000, `att` = 0.300000000000 |
| FIX-B | **no fallback**; `ame` = 0.250000000000, `att` = 0.620637773194 | fallback still FIRES; unchanged |

So the single-level `i.null` in `absorb()` is confirmed as the fallback's cause, the fix removes it
where `untreatedZero` is concerned, and `basic` is left exactly as it was. Wall-clock for the
whole design-7 leg went 61s (unfixed) -> 11s (FIX-B) -> 10s (FIX-A), but the machine was under
varying load from other agents and the `basic` half still falls back in every build, so treat the
timing as suggestive and the sentinel as the evidence.

Whichever fix ships, it should use the SAME predicate as `ado:1313` (that is what FIX-A and FIX-B
both do), so the fix and the branch can never disagree about what counts as `untreatedZero`.

---

## 6. What this harness does NOT cover (be honest about the gaps)

* **Standard errors, p-values and intervals are not checked at all.** The harness reads column 1
  of `e(wdidmainresults)` only. The defect moves point estimates; whether the SEs on the now-correct
  broken paths are also right is a separate question, and the exact-fit design cannot answer it
  (zero noise means every analytic SE is degenerate). A coverage simulation is the instrument for
  that and it is not here.
* **`subgroups()`, `interactivecontrols()`, `controls()`, `timetrends`, `coarsecohortcontrols()`,
  `ccc_absorb()`, `oostreatedcontrols()`, `anticipation()`, `esfixedbaseperiod`, `contreatpoly()>1`,
  `semielasticity`, `contreatwithin`, `impvalues`, `bootstrap()`, `permute()`** are all absent.
  Several of them add their own terms to the same regression and could interact with this block.
  `contreatpoly()>1` in particular puts more than one ringer column in `ringerpoly` and deserves
  its own leg.
* **The cfx materials** (`makecfxplots`, `e(cfxame_*)`) and `e(b)`/`e(V)` are not read. If a fix
  changed column names, the point-estimate matrix would not notice.
* **Design 5 is baselined under `[pw=w]`** and is in the merged matrix above. It is worth saying
  why it changed: as first written it used `[fw=w]` and got `r(198)`, because `ado:628-629`
  refuses fweights from users by design and demotes them to aweights only inside a bootstrap
  replication. Correct behaviour, not a finding. Under `[pw=w]`, `customvce(ols)` returns `r(198)`
  — `reg [pw=w], vce(ols)` is not a legal combination — so design 5 exercises six engines, not
  seven.
* **Only two cohort-count/period-count shapes** (30 units/8 periods, and its ragged variant) were
  used. The dropped cell was always the last in factor order; a design where the last cell is
  empty, or where `emptycellsoverride` is needed, is untested.
* **One seed.** The designs are deterministic given `set seed 20260918`; the oracle is exact so
  re-seeding should change nothing structural, but it was not tried.

* **Exact-fit designs cannot test standard errors, and the harness would be misread if it tried.**
  A side observation from `spell.do`: `wooldid yp i t ttre, att poisson` on a zero-noise design
  returns an all-zero VCV and wooldid's degenerate-covariance caution. That is correct behaviour on
  noiseless data, not a defect — but it is why the oracle leg reads column 1 only, and why a fix's
  effect on inference needs a separate, noisy simulation.
* **Design 5 under `[pw=w]`: `customvce(ols)` returns r(198)**, because `reg [pw=w], vce(ols)` is
  not a legal combination. Engine/weight incompatibility, not a defect; the other six engines show
  the usual partition and the usual 0.72 shift.

---

## 7. HOW TO USE THIS AS THE SHIP GATE

A candidate patch is accepted iff, over the ten designs:

1. `(O)` no engine fails any row. Linear rtol 1e-9, poisson 1e-7.
2. `(X)` `uz_compare.py` prints *"no split anywhere"*.
3. `(N)` `moved, MUST NOT MOVE : 0` and `present in base, missing : 0`.
4. `(N)` `fit-invariants (N, df_r, r2/ll) that moved : 0`.
5. every allowed move lands on the truth (`of the allowed moves, k/k now hit the design truth`).
6. the ame-family moves equal one constant per design, and that constant matches a design cell's
   `beta` to ~1e-12 (this is what distinguishes "the fix restored the right slope" from "the fix
   changed the answer").
7. `spell.do` part 2 shows the rc-3200 caution gone under `untreatedZero` and still present under
   `basic`.

FIX-B meets 1-7. FIX-A meets 1-6 on the four designs it was run on and is bit-identical to FIX-B
there.

## 8. Files a reviewer should open first

| what | path |
|---|---|
| the full unfixed readout (211 lines, the failing matrix) | `probes/uz-verify/out/READOUT_unfixed.txt` |
| the fixB readout, including the three identity blocks | `probes/uz-verify/out/READOUT_fixB.txt` |
| merged CSVs | `probes/uz-verify/out/uz_unfixed_ALL.csv`, `uz_fixB_ALL.csv` |
| per-design CSVs and Stata logs | `probes/uz-verify/out/uz_{base,fixB,fixA}_d<k>.csv`, `probes/uz-verify/{base,fixB,fixA}_d<k>.log` |
| spelling + ppmlhdfe sentinel, unfixed and fixB | `probes/uz-verify/spell.log`, `probes/uz-verify/spellB.log` |
| running state and every intermediate result | `probes/uz-verify/NOTES.md` |
