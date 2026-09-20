# uz-minfix — the narrowest fix for CRIT-2 / CRIT-3, written and tested

**Agent:** uz-minfix (angle: design and test the NARROWEST fix).
**Workdir:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uz-minfix`
**Patched build (mine, never the repo):** `.../probes/uz-minfix/buildfix/wooldid.ado`
**The repo was not written to.** Everything below was measured on my own copy.

## VERDICT IN ONE PARAGRAPH

The narrow fix — excluding `untreatedZero` from the `ado:1234` gate — **works, and it is
provably inert everywhere else**. It restores the design truth on `customvce(ols|robust|cluster)`,
`unconditionalse` and all five poisson arms, on three designs (two exact-fit, one noisy); it
leaves the default reghdfe path and every other `contreatcontrols()` value **bit-identical** in
a 74-config battery whose pinned-twice control is itself byte-identical; and as a free second
benefit the ppmlhdfe rc-3200 fallback **stops firing entirely** under `untreatedZero`
(2 occurrences pinned, 0 patched), with one r(3201) hard failure turning into a clean run.
**It is not sufficient on its own**: under `contreatpoly(>1)` the same collinear column arrives
through `contreatpolypart`'s `##` instead of the gate, and `att` (plus every att-family event-study
horizon, plus the pre-treatment `ame` row) stays at `contreatcontrols(basic)`'s answer. §7 gives a
second, equally narrow edit (form C) that closes that, tested, and measured inert at the default
`contreatpoly(1)`. **Ship A and C together.**

| what | evidence |
|---|---|
| truth restored on every repaired arm | §3, relerr <= 1.3e-15 on exact-fit designs |
| default path unmoved | §3, §4 — bit-identical, 2 designs, 4 engines |
| every other `contreatcontrols()` unmoved | §4 — 67 of 74 configs bit-identical |
| ppmlhdfe fallback no longer needed | §5 — 2 -> 0 occurrences |
| narrow fix NOT enough | §7 — `contreatpoly(>1)`; form C fixes it |
| form C validated against real dose curvature | §7b — `ame` 0.000129 -> 0.477908937767 (truth) |

---

## 1. The edit

One logical line, at `wooldidinstall/wooldid.ado:1234`.

BEFORE:
```stata
    if strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
```

AFTER:
```stata
    * UZ-FIX: exclude untreatedZero. Under untreatedZero `contreatinj' (ado:1313-1315) carries
    * NO margins-visible base dose term, so there is nothing for the drop-order manipulation
    * described above to protect; and the dose is exactly 0 off treated cohort-periods, so the
    * pooled i.`null'##c.(ringer) column is EXACTLY the sum of the cell#dose columns. Inserting
    * it makes the design rank-deficient by one on every explicit-regressor engine, the sweep
    * drops one cell dose slope, and every reported quantity returns as truth - beta(dropped).
    * Same predicate as the untreatedZero branch at ado:1313; wooldid_validate (ado:4797)
    * guarantees untreatedZero is the only word in contreatcontrols(), so it is a total test.
    if (strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson") & "`contreatcontrols'" != "untreatedZero" {
```

Two deliberate choices:
- **Parenthesise the existing disjunction and AND a single new conjunct.** Nothing inside the
  old condition is touched, so a reader can see the whole change is "…and not untreatedZero".
- **Use `== "untreatedZero"`, the same predicate as ado:1313**, not `strpos(...)`. The two sites
  must never disagree about what "this run is untreatedZero" means: if they ever did, the run
  would get the untreatedZero `contreatinj` without the matching control decision, or vice
  versa. `wooldid_validate` (ado:4792-4800) refuses `untreatedZero` with any other word
  (`wordcount > 1` -> r198), so exact equality is a total test of the case.

**Optional belt-and-braces (form B), not tested and not required.** Inside the untreatedZero
branch at ado:1313-1315, additionally clear the two macros:
```stata
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
      local contreatcontrolinserterFE
      local contreatcontrolinserterABS
    }
```
On this build form B is behaviourally identical to form A (I verified by reading that nothing
between 1234 and 1313 can populate those macros under untreatedZero other than the 1234 gate:
`ccringertail` is fed only by byTreated/byTtre/byCohort/byTime, all of which validate forbids
here). Its only value is that a future edit to the gate cannot re-open the hole. I did not
test it because testing two behaviourally identical patches would muddy the diff; if the
author wants both, A is the one that must be kept, B is free insurance.

---

## 2. Why this is the narrowest possible fix — the structural argument

`strpos("untreatedZero","basic") == 0`. So under `contreatcontrols(untreatedZero)` the gate at
1234 can fire **only** through the three engine disjuncts. Therefore:

> The blast radius of the edit is EXACTLY
> `contreatcontrols(untreatedZero)` x `{unconditionalse | customvce() | poisson}`
> and nothing else in the program can move.

In particular the **default reghdfe untreatedZero path cannot move**, because there the gate is
already false and `contreatcontrolinserterFE` / `...ABS` are already empty. And **no other
`contreatcontrols()` value can move**, because for them the new conjunct is identically true.

That is an argument, not a measurement; §4 is the measurement.

---

## 3. Truth recovery (design A, exact fit, 378 obs)

Design: 6 groups x 7 units x 9 periods; onsets 3/5/7 for groups 1/2/3, groups 4-6 never
treated; `dose = exp(N(0.2,0.55))` on treated cohort-periods and **exactly 0** elsewhere;
`beta_it = 1.45 - 0.31*i + 0.22*r - 0.06*r^2`; `y = a_i + g_t + kappa_it + beta_it*dose` with
no noise, so a correct run must return the design truth to machine precision.
(This is the `oob-untreatedZero` design; I reused it deliberately, because its pinned-build
numbers are already in the record and the patched numbers can be laid straight against them.)

Truths: `att` 2.117554176197, `ame` 0.907333333333. Poisson twin (exact-fit log-linear,
heterogeneous semi-elasticity `gam = beta/4`): `ame` 0.226833333333.
Broken-value predictions: `ame` 0.187333333333 = truth - beta(i=3,t=9) = truth - 0.72;
poisson `ame` 0.046833333333 = truth - gam(3,9) = truth - 0.18.

| arm | PINNED `att` | PINNED `ame` | FIXED `att` | FIXED `ame` |
|---|---|---|---|---|
| default (reghdfe)          | 2.117554176197 | 0.907333333333 | **bit-identical** | **bit-identical** |
| `robust`                   | 2.117554176197 | 0.907333333333 | **bit-identical** | **bit-identical** |
| `cluster(unit)`            | 2.117554176197 | 0.907333333333 | **bit-identical** | **bit-identical** |
| `ccnoabsorb`               | 2.117554176197 | 0.907333333333 | **bit-identical** | **bit-identical** |
| `ccnoabsorb robust`        | 2.117554176197 | 0.907333333333 | **bit-identical** | **bit-identical** |
| `customvce(ols)`           | 1.148673365330 | 0.187333333333 | 2.117554176197 | 0.907333333333 |
| `customvce(robust)`        | 1.148673365330 | 0.187333333333 | 2.117554176197 | 0.907333333333 |
| `customvce(cluster unit)`  | 1.148673365330 | 0.187333333333 | 2.117554176197 | 0.907333333333 |
| `unconditionalse robust`   | 1.148673365330 | 0.187333333333 | 2.117554176197 | 0.907333333333 |
| `poisson`                  | ame 0.046833333333 | | ame **0.226833333333** | |
| `poisson robust`           | ame 0.046833333333 | | ame **0.226833333333** | |
| `poisson ccnoabsorb`       | ame 0.046833333333 | | ame **0.226833333333** | |
| `poisson customvce(robust)`| ame 0.046833379047 | | ame **0.226833379057** | |
| `poisson unconditionalse robust` | ame 0.046833379047 | | ame **0.226833379057** | |

Relative error against truth on the repaired arms is <= 1.3e-15 on the linear arms. The two
`poisson customvce/unconditionalse` arms sit 2.0e-7 off truth on BOTH builds — that is the
pre-existing `poisson` vs `ppmlhdfe` convergence-tolerance gap at `regtol(12)`, not the fix
(the same 4.57e-8 absolute offset is present pre-fix, displaced by the same 0.18).

**Every reported quantity moves, not just `ame`.** The event-study rows move with it: on the
pinned build 18 of 25 point-estimate cells (`e(wdidmainresults)` col 1, `e(b_es_att)`,
`e(b_es_ame)`, `e(N)`) differ between `customvce(robust)` and the default arm; on the patched
build **0 of 25** differ beyond 1e-9 absolute.

**`e(N)`, `e(df_r)`, `e(r2)` do not move, before or after** (378 / 329 / 1 on every linear arm
of both builds). That is exactly why the defect is silent: the pooled ringer takes the rank
that the dropped cell slope gave up, so every diagnostic a user would check is unchanged.

---

## 4. Inertness: the blast-radius battery, measured

`p2_identity.do`: 7 `contreatcontrols()` values x 10 engines + 4 binary anchors = 74 runs per
build, on a design carrying **real noise** (`eps ~ N(0,0.3)`), dumping `e(wdidmainresults)`
(every cell) plus `e(N)`, `e(df_r)`, `e(r2)`, `e(rmse)`, `e(ll)` in `%21x`.

**Noise control first.** The pinned build was run TWICE. `dump_p2_pin1.txt` and
`dump_p2_pin2.txt` are **byte-identical**, so on this design the envelope the brief warns
about for clustered explicit paths is exactly zero and every pinned-vs-patched difference
below is the patch, not run-to-run drift.

| | result |
|---|---|
| **MOVED (6 of 74)** | `untreatedZero` x `customvce(robust)`, `customvce(cluster unit)`, `unconditionalse robust`, `poisson`, `poisson ccnoabsorb`, `poisson customvce(robust)` |
| **BIT-IDENTICAL (67 of 74)** | every `basic` / `byTreated` / `byTtre` / `byCohort` / `byTime` / `basic byTime` config on all ten engines; the binary (no-`contreat()`) anchor; and `untreatedZero` x `{default, robust, cluster(unit), ccnoabsorb}` |

`byCohort` exits r(459) on this design, identically on both builds (10 configs); that is the
design, not the patch. **The predicted blast radius is exactly the measured one.**

Magnitude on the noisy design, the `ame` row:

| arm | PINNED b / se / p | PATCHED b / se / p |
|---|---|---|
| default, robust, cluster(unit), ccnoabsorb | 0.914042428 / .0585 / <.0001 | unchanged, bit for bit |
| `customvce(robust)` | 0.604743045 / .335374 / **0.0723** | 0.914042428 / .055504 / **<0.0001** |
| `customvce(cluster unit)` | 0.604743045 / .338909 / 0.0818 | 0.914042428 / .057355 / <0.0001 |
| `unconditionalse robust` | 0.604743045 / .337438 / 0.0740 | 0.914042428 / .066852 / <0.0001 |
| `poisson` / `poisson ccnoabsorb` | 0.091872539 / .0513 / 0.0732 | 0.227800975 / .00754 / <0.0001 |

Two things worth noticing there. The pinned `customvce` arms are **not significant at 5%**
where the truth is overwhelming, so the defect flips conclusions and not only digits. And the
patched `customvce(robust)` se (.055504) equals the `robust` arm's se **exactly**, and the
patched `customvce(cluster unit)` se (.057355) equals the `cluster(unit)` arm's se exactly —
which is what those options are supposed to reproduce, and independent confirmation that the
repaired explicit paths are now estimating the same model as the default path.

---

## 5. Second benefit: the ppmlhdfe rc-3200 fallback stops firing

`i.\`null'` is a single-level absorbed factor, which is precisely what ppmlhdfe refuses; the
`_rc == 3200` branch at ado:1463-1469 then re-runs with the whole control set as explicit
regressors. In my p1 run the visible caution

> `Caution: initial estimation of ppmlhdfe failed, possibly due to inclusion of an absorbed
> categorical fixed effect with only 1 level.`

appears **twice on the pinned build and zero times on the patched build** (the `poisson` and
`poisson robust` arms; `poisson ccnoabsorb` never took the fallback because `ccnoabsorb` had
already emptied `absorb()`, and it was broken by the explicit-regressor route instead).
So under `untreatedZero` the fix removes the fallback entirely: the first ppmlhdfe attempt now
succeeds, which is both correct and faster. The fallback of course remains for every other
specification that needs it.

---

## 6. What a user sees — before and after

Unchanged, and this is worth stating plainly for the write-up: the generic caution

> `Caution: a treatment period coefficient interacted with something ... was dropped for
> multicollinearity.` / `Note: the multicollinearity caution above was triggered N times total.`

**still fires after the fix** (N goes 6 -> 5 on the repaired `customvce(robust)` arm). The
remaining drops are the structurally all-zero dose columns of untreated cells, which
untreatedZero always produces. So this caution was never a usable signal of the defect and
still is not one; do not propose it as the detector.

---

## 7. Where the narrow fix is NOT enough — `contreatpoly(>1)`, and the second edit that fixes it

The narrow fix repairs `lattice()` arms completely (`lattice(all 3)` and
`lattice(baseTreatedPost 3)`, both engines, both designs: `att` and `ame` back to truth).
It does **not** fully repair `contreatpoly(>1)`.

Measured (design B, 6 groups x 20 units x 9 periods, exact fit; truths `att` 2.267875862773,
`ame` 0.907333333333), `contreatcontrols(untreatedZero) contreatpoly(2)`:

| arm | PINNED att / ame | NARROW-FIX att / ame |
|---|---|---|
| default (reghdfe) | 2.267875862773 / 0.907333333333 | unchanged (correct) |
| `robust` | 2.267875862773 / 0.907333333333 | unchanged (correct) |
| `poisson` (ppmlhdfe) | broken / broken | repaired on THIS design — but see §7b, where with real curvature its `att` is repaired and `poisson customvce()` / `poisson unconditionalse` are NOT |
| `customvce(robust)` | 1.209509094751 / 0.187333333332 | 1.209509094752 / **0.907333333333** |
| `unconditionalse robust` | 1.209509094751 / 0.187333333332 | 1.209509094752 / **0.907333333333** |

`ame` is repaired; **`att` does not move at all**, and 1.20950909475 is exactly
`contreatcontrols(basic)`'s `att` on the same data. Per horizon on the patched build, the
`customvce(robust)` arm against its own default arm:

```
b_es_att post0   1.768368320  ->  0.758076161     WRONG
b_es_att post1   2.186324924  ->  1.103092788     WRONG
b_es_ame PRE     0.000000000  ->  0.720000000     WRONG  (= beta(i=3,t=9), the dropped slope)
b_es_ame post0/1 0.830 / 0.990                    correct
```
So after the narrow fix alone, a `contreatpoly(2)` run on those engines returns a **hybrid**:
the `ame` post rows from the right model and the whole `att` family plus the pre-treatment
`ame` row from `basic`. Reproduced on a second design (`p4_edge.do`, 7 units per cell:
patched `att` 1.148673365335 vs truth 2.117554176197).

**Mechanism.** It is the same collinear column arriving by a different route.
`contreatpolypart` is built with `##` at ado:1110, so the untreatedZero injection at ado:1314,
`` `latticepart1hash'#`contreatpolypart' ``, expands to `cell#dose` + **pooled `c.dose`** +
`cell#dose#dose`. That pooled column is exactly what the 1234 gate used to insert as the
ringer; under `untreatedZero` it is again identically the sum of the cell dose columns.
`margins ... dydx(contreat)` still sees it (so `ame` survives), but the `att` family's
`targetflag` contrast does not (so `att` loses `beta_dropped * mean(dose)`). Confirmed with a
direct `reg` replication: `_b[dose]` resolves, and the only omitted dose columns are the
structurally empty untreated-cell ones. `reghdfe` drops the pooled column instead of a cell slope, which is why the
**default** path is unharmed — measured with real curvature in §7b below, relerr 9.3e-16.
`ppmlhdfe` does the same, so the `poisson` (ppmlhdfe) arm's `ame` is repaired by form A alone;
but `poisson customvce()` and `poisson unconditionalse` go through plain `poisson` at ado:1485,
the same explicit sweep as `reg`, and they are affected exactly as the `reg` arms are.

### 7b. Validating form C needs CURVATURE — the test that actually bites (`p5_curv.do`)

My first poly probes used a DGP linear in the dose, so every true `cell#dose#dose` coefficient
was 0 and `dydx(contreat)` came back as `mean(beta)` however margins differentiated the
quadratic. That cannot validate a change to the spelling margins reads. `p5_curv.do` adds real
curvature: cell-specific `b1` **and** `b2` (`mean(b2) = -0.1613`),
`y = a_i + g_t + kappa + b1*dose + b2*dose^2`, exact fit, so the truth is
`ame = mean(b1 + 2*b2*dose | treated) = 0.477908937767` and `att = 1.740277602211`.
The poisson twin is `exp(... + g1*dose + g2*dose^2)`, `g = b/4`, truth `ame 0.119477234442`.

`contreatcontrols(untreatedZero) contreatpoly(2)`, all three builds, same data:

| engine | PINNED att / ame | A only (`buildfix`) att / ame | A+C (`buildfix2`) att / ame |
|---|---|---|---|
| default (reghdfe) | 1.740277602211 / 0.477908937767 | same | same |
| `robust` | 1.740277602211 / 0.477908937767 | same | same |
| `customvce(robust)` | 0.986295713439 / **0.000129140480** | 0.771396791332 / 0.477908937767 | **1.740277602211 / 0.477908937767** |
| `unconditionalse robust` | 0.986295713439 / **0.000129140480** | 0.771396791332 / 0.477908937767 | **1.740277602211 / 0.477908937767** |
| `poisson` (ppmlhdfe) | 0.312907261693 / **0.000032285121** | 0.501402733886 / 0.119477234442 | same as A |
| `poisson customvce(robust)` | 0.312907261944 / **0.000032285133** | 0.259182531890 / 0.119477234442 | **0.501402734183** / 0.119477234442 |
| `poisson unconditionalse robust` | 0.312907261944 / **0.000032285133** | 0.259182531890 / 0.119477234442 | **0.501402734183** / 0.119477234442 |

Three things this settles that the linear designs could not.
1. **Form C is validated.** With `b2 != 0` the `#`-spelled quadratic reproduces the truth to
   2.8e-15 on the `reg` engines. If the pure-`#` spelling had broken margins' derivative, this
   is the test that would have shown it.
2. **The curvature makes the pre-fix damage far worse**: `ame` comes back as **0.000129** where
   the truth is **0.478** — not a shifted estimate, an annihilated one, and a user reading it
   would conclude the dose does nothing.
3. **`poisson customvce()` / `poisson unconditionalse` ARE affected**, and form A alone leaves
   their `att` at 0.259182531890 against the ppmlhdfe arm's 0.501402733886. A+C brings them to
   0.501402734183, agreeing with the ppmlhdfe arm to 6e-10 (the pre-existing `poisson` vs
   `ppmlhdfe` tolerance gap, present on every build). I had earlier written "poisson is not
   affected" from the ppmlhdfe arm alone; that was wrong, and this table is the correction.

At `contreatpoly(1)` on the same data all three builds agree arm for arm (`ame` 0.424308587232
on every `reg`-family engine of both patched builds; the 0.112 relative distance from the
quadratic truth is a misspecified linear model, not a defect), so C is again visibly inert.

### Form C — the second edit, written and TESTED (`buildfix2`)

At ado:1313-1315, build the untreatedZero injection from pure `#` powers:

```stata
    if "`contreatcontrols'" == "untreatedZero" {
      * UZ-POLY FIX: `contreatpolypart' is spelled with ## (ado:1110) so that on the DEFAULT
      * injection (ado:1117) the pooled base dose slope is a term the `basic' ringer can then
      * partial out. untreatedZero has no ringer, so under contreatpoly(>1) the ## leaves a
      * POOLED c.`contreat' column margins-visible; it is exactly collinear with the sum of the
      * cell#dose columns (the dose is 0 off treated cohort-periods) and on the reg engines the
      * sweep hands it one cell's slope, so every att-family row and the pretreatment ame row
      * come back as contreatcontrols(basic)'s answer. Spell the powers as pure # interactions
      * instead. At contreatpoly(1) - the default - the string built here is
      * character-for-character the old one, so this branch is inert unless a polynomial is
      * asked for.
      local contreatinj   `latticepart1hash'#c.`contreat'
      local uzpow c.`contreat'
      forval uzp = 2(1)`contreatpoly' {
        local uzpow `uzpow'#c.`contreat'
        local contreatinj   `contreatinj' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'#`uzpow'
      }
      local contreatinj   `contreatinj' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
    }
```

Results with A **and** C applied (`buildfix2`):
- `p2_identity.do` (74 configs, all at the default `contreatpoly(1)`): **byte-identical to the
  A-only build**. The `forval` runs zero times and the string it builds is character-for-
  character the old one, so C is provably inert at the default.
- `p3_latpoly.do`: exactly **4 configs move** relative to A-only —
  `untreatedZero x {poly2, poly2lat2} x {customvce(robust), unconditionalse robust}` — and each
  lands on its default arm's numbers with **0 of 16 point cells differing, max abs 0.000e+00**
  (bit-identical `att`, `ame`, every ES horizon).
- `p4_edge.do` second design: `contreatpoly(2) customvce(robust)` goes
  1.148673365335 -> **2.117554176197** (= truth) with `ame` still 0.907333333333.
- every `basic` config in `p3` unchanged from A-only.
- **C does not rewrite the default path's numbers even at `contreatpoly(>1)`.** It does rewrite
  the string reghdfe is handed there, and a reviewer will ask; the answer is measured twice:
  `untreatedZero_poly2_def` and `untreatedZero_poly2lat2_def` are byte-identical between
  `buildfix` and `buildfix2` in the `p3` dumps, and the `p5` default row is identical across all
  three builds. reghdfe was already dropping the pooled column, so removing it from the spelling
  changes the basis it never used.
- validated against real dose curvature in §7b, which is the only test that can validate C.

**Recommendation: ship A and C together.** A alone leaves a silently wrong `att` on a
documented option pair; C alone does not touch the plain (non-polynomial) case at all. Neither
is a substitute for the other, and C's inertness at `contreatpoly(1)` is measured, not argued.

---

## 8. Regression risk, and how to prove it did not happen

**What could regress, honestly stated.**
1. *Nothing outside `untreatedZero`.* Structural (the new conjunct is identically true for
   every other `contreatcontrols()` value) and measured (67 of 74 configs bit-identical, plus
   20 `basic` lattice/poly configs bit-identical).
2. *Nothing on `untreatedZero`'s default path.* The gate was already false there; measured
   bit-identical on four engines and two designs.
3. *Published numbers WILL move* for anyone who ran `untreatedZero` with `customvce()`,
   `unconditionalse` or `poisson`. That is the point — they were `basic`'s numbers — but by the
   project's own convention this is a **number-mover** and needs a cert leg and a ledger entry,
   not a silent bug-fix commit. On the exact-fit designs the moved numbers are provably the
   truth, so the leg can assert equality with a design value rather than re-blessing a baseline.
4. *The one risk I could not rule out by structure* was that removing the explicit
   `i.`null'##c.(ringer)` would change the collinearity DROP ORDER on the reg/poisson engines in
   some way the long comment at ado:1225-1233 was protecting against. It does not, because under
   `untreatedZero` there is no margins-visible base dose term to protect: `contreatinj` at
   ado:1314 has no `##c.`contreat'` main effect (at `contreatpoly(1)`; at `contreatpoly(>1)` it
   does, which is exactly finding §7, and form C removes it).
5. *`e(N)`, `e(df_r)`, `e(r2)` do not move* on the repaired arms (378 / 329 / 1 before and
   after). The rank of the design is unchanged; only the basis is. Nobody's degrees of freedom
   change.

**The test that would prove the fix correct (what I would put in the suite).**
- An exact-fit continuous-treatment design whose dose is identically 0 off treated
  cohort-periods, with a known per-cell slope. Assert `ame == mean(beta)` and
  `att == mean(kappa + beta*dose)` to 1e-10 on **every** engine:
  default, `robust`, `cluster()`, `ccnoabsorb`, `customvce(ols|robust|cluster)`,
  `unconditionalse robust`, and the five poisson arms. Pre-fix six of those fail; post-fix none.
- A cross-engine identity assertion that needs no truth at all and would have caught this in
  2023: under `contreatcontrols(untreatedZero)`, **every** engine's point-estimate vector
  (`e(wdidmainresults)` col 1, `e(b_es_att)`, `e(b_es_ame)`) must agree with the default
  reghdfe path's. It is cheap, it is oracle-free, and it generalises to every other
  `contreatcontrols()` value.
- A guard assertion that `untreatedZero` inserts no controls: after estimation under
  `verbose`, the echoed command must contain neither `contreatcontrolinserterFE` content nor a
  pooled `c.`contreat'` term.
- The inertness leg: the 74-config `contreatcontrols x engine` battery, pinned-vs-patched,
  demanding bit-identity everywhere except the six `untreatedZero` explicit-engine cells.

**One more leg the §7b result adds.** Any suite case for this must carry **real dose curvature**
(`b2 != 0`) under `contreatpoly(2)`, or it cannot see form C at all: with a dose-linear DGP the
quadratic coefficients are 0 and `ame` comes back right under either spelling. A linear
`contreatpoly(2)` test is a test that cannot fail.

**Files.** All probes and logs are in
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uz-minfix`:
`p1_truth.do` (truth recovery + whole-vector agreement), `p2_identity.do` (blast radius),
`p3_latpoly.do` (lattice and poly), `p4_edge.do` (verbose + direct `reg` replication +
bootstrap smoke), `p5_curv.do` (real dose curvature; the test that validates form C),
`cmp.py` / `agree3.py` / `sc.py` (the comparators), `NOTES.md`,
and the patched builds `buildfix/wooldid.ado` (form A) and `buildfix2/wooldid.ado` (A + C).
Dumps: `dump_p{1,2,3}_*.txt`; logs `p{1,2,3,4}_*.txt`.

## 9. Things I did NOT test
- `bootstrap()` under `untreatedZero x customvce(robust)` is refused **r(198)** on both builds
  (I did not chase which screen refuses it), so I could not check the bootstrap frame. `permute()`
  is refused by design (ado:5841). The bootstrap/permute re-entry paths pass
  `contreatcontrols()` through (`wdb_ro`/`wdp_ro` at ado:141/337), so they should inherit the fix,
  but I did not demonstrate it.
- `contreatpoly(3)` and above. The `forval` in form C generalises, but only `contreatpoly(2)`
  was run.
- `weights`, `subgroups()`, `interactivecontrols()`, `contreatwithin`, `ccc_absorb()` and
  `timetrends` combined with `untreatedZero` on the explicit engines. None of them touch the
  1234 gate, so I expect nothing, but they are untested.
- Form B (clearing the macros inside the 1313 branch) was reasoned about, not run.
- The `poisson customvce/unconditionalse` arms sit 2.0e-7 from truth on **both** builds. I
  attribute that to the `poisson` vs `ppmlhdfe` convergence tolerance at `regtol(12)` because the
  same offset is present pre-fix; I did not confirm it by varying `regtol`.
- Whether the project's own suite has any `untreatedZero x customvce/unconditionalse/poisson`
  config at all. The suite is off-limits to me; if it does not, that absence is itself worth a
  line in the write-up.
