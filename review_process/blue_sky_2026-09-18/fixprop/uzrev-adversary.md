# uzrev-adversary — adversarial review of the converged `untreatedZero` fix

**Mandate: break it.** Workdir `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uzrev-adversary`
(five ado builds, nine datasets, four probe grids, 368 config x build runs, every dump and log kept).

## VERDICT IN ONE LINE

**SHIP, WITH CHANGES.** The primary (either form) survived everything I threw at it: I could not
find a single number outside `untreatedZero x {customvce()|unconditionalse|poisson}` that moves,
and Form A and Form B are byte-identical on 368 configurations. But **four of the synthesis's
supporting claims are false as stated** — "0 return-code changes" (there are five), "Form C moves
ONE default-path number" (six), Form C's "the number it moves is already the as-if-`basic` answer"
(not on my curvature design), and the proposers' own suggested regression assertion (cross-engine
identity) is unsafe as a suite leg. And I found **one wrong number the fix does not fix**, on
entirely ordinary data, carrying nothing but wooldid's existing generic multicollinearity caution.

---

## 0. What I built

Five ado copies from the pinned `c1d8979` build:
`base`, `fa` (Form A, gate at ado:1234), `fb` (Form B, three `local` clears after ado:1314),
`fac` (A+C), `fbc` (B+C). Diffs are one hunk each; `fa`/`fb` verified byte-exact against the
proposal text.

Nine designs, all mine, all exact-fit with `e(r2)` asserted `== 1` (the coordinator's hazard (e)):
`i` is the COHORT and there is no unit-level heterogeneity, so every term lies in the model's span.

| data | shape |
|---|---|
| **d1** | 6 cohorts x 8 units x 9 periods; onsets 3/5/7; dose **exactly 0** off treated cohort-periods. Truth `att` 1.188350230862, `ame` 0.795333333333 |
| **d2** | d1 but dose drawn for **every** row plus a real nuisance slope 0.83 — `untreatedZero`'s precondition VIOLATED |
| **d3 / d4** | d1 / d2 plus real dose **curvature** (cell-specific `b2`) |
| **d6** | d1 at **40 units per cell** (2160 obs) + a clean 3-level premade lattice. Truth `att` 1.314259241731, `ame` 0.795333333333 |
| **d7** | d1 with a tiny measurement noise on the off-treatment dose (near-violation) |
| **d8** | dose is a **cohort-by-period** variable (constant within every treated cell) — an ordinary policy-intensity panel |
| **d9** | d6 with **one** treated cell (cohort 2, onset period, 40 obs) given a constant dose |

Grids: 221 configs (engines x lattice types x poly x weights x `if` x `timetrends` x
`interactivecontrols` x `ccc_absorb` x `controls` x `contreatwithin` x `subgroups` x
`semielasticity` x `contreatelasticitytype` x `poisexpresults` x `fastnose` x `nofastcollapse` x
`jointtests` x `impvalues` x `esfixedbaseperiod` x `anticipation` x `oostreatedcontrols` x
`emptycellsoverride` x `safety(off)` x `regtol()` x `makecfxplots` x `gteffects`), then 94 lattice /
near-violation configs, 31 bootstrap / `fe()` / `subgroups` / `customweightmultiplier` /
`adversarial` / `pretrendpower` / poly-3-4 configs, and 22 collinear-dose configs. Every `e()`
scalar and matrix I could enumerate is dumped in `%21x` hex and compared as **strings** (exact
bit identity), not as decimals.

---

## 1. WHAT SURVIVED (the fix's core claims hold, on a much wider grid)

### 1.1 Form A and Form B are byte-identical — verified, not asserted
`0` differing keys, `base`-relative and against each other, over **368 config x build runs** in four
independent grids, and the same for `A+C` vs `B+C`. This includes every configuration the
proposers did not run: `lattice(byCohort|byTtre|byTreated|byTime|baseTreatedRef|baseTreatedPost)`,
a premade lattice variable, `aw`/`pw`/`fw`, `if`, `timetrends`, `interactivecontrols()`,
`ccc_absorb()`, `controls()`, `contreatwithin`, `subgroups()`, `fe()`,
`customweightmultiplier()`, `semielasticity`, `contreatelasticitytype()`, `poisexpresults`,
`fastnose`, `nofastcollapse`, `impvalues`, `adversarial`, `bootstrap(standard)`, `contreatpoly(2/3/4)`,
`safety(off)`, `regtol(9|12|14)`, `emptycellsoverride`.

The only structural difference I can find by reading is that Form B still runs `tempvar null` +
`qui gen \`null' = 1` and still builds the dead `ccringerhead`/`ccringerheadabs` locals. I checked
the one place that could have cared — the SI-8 preflight at `ado:1352-1378` — and it counts cells
and compares against `c(maxvar)`, never `c(k)`, so the extra constant variable cannot flip it.
**Verdict: the byte-identity claim is CONFIRMED. Either form is safe. I prefer B** because it
cannot be reopened by a later edit to the 1234 gate, and because it leaves the DF-1 comment block
at 1225-1300 untouched; A's only advantage is that it never creates the dead variable.

### 1.2 Blast radius — I could not move anything outside the allowed set
Across all four grids, **zero** non-`untreatedZero` tags moved: every `contreatcontrols(basic)`
arm, every binary-mode arm, on every engine, bit-identical. And every `untreatedZero` arm on the
default `reghdfe` path — plain, `robust`, `cluster(i)`, `cluster(u)`, `ccnoabsorb`,
`ccnoabsorb robust`, and all of the feature crosses above — is bit-identical.

I also checked the structural claim that the default `reghdfe` path *cannot* move: `strpos` over
the ado shows no internal assignment to `customvce`/`unconditionalse`/`poisson`, and `ado:1440`
routes to `reg`/`poisson` whenever any of them is set, so under `untreatedZero` the 1234 gate is
reachable only through the three engine disjuncts. **Structure and measurement agree.**

### 1.3 Truth recovery is exact, on every repaired arm
On d1 and d6 every repaired arm lands on the design truth to the last bit — `att`
1.188350230862 / 1.314259241731 and `ame` 0.795333333333 — under `customvce(ols)`,
`customvce(robust)`, `customvce(cluster u)`, `unconditionalse robust`, `ccnoabsorb customvce()`,
`ccnoabsorb unconditionalse`, and in combination with `fe()`, `subgroups()`,
`customweightmultiplier()`, `adversarial`, every lattice type, weights and `if`.
The five poisson arms agree with each other to 2e-10 (the pre-existing `poisson`-vs-`ppmlhdfe`
tolerance gap, present on `base` too).

`bootstrap(standard) reps(19)` inherits the fix correctly: `untreatedZero poisson bootstrap()`
moves with the fix and lands on the ppmlhdfe answer, while `untreatedZero bootstrap()` and
`untreatedZero cluster(i) bootstrap()` on the default path are bit-identical to `base`.
(`untreatedZero customvce() bootstrap()` is r(198) on all builds — the stray-companion screen.)

### 1.4 e(N), e(df_r), e(r2) never move
Confirmed across every moved tag. `e(ll)`, `e(rss)`, `e(rmse)` move by 1 ulp on the poisson arms
only (the ppmlhdfe path stops taking the rc-3200 fallback), which is benign but is a literal
exception to the "0 fit invariants moving" phrasing.

---

## 2. WHAT I BROKE — five defects in the proposal as written

### DEF-1 (moderate) — "0 return-code changes" is false: there are FIVE
Four are caused by the **primary** and one by **Form C**. All are error -> success, so the
direction is benign, but each one means a user who saw an error yesterday sees a number tomorrow,
and one of them creates a new cross-engine asymmetry (DEF-2).

| config | base | patched |
|---|---|---|
| `untreatedZero lattice(all 3) contreatpoly(2) poisson` (d1) | **r(3201)** | rc 0, `att` 0.528201781523 |
| `untreatedZero customvce(ols)` (d8) | **r(322)** | rc 0, `att` 2.254205175163 = truth, `ame` **missing** |
| `untreatedZero customvce(robust)` (d8) | **r(322)** | rc 0, same |
| `untreatedZero unconditionalse robust` (d8) | **r(322)** | rc 0, same |
| `untreatedZero contreatpoly(2) customvce(robust)` (d8) | r(322) on base **and on A-only** | rc 0 **only under Form C** |

r(322) here is wooldid's own informative guard —
*"Error: a treated cohort-period dummy coefficient was dropped for multicollinearity … i == 1,
relative-t == 0"*. That guard is the good design; the point is that the fix silences it on three
engines.

### DEF-2 (high) — after the fix, the default engine refuses a specification that the explicit engines now answer
On **d8** (dose is a cohort-by-period intensity variable — minimum wage, subsidy rate, capacity
cap: an entirely ordinary shape), post-fix:

| engine | result |
|---|---|
| default `reghdfe`, `robust`, `ccnoabsorb`, `poisson` | **r(322)**, the informative guard, on every build |
| `customvce(ols)`, `customvce(robust)`, `unconditionalse robust` | **rc 0**, `att` = 2.254205175163 (the truth), `ame` = `.` |

Before the fix all seven refused. The `att` the explicit engines now report is *correct*, the
missing `ame` is *honest*, and the run does carry the generic "a treatment period coefficient
interacted with something ... was dropped for multicollinearity" caution (verified in
`run4_fa.log`), so the new behaviour is defensible and is arguably better than the hard error. But
it is **new**, it is engine-dependent in the opposite direction from CRIT-3, and it is not inside
the proposal's stated envelope. A user refused by the default engine who adds `customvce(robust)`
on a referee's advice now gets a table with a caution above it rather than an error.

### DEF-3 (MEDIUM — a wrong `ame` the fix does not fix; cautioned, but only generically)
On **d9** — `d6` with **one** of fifteen treated cells given a constant dose (40 obs; think "every
unit in the 2010 cohort got the same first-year subsidy") — post-fix:

| engine | `att` | `ame` |
|---|---|---|
| truth | 1.314259241731 | **0.795333333333** |
| default `reghdfe` / `robust` / `ccnoabsorb` | **r(322)** (guard fires) | — |
| `customvce(ols)` / `customvce(robust)` / `unconditionalse robust`, **patched** | 1.314259241731 (exact) | **0.729333333333** |
| the same arms on `base` | 0.401894023221 (= `basic`'s answer) | 0.132 |

The deficit is **0.066000000000**, and it is exactly `beta(dropped cell) x cell share` =
`0.99 x 40/600`. That is the **same signature as CRIT-2** — *truth minus one dropped cell slope* —
arriving through a second route the fix does not close: when the dose is constant inside a treated
cell, `i.cell#c.contreat` is proportional to `i.cell`, one of them is swept out, and on the
explicit engines it is the **dose slope** that goes, so the hard `r(322)` guard — which is keyed on
a treated cohort-period **dummy** being dropped — does not fire.

**What the user does see (I checked the logs, at the advisor's insistence, rather than assuming).**
Every one of these rc-0 arms on d8 and d9, on `base`, `fa` and `fac` alike, prints wooldid's
existing generic caution: *"Caution: a treatment period coefficient interacted with something
(likely, a term implementing a continuous treatment variable or an interactive control variable)
was dropped for multicollinearity. This could imply your results are not identified... Re-estimate
with option verbose for details."* So this is **not silent** — it is the observable-failure class,
and it is the same generic signal CRIT-3 already names as the only warning a user gets. That is why
I rate it medium and not high. What is still true and still worth the ledger line: the number is
wrong by a specific, identifiable amount; the caution does not say which quantity is affected or
by how much; the `att` beside it is exactly right, which makes the caution easy to dismiss; and
`e(r2) = 1`.

The fix did not create this (base is wrong there too, differently), and it is its own defect rather
than the fix's. It matters here only because it shows that the fix's own success criterion —
"under `untreatedZero` every engine agrees with the default path" — is still violated on ordinary
data. The natural remedy, if the author wants one, is to make the drop check fire when a treated
cell's `contreat` slope is dropped and not only when the cell dummy is; that is a new feature,
outside this fix's scope, and I list it as a recommendation rather than a ship gate.

### DEF-4 (moderate) — Form C moves SIX default-path numbers, not one, and moves them by O(1)
The synthesis says C "moves ONE default-path number (att 2.444765194869 -> 4.298562435389)". On my
designs the primary-plus-C build moves six default-engine (`reghdfe`) results materially:

| config | A-only | A+C |
|---|---|---|
| `untreatedZero contreatpoly(2)` (d2, default) | 1.372267694431 | **2.578226292560** |
| `untreatedZero contreatpoly(2)` (d4, default) | -0.076915721698 | **1.131312999789** |
| `untreatedZero contreatpoly(2) robust` (d4) | -0.076915721698 | **1.131312999789** |
| `untreatedZero contreatpoly(3)` (d4, default) | 0.237498926680 | **1.131312999789** |
| `untreatedZero contreatpoly(3) robust` (d4) | 0.237498926680 | **1.131312999789** |
| `untreatedZero contreatpoly(4)` (d4, default) | 0.237498926680 | **1.131312999789** |

plus six more on d7 that move only at 1e-15. `e(ll)`, `e(rss)`, `e(rmse)` move on these too
(`e(ll)` 14059.16 -> 13616.99), which the "0 fit invariants" line does not allow for. On
precondition-**satisfying** data C is inert on the default path at every `contreatpoly()` I tried
(1,2,3,4), confirming the proposers. The movement is confined to violated / near-violated data —
but that is the population that will actually hit it, because nothing checks the precondition.

### DEF-5 (moderate) — Form C's justification does not reproduce
The stated defence is that the pre-C default number is *"already incoherent … the as-if-`basic`
answer on an `untreatedZero` request"*. Tested directly:

| design | pre-C `untreatedZero` default `att` | `basic`'s `att`, same data/engine | post-C |
|---|---|---|---|
| d2 (violating, linear) poly 2 | 1.372267694431 | **1.372267694431** — identical | 2.578226292560 |
| **d4 (violating + curvature) poly 2** | **-0.076915721698** | **0.237498926680** — *different* | 1.131312999789 |
| d4 poly 3 | 0.237498926680 | 0.237498926680 — identical | 1.131312999789 |
| d4 poly 4 | 0.237498926680 | 0.237498926680 — identical | 1.131312999789 |

So at `contreatpoly(2)` with real dose curvature the number Form C moves belongs to **neither**
model. It is whatever the collinearity sweep produced, and it changes with `contreatpoly()` on the
same data. **This strengthens the case for C** — the pre-C default-path answer is solver-dependent,
which is a worse thing to be than "as-if-`basic`" — but the ledger entry and changelog must not
repeat the "it is the `basic` answer" sentence, because it is only sometimes true.

---

## 3. THE FIX'S OWN INCOMPLETENESS, CONFIRMED HARDER THAN THE PROPOSERS STATED

**Form A alone is not shippable.** Independently reproduced and extended:

- `untreatedZero contreatpoly(2) customvce(robust)` (d1): A-only `att` 0.380010603028, truth
  1.188350230862. `ame` repaired, `att` not. (= minfix §7.)
- **`contreatpoly(3)` and `contreatpoly(4)`**, which nobody tested: A-only `att` 0.401894023222 /
  0.401894022058 on d6 against a truth of 1.314259241731 — i.e. still exactly `basic`'s answer.
  A+C gives 1.314259241731 exactly at both orders.
- `contreatpoly(2) lattice(all 3)` and `contreatpoly(2) lattice(<premade>)`: A-only 0.401894023212,
  A+C 1.314259241731 = truth.
- **The ppmlhdfe DEFAULT poisson engine is affected too**, which the proposal explicitly denied
  ("the poisson (ppmlhdfe) arm's `ame` is repaired by form A alone"). On d7,
  `untreatedZero contreatpoly(2) poisson`: base `att` 4.903282442544, **A-only 4.878692729346**,
  A+C 0.579867572686, against the poly-1 answer 0.581271768823. A-only leaves the default poisson
  path wrong by a factor of eight.

**Form C is inert at `contreatpoly(1)`** — confirmed by string construction (the `forval` runs zero
times and the assembled macro is character-for-character the old one) and by measurement: of the
368 runs, the only `fa`-vs-`fac` differences are `contreatpoly(>1)` configs. **Form C's derivative
is correct at poly 2, 3 and 4** on real curvature (d3: `ame` -0.301153004274 and `att`
0.272703537408 = truth to 3e-15 on every engine at poly 2 and poly 3).

---

## 4. THINGS I CHASED THAT DID **NOT** BREAK IT (negative results, stated so nobody repeats them)

- **`contreatcontrols( untreatedZero )` with padding.** Stata strips the blanks; the `==` predicate
  both forms rely on is total. `wooldid_validate` (`ado:4792-4800`) forbids `untreatedZero` with any
  other word, so `ccringertail` is provably empty under it and A and B must agree. Verified.
- **Precedence.** Form A's parenthesisation is necessary (`&` binds tighter than `|`) and correct;
  measured identical to B, which sidesteps the question entirely.
- **`ccnoabsorb`** (which moves `...FE` into the explicit regressor slot at `ado:1318-1323`, *after*
  Form B's clears): both forms leave it empty, verified on `ccnoabsorb customvce(robust)`,
  `ccnoabsorb unconditionalse robust`, `ccnoabsorb poisson customvce(robust)` and
  `lattice() ccnoabsorb customvce(robust)`.
- **The ppmlhdfe rc-3200 fallback** (`ado:1468`, the only line that reads `...FE` under poisson):
  post-fix it is empty; the fallback stops being taken at all on the arms I measured.
- **`interactivecontrols`, `coarsecohortcontrols`, `ccc_absorb`, `controls`, `timetrends`** carry no
  dose term (read at `ado:1136-1177`), and none of them moved a non-`untreatedZero` number.
- **`subgroups()`** enters only as a `margins ... over()`, never the model.
- **`permute()`** is refused under `untreatedZero` at `ado:5841`, so there is no permutation surface.
- **`regtol(9|12|14)` on near-violating data**: the pre-C default path is *not* tolerance-sensitive
  on d7; C does not destabilise it (movement 1e-15).
- **A lattice residual I thought I had.** On the thin d1 (8 obs per cell) the post-fix `ame` still
  differs between the default and `customvce()` engines by 3-9% for `lattice(byCohort 2)`,
  `lattice(baseTreatedPost 3)` and a premade lattice. On d6 (40 obs per cell) **every** lattice
  arm agrees to the last bit and equals the truth. So it is thin-cell non-identification of the
  per-(cell x lattice) slope, pre-existing and not the fix's — but it has one consequence for the
  fix: **the proposers' recommended suite leg ("under `untreatedZero`, every engine's
  point-estimate vector must agree with the default `reghdfe` path") would be FLAKY.** It needs a
  minimum-obs-per-(cell x lattice) condition, or it will fail on a legitimate build.

---

## 5. REQUIRED CHANGES BEFORE THIS SHIPS

1. **Ship the primary and Form C together**, as the synthesis says. My evidence is stronger than
   theirs: A alone leaves `att` wrong at `contreatpoly(2/3/4)` on `customvce()`/`unconditionalse`
   **and on the default `poisson`/ppmlhdfe engine** (4.879 vs 0.581).
2. **Correct the ledger text.** Strike "0 return-code changes" (five, listed in DEF-1) and
   "Form C moves ONE default-path number" (six, listed in DEF-4); strike the claim that the moved
   number "is the as-if-`basic` answer" (DEF-5) and replace it with "is solver-dependent and
   belongs to no model — it changes with `contreatpoly()` on the same data". Note that `e(ll)`,
   `e(rss)` and `e(rmse)` do move (1 ulp on the poisson arms under the primary; materially on the
   six default-path poly arms under C), even though `e(N)`, `e(df_r)` and `e(r2)` do not.
3. **Record the d9 residual in the same ledger entry** (DEF-3): with a constant dose inside a
   treated cell, `customvce()`/`unconditionalse` return `att` exactly right and `ame` short by
   `beta(dropped cell) x weight(cell)` under the generic multicollinearity caution, while the
   default engine refuses with `r(322)`. Pre-existing, not caused by the fix, but it is the reason
   the fix does not achieve full cross-engine identity and it should not be discovered later as a
   regression.

## 5b. RECOMMENDATIONS (not ship gates)

- **Prefer Form B** (or ship both): measured identical on 368 runs, and it is the one a future edit
  to the 1234 gate cannot reopen. If Form A is chosen, keep its comment pointing at `ado:1313`,
  because the two `== "untreatedZero"` predicates must never drift apart.
- **Do not adopt the proposed cross-engine-identity suite leg unguarded** (section 4, last bullet):
  add a minimum-obs-per-(cell x lattice) condition, or restrict the assertion to lattice-free runs,
  or it will fail on a legitimate build.
- **Consider extending the drop check to the dose slope** (DEF-3). Today `r(322)` fires only when a
  treated cohort-period **dummy** is dropped; on the explicit engines the drop lands on
  `i.cell#c.contreat` and only the generic caution fires. A new feature, outside this fix's scope.

## 6. WHAT I DID NOT TEST
`makeplots`/`saveplots` output files; `permute()` (refused by design); `boottstat`/`bootdisplay`
families; `honestdid`/`pretrends` downstream consumers; `contreatpoly() > 4`; panels large enough
for the cell collapse to matter; `oldsyntax`; multi-way `cluster()`; and any design with more than
six cohorts. Nothing I read suggests any of them touches the 1234 gate, but they are untested.

## 7. FILES
All in `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uzrev-adversary`:
`base|fa|fb|fac|fbc/wooldid.ado` (the five builds), `mkdata.do`/`mkdata2.do`/`mkdata3.do` (d1-d9),
`harness.do`/`harness3.do` + `cfgs.do`/`cfgs2.do`/`cfgs3.do`/`cfgs4.do` (the four grids),
`runbuild*.sh`, `out*_<build>.txt` (hex dumps), `run*_<build>.log`, `cmp.py` (the comparer),
`tables.txt`, `NOTES.md`.
