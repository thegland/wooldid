# FINDINGS -- hv-lattice-maxvar (independent verification)

## VERDICT ON THE SEEDED CLAIM: **TRUE on the mechanism and the size of the miss;
## TRUE-BUT-DIFFERENT on "no message".**

The preflight really does omit `lattice()`, the miss is MULTIPLICATIVE (~1+L) rather
than additive, and a perfectly ordinary continuous-treatment specification at Stata's
**default** `maxvar` sails past the guard and dies r(907).  But the death is not
message-less: Stata's own r(907) text prints and it names the remedy (`set maxvar`).
What is lost is wooldid's *diagnosis* -- which option blew the budget -- and the
documented guarantee that this failure happens **before estimation**.

---

## F1 (CONFIRMED, self-verified from scratch) -- `si8demand` omits `lattice()` and `contreatpoly()`, so the preflight under-predicts by a factor of (1+L)(1+p)/2

### Mechanism (read from the pinned build, then measured)
- `ado:1476` regressor block: `` i.`_i'#i.`_t'#i.`targetflag'`contreatinj' ``
- `ado:1117` `` contreatinj = `latticepart2hash'##`contreatpolypart' ``
- `ado:1087` `` latticepart2hash = ##i.`latticevar' ``

so the fitted term is the FULL FACTORIAL `cell ## lattice ## c.dose(^1..^p)`.
`ado:1371` computes `si8demand = si8cells*(1 + (contreat!="")) + si8ilev + si8tlev + 200`.
The `*2` models exactly one dose slope per cell -- i.e. L=1, p=1.  Every additional
lattice level and every additional polynomial degree multiplies the block; the
preflight never sees either.

### Measurement 1 -- structural column counts (`fvexpand`, `set emptycells drop`)
Design A: 25 states (= `i`), 18 treated with staggered 2005-2012 onsets, 20 counties
per state, 20 years, 10,000 rows, continuous time-invariant county EXPOSURE.
si8cells = 213, si8ilev = 25, si8tlev = 20 -> **si8demand = 671**.

| wooldid's actual regressor block            | columns | vs si8demand 671 |
|---------------------------------------------|---------|------------------|
| `cell` only                                  |   214   |                  |
| `cell##c.exposure`      (what si8 models)    |   429   | 0.64x            |
| `cell##i.lat3##c.exposure`                   |  1719   | 2.6x             |
| `cell##i.lat5##c.exposure`                   |  2579   | **3.8x**         |
| `cell##i.lat3##c.exp##c.exp##c.exp` (poly 3) |  3439   | **5.1x**         |
| `cell##i.lat5##c.exp##c.exp`       (poly 2)  |  3869   | **5.8x**         |

429 -> 2579 is exactly x6.01 for L=5, i.e. the (1+L) factor.

### Measurement 2 -- the reconstruction is EXACT
Design B (120 states, 110 treated, 4 counties, 20 years) is large enough to trip the
preflight at maxvar 2048.  My hand reconstruction predicted cells=1271, si8demand=2882.
wooldid printed, verbatim:

    Error: this specification requires at least 2882 variables (about 1271 treated
    cohort-by-period interaction terms) but maxvar is currently 2048, so estimation
    would fail with r(907).

Exact match on both numbers.  The formula reconstruction is therefore validated, not assumed.

### Measurement 3 -- the bisection (`set maxvar`, identical command, no `verbose`)
Design A, `wooldid y state t ttre, att contreat(exposure) contreatcontrols(basic) [lattice(all 5)]`:

| `set maxvar` | no lattice | `lattice(all 5)` |
|--------------|-----------|------------------|
| 2048         | **rc 0**  | **rc 907**       |
| 2500         | -         | rc 907           |
| 3000         | -         | rc 907           |
| 3500         | -         | **rc 0**         |
| 12000        | -         | rc 0, 162 s, att = 1.0478 (design truth 0.8*E[exposure] ~ 1.02) |

Preflight predicted 671.  True demand is in (3000, 3500].  **Miss factor 4.5x-5.2x,
measured by bisection** (the structural column count is 2579; the extra ~500-900 is
dataset variables plus wooldid's ringer copies, reghdfe's and margins' own tempvars).
2048 is not an arbitrary choice: it is Stata/BE's *fixed, unraisable* ceiling, and it
is also Stata/MP's `set maxvar` minimum, so the emulation is exact.

`contreatpoly()` shares the omission and compounds it, same data, same maxvar 2048,
same si8demand of 671 for both rows:

| spec (same data, same si8demand = 671)  | fvexpand cols | rc at maxvar 2048 |
|-----------------------------------------|---------------|-------------------|
| `lattice(all 3)`                        | 1719          | **0**             |
| `lattice(all 3) contreatpoly(3)`        | 3439          | **907**           |

### Measurement 4 -- it fires at Stata's DEFAULT maxvar on a textbook design
Design p2: 50 states, 45 with staggered 2005-2012 onsets, 15 counties each, 20 years
(15,000 rows), continuous county exposure, `lattice(all 5)`.  Stata 19 MP, **default
`c(maxvar) = 5000`, nothing set by hand**.
si8cells = 525, si8ilev = 50, si8tlev = 20 -> **si8demand = 1320**, a quarter of
maxvar, so the guard is silent.

| run at the UNTOUCHED default maxvar 5000            | fvexpand cols | rc      | time |
|------------------------------------------------------|---------------|---------|------|
| `att contreat(exposure) contreatcontrols(basic)`     | 1053          | **0**   | 42 s |
| ... same + `lattice(all 5)`                          | **6115**      | **907** | -    |

`1053 -> 6115` is x5.81 for L=5. Preflight predicted **1320**; the run needed **6115**.
The bound is short by **4.63x**, and the gap straddles the default limit exactly:
1320 < 5000 < 6115.  Corroborating datum: a bare `fvexpand` of that same block at
maxvar 5000 itself returns r(907) (see p2_default_maxvar.log, the aborted first run).

The no-lattice arm of this same panel runs in **42 seconds**, so this is not a design
where runtime is the binding constraint - which matters, see below.

### Why this is realistic, and why it is worst exactly where the help sends you
- The shape is a Bartik / shift-share exposure DiD on a county panel with a state-level
  staggered policy -- among the most common continuous-treatment DiD designs there is.
- `i` = state = 50 distinct values, so the **>100-distinct-`i` caution never fires**.
  The user gets no warning of any kind.
- **The blind spot is exactly the shape wooldid recommends.**  Measured on one panel
  (25 states x 20 counties x 20 years), changing only which variable is `i`
  (`p4b_thinvsthick.do`; the identity `cols = 2*(n_cell + n_cell_lat) + 2L + 1`
  reproduces fvexpand's 2579 exactly, so it is calibrated):

  | `i`      | cells | si8demand | actual cols, lattice(all 5) | ratio | guard fires at maxvar 5000? |
  |----------|-------|-----------|------------------------------|-------|------------------------------|
  | state (20 rows/cell) | 213  | 671  | 2579  | **3.84x** | **no** - and the run needs 3000-3500 |
  | county (1 row/cell)  | 4260 | 9240 | 17063 | 1.85x     | yes (9240 > 5000)            |

  At unit-level `i` the cells are thin, `emptycells drop` collapses most of the lattice
  interaction, AND the bound is already above maxvar, so the guard does its job.  At
  coarse `i` with thick cells the ratio is twice as bad AND the bound sits far below
  maxvar, so the guard is silent.  `ado:1371` and `sthlp:1163` both tell the user, when
  a specification gets big, to "specify unit-level fixed effects via fe() and swap `i`
  to a coarser level of aggregation" -- which is the move from row 2 to row 1 of that
  table.  **Following wooldid's own advice walks the user into the blind spot.**

### What is and is NOT contradicted (this matters for the ruling)
The code comment at `ado:1342-1348` claims si8demand is "a lower bound on the
regression's variable demand, so an early refusal can only fire when the late r(907)
was certain anyway."  That **zero-false-positive invariant is preserved**: `lattice()`
only ADDS columns, so the bound stays below the truth.  The retracted 2026-09-16 defect
("a lower bound can be exceeded") therefore stays retracted, and this finding is not a
rediscovery of it.  What is broken is different: on the whole continuous-treatment
family the bound is not merely loose, it is **vacuous** -- short by a multiplicative
factor, so the guard cannot fire on the specifications most likely to need it.  And the
help makes a stronger promise than the code comment does.

### Documented-vs-actual divergence (brief category 5)
`sthlp:1167`: "if the count of variables the specification requires exceeds the current
maxvar setting, wooldid halts with an informative error before estimation rather than
failing mid-run with r(907)."  Unqualified.  `sthlp:1207` lists "the variable-count
preflight that would otherwise fail mid-run at r(907)" among what default safety buys.
Under `lattice()` the count computed is not the count the specification requires, and
the run does fail mid-run at r(907).  Neither sentence carries a lattice carve-out.

---

## F2 (the UX half of the claim -- REFINED, the claim overstates it)

What the user actually sees at maxvar 2048, `lattice(all 5)`, no `verbose` (verbatim):

    Caution: the default cell collapse does not apply with: lattice(). ...
    maxvar too small
        You have attempted to use an interaction with too many levels or
        attempted to fit a model with too many variables.  You need to increase
        maxvar; it is currently 2048.  Use set maxvar; see help maxvar.
        ... try set emptycells drop ...
        ... you might have accidentally treated a continuous variable as a
        categorical ... Use the c. operator on such variables.

So the claim's "no wooldid message" is literally true, but "uninformative" is too strong:
Stata's canned text names `set maxvar`.  Three things are nevertheless lost:

1. **The diagnosis.** Stata's text talks about "an interaction with too many levels" and
   about applying `c.` -- but the offending interaction is wooldid's own tempvar block
   (`__000006#__000007#__000008##__000009##c.exposure`), which the user never wrote and
   cannot edit.  Nothing connects the failure to the word `lattice`.  wooldid's own
   message *would* have named "about N treated cohort-by-period interaction terms".
2. **The remedy that actually applies on Stata/BE**, where `set maxvar` is impossible
   (2048 is fixed).  There wooldid's advice -- coarsen `i`, use `fe()`, or lower
   `latticenum` -- is the only way out, and it never prints.
3. **The pre-estimation guarantee** the help gives.

**Counter-observation worth recording for the fix:** when the preflight *does* fire it
raises `error 907`, so Stata's canned "Use set maxvar" text prints after wooldid's two
lines.  The `set maxvar` remedy therefore reaches the user on BOTH paths.  Conversely,
wooldid's own preflight message (`ado:1372-1373`) never mentions `set maxvar` itself --
it only offers `fe()` + a coarser `i`.  The remedy split is: Stata supplies
`set maxvar`, wooldid supplies the redesign advice, and under `lattice()` only the
former survives.

---

## F3 (negative result, reported as such) -- P6 holds: the preflight DOES still run under `contreatwithin`
`contreatwithin` rewrites the working `safety` to off; the P6 ruling (2026-09-07) says
the preflight must gate on the user's `safety()` instead.  Verified: design B at
maxvar 2048 with `contreatwithin` and no `lattice()` still produced wooldid's own
"requires at least 2882 variables" error and rc 907.  **No contradiction.** Clean.

---

---

## F4 (SELF-FOUND, same mechanism, different surface) -- the omission is a FAMILY, not a lattice special case

`si8demand` models exactly one term family: `cells x (1 or 2)`.  Four constructs in the
worker multiply the cell block and none of them is in the formula:

| construct              | ado line | what it adds to the regressor list                     | in si8demand? |
|------------------------|----------|--------------------------------------------------------|---------------|
| `lattice()`            | 1087/1117| `## i.latticevar` (full factorial)                     | **no**        |
| `contreatpoly(p)`      | 1107-1116| `## c.dose ## c.dose ...`                              | **no**        |
| `coarsecohortcontrols()`| 1138    | `i._i#i._t#(vars)`                                     | **no**        |
| `interactivecontrols()`| 1171     | `i._i#i._t#1.targetflag#c.(vars)` + `i.t#c.(vars)`     | **no**        |

`p3_family.do` measures the last two with **no lattice at all**, same design A, same
si8demand of **671** on every row, same `set maxvar 2048`:

| spec (contreat, no lattice)                      | fvexpand cols | rc      |
|---------------------------------------------------|---------------|---------|
| baseline                                          | 429           | **0**   |
| `interactivecontrols(x1 x2 x3 x4 x5 x6 x7 x8)`    | 2293          | **907** |
| `coarsecohortcontrols(i.reg)` (12 regions)        | 2997          | **907** |

Eight interactive controls is an ordinary demographic vector; a 12-level region or
industry factor in `coarsecohortcontrols()` is ordinary too.  Neither involves
`lattice()` and neither moves `si8demand` by a single variable.

This matters for the FIX: patching `si8demand` for `lattice()` alone leaves three
other multipliers uncovered.  The structurally right fix is to compute the bound from
the assembled macros (`contreatinj`, `interinject`, `coarsecontrolinserter`) rather
than from `si8cells * (1 + (contreat != ""))`, or simply to wrap the estimation command
in `capture` and translate rc 907 into wooldid's own message.

---

## F5 (low, confidence PROBABLE not certain) -- the granular-`i` caution's framing points the user away from the constraint that actually bit

`ado:1371` / `sthlp:1163` tell the user that under `lattice()` "runtime grows ... about
its 1.9th power ... so a granular `i` can dominate runtime long before it approaches
Stata's variable limits."  On **both** designs I measured, the variable limit bit first
at runtimes nobody would call prohibitive:

| design | `i` levels | cells | lattice | runtime of the comparable completed run | variable limit breached? |
|--------|-----------|-------|---------|------------------------------------------|--------------------------|
| A, 10k rows  | 25 | 213 | all 5 | 162 s at maxvar 12000                | yes, above 3000          |
| p2, 15k rows | 50 | 525 | all 5 | 42 s for the NO-lattice arm          | yes, above 5000          |

**Scope limit, stated honestly:** the caution's sentence is an asymptotic claim about
large cell counts K, and I measured two points, neither of them large.  What I can say
is narrower than the sentence: at the K where the variable ceiling is first breached
under `lattice()`, runtime is still in the tens-to-low-hundreds of seconds, so a user
reading that sentence has no reason to suspect `maxvar`.  I did not measure where the
two curves actually cross.


---

## Files (all absolute paths under
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-lattice-maxvar/`)

| file | what it establishes |
|------|---------------------|
| `repro_lattice_maxvar.do` / `.log` | **the self-contained repro for F1/F2**, RUN END-TO-END: prints `si8demand = 1320`, `ACTUAL columns no lattice = 1053`, `ACTUAL columns lattice(all 5) = 6115`, then `RUN1_RC = 0`, `RUN2_RC = 907`, `CONFIRMED`. It counts columns with the nonempty-cell identity rather than `fvexpand`, because `fvexpand` materialises the full factor grid before applying `set emptycells drop` and so hits r(907) itself at the default maxvar |
| `p2b_default_maxvar.do` / `.log` | the same at greater length. Its final OPTIONAL arm (estimability at maxvar 20000) was killed at its deadline inside the margins stage and is not load-bearing -- design A's lattice run at maxvar 12000 (rc 0, 162 s, att 1.0478 vs truth ~1.02) already establishes the spec class is estimable |
| `p1_bisect.do` / `.log` | design A, fvexpand counts, the `set maxvar` bisection, the verbatim user-visible output at 2048 |
| `p1b_extend.do` / `.log` | validation of the si8demand reconstruction against wooldid's own printed number, `contreatpoly` arms, the `contreatwithin` (P6) check, the verbose death-location run |
| `p3_family.do` / `.log` | F4: `interactivecontrols()` and `coarsecohortcontrols()` break the same bound with no lattice at all |
| `p4b_thinvsthick.do` / `.log` | F1's blind-spot table (coarse `i` vs unit `i`) |
| `p4_thinvsthick.do` / `.log` | the first attempt; kept because it shows `fvexpand` itself builds the full grid before `emptycells drop` and so cannot be used at unit-level `i` |
| `p0_orient.do` / `.log` | orientation: flavor/maxvar, the verbose regression string |
| `p2_default_maxvar.do` / `.log` | superseded first cut of p2b; kept because its abort IS a datum (a bare `fvexpand` of the lattice block returns r(907) at the default maxvar) |


---

## Harness facts a re-runner needs
- Stata 19.0, banner "MP-Parallel Edition", single-user 8-core.  `c(flavor)` reports
  `IC` here, which is a Stata 19 reporting quirk, not the licence: the operative fact
  is `c(maxvar) = 5000` by default, with a note in the startup banner saying it can be
  increased.
- `c(max_k_theory)` tracks `c(maxvar)`: 5000 by default, 20000 after `set maxvar 20000`.
  So `maxvar` is the single operative ceiling; there is no separate model-size limit
  firing here.
- `set maxvar 2048` is an exact emulation of Stata/BE's fixed, unraisable ceiling, and
  2048 is also Stata/MP's `set maxvar` minimum.
- The dataset in memory survives the r(907): `desc, short` after the failure shows the
  panel intact, and `WOOLDID_LASTRC` is set so a bare `wooldid` replay is refused.
