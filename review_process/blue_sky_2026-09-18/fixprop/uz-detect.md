# uz-detect — can wooldid detect this class generally?

**Agent:** uz-detect. **Angle:** step back from `untreatedZero`; the general failure is "a
treated-cell dose term is dropped for collinearity and the reported average silently becomes
`mean(beta) - beta_dropped`, or averages a structural zero in".
**Workdir:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uz-detect` (NOTES.md, seven probe do-files and logs).
**The repo was not touched.** Everything below was tested on an instrumented private copy at
`.../probes/uz-detect/adocopy/wooldid.ado`.

---

## HEADLINE (four things the synthesis needs from this probe)

1. **Yes, wooldid can detect this class generally, and cheaply.** A two-stage screen —
   colnames first, one `gegen` pass only if the colnames flagged something — separated the two
   classes perfectly. Predicate design was measured on 41 runs' colnames (probes C and F) and the
   working prototype was then run end-to-end on 28 arms (probe G, section 5): **13 correct arms
   clear and silent, 15 wrong arms caught, zero false positives and zero misses** — lattice,
   `lattice(all K)`, `contreatpoly(2)`, `byTime`, `byCohort`, `interactivecontrols`, `ccnoabsorb`,
   `esfixedbaseperiod` and event-study pre-windows all included. Proposal **FP-D1**, section 4.
2. **The naive widening the brief suggests (`amedrop > 0`) must not ship.** `lattice(all 5)`
   legitimately drops 80 treated-cell dose terms on a run whose ame is right to 1e-11.
3. **The class is bigger than `untreatedZero` and it reaches the DEFAULT reghdfe engine.**
   `contreatcontrols(basic)` with one treated cell that has no within-cell dose variation returns
   `ame = 0.6222222222` against a truth of 0.7 — `0.7 * 24/27`, three of twenty-seven cells
   averaged in as structural zeros — on both engines, with a full standard error and p < 0.001.
   The narrow `untreatedZero` fix does not touch this. Section 2a.
4. **Two contradictions of the standing record.** (a) Under `untreatedZero` + `lattice()` there is
   **no multicollinearity caution at all**, not the "generic caution" the write-up credits —
   ado:11897 skips the branch when a lattice is active and `verbose` is off. (b) The block comment
   at ado:11836-11844 says the high-cell-count failure is `r(1000)` from `local bvars: colnames e(b)`;
   measured, it is **`r(920)` from `strpos("`bvars'", "o.")` at ado:11851**, and `: colnames`
   succeeds. Sections 1 and 4/FP-D2.

---

## 0. Correction to the task brief's premise (read this first)

The brief asks what `wooldid_amescreen` mode(detect) (ado:12861, called ado:1542) detects "on
`e(b)` colname markers". **It does not read colnames at all.** mode(detect) is the AMESE R1
*data-side* screen: for each cluster variable it `gegen ... nunique()`s the cluster count inside
each treated `(icell,tcell)` cell and fires when any treated cell lies wholly inside one cluster.
It prints a caution and returns `r(fire)`/`r(ncells)`; the caller uses `r(fire)` to scrub the ame
`se/t/p/ci`. It has nothing to do with omitted columns.

The `o.`/`co.` colname screen the angle is about is **`wooldid_ccscreen`** (ado:11816-11982),
which contains two separate things:

1. **the main loop** (ado:11851-11945) — the dropped-treated-dummy error 322 branch
   (`wordcount<=3`) and the **generic interaction caution** (`wordcount>3`, ado:11926-11941);
2. **the SI-A4/R6 block** (ado:11958-11976) — the ame void, predicate
   **`amelive == 0 & amedrop > 0`**.

Everything in this file is about those two.

---

## 1. What the current screens actually do — established exactly

**R6's token test** (ado:11964). A colname token counts as a treated-cell dose term iff it
contains `_i`, contains `_t`, contains `1.targetflag` or `1o.targetflag`, and contains
`c.<contreat>#` or `co.<contreat>#` (the token has a `#` appended so the test can match at
end-of-token). It is **live** iff the token contains no `o.` anywhere.

**R6's response** (ado:11969-11975): only when `amelive == 0 & amedrop > 0`. Then a caution, and
`c_local amevoid 1` / `amescrub 1` / `amescrublo 1` — the ame block is voided **including the
point estimates**.

**Gates.** The whole `wooldid_ccscreen` call is gated `if "safety" != "off"` (ado:1706). The R6
block additionally needs `contreat != ""`, `fastnose == ""`, and `strpos(bvars,"o.") != 0`.

**MEASURED: the R6 predicate fires on nothing that matters.** Across 18 arms (probe C,
`pC_colnames.txt`, analysed by `analyze2.py`) R6 returned 0 on every single one, including seven
arms whose reported `ame` was wrong — one of them wrong by 100%, one by sign.

**The generic caution is weaker than the standing write-up says.** The write-up (CRIT-2/CRIT-3)
records "the user's only signal is a GENERIC multicollinearity caution". Measured (probe E,
`pE.log`): under **`untreatedZero` + `lattice()` + `customvce(robust)` there is NO caution at
all.** The wordcount>3 branch is gated `if wordcountproblem <= 3 | "lattice" == "" | verbose != ""`
(ado:11897), so with a lattice active and `verbose` off it never executes. That run printed
`att = 8.349e-14, se 6.108e-14, t 1.37, p 0.173` against a truth of 0.8066, and said nothing.
**This is a contradiction of the existing write-up and it widens the defect.**

---

## 2. New evidence: the class is bigger than `untreatedZero`, and it reaches the DEFAULT engine

### 2a. The K-class — `contreatcontrols(basic)`, default reghdfe, one cell with no within-cell dose variation

Design: 12 cohort groups x 5 units x 6 periods, staggered onsets 3/4/5, dose `x` unit-level and
nonzero everywhere, `y = gfe + tfe + 0.7*x*post` (exact fit). Then `replace x = 1.0 if i<=3 & t==3`
— three treated cohort-period cells lose their within-cell dose variation.

| run | engine | reported ame | truth |
|---|---|---|---|
| `contreatcontrols(basic)` | **default reghdfe** | **0.622222222222222** | 0.7 |
| `contreatcontrols(basic) customvce(robust)` | reg | **0.622222222222219** | 0.7 |

`0.7 * 24/27 = 0.6222222...` to 1e-15: three of twenty-seven treated cells lost their dose slope
and were averaged in as structural zeros. `att` also moves (1.2488 vs 1.1488). R6 does not fire
(`amelive = 24 > 0`). On the default path the user gets only the generic caution, a full standard
error and `p < 0.001` on the wrong number. Repro: `pC.do` arms C_K1/C_K2, `pE.do` part 2.

**This matters for the synthesis: the narrow `untreatedZero` fix does not touch this.** It is the
same class on the engine the brief calls "correct".

### 2b. The shift algebra, confirmed on a heterogeneous design

Same panel, `untreatedZero`, cell slopes `0.2 + 0.15*i` (so they differ across cohorts).

| run | reported ame | truth |
|---|---|---|
| default reghdfe | 0.850000000000000 | 0.85 |
| `customvce(robust)` | **-0.699999999999059** | 0.85 |

`0.85 - 1.55 = -0.70`, and `1.55 = 0.2 + 0.15*9` is exactly the slope of cohort 9 — the cohort the
dropped token names (`9o.__000006#2o.__000007#1o.__000008#co.x`). Agreement to 1e-12. This is the
first heterogeneous-slope confirmation of `truth - beta_dropped`; the earlier exact-fit designs
could not distinguish "shifted by beta_dropped" from "collapsed to zero".

### 2c. Why the cell term and not the ringer

On the plain `reg`/`poisson` path (ado:1485) the order is
`... trendinjector contreatcontrolinserter contreatcontrolinserterFE ccc_absorbinserter i._i#i._t#i.targetflag<contreatinj> ...`
— the ringer set is **before** the treated-cell dose terms, and Stata drops the later collinear
column. Algebra: with `R = sum_c D_c` exactly, fitting `a*R + sum_{c>=2} b_c D_c` reproduces
`beta_1*D_1 + sum_c (a+b_c) D_c`, so `a = beta_1` and every reported `b_c = beta_c - beta_1`,
with the dropped cell reporting `0 = beta_1 - beta_1`. Every cell is shifted, which is why `att`
moves and not only `ame`.

---

## 3. Designing the widened predicate — what works and what does not

Four candidate predicates, evaluated on the colnames of 18 + 23 runs (probes C and F):

- **P0** = current R6: `amelive == 0 & amedrop > 0`
- **P1** = the naive widening the brief suggests: `amedrop > 0`
- **P3** = *subset-aware*: a dropped dose token is a hit iff **no live dose token has a key-set
  that is a subset of its own** (the key-set is the token's factor components with the `c.x`
  parts removed). Rationale: under a `##` expansion a coarser live term still supplies the slope.
- **P4** = *cell-keyed*: a treated `(_i,_t)` cell is a hit iff it has a dropped dose term and
  **no live dose term at any key**.

### 3a. P1 is unusable — lattice designs drop dozens of dose terms legitimately

| arm | spec | ame | truth | P0 | P1 | P3 | P4 |
|---|---|---|---|---|---|---|---|
| C_B4 | basic + `lattice(lat2)` + customvce(robust) | 0.700000000000 | 0.7 | 0 | **1 FP** | 0 | 0 |
| C_B5 | basic + `lattice(all 3)` + customvce(robust) | 0.700000000001 | 0.7 | 0 | **1 FP** | 0 | 0 |
| C_B11 | basic + `lattice(all 3)`, reghdfe | 0.699999999999 | 0.7 | 0 | **1 FP** | 0 | 0 |
| C_B18 | basic + `lattice(all 5)` + customvce(robust) | 0.699999999994 | 0.7 | 0 | **1 FP** | 0 | 0 |

C_B18 drops **80** dose terms on a run whose answer is right to 1e-11. An `amedrop > 0` screen
would caution on every lattice run wooldid does. It must not ship.

### 3b. P3 and P4 separate the classes cleanly

On all 18 probe-C arms: **P3 and P4 = 0 on all nine correct arms** (every lattice, `contreatpoly(2)`,
`byTime`, `byCohort` and `interactivecontrols` arm included) and **>= 1 on all seven wrong arms**
(C_U2, C_U4, C_U5, C_U8, C_K1, C_K2, C_H2). Full table in `NOTES.md`.

Five-seed check on the hardest false-positive surface (probe F): ordinary design +
`lattice(all 4)` + `customvce(robust)`, seeds 700001-700005 — P3 = P4 = 0 in all five, ame correct
to 1e-9. And on `untreatedZero` + `lattice(all 3)`/`lattice(all 4)`, same five seeds — P3 = 6..17,
P4 = 1 or 2, so P4 never missed, but in one arm it was one hit away from missing.

### 3c. A colname-only predicate is NOT sufficient — the pre-period carve-out

**This is the trap that would have shipped a false-positive screen.** Under `untreatedZero` the
dose is exactly 0 at pre-treatment relative periods, so with an event-study pre-window those cells'
dose columns are identically zero and are omitted — on the **correct default reghdfe path**.

Measured (probe F, `F_U1p`): `untreatedZero esprelength(2) espostlength(2)` on the default engine
gives **15 dropped dose tokens with no live twin**; P3 and P4 both fire on a run whose answers are
right. `contreatcontrols(basic)` with the same window drops nothing (`F_B1p`, drop = 0), so the
carve-out is about *structurally-zero cells*, not about pre-periods as such. This is exactly the
exemption `sthlp:903` already claims to make ("expected zero doses in untreated periods, as under
`contreatcontrols(untreatedZero)`, do not trigger this") and which R6 only achieves by accident,
because it needs ALL dose terms dropped.

**So the shipping predicate is two-stage: colnames first, then one data pass to confirm.**

---

## 4. FIX PROPOSALS

### FP-D1 (recommended, default on, moves no number) — per-cell dose-slope screen with a data confirmation

**Which predicate: ship P4.** One sentence so nothing below hedges differently — **ship P4** (a
treated `(_i,_t)` cell with a dropped dose term and no live dose term at any key), because it is
the cheap one, it is what section 5's prototype actually runs, and it caught every one of the 15
wrong arms. **Its one residual miss mode**, which the synthesis should record: under a lattice,
a cell where one bin's dose term dies while another bin's stays live — that bin's ame is a
structural zero and P4 stays silent, because the cell still has a live dose term. It did not
occur in 10 lattice arms across 5 seeds (P4 returned 1 or 2 hits, never 0), but it came within
one hit of occurring. **P3 closes it** (section 6d) at a cost I did not measure at scale, so P3
is an upgrade to put to the author with that measurement attached, not the shipping default.

**Site.** `wooldidinstall/wooldid.ado:11958-11976`, the SI-A4/R6 block inside `wooldid_ccscreen`,
plus one new `mode(dosecells)` branch in `wooldid_amescreen` (ado:12861) holding the Mata helper.

**What it becomes.** Keep the existing R6 block byte-for-byte — its all-dropped void is a
different, harder response and it is already ruled. **Add** a second, weaker screen after it:

*Stage 1 (colnames only, no data touched).* Walk the columns with the existing R6 token test.
For each dose token parse the `_i` and `_t` factor levels out of its own name and record the cell
`(_i,_t)` together with whether the token is omitted. A cell is a **candidate** iff it carries at
least one omitted dose token and **no live dose token at any key** (P4). If no cell is a
candidate, stop — nothing further is computed and nothing is printed.

*Stage 2 (runs only if stage 1 produced candidates).* One `gegen` pass on the estimation sample,
`max(abs(contreat))` by `(_i,_t)` over `targetflag == 1` rows, and drop every candidate whose cell
has a dose that is identically zero. That is the `sthlp:903` carve-out, made explicit rather than
accidental, and it is what clears `untreatedZero`'s pre-periods.

*Response.* A caution, **and no change to any number**:

> `Caution: <n> treated cohort-period cell(s) have a nonzero <contreat> but no estimable dose-slope coefficient: the coefficient was dropped for multicollinearity, so the average marginal effects reported below average a structural zero over those cells (and, when a pooled dose control absorbed the dropped column, every other cell's slope is shifted by it). Affected cell(s) (<i>: relative period): <list, first 5>. Check contreatcontrols(), any control set collinear with the treatment cells, and whether those cells contain within-cell variation in <contreat>.`

Cell naming is free: the `wordcount>3` branch at ado:11930-11939 already parses `problem_i` and
`problem_t` out of the colname and already converts to relative-horizon units with
`problem_t = problem_t - min_t - 1`. Reuse that identity.

**Evidence.** Prototype implemented at `.../probes/uz-detect/wdpred.do` and run over 28 arms in
`pG.do`; the result table is section 5 below.

**Blast radius.** Display only. No `e()` scalar or matrix changes, no estimate changes. It can
only add a line of output. It is inside the `safety != off` gate, so `safety(off)` still silences
it.

**What it costs, honestly.** Three tiers. A run with **no dropped dose term** pays *nothing* — the
call is gated on `amedrop > 0` and never happens. A run **with dropped dose terms but no candidate
cell** — which is every lattice run, `lattice(all 5)` dropping 80 of them (C_B18) — pays one Mata
walk of the `e(b)` colstripe on top of R6's existing interpreted loop: **0.130 s at 20,000 columns,
0.454 s at 50,000** (probe E). A run where stage 1 **does** flag a cell additionally pays one
`gegen`/`levelsof` pass over the estimation sample. The data never gets touched on a clean run,
which is the property the author's "no expensive guards on large datasets" rule is about.

**Risks and how to prove they did not fire.**
- *False positive on a design the author considers fine.* Prove with the 28-arm matrix in section
  5 plus the two five-seed lattice sweeps (`pF.do`). Add `contreatwithin`, weighted and
  unbalanced-panel arms, which I did not reach (section 7).
- *A cell whose dose is nonzero but constant (the K-class).* It fires there, and it should: that
  ame is wrong by 11% in the measured case.
- *Baseline movement.* None possible — no number changes. A blessed-baseline rerun only needs to
  confirm no new line appears on suite configs; if one does, it has found a real K-class member.

### FP-D2 (recommended alongside D1) — do the colname scan in Mata over `st_matrixcolstripe("e(b)")`

**Site.** `wooldid.ado:11846` (`local bvars: colnames e(b)`) and everything downstream of it in
`wooldid_ccscreen`; minimally, the new D1 screen should be written this way from the start.

**Why.** Two measured reasons.

*Cost.* The existing interpreted R6 loop (probe D, `pD.log`):

| columns | macro length | interpreted R6 loop |
|---|---|---|
| 500 | 19,491 | 0.121 s |
| 2,000 | 79,292 | 0.474 s |
| 8,000 | 320,492 | 1.221 s |
| 16,000 | 645,200 | **dies r(920)** |

The same predicate in Mata reading the colstripe (probe E, `pE.log`): **0.130 s at 20,000 columns
and 0.454 s at 50,000** — sub-second at more than three times the column count at which the
interpreted route dies. (The 2,000 and 8,000 Mata timings, 0.283 s and 0.719 s, are contaminated
by first-call overhead and by six other Stata processes on the machine; the 50,000 figure is the
one to quote.)

*A live ceiling defect, and the comment that describes it is wrong.* The block comment at
ado:11836-11844 says the failure at high cell counts is `r(1000)` out of `local bvars: colnames e(b)`.
**Measured: the first failure is `r(920)` — "macro substitution results in line that is too long" —
raised by `if strpos("`bvars'", "o.") != 0` at ado:11851**, because the *line* limit is
`macro limit + 16` and a 645,200-character macro inside a `strpos()` exceeds it. `: colnames`
itself succeeded. A maintainer following the comment would look at the wrong statement and the
wrong return code. Reading the colstripe in Mata removes the ceiling entirely, since nothing
passes through a macro.

**Blast radius.** `wooldid_ccscreen` is deliberately plain (not rclass, not eclass) and the Mata
call must stay that way — `st_matrixcolstripe` reads, it does not post. The one behaviour change
is that runs that used to die `r(920)`/`r(1000)` now complete. Prove with a large-cell-count run
before and after.

**Guard the Mata (advisor's catch, implemented in the prototype).** A token whose `_i` or `_t`
level fails to parse yields a missing cell code, and `anyof(live, .)` would then match another
missing on the live side and *silently clear* a real hit. The prototype `continue`s on an
unparsed level rather than recording it; a shipped version should count such tokens and, if any
exist, fall back to the conservative branch instead of clearing.

### FP-D3 (recommended, small) — the new screen must not inherit the lattice/verbose suppression

**Site.** ado:11897, `if wordcountproblem <= 3 | "lattice" == "" | verbose != ""`.

That gate is why `untreatedZero` + `lattice()` prints nothing at all (section 1). D1 as specified
sits **outside** the main loop and so is unaffected — this item is only to make sure the
synthesis does not "reuse the existing caution" by hanging D1 inside that branch, which would
reintroduce the silent-under-lattice hole. If instead the author wants the existing generic
caution to fire under a lattice, that is a separate change with a real false-positive cost
(C_B18 would print 80 times), and D1 is the cheaper way to get the same coverage.

### FP-D4 (number-mover — for ruling, not recommended by default) — void the affected cells

Same detection, stronger response: set `amescrub`/`amevoid` for the affected cells, or for the
whole ame block, as R6 does today. This **moves numbers** (an ame that currently prints would
print missing), so by the author's standing rule it is held for ruling, not shipped. Two sub-options:

- **D4a — void the ame block** whenever D1 fires. Clean and consistent with R6, but it turns the
  K-class from "11% wrong" into "no ame at all", which on a large panel with one bad cell is a
  heavy price.
- **D4b — refuse (error 459) only when the affected cells carry more than a stated share of the
  ame weight.** More proportionate, but needs the ame weights, which are built *after* this
  screen runs (ado:1570 onward), so it is a bigger structural change.

I would put D4a to the author as an opt-in option (`amestrict`, say) rather than a default.

### FP-D5 — relationship to the narrow `untreatedZero` fix

D1 is a **complement, not a replacement**, and the synthesis should say so explicitly:

- The narrow fix (exclude `untreatedZero` at the ado:1234 gate, or clear
  `contreatcontrolinserterFE`/`ABS` in the ado:1313 branch) removes the *cause* on that one
  option. It is the right fix and D1 does not substitute for it: a caution on a wrong number is
  worse than a right number.
- D1 catches what the narrow fix cannot: the K-class on the default engine (section 2a), any
  future control set that eats a treatment column, and any regression of the narrow fix itself.
- After the narrow fix lands, D1's *live* job on `untreatedZero` should be **zero hits** — which
  is itself the cheapest regression test for the narrow fix. Run the seven wrong arms from section
  3b against the patched build and require D1 to fall silent on all of them.

---

## 5. Prototype results — the proposed screen run over 28 arms

Prototype at `.../probes/uz-detect/wdpred.do`, driver `pG.do`, log `pG.log` (GRES lines).
`stage1` = candidate cells from the colname walk alone; `stage2` = candidates surviving the
one-pass zero-dose exemption, i.e. what the caution would actually report. The instrumentation
that exposes `e(b)` and the internal cell variables adds four variables to the data and changed
**no reported number** (G_B1's att is 1.148830934096, byte-for-byte the uninstrumented C_B1).

| arm | specification | stage1 | stage2 | cells reported | reported ame | truth | verdict |
|---|---|---|---|---|---|---|---|
| G_B1 | `contreatcontrols(basic)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B2 | `contreatcontrols(basic) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B4 | `contreatcontrols(basic) lattice(lat2) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B5 | `contreatcontrols(basic) lattice(all 3) customvce(robust)` | 0 | 0 | - | 0.700000000001 | 0.7 | clear (correct, silent) |
| G_B18 | `contreatcontrols(basic) lattice(all 5) customvce(robust)` | 0 | 0 | - | 0.699999999994 | 0.7 | clear (correct, silent) |
| G_B9 | `contreatcontrols(basic) contreatpoly(2) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B6 | `contreatcontrols(byTime) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B14 | `contreatcontrols(byCohort) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B16 | `contreatcontrols(basic) interactivecontrols(ctl) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B1p | `contreatcontrols(basic) esprelength(2) espostlength(2)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B2p | `contreatcontrols(basic) esprelength(2) espostlength(2) customvce(robust)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_B5p | `contreatcontrols(basic) esprelength(2) espostlength(2) lattice(all 3) customvce(robust)` | 0 | 0 | - | 0.699999999996 | 0.7 | clear (correct, silent) |
| G_B1pf | `contreatcontrols(basic) esprelength(2) espostlength(2) esfixedbaseperiod` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_U1 | `contreatcontrols(untreatedZero)` | 0 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_U2 | `contreatcontrols(untreatedZero) customvce(robust)` | 1 | 1 | 9:1 | 0.000000000001 | 0.7 | **TRUE POSITIVE** |
| G_U4 | `contreatcontrols(untreatedZero) lattice(lat2) customvce(robust)` | 2 | 2 | 8:1 9:1 | 0.000000000000 | 0.7 | **TRUE POSITIVE** |
| G_U5 | `contreatcontrols(untreatedZero) lattice(all 3) customvce(robust)` | 1 | 1 | 9:1 | 0.000000000000 | 0.7 | **TRUE POSITIVE** |
| G_U8 | `contreatcontrols(untreatedZero) contreatpoly(2) customvce(robust)` | 1 | 1 | 9:1 | -0.000000000016 | 0.7 | **TRUE POSITIVE** |
| G_U6 | `contreatcontrols(untreatedZero) ccnoabsorb customvce(robust)` | 1 | 1 | 9:1 | 0.000000000001 | 0.7 | **TRUE POSITIVE** |
| G_U1p | `contreatcontrols(untreatedZero) esprelength(2) espostlength(2)` | 15 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_U1pf | `contreatcontrols(untreatedZero) esprelength(2) espostlength(2) esfixedbaseperiod` | 18 | 0 | - | 0.700000000000 | 0.7 | clear (correct, silent) |
| G_U2p | `contreatcontrols(untreatedZero) esprelength(2) espostlength(2) customvce(robust)` | 16 | 1 | 9:1 | 0.000000000001 | 0.7 | **TRUE POSITIVE** |
| G_U4p | `contreatcontrols(untreatedZero) esprelength(2) espostlength(2) lattice(lat2) customvce(robust)` | 17 | 2 | 8:1 9:1 | 0.000000000000 | 0.7 | **TRUE POSITIVE** |
| G_K1 | `contreatcontrols(basic) customvce(robust)` | 3 | 3 | 1:0 2:0 3:0 | 0.622222222222 | 0.7 | **TRUE POSITIVE** |
| G_K2 | `contreatcontrols(basic)` | 3 | 3 | 1:0 2:0 3:0 | 0.622222222222 | 0.7 | **TRUE POSITIVE** |
| G_K3 | `contreatcontrols(basic) esprelength(2) espostlength(2)` | 3 | 3 | 1:0 2:0 3:0 | 0.622222222222 | 0.7 | **TRUE POSITIVE** |
| G_H1 | `contreatcontrols(untreatedZero)` | 0 | 0 | - | 0.850000000000 | 0.85 | clear (correct, silent) |
| G_H2 | `contreatcontrols(untreatedZero) customvce(robust)` | 1 | 1 | 9:1 | -0.700000000000 | 0.85 | **TRUE POSITIVE** |
| G_NC | `contreatcontrols(basic) verbose` | 3 | 3 | 107:0 114:0 121:0 | 0.622222222222 | 0.7 | **TRUE POSITIVE** |

Notes on the table. The four `G_K*`/`G_NC` arms are the K-class design (`replace x = 1.0 if i<=3 & t==3`);
their DGP dose slope is 0.7 in every cell, and the 0.6222222222 wooldid reports is `0.7 * 24/27`, the
defect, so they are true positives. `G_NC` is the same design with non-contiguous cohort codes
`ibig = 100 + 7*i`: the screen named cells `107:0 114:0 121:0`, which are exactly `ibig` for cohorts
1, 2 and 3, and `verbose` named the same three — so **the `_i` colname level is the user's cohort
value, and both the existing verbose message and the proposed caution label the right cohort.**
(That was on my untested list; it is now measured and clean.)

**Read the `G_U1p` / `G_U1pf` / `G_U2p` rows together — they are the whole argument for the
two-stage design.** On the correct default-engine runs with an event-study pre-window the
colname walk alone flags 15 and 18 cells; the zero-dose pass clears every one of them and the
screen stays silent. On the *broken* `customvce(robust)` twin of the same specification it clears
15 of 16 and reports the one cell that is genuinely wrong. A colname-only screen would have
cautioned on two correct runs and buried the real hit in fifteen spurious ones.

## 6. Exact edit sketch (so the synthesis can hand it to a reviewer)

Read on the pinned build `c1d8979`. Nothing below was applied to the repo.

**6a. New branch in `wooldid_amescreen` (ado:12861, `syntax` line 12864).** Add `dosecells` to the
mode set and, alongside the existing `detect`/`bycoh`/`sescan`/`vzero` branches:

```stata
  if "`mode'" == "dosecells" {
    * stage 1: colnames only. reads the stripe in Mata - no macro, so no r(920)
    * ceiling and no interpreted pass over the column list
    mata: st_local("wdcand", wdid_dosecell_cands("e(b)", "`icell'", "`tcell'", "`tflag'", "`contreat'"))
    local nfire 0
    local celllist
    if "`wdcand'" != "" {
      * stage 2: ONE gtools pass, reached only when stage 1 found a candidate.
      * exempts cells whose dose is identically zero on the estimation sample -
      * the sthlp:903 carve-out, which is what clears untreatedZero's pre-periods
      * the `w' > 0 qualifier mirrors ado:11898 and is load-bearing, not cosmetic: under
      * pw/aw a cell whose only nonzero-dose rows carry zero weight contributes nothing to
      * the ame, so flagging it would be a false positive
      tempvar dcmx dctg dccc
      qui gegen double `dcmx' = max(abs(`contreat')) if `esamp' & `tflag' == 1 & `w' > 0, by(`icell' `tcell')
      qui gegen byte   `dctg' = tag(`icell' `tcell')  if `esamp' & `tflag' == 1 & `w' > 0
      qui gen double `dccc' = `icell'*1000000 + `tcell'
      qui levelsof `dccc' if `dctg' == 1 & `dcmx' > 1e-12 & `dcmx' < ., local(dcnz)
      foreach c of local wdcand {
        local hit 0
        foreach z of local dcnz {
          if `c' == `z' local hit 1
        }
        if `hit' {
          local ++nfire
          local ci = floor(`c'/1000000)
          local cr = `c' - `ci'*1000000 - `mint' - 1
          if `nfire' <= 5 local celllist `celllist' `ci':`cr'
        }
      }
    }
    if `nfire' > 0 {
      dis as text "Caution: `nfire' treated cohort-period cell(s) have a nonzero `contreat' but no estimable dose-slope coefficient - it was dropped for multicollinearity. The average marginal effects below average a structural zero over those cells, and where a pooled dose control absorbed the dropped column every other cell's slope is shifted by it, which moves the att results too. Affected (cohort:relative period): `celllist'`=cond(`nfire'>5, " ...", "")'. Check contreatcontrols(), any control set collinear with the treatment cells, and whether those cells contain within-cell variation in `contreat'."
    }
    return scalar nfire = `nfire'
    return local  cells "`celllist'"
    exit
  }
```

**6b. The Mata helper** (working code, run in `wdpred.do`; guards the unparsed-level case):

```mata
real scalar wdid_dc_lvl(string scalar p)
{
  real scalar d
  string scalar h
  d = strpos(p, ".")
  if (d == 0) return(.)
  h = substr(p, 1, d-1)
  return(strtoreal(subinstr(subinstr(h, "o", ""), "b", "")))
}
string scalar wdid_dosecell_cands(string scalar mname, string scalar I,
                                  string scalar T, string scalar F, string scalar CT)
{
  real scalar k, j, n, isom, ii, tt, m
  string matrix cs
  string rowvector pc, out
  string scalar tok, p
  real colvector cell, omv, live, drp
  cs = st_matrixcolstripe(mname)
  n = rows(cs); cell = J(n,1,.); omv = J(n,1,0); m = 0
  for (k=1; k<=n; k++) {
    tok = cs[k,2]
    if (!strpos(tok, I)) continue
    if (!strpos(tok, T)) continue
    if (!(strpos(tok, "1."+F) | strpos(tok, "1o."+F))) continue
    if (!(strpos(tok+"#", "c."+CT+"#") | strpos(tok+"#", "co."+CT+"#"))) continue
    isom = (strpos(tok, "o.") > 0)
    pc = tokens(subinstr(tok, "#", " "))
    ii = .; tt = .
    for (j=1; j<=cols(pc); j++) {
      p = pc[j]
      if (strpos(p, "."+I)) ii = wdid_dc_lvl(p)
      if (strpos(p, "."+T)) tt = wdid_dc_lvl(p)
    }
    if (ii >= . | tt >= .) continue          /* never let a missing key clear a hit */
    m++; cell[m] = ii*1000000 + tt; omv[m] = isom
  }
  if (m == 0) return("")
  cell = cell[1..m]; omv = omv[1..m]
  live = uniqrows(select(cell, omv :== 0))
  drp  = uniqrows(select(cell, omv :== 1))
  if (rows(drp) == 0) return("")
  out = J(1,0,"")
  for (k=1; k<=rows(drp); k++)
    if (rows(live)==0 | !anyof(live, drp[k])) out = out, strofreal(drp[k], "%21.0f")
  return(invtokens(out))
}
```

**6c. Call site.** One line at the end of the R6 block, `wooldid.ado:11976` (after the existing
`if amelive == 0 & amedrop > 0 {...}` closes, inside the same `contreat`/`fastnose`/`o.` gate) —
or, better, moved next to the ado:1706 `wooldid_ccscreen` call so it is not skipped when `bvars`
is huge:

```stata
    if `amelive' > 0 & `amedrop' > 0 wooldid_amescreen , mode(dosecells) icell(`_i') tcell(`_t') tflag(`targetflag') esamp(`esamp') contreat(`contreat') mint(`min_t')
```

The `amelive > 0` qualifier keeps the two screens disjoint: the all-dropped case stays R6's (and
keeps its void), the partial case is the new one. Everything the new branch needs is already in
`wooldid_ccscreen`'s scope: `_i`/`_t`/`targetflag` from `iv()`/`tv()`/`tf()`, `esamp` from `es()`,
`contreat` from `ct()`, `min_t` from `mint(real)` and the weight variable from `wt()`. Only
`wooldid_amescreen`'s own `syntax` line (ado:12864) needs two new optionals, `mint(string)` and
`wt(string)`.

**6d. The P3 upgrade** (see the "Which predicate" paragraph in FP-D1: P4 ships, P3 is the
upgrade offered with a cost measurement attached). The only change is inside
the Mata helper: record each dose token's full component key-set rather than just `(_i,_t)`, and
mark a dropped token a candidate iff no live token's key-set is a subset of it. Same single pass
over the stripe, plus a bounded subset test per dropped token. P3 fired on every arm where P4 did
and on more of them (6..17 vs 1..2 on the lattice-untreatedZero arms), with no extra false
positives across 10 lattice arms and five seeds. It costs a nested comparison whose bound I did
**not** measure at scale — see section 7.

---

## 7. What I did not test (honest list)

- **P3 at scale.** I measured P4 in Mata at 50,000 columns (0.454 s). The subset test P3 adds is
  `O(ndrop x nlive)` in the worst case and I did not time it on a large model. If the synthesis
  prefers P3, that measurement is required before it ships.
- **`contreatwithin`, weights (`aw`/`fw`/`pw`), unbalanced panels, `subgroups()`, `impvalues`,
  `bootstrap()`/`permute()` replications, and `ccc_absorb()`** were not run against the prototype.
  **No weighted arm was run at all.** The spec in 6a carries the `& `w' > 0` qualifier on both
  stage-2 `gegen` lines for the reason given there (a cell whose only nonzero-dose rows carry zero
  weight contributes nothing to the ame and must not be flagged), mirroring ado:11898 — but the
  prototype in `wdpred.do` was run unweighted, so that qualifier is reasoned, not measured. A
  weighted arm is the first thing an adversarial reviewer should ask for.
- **`poisson`.** Probe B confirmed the defect fires there (`att = 0.0000000000` on the default
  ppmlhdfe engine after the rc-3200 fallback, dropped token
  `9o.__000005#2o.__000006#1o.__000007#co.x`), but I did not run the prototype screen on a poisson
  arm. The predicate is engine-agnostic by construction; confirm it anyway.
- ~~Whether the `_i` colname level equals the user's cohort value.~~ **Measured and clean** — see
  the `G_NC` row and the note under section 5's table. Non-contiguous cohort codes
  `ibig = 100 + 7*i` were named correctly by both the existing verbose message and the prototype.
  Not tested: a *string* cohort variable, or cohort codes that collide after wooldid's recode.
- **Real-scale timing of the whole screen inside a live wooldid run.** The 0.454 s figure is the
  Mata helper on a synthetic 50,000-column stripe, not an end-to-end run at that cell count (which
  cannot currently complete: the interpreted screens die first).
- **The `e(V)`/`e(b)` consequences.** I checked reported `att`/`ame` only. Whether the event-study
  horizon matrices and the `e(gteffects)` cells carry the same shift on the K-class was not
  measured (the write-up establishes it for `untreatedZero`).

