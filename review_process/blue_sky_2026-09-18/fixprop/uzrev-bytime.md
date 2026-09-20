# uzrev-bytime — the study FIX-2 needs before it can ship

**Agent angle:** measure the FIX-2 reference-pool screen's FIRE RATE on realistic staggered
panels, per `contreatcontrols` type, and decide TRUE POSITIVE vs FALSE POSITIVE on every fire.
Pinned build `c1d8979`. Workdir
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uzrev-bytime`
(every do-file, log and CSV is there; `NOTES.md` carries the running record).

---

## 0. VERDICT IN ONE PAGE

**SHIP FIX-2 as a hard `r(459)` refusal for all four types, `byTime` included — but not as
written. Four things in the proposal are wrong or missing, and one of them (the
`contreatwithin` carve-out) would leave a silently sign-flipped answer in place.**

The author's own open question — *would the byTime screen fire on ordinary staggered panels at a
rate the author would accept?* — now has a measured answer in two parts:

1. **The rate is not small. On the panel families where it fires, it fires on 100% of draws**
   (20/20 every time), on 1 to 16 calendar periods out of T. It is not a tail event: the fire is
   *structural*, it is a deterministic function of the panel geometry, and on
   dose-zero-for-never-treated data it is guaranteed for every calendar period at or after the
   last cohort's onset.
2. **Every single RUN whose screen fired was an actually-wrong run.** Across 48
   (dose-shape x type) cells at `e(r2) == 1`, plus 16 further boundary designs:
   **0 run-level false positives, 0 false negatives.** Read §4.1 before quoting that: the
   predicate detects *rank deficiency* exactly (corroborated 6/6 against wooldid's own dropped-
   coefficient count), but whether a given deficiency moves `main:ame`, moves the *pre-treatment*
   rows, or moves nothing **depends on the arbitrary sort order of the cohort identifier** — I
   reproduced the label-dependence the coordinator could not, and it changes the argument in
   FIX-2's favour rather than against it.

So the byTime screen is not "a refusal that fires on normal data". It is a refusal that fires on
normal data **that is already producing a wrong number today**, and today the user's only signal
is a three-line generic multicollinearity caution — or, under `contreatwithin`, nothing at all.

**Required changes before shipping (detail in §5):**

| # | change | why |
|---|---|---|
| **R1** | **The screens MUST fire under `contreatwithin`.** Reverse the proposer's recommendation. | Measured: `contreatcontrols(basic) contreatwithin` on a within-treated-variation design returns **-0.2389 against a truth of 1.4111 — a sign flip — with ZERO cautions printed**. And the *shipped* `byCohort` screen already fires under `contreatwithin` today, so the carve-out would make the new screens inconsistent with their own template. |
| **R2** | **The hard leg must be a TOLERANCE test, not `sd == 0` — and the near-tier's TEXT is wrong whatever threshold is picked.** | Measured miss band, bracketed: with reference-pool `sd / sd(dose)` at or below **1.1e-5** the reported `ame` is the exactly-collinear wrong answer; by **3.2e-5** it is correct to 5e-11. Across that whole band FIX-2 as written prints only "nearly unidentified ... likely computationally fragile: consider increasing `regtol()`" — false there: the answer is deterministic, not fragile, and `regtol()` does not move it. |
| **R3** | **Run every requested type's screen, including on compound spellings.** | `contreatcontrols(basic byTime)` returns the *same* wrong number as plain `byTime` (0.752778 vs truth 1.330556) while plain `basic` on the same data is exact. A screen that only looked at "the" type would pass it. |
| **R4** | **The `byTime` message must name the offending calendar periods and the two remedies that actually work** — and must not say "your dose has no variation", because on the `byTime` fires the dose usually *does* vary, just not inside the period's own reference pool. | `byTime` fires on panels where `basic` is exactly right, so the user's first reaction will be "my data are fine". Measured working remedies: (a) keep in-sample the units treated *after* the window — adding one out-of-window-treated cohort took `byTime` from wrong (0.752778) to exact (1.330556) and silenced the predicate 20/20 (probe E5, family P1); (b) drop `byTime` from the spec. |

---

## 1. METHOD, AND WHY THE PREDICATE IS EXACT RATHER THAN HEURISTIC

### 1.1 The geometric fact that decides `byTime`

Onset is constant within `i` (wooldid enforces it, `ado:772-786`). Therefore **every treated cell
`(_i,_t)` lies entirely inside one calendar period.** Consequently, for calendar period `t`, the
`byTime` ringer column `i.t#c.(ringer)` decomposes as

```
dose x 1{t}  =  dose x 1{t, treated}            +  dose x 1{t, reference}
             =  SUM of that period's cell#dose  +  dose x 1{t, _t==0}
                 (margins-visible, in span)
```

and `1{t, _t==0} = 1{t} - SUM(cell dummies in t)`, with `i.t` absorbed. So:

> **If the dose is CONSTANT (including 0) over period `t`'s `_t == 0` rows, the `byTime` ringer
> column lies EXACTLY in the model's span.** Not approximately — exactly.

The sweep then keeps the ringer (it is inserted first) and drops a margins-visible column. With
`K` cells in that period and columns `x_1..x_K` plus their sum present, the fit becomes
`SUM_{k<K}(gamma + beta_k) x_k + gamma x_K`, so every surviving cell in that period is reported as
`beta_k - beta_K`: **shifted by the dropped cell's true slope.** That is the blast report's
"truth minus the slope of whichever cell sorts last", derived rather than observed — and it is why
the only possible false positive is `beta_K == 0`.

The same decomposition explains why `byTreated`/`byTtre`/`byCohort` behave differently: their group
indicators are unions of `i` levels and therefore already inside `absorb(i.i)`, so a group whose
ringer column is *identically zero* collapses harmlessly. That is the design-Q false positive the
live-group refinement exists to kill, and §3 measures how much it buys.

### 1.2 What I built

- **Own generator, no project suite.** `probes/uzrev-bytime/genlib.do`: `wdgen` (panel shapes),
  `wdrelt`/`wdrelt3` (a faithful re-implementation of wooldid's default-regime `_t`), `wdpred`
  (the FIX-2 predicate for an arbitrary grouping variable), `bdesign2` (exact-fit truth designs).
- **`_t` replication** (default pooled-reference regime, `esprelength` default 0, no
  `timetrends`, no `anticipation`): `_t = t - ttre`; `_t = .` if `_t < -esprelength` or
  `t == mint(i)`; `min_t = max(0,-min(_t))`; `_t = _t + min_t + 1`; `_t = 0` if missing. So the
  `_t == 0` pool is **all never-treated rows + all pre-onset rows of treated units**, and it
  *shrinks* as `esprelength()` grows.
- **Truth by exact fit.** `y = a_i + g_t + g0*dose + tau_{cell}*dose`, every term inside the
  model's span (`i` = the group actually passed to wooldid, with `nrep` units inside each
  group-period cell so cells carry within-cell dose variation). `e(r2) = 1.000000` on **every**
  run reported below — checked, per the brief's warning.
  Truth = `mean(tau)` over treated-post rows = the att-weighted AME.

### 1.3 Independent corroboration of the predicate (this is the part to trust)

The predicate's **count of dangerous groups equals wooldid's own count of dropped
treatment-interaction coefficients, exactly, in every case checked**:

| design x type | predicate `ndang` | wooldid's *"the multicollinearity caution above was triggered N times"* |
|---|---|---|
| D5 x `basic` | 1 | 1 |
| D5 x `byTime` | 7 | 7 |
| D5 x `byTreated` | 1 | 1 |
| D5 x `byTtre` | 4 | 4 |
| D9 x `byTime` | 2 | 2 |
| D4 x `byTime` | 2 | 2 |
| every `ndang == 0` cell | 0 | *no caution printed at all* |

Two different pieces of code, written from opposite ends (mine from the data, wooldid's from the
`e(b)` name list), agreeing on the count. That is what makes the fire table below evidence rather
than a model.

---

## 2. TABLE 1 — FIRE RATE, type x panel shape (20 independent draws each)

Pure-data measurement; `i` = unit-level id (the new screens' predicate is **independent of the
user's choice of `i`** — only `byCohort`'s depends on it). `stageA.do` / `stageF.do`,
raw in `stageA_fire.csv` / `stageF_fire.csv`. Entry = fires / 20 draws, with the mean number of
dangerous groups in brackets.

| fam | panel shape | dose shape | `basic` | `byTime` | `byTreated` | `byTtre` | `byCohort` (shipped) |
|---|---|---|---|---|---|---|---|
| A | G40 T10, 4 cohorts (5/6/7/5), NT 17 | everywhere, time-varying | 0 | **0** | 0 | 0 | 0 |
| B | " | everywhere, unit-constant | 0 | **0** | 0 | 0 | 20 [23] |
| C | " | everywhere, 40% spike-at-0, heavy skew | 0 | **0** | 0 | 0 | 5 [0.35] |
| D | " | treated-only, unit-constant, NT dose 0 | 0 | **20 [2 of 10]** | 0 | 0 | 20 [23] |
| E | " | 0 pre-onset / `D_u` post, NT 0 (*untreatedZero shape*) | 20 [1] | **20 [7 of 10]** | 20 [1] | 20 [4] | 20 [23] |
| F | G10 T5, 2 cohorts, NT 5 | everywhere, time-varying | 0 | **0** | 0 | 0 | 0 |
| G | G200 T8, 6 cohorts (8..30), NT 90 | everywhere, spike-at-0 | 0 | **0** | 0 | 0 | 18 [2.6] |
| H | G200 T5, 3 cohorts, NT 120 | treated-only, unit-constant | 0 | **20 [1 of 5]** | 0 | 0 | 20 [80] |
| I | G40 T40, 5 cohorts, NT 11, **10% of rows dropped** | everywhere, time-varying | 0 | **0** | 0 | 0 | 0 |
| J | G40 T10, **15% of rows dropped**, unequal cohorts | treated-only, time-varying post, 0 pre | 20 [1] | **20 [7.95 of 10]** | 20 [1] | 20 [4] | 20 [22.9] |
| K | G40 T10, NT 17 | treated dosed + **20% of NT dosed** (design-P shape) | 0 | **0** | 0 | 0 | 20 [23] |
| L | G40 T10, NT 17 | everywhere, but **NT pool dose-CONSTANT** | 0 | **20 [2 of 10]** | 0 | 0 | 0 |
| M1 | G40 T10, NT 17 | everywhere, time-varying | 0 | **0** | 0 | 0 | 0 |
| M2 | G40 T10, NT 13 | everywhere, time-varying | 0 | **0** | 0 | 0 | 0 |
| N1 | G40 T10, NT 17, **`esprelength(4)`** | everywhere, time-varying | 0 | **0** | 0 | 0 | **20 [11]** |
| N2 | G40 T10, NT 17, **`esprelength(4)`** | treated-only, unit-constant | 0 | **20 [6 of 10]** | 0 | 0 | 20 [23] |
| O1 | G40 T10, NT 17, **`esprelength(2)`** | treated-only, unit-constant | 0 | **20 [4 of 10]** | 0 | 0 | 20 [23] |
| P1 | G40 T10, NT 11 **+ one cohort treated OUT of window (onset 14)** | treated-only, unit-constant | 0 | **0** | 0 | 0 | 20 [23] |
| Q1 | G40 T40, 5 cohorts, **10% rows dropped** | treated-only, unit-constant | 0 | **20 [9 of 40]** | 0 | 0 | 20 [29] |
| R1 | G200 T8, 6 cohorts, NT 90 | treated-only, unit-constant | 0 | **20 [1 of 8]** | 0 | 0 | 20 [110] |

Plus two single-design probes that matter more than their row count:

| probe | shape | `byTime` |
|---|---|---|
| **E3** | G40-style, 4 cohorts, **exactly ONE never-treated unit**, dose defined and varying EVERYWHERE | **fires, 2 of 10 periods** |
| C1a | 3 cohorts with **compressed onsets (3/4/5) and a long post window (T=20)**, NT dose 0 | **fires, 16 of 20 periods** |

### What Table 1 says

1. **`basic`, `byTreated`, `byTtre`: 0 fires in 17 of 20 families, 400+ ordinary draws.** They fire
   only on the two families whose dose is identically zero outside the treated cohort-periods —
   i.e. exactly the shape `contreatcontrols(untreatedZero)` is documented for, where the answer
   really is wrong (§3). **Their false-positive exposure on ordinary panels is nil.**
2. **`byTime` is the exposed one, and the exposure is entirely driven by ONE property:** does the
   comparison pool *present in that calendar period* carry dose variation? Families D, H, L, N2,
   O1, Q1, R1 and probe E3 are ordinary staggered panels — `basic` is exactly right on all of
   them — and `byTime` still fires, because the tail periods' reference pool is either
   dose-zero (D/H/Q1/R1), dose-constant (L), or a single row (E3).
3. **The fire is structural, not stochastic.** On dose-zero-for-never-treated panels the number of
   dangerous periods is exactly `T - max(onset) + 1` (D: 10-9+1 = 2 ✓; H: 5-5+1 = 1 ✓;
   R1: 8-8+1 = 1 ✓; Q1: 40-32+1 = 9 ✓). Compressed onsets with a long post window are the worst
   case: C1a's 16 of 20.
4. **`esprelength()` amplifies the fire COUNT, and the extra deficiencies are real but their
   damage is label-dependent.** Event-study windows eat the reference pool: family D -> N2/O1,
   dangerous periods go 2 -> 4 (`esprelength(2)`) -> 6 (`esprelength(4)`), and **wooldid's own
   dropped-coefficient counter agrees at every step** (2,3,4,5,6 for `epre` 0..4 — so all six
   deficiencies are real, not predicate artefacts). On the *as-generated* cohort labelling those
   four extra periods happened to displace pre-horizon cells whose true slope is 0, so `main:ame`
   stayed at 0.752777777778 for every window width. **Relabel the cohort identifier and that
   stops being true** — see §4.1. On ordinary dose (N1) the count stays 0 at every window width
   (measured `epre` 0,1,2,3,4 -> 0,0,0,0,0, `ame` exact at all five), so `esprelength()` creates
   no new fires on data that was fine.
5. **A calibration datum the author should weigh.** Family N1 shows the **shipped** `byCohort`
   screen firing 20/20 on a perfectly ordinary panel once `i` is unit-level and
   `esprelength(4)` is requested. Verified against the real build, not my replication
   (`stageG.log`): `wooldid y unit t ttre, att contreat(dose) contreatcontrols(byCohort)
   esprelength(4)` on a plain noisy panel with dose defined everywhere exits **`r(459)`**, while
   the same call at `esprelength(0)` runs. So an unconditional `r(459)` that fires on ordinary
   data is **already shipped precedent**, and `byCohort`'s screen has **no `safety()` gate at all**
   (`ado:1383`, "always-on").

---

## 3. TABLE 2 — THE COLUMN THAT DECIDES IT: on the fires, is the run actually wrong?

Exact-fit designs, `e(r2) = 1.000000` on every row, truth = 1.330555555556 unless noted.
`i` = cohort group, 10 units per group-period, T = 10, onsets 4/6/7/9, 3 never-treated groups,
default reghdfe path + `robust`. `stageB.do` / `stageB_truth.csv`.

| dose shape | type | predicate | reported `ame` | truth | verdict |
|---|---|---|---|---|---|
| D1 everywhere, time-varying | all 5 | no fire | 1.330555556 (all) | 1.330555556 | correct, silent — **OK** |
| D2 everywhere, unit-constant | all 5 | no fire | 1.330555556 (all) | 1.330555556 | **OK** |
| D3 everywhere, 40% zeros | all 5 | no fire | 1.330555556 (all) | 1.330555556 | **OK** |
| D8 treated + 20% of NT dosed | all 5 | no fire | 1.330555556 (all) | 1.330555556 | **OK** |
| **D4 treated-only, NT dose 0** | `basic` | no fire | 1.330555556 | 1.330555556 | **OK — correctly NOT refused** |
| | `byTreated` | no fire (raw predicate WOULD have fired) | 1.330555556 | | **OK — the live-group refinement earns its keep** |
| | `byTtre` | no fire (raw would fire) | 1.330555556 | | **OK** |
| | `byCohort` | no fire | 1.330555556 | | **OK** |
| | **`byTime`** | **fires, 2 periods** | **0.752777778** | 1.330555556 | **WRONG by -43.4% — TRUE POSITIVE** |
| **D9 everywhere, NT pool dose-CONSTANT** | `basic`/`byTreated`/`byTtre`/`byCohort` | no fire | 1.330555556 | | **OK** |
| | **`byTime`** | **fires, 2 periods** | **0.752777778** | | **WRONG -43.4% — TRUE POSITIVE** |
| **D5 zero-pre (untreatedZero shape)** | `basic` | fires | **-0.069444444** | 1.330555556 | **WRONG, SIGN FLIPPED — TRUE POSITIVE** |
| | `byTime` | fires, 7 periods | **0.175000000** | | **WRONG -86.8%** |
| | `byTreated` | fires | **-0.069444444** | | **WRONG, sign flipped** |
| | `byTtre` | fires, 4 groups | **-0.422222222** | | **WRONG, sign flipped** |
| | `byCohort` | fires | **r(459) today** | | already refused |
| | `untreatedZero` | (no screen) | **1.330555556** | | correct — the spec the data wants |
| **D6 zero-pre, time-varying post** | identical pattern to D5 | | | | |
| every ordinary shape | `untreatedZero` | (no screen) | 1.930555556 = truth + `g0` | | its own (violated) assumption, not a screen matter |

**48 cells. 0 false positives. 0 false negatives.**

### 3.1 Boundary probes — the whole-panel questions

| probe | design | predicate | `ame` | truth | reading |
|---|---|---|---|---|---|
| **C1a** | compressed onsets 3/4/5, **T = 20**, NT dose 0 | `byTime` fires **16 of 20** periods | **0.098039216** | 2.498039216 | **-96.1%.** Long post windows are the catastrophic case |
| **C1b** | onsets 4/6/8/10, T=10 — **exactly ONE dangerous period** | `byTime` fires 1 of 10 | **1.006250000** | 1.306250000 | **-23.0%. There is no "small contamination" escape hatch**: one bad period out of ten moves the headline by a quarter |
| **C2** | D4 shape, `contreatcontrols(basic byTime)` | fires (byTime leg) | **0.752777778** | 1.330555556 | compound spellings carry the full damage (R3) |
| **E3** | **ONE never-treated unit**, dose everywhere and varying | `byTime` fires 2 of 10 (the `n < 2` leg) | **0.752777778** | 1.330555556 | **the `n < 2` leg is a TRUE positive too**, and it fires on ordinary dose data |
| **E5 / P1** | D4 shape **+ one cohort treated at t=14, window T=10** | `byTime` **silent** | **1.330555556** | 1.330555556 | **the remedy works and must be in the message** |
| **E1** | D4 shape with `tau == 0` in the dangerous periods (constructed) | `byTime` fires 2 of 10 | **0.630555556** | 0.630555556 | **the only class in which a fire leaves `main:ame` intact** — exactly as the algebra predicts (`beta_dropped == 0`). Constructed here, but §4.1 shows it also arises *naturally* under `esprelength()>0`, and which periods land in it is decided by the label order, not the data |

### 3.2 What the user sees today, on the fires

- `byTime` / `byTreated` / `byTtre` / `basic` on the reghdfe default path: the generic
  *"a treatment period coefficient interacted with something ... was dropped for
  multicollinearity ... This could imply your results are not identified, depending on exactly
  what was dropped"* plus *"triggered N times"*. Nothing names the periods, nothing states that
  every affected `ame` row is shifted by a determinate amount, and the wording invites the reader
  to conclude it is innocuous.
- **Under `contreatwithin`: nothing.** See R1 below.

---

## 4. ANSWER TO THE AUTHOR'S OPEN QUESTION

> *Would the `byTime` screen fire on ordinary staggered panels at a rate the author would accept?*

**It fires on ordinary staggered panels, deterministically, whenever the panel's comparison pool
stops carrying dose variation inside a calendar period — and every such fire is a wrong number.**
Concretely, `byTime` is unidentified for calendar period `t` unless period `t`'s own reference rows
(never-treated units, plus units not yet treated at `t`) carry at least two values of the dose.
That fails:

- in **every period at or after the last onset** whenever the dose is zero (or constant) for
  never-treated units — an extremely common encoding of "treatment intensity";
- in **every period** where the reference pool has shrunk to a single row (one never-treated unit,
  or a ragged panel's hole, probe E3);
- **earlier and earlier as `esprelength()` grows**, because ES pre-horizons are treated cells, not
  reference rows.

That is a *high* rate on a recognisable family of designs, and a *zero* rate on the family where
the dose is measured for the comparison group too (A/B/C/F/G/I/K/M1/M2/N1: 0 fires in 200 draws).
**The fire rate is not the right question — the false-positive rate is, and it is zero.**

### 4.1 THE FIRE IS A PROPERTY OF THE DATA; THE DAMAGE IS A PROPERTY OF THE LABELS

This is the single most important thing I found after the fire table, and it **contradicts two
statements in `critical.md`**. Design: family-D shape (onsets 4/6/7/9, T=10, NT dose 0),
`contreatcontrols(byTime)`, default reghdfe path, truth `main:ame` = 1.330555555556 and truth for
every pre-treatment row = 0. `grp2 = 8 - grp` is a **pure relabel**: identical partition,
identical rows, identical onsets. (`stageH.do` / `stageH.log`.)

| `i` passed | `esprelength` | `main:ame` | `pretreatment:ame` | wooldid's own drop count |
|---|---|---|---|---|
| `grp` | 0 | 0.752777777778 | (no pre rows) | 2 |
| `grp` | 4 | 0.752777777778 | **-0.000000000000** (= truth) | 6 |
| **`grp2 = 8 - grp`** | 0 | **0.508333333333** | (no pre rows) | 2 |
| **`grp2 = 8 - grp`** | 4 | **-0.155555555556** | **-0.575000000000** (truth 0) | 6 |

Four consequences.

1. **Label dependence REPRODUCES.** `critical.md` CRIT-5 records the probe's label-dependence
   claim as *"NOT reproduced by the coordinator ... it should not be repeated without its own
   repro"*. Here is a repro: same data, same rows, a monotone relabel of the cohort id moves the
   headline `ame` from 0.752777777778 to 0.508333333333 (`esprelength(0)`) and to
   **-0.155555555556** (`esprelength(4)`). The coordinator's null was presumably on a design where
   the relabel did not change which cell sorted last; it is not a refutation.
2. **"Every pre-treatment row is bit-identical between the correct and the broken run" is NOT
   universal.** CRIT-5's fact 1 and the blast report's headline #4 both state it. Under
   `grp2` + `esprelength(4)` the pre-treatment `ame` comes back **-0.575 against a truth of 0**.
   So a pre-trend check *can* catch this — on some labellings, which is worse than never, because
   it means the pre-trend evidence is itself label-dependent.
3. **The `esprelength` puzzle is resolved, and it cuts FIX-2's way.** The reason `main:ame` sat at
   0.752777777778 for `epre` 0..4 under `grp` is that the four extra rank deficiencies displaced
   *pre-horizon* cells whose true dose slope is 0 — E1's class, arising naturally under an
   event-study window, where "the pre-period dose slope is 0" is just the parallel-trends
   assumption the user is already making. Under `grp2` the same four deficiencies displace post
   cells and the number moves.
4. **Therefore a screen cannot be narrowed to "only the deficiencies that do damage".** Whether a
   dangerous period damages `main:ame`, damages the pre rows, or damages nothing is not a function
   of the data — it is a function of Stata's factor sort order. "Per-period false positive" is not
   a stable category, so the run-level statement is the only defensible one: **every run whose
   screen fired was wrong somewhere, and no run whose screen stayed silent was wrong anywhere.**

---

## 5. THE FOUR REQUIRED CHANGES

### R1 (CRITICAL) — the screens MUST fire under `contreatwithin`. The proposal has this backwards.

The proposal recommends *"the new screens do NOT fire under `contreatwithin`"*, reasoning from
`ado:11918`'s existing steer. Measured (`stageD.log`, designs D2a/D2b: no never-treated units at
all, onsets 3/5/7, T=10, truth 1.411111111111):

| design | spec | reported `ame` | cautions printed |
|---|---|---|---|
| dose nonzero pre-onset (D4) | `basic contreatwithin` | 1.411111111 ✓ | — (predicate silent ✓) |
| | `byTreated` / `byTtre` / `byCohort` + `contreatwithin` | 1.411111111 ✓ | — (predicate silent ✓) |
| | `byTime contreatwithin` | — | wooldid **refuses r(198)** at parse *under default safety*: *"contreatwithin is specified in conjunction with the byTime contreatcontrol type ... To proceed anyway in a new estimation run, manually set safety(off)"* |
| **dose zero pre-onset (D5)** | **`basic contreatwithin`** | **-0.238888889** | **NONE. Zero cautions.** `se` 1.6e-15, `t` -1.5e14, `p` 0.000 |
| | `byTreated contreatwithin` | **-0.238888889** | none |
| | `byTtre contreatwithin` | **-0.544444444** | none |
| | `byCohort contreatwithin` | **r(459)** | **the shipped screen ALREADY fires under `contreatwithin`** |
| | `untreatedZero contreatwithin` | **1.411111111** ✓ | — |

Three things follow.

1. `contreatwithin` **does not rescue** the misspecification — it produces a sign-flipped point
   estimate printed as overwhelmingly significant, with **no caution of any kind**. This is the
   worst user-visible state I found anywhere in this study.
2. The carve-out would make the new screens **inconsistent with their own template**: `byCohort`'s
   screen is unconditional (`ado:1383` has no `safety()` and no `contreatwithin` gate) and fires
   here today.
3. `byTime` + `contreatwithin` is **already refused at parse under default safety**, so the
   proposal's worry about `byTime` under `contreatwithin` is largely moot — but not entirely: that
   refusal's own text offers `safety(off)` as the way through, and the D2a predicate flags 4 of 10
   periods there. The main exposure is on `basic`/`byTreated`/`byTtre`.

The `ado:11918` steer the proposal cites is for a different failure (a dropped treated *dummy*
when the dose is 0 in the untreated period). The right reconciliation is to have the new screens
fire under `contreatwithin` with a `contreatwithin`-aware remedy line pointing at
`contreatcontrols(untreatedZero)` — which returns the exact truth on precisely this design.

### R2 (should-fix) — `sd == 0` misses a band in which the answer is fully wrong

Jitter sweep (`stageD.log`, D1 block): never-treated reference dose `= 1 + j*U(0,1)`, truth
1.330555555556, `byTime`.

| jitter `j` | ref-pool `sd` | hard predicate | near-tier | reported `ame` |
|---|---|---|---|---|
| 0 | 0 | **FIRES** | – | 0.752777777778 |
| 1e-12 | ~3e-13 | silent | fires | **0.752777777778** |
| 1e-10 | ~3e-11 | silent | fires | **0.752777777778** |
| 1e-9 | ~3e-10 | silent | fires | **0.752777777778** |
| 1e-8 | ~3e-9 | silent | fires | **0.752777777778** |
| 1e-6 | ~3e-7 | silent | fires | **0.752777777778** |
| 1e-5 | 2.82e-6 (rel 2.9e-6) | silent | fires | **0.752777777862** |
| 3e-5 | 8.96e-6 (rel 1.1e-5) | silent | fires | **0.752777777997** |
| 1e-4 | 2.84e-5 (rel 3.2e-5) | silent | fires | 1.330555555509 (correct, 5e-11) |
| 1e-2 | ~3e-3 | silent | silent | 1.330555555555 ✓ |
| 1e-1, 1 | – | silent | silent | 1.330555555556 ✓ |

Up to a relative reference-pool sd of **1.1e-5** the answer agrees with the exactly-collinear one
to 9 decimals; by **3.2e-5** it is exact. The flip therefore sits three to four orders of magnitude
ABOVE `regtol(9)`'s 1e-9, so "tie the hard leg to `regtol()`" would not close the gap either. Two
things follow, and I deliberately do not name a constant I have bracketed on only one design:

- **Text fix (cheap, do it regardless).** The `1e-3` near-tier's message must stop calling the
  estimate *computationally fragile* and must say it **may be exactly the collinear answer** —
  deterministic, not noisy, and unmoved by `regtol()`. The band from ~1.1e-5 to 1e-3 is a genuine
  caution band (at rel 3.2e-5 the answer is right to 5e-11), so the tier should stay a caution;
  only its wording is wrong.
- **Threshold, if the author wants one.** The hard leg would need `sd <= tau * sd(dose|sample)`
  with `tau` around 1e-5. Bracket it on more than one design before fixing the constant.

**This is inherited from the shipped `bycoh` screen**, so it is a pre-existing gap — but FIX-2 is
the moment to close it, and the `bycoh` call site should move with it.

*Honest caveat:* a reference pool with 1e-5 relative dose variation is not an everyday dataset. I
rank this "should-fix", not "must-fix"; its main value is the text fix.

### R3 — screen every requested type, on compound spellings

`contreatcontrols(basic byTime)` returns 0.752777777778 (= plain `byTime`) where plain `basic`
returns the exact 1.330555555556. The proposal's loop-over-requested-types does this correctly as
long as it uses `strpos()` per type rather than an exact match — say so explicitly, because
`byTime` is documented as *"in addition to, not in lieu of the others"* (`ado:1254`).

### R4 — the `byTime` message

`byTime` fires on panels where `basic` is exactly right. The message must therefore not read as
"your dose has no variation". It should:

- name the **calendar periods** (not cohorts) that are unidentified, and how many of T;
- state the consequence, not the possibility: every `ame` row drawing on those periods is shifted
  by minus the dropped cell's slope, `att` rows by that slope times the cell's mean dose, and
  **pre-treatment rows are unaffected, so a pre-trend or placebo check cannot detect it**;
- give the two measured remedies: **(a)** keep in sample the units treated *after* the window, so
  they serve as that period's dosed reference pool — this is verified to restore the exact truth
  (E5/P1); **(b)** drop `byTime` (on every one of my designs, plain `basic` is exact wherever
  `byTime` fires and `basic` does not).

---

## 6. RECOMMENDATION PER TYPE

| type | recommendation | evidence |
|---|---|---|
| **`basic`** | **hard `r(459)`**, unconditional (no `safety()` gate), matching `byCohort` | 0 fires in 400+ ordinary draws (17 families); 2/2 fires wrong, one of them sign-flipped; the remedy (`untreatedZero`) returns the exact truth |
| **`byTreated`** | **hard `r(459)`**, unconditional | identical: 0 ordinary fires, all fires wrong (sign-flipped on D5/D6) |
| **`byTtre`** | **hard `r(459)`**, unconditional | identical; largest measured error of the three (-0.4222 vs +1.3306) |
| **`byTime`** | **hard `r(459)`, unconditional — with the R4 message** | 9 of 20 families fire; 0 run-level false positives across 48 truth cells + 8 boundary designs; one dangerous period out of ten already costs 23%; and §4.1 shows the damage is label-dependent, so there is no stable target to narrow the screen to |
| **`byCohort`** | keep exactly as shipped | but apply R2's tolerance leg to it as well |

### Evaluating the proposer's SI-A1 suggestion for `byTime` (error under default safety, caution under `safety(off)`)

**My evidence does not support making `byTime` the odd one out, and I recommend against it —
with one condition.**

Arguments from measurement:
- The SI-A1 pattern exists for a screen about **inference** (oversized joint tests), where
  proceeding yields a usable point estimate with a bad p-value. This screen is about **point
  identification**: proceeding yields a *wrong number*. The `bycoh` template already draws that
  distinction in its own comment — *"Runs regardless of VCE/fastnose/impvalues: point
  identification is at stake, not variances"* — and its call site carries no `safety()` gate at
  all.
- A `safety(off)` escape hatch is only worth having if some user can legitimately want the number.
  I could not construct such a case: the damage is not a variance inflation and not a bias that
  shrinks with N; it is `truth - beta(dropped cell)`, fixed by the design.
- The run-level false-positive rate that would justify an escape hatch is zero in 48 truth cells.
- **The §4.1 label-dependence result is decisive against softening.** A caution invites the user to
  read the number anyway. But the number is not merely wrong by a determinate amount: under a pure
  relabel of the cohort identifier the *same data* return 0.752777777778, 0.508333333333 or
  -0.155555555556, and the pre-treatment row goes from 0 to -0.575. There is no defensible reading
  of an estimate that is a function of Stata's factor sort order, so "proceed with a warning" has
  no coherent meaning here.

**Condition, if the author wants the hatch anyway** (and there is a fair argument for symmetry
with the rest of the package): under `safety(off)` the caution must state the *determinate* error
— "every reported `ame` row that draws on period(s) X is shifted by minus the dose slope of one
cell in that period; the shift does not shrink with sample size; pre-treatment rows are unaffected
so a pre-trend check will not reveal it" — and **not** the current
`ado:11935` language ("could imply your results are not identified, depending on exactly what was
dropped"), which understates it by a wide margin. If the author prefers the hatch, apply it
uniformly to all four new screens rather than singling out `byTime`, since `basic`/`byTreated`/
`byTtre` have the same character and a *lower* fire rate.

### On the live-group refinement

**Keep it; it is load-bearing and I measured what it buys.** On families D, H, L and R1 the raw
`bycoh` predicate (no live restriction) fires 20/20 for `byTreated` and `byTtre` — and Stage B
shows those exact specifications returning the **exact truth** on that dose shape. Without the
refinement FIX-2 would refuse four correct specifications on ordinary panels. The refinement also
correctly leaves `byCohort` silent on family L (raw would fire 20/20), matching the shipped
screen's treated-cohorts-only loop.

**One narrow gap in the refinement, for the record.** "Live" is defined as *"contains at least one
treated cohort-period row with a nonzero `contreat`"*. The geometrically exact condition is
*"the group's ringer column is not identically zero"*. They differ only for a group whose dose is
nonzero **only on reference rows** — then the column is live but confined to `_t == 0`, and for
`byTime` (whose group indicator is not a union of `i` levels) it is still exactly in span and will
displace a cell **dummy**, damaging `att` rather than `ame`. I did not construct this case and I
do not think it is worth chasing; noting it so the author knows the refinement is calibrated to
the `ame` side.

---

## 7. WHAT I DID NOT TEST (honest list)

- **Engines other than the reghdfe default path with `robust`.** The blast report establishes that
  the `basic`/`byX` defect is engine-invariant; I did not re-measure `poisson`, `customvce()`,
  `unconditionalse` or `ccnoabsorb` here. If FIX-1 changes the drop order on those paths the fire
  table is unchanged (it is a data predicate) but the *wrongness* column should be spot-checked.
- **`lattice()`, `contreatpoly(>1)`, `subgroups()`, weights, `timetrends`, `anticipation()`.**
  `lattice()` in particular multiplies the ringer set (`latticepart2hash` twins at `ado:1240-1256`)
  and the blast report already flags `lattice(all 2)` as making things worse; the screen's
  grouping would probably have to become (type x lattice level).
- **`esfixedbaseperiod`.** My `_t` replication covers the pooled-reference default only. Under
  `esfixedbaseperiod` the `_t == 0` pool is a *single* relative period per cohort, which will make
  every reference pool much thinner — **I expect the fire rate to be markedly higher there and
  this is the largest untested amplifier.** It should be measured before shipping.
- **`permute()` / `bootstrap()` replications.** A screen that errors must not error inside a
  replication; the `wdrepmode` convention used by the F13 and si8 blocks (`ado:772`, `ado:1355`)
  is presumably the right pattern, but I did not test it.
- **Cost.** I did not time the screen. The proposal's O(N) claim is consistent with the shape of
  the code, but `guard_cost_2026-09-16.md` measured the existing F13+si8 blocks at 5.3% of runtime,
  so five `gegen sd(), by()` passes are not free and should be measured on a large panel.
- **Real data.** Everything here is synthetic and exact-fit by construction.
- The **`att`-side** consequences of the fires (I checked `main:ame` only; the blast report's
  identity `shift_att = beta_dropped x mean(dose)` predicts them and is worth asserting in a test).

---

## 8. REPRO

All self-contained; none touches the project test suite or the repo.

```sh
cd /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/uzrev-bytime
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageA.do   # fire rate, 12 families x 20 draws  -> stageA_fire.csv
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageB.do   # 8 dose shapes x 6 types vs closed-form truth -> stageB_truth.csv
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageC.do   # long panel / one-bad-period / basic+byTime / contreatwithin (first pass)
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageD.do   # jitter band, contreatwithin (corrected), esprelength sweep
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageE.do   # constructed false positive, thin NT pool, out-of-window rescue
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageF.do   # 8 further families -> stageF_fire.csv
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageG.do   # SHIPPED byCohort r(459) on an ordinary panel at esprelength(4)
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do stageH.do   # label-dependence repro (SS4.1) + the R2 tolerance bracket
```

`genlib.do` holds the generator, the `_t` replication and the predicate.

**One harness trap worth passing on:** `syntax , [nocc]` does **not** create a local named `nocc` —
Stata parses a `no`-prefixed option name as the negation of an option called `cc`, so the local is
`cc` and the test `if "\`nocc'" == ""` is always true. It silently cost me one `contreatwithin`
measurement (stageC) before I caught it in stageD. Never name a `syntax` option `no<something>`
unless you mean the negation.
