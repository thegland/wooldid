# SKEPTIC report — fp-subgroups (fpsg-01, fpsg-02) + three new findings

Build asserted in every log (`which wooldid`):
`/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`.
Every repro below generates its own data. None of them runs the original prober's script,
and none touches `test/`.

---

## VERDICT fpsg-01 — **CONFIRMED** (severity kept at HIGH; at the top of high)

**My repro:** `r01_sgmissing.do`, corroborated by `r02_mechanism.do`, `r04_realistic.do`,
`r05_contrast.do`.

### 1. Independent reproduction on my own design (`r01`)
Five units in subgroup 1, five DISTINCT effects (1,2,3,4,5) and five DISTINCT observed
shares (0.125, 0.25, 0.5, 1, 1), so the two candidate answers are far apart and the
share-weighted prediction is one exact number that cannot be hit by accident.

| | value |
|---|---|
| RUN A, clean sg: `att_i` = `att_itime` | **3.0000000000** (equal-unit target 3; dev −7.1e-15) |
| RUN B, sg missing on some post rows: `att_i` = `att_itime` | **3.8695652174** |
| share-weighted prediction 11.125/2.875 | 3.8695652174 (dev **−3.1e-15**) |
| equal-unit target | 3.0 → error **+0.8696, +29 %** |
| RUN C, same rows but sg CHANGES 1→0 | **r(459)** refusal |

No caution of any kind in RUN B. The header prints `N = 240 (0 obs dropped from initial
sample)`; `e(N_dropped) = 0`; no e() object records the lost weight mass (full
`ereturn list` in `r05_contrast.log`).

### 2. The mechanism is *identified*, not merely the error (`r02` PART 1)
`r01` (like the original prober's design) gives every unit the same number of post rows, so
"share of the unit's own post weight" and "row count" coincide. An UNBALANCED panel pulls
the three candidates apart:

| candidate weighting | target | reported |
|---|---|---|
| equal-unit (what `att_i` is documented to be) | 5.0000000000 | — |
| **share of own post weight** | **3.6666666667** | `att_i` = `att_itime` = 3.6666666667 (dev 1.3e-15) |
| row-weighted | 3.0000000000 | `att` = `att_it` = 3.0 |

So it is specifically `w / (unit's TOTAL post weight over the whole estimation sample)`,
exactly as claimed. Error on `att_i`: **−1.333, −27 %**.

Code corroboration (read-only): `w_per_i` at `ado:1621` is
`gegen total(w*esamp), by(posttreatment i)` — no subgroup term; `subgroups` is deliberately
omitted from the must-be-non-missing markout list at `ado:~703` ("normally, it is ok to have
missing values in the sg variable"); the F13 split guard at `ado:12266-12280` folds missing
through `cond(sg>=., ., sg)` with the comment "missing sg rows belong to no subgroup by
design"; the subgroup margins at `ado:2819` run `over(subgroups)` on the full post subpop.
**So this is not a guard oversight — the author ruled the guard's DOMAIN without addressing
the estimand consequence.** Not covered by `known_bydesign.md` §3, which is the
within-cell straddle (refused) only.

### 3. Magnitude on realistic data (`r04`, `r05`) — this is what decides severity
`r04`: 90 units, staggered onsets 5/7/9, unbalanced entry and exit, N(0,0.6) noise.
Oracle built from `e(gteffects)` (NOT `e(histogramestimates)` — the margins call that builds
the histogram, `ado:12405`, carries no subgroup term in its subpop, so it averages each unit
over ALL its post cells, labelled or not).

| arm | reported `att_i`(sg1) | share-wtd oracle | equal-unit oracle | gap in SE | per-unit weight mass |
|---|---|---|---|---|---|
| clean sg | 2.1515521370 | 2.1515521370 | 2.1515521370 | 0.00 | min 1.00, max 1.00 |
| MCAR 15 % of post rows | 2.1337963363 | 2.1337963363 (gap **0.0e+00**) | 2.1438244188 | −0.115 | 0.60 – 1.00, mean 0.847 |
| label collected only in early post years | 2.0375397973 | 2.0375397973 (gap 4.4e-16) | 2.1943162295 | **−1.468** | 0.20 – 1.00, mean 0.292 |

`r05` is the decisive one because it has a TRUTH. sg1 unit effects spread 0.5…5.5
(equal-unit truth 3), sg0 all 1, and the ten LARGEST-effect sg1 units lose the subgroup
label from relative year 2 on — non-ignorable item nonresponse, the realistic case.
**58 of 1260 rows (4.6 %) lose their label; nothing about the outcome data changes.**

| | clean | labels lost | truth |
|---|---|---|---|
| `att_i`(sg1) | 3.018843 | **2.409209** (−20 %) | 3 |
| `att_i`(sg0) | 0.922207 | 0.922207 (identical — internal control) | 1 |
| `e(pwsgs)` `att_i:1vs0` | 2.09663544 | **1.48700193** | 2 |
| its 95 % CI | [1.95302654, 2.24024435] — covers | **[1.33407235, 1.63993151] — excludes the truth by 8.3 SE** | |

No caution; `N = 1260 (0 obs dropped from initial sample)`.

### 4. Corrections to the finding as filed
- Their note *"att and att_it also read 1.16 in run B and that is CORRECT for them"* is true
  only because their design (and my `r01`) has exactly one row per (i,t) cell. It is **false
  in general** — see NEW finding **skp-01** below.
- Their noisy run D compared against `e(histogramestimates)`, which is the wrong oracle
  (see above); that is why their 0.9545 did not match their own 0.9486. `e(gteffects)`,
  which carries the `w_att_i` column directly, is the right one and matches to 4e-16.
- Everything else in the claim reproduces exactly.

---

## VERDICT fpsg-02 — **CONFIRMED** (severity kept at MEDIUM)

**My repro:** `r03_names.do`. Subgroup levels 1, 3, 7, 99; `e(b)` has 80 columns.

| name a user types after reading sthlp:1307 | `lincom` result |
|---|---|
| `main_att_s1`, `main_att_s2`, `main_att_it_s1`, `pretreatment_att_s1`, `ES_att_post0_s1`, `main_att_s99` | **r(111)** — all six |
| `main_SG1_att`, `main_SG3_att`, `main_SG99_att`, `ES_att_SG1_post0`, `pretreatment_SG1_att` | rc 0 — all five |

`e(b_es_att_s1)`, `e(b_es_att_s99)`, `e(escomp_s3)`, `e(numpre_att_s99)` all exist: the
`_s<value>` token is real, but it belongs to the e() NAME family, not to `e(b)`/`e(V)`
columns. sthlp:1307 is the one paragraph that tells `lincom`/`test`/`esttab`/`coefplot`
users how to address a coefficient, and it is wrong on both position and spelling.

**One correction to the claim:** "carries the subgroup VALUE, not a sequence number" is not
itself a defect. `#` is Stata help idiom for "a number", so `_s99` satisfies `_s#`. The
defect is POSITION (infix `_SG1_`, not a suffix) and SPELLING (`SG`, not `s`).

**Duplicate notice:** agent `hv-sg-names` filed the same divergence independently as
HV-SGNAME-01 with broader surface coverage (including `grep -c '_SG' wooldid.sthlp == 0`).
Verdict unchanged; the coordinator should merge.

---

## NEW skp-01 — `att_it` is cell-SHARE-weighted too when a cell is only partially labelled [HIGH / certain]

**Repro:** `r06_attit.do` (self-contained, noiseless, every target exact).

**What this does and does NOT contradict.** `sthlp:887` and wooldid's r(459) text say
"by-subgroup `att` and `att_it` are unaffected" — that sentence is about the cross-period
SPLIT under `safety(off)`, where the F13 within-cell guard has already guaranteed each cell
sits wholly in one subgroup, so `att_it` genuinely is unaffected there whatever the rows per
cell. **The help is not wrong about the split case.** What r06 contradicts is
(a) the prober's generalization of it — "att and att_it ... that is CORRECT for them" is
false as soon as a cell is partially labelled — and (b) the help's SILENCE: its one sentence
on missing subgroup values, "Observations with missing subgroup values remain allowed and
belong to no subgroup", says nothing about the cell those observations sit in being
down-weighted to its labelled share.

**When this fires.** The cell must hold several rows and only some of them carry the label,
with the labelled rows agreeing (otherwise F13 refuses outright). So the trigger is not
individual-level item nonresponse on an individual-level characteristic — that varies within
the cell and is refused — but a CELL-LEVEL classification carried on individual rows where
the row-level linkage or recording fails for some rows: employer sector on worker-year
records, county type on establishment records, plan type on claim lines. Administrative
linkage failure of exactly this shape is common. The geometry itself (several rows per (i,t)
cell) is the one wooldid's own help RECOMMENDS — a coarser `i` with unit fixed effects
through `fe()`; `known_bydesign` §5 notes wooldid steers users there above 100 distinct `i`.
The F13 within-cell guard folds missing through `cond(sg>=., ., sg)`, so the run passes.

Design: i=1 (sg1, effect 10) has 4 post cells of 5 rows, all labelled; i=2 (sg1, effect 2)
has post cells of 5,5,2,2 rows labelled on 1,1,2,2 rows → cell shares .2 .2 1 1.

| estimand | clean | partially labelled | documented target | share-weighted target |
|---|---|---|---|---|
| `att_it` | 6.0000000000 | **7.0000000000** | 6.0 (equal-cell) | 7.0 (dev **−1.3e-14**) |
| `att_i` | 6.0000000000 | 7.6000000000 | 6.0 | 7.6 (dev −9.8e-15) |
| `att_itime` | 6.0000000000 | 7.0000000000 | 6.0 | 7.0 (dev −1.3e-14) |
| `att` | 6.7058823529 | 8.1538461538 | (row-weighted) | 8.1538461538 (dev −1.1e-14) |

`att_it` error **+16.7 %**. 8 of 234 rows lost their label; no caution; `0 obs dropped`.

---

## NEW skp-02 — the PRE-treatment rows are share-weighted while `main` is exact, so one unit is weighted 1.00 and 0.25 in two rows of the SAME printed table [MEDIUM / certain]

**Repro:** `r02_mechanism.do` PART 2 (self-contained).

`w_per_i` is grouped by `(posttreatment, i)`, so the pre-treatment bucket has its own
denominator. Missingness confined to a unit's PRE rows therefore leaves `main:att_i` exact
and silently share-weights `pretreatment:att_i` / `att_itime`.

Three sg1 units with pre-window deviations 0.2, 1.0, 3.0 at relative years −3 and −2;
unit 3's subgroup label is unknown before t = 8.

| row | clean | label missing on unit 3's pre rows |
|---|---|---|
| `main:att_i` | 5.0000000000 | **5.0000000000** (unaffected) |
| `pretreatment:att_i` | 0.7000000000 | **0.2666666667** |
| equal-unit pre target | | 0.2000000000 |
| share-weighted pre target 0.6/2.25 | | 0.2666666667 (exact) |
| `ES_att:pre3` / `pre2` | 1.4 | 0.6 |

The printed subgroup table therefore describes two different sets of units under one
heading, and this is the surface that feeds `jointtests` pre-trend tests, `e(escomp_s#)`,
`e(numpre_att_s#)` and the honestdid handoff. Same root cause as fpsg-01, different reported
rows and different consumers; the original prober did not touch the pre-treatment side.

---

## SEVERITY EVIDENCE for fpsg-01 (NOT filed separately — same root cause) — the distortion reaches `e(pwsgs)` and its interval stops covering

**Repro:** `r05_contrast.do` (self-contained, truth known by construction).

Measured above under fpsg-01 §3. 4.6 % item nonresponse on the subgroup variable, with no
change whatsoever to the outcome data, moves the reported `att_i` subgroup contrast from
2.0966 (truth 2, CI covers) to **1.4870 with a 95 % CI of [1.334, 1.640] that excludes the
truth by 8.3 standard errors**, under a header reading `0 obs dropped` and with no caution.
The untouched subgroup's `att_i` is bit-identical across the two arms, which rules out any
estimation-side explanation. This is NOT filed as a separate finding — it is fpsg-01's own
mechanism measured on the downstream quantity a heterogeneity analysis actually reports, and
it is the measurement that fixes fpsg-01's severity.

---

## NEGATIVE RESULTS (attacked, clean)

- **Clean subgroups on an unbalanced panel.** `r02` PART 1 RUN A: a unit with 2 post rows and
  a unit with 8 both get weight exactly 1.00; `att_i` = 5.0000000000, dev 0.0e+00 from the
  equal-unit target. The clean path is exact, so fpsg-01 is not a general `att_i` defect.
- **The untouched subgroup never moves.** `r05`: `att_i(sg0)` = 0.922207 in both arms,
  bit-identical. The distortion is confined to the subgroup whose members lost labels.
- **MCAR missingness is nearly benign.** `r04` arm 2: 15 % MCAR gives −0.115 SE, −0.47 %.
  Anyone worried this fires on any missingness at all is wrong; it needs the shares to
  correlate with the effects, which is exactly what non-ignorable item nonresponse does.
- **`e(gteffects)` is honest.** Its `w_att_i` column reproduces the share structure exactly,
  so a user who digs that far CAN see the problem. Nothing in the printed output or in the
  headline e() scalars discloses it.
- **The `_SG` infix is internally consistent.** `r03`: every infix name resolves, every
  suffix name fails. There is no third scheme and no collision.

## Incidental, not filed
`r02` PART 2: `esprelength(4)` with `subgroups()` prints
`Caution: espostlength set to 1 due to use of subgroups in conjunction with esprelength`.
Disclosed by a caution, so not a defect — noted only because it silently narrows the
event-study window a user asked for.
