# Author rulings on the blue-sky slate — the working queue

Reconstructed 2026-09-20 from the ruling session of 2026-09-19/20 (the multiple-choice batches and
the free-text rulings that followed). `fix_proposal_status.md` says what adversarial review
concluded; this file says what the author DECIDED, and is the queue of record.

Status: `DONE` = built and verified; `TODO` = ruled, not yet built; `NONE` = ruled to take no action.

## Built

| # | item | ruling | state |
|---|---|---|---|
| R1 | `impvalues` IM geometry | escalate the caution to `r(459)` | **DONE** 2026-09-20 |
| R2 | `poisexpresults` / `semielasticity` / `contreatelasticitytype(eydx\|eyex)` | remove outright, not repair | **DONE** 2026-09-20 |
| R3 | replay of a pre-removal stored estimate | accept the mislabelled header; unsupported | **DONE** (no code) |

## Queue

| # | item | ruling | spec |
|---|---|---|---|
| ~~R4~~ | `BAND-1` | ship the built fix, stop there | **DONE** 2026-09-20. What shipped is **BAND-2b**, not BAND-1: BAND-1 substituted a floor value and the author rejected it outright ("absolutely never ever heuristically edit a SE"), and a package-wide audit found it would have been the only such site in 17,021 lines. BAND-2b is the replacement that was already built and verified at `probes/rev6-band2fix/` — one viability verdict per replication, non-viable draws withheld from the inference kernels by masking a COPY, the sampled values themselves untouched. Applied here as 24 hunks, reconstructed against the moved build by resolving the diff's append-anchors from pinned `c1d8979` and matching each change-hunk verbatim (all 22 change-hunks verified against pinned before application; the two that matched twice are the physical and frequency replication legs, which take the identical edit). **Its own honest limit, from `BAND_STATUS.md` and repeated here:** this removes the numerical absurdity (3.85e14 → 155) and does NOT make the band usable on a near-degenerate design — design H sits at 29.07 against 3.4-4.5 healthy with ZERO non-viable draws, which no per-replication verdict can touch. Smoke-verified on the current build: `boottstat` + `simultaneousbands` and a plain percentile run both rc 0, `e(boot_novar)` 0 on a healthy design |
| ~~R5~~ | `untreatedZero` | accept; re-bless the six. Engines **align down** | **DONE** 2026-09-20. Form B (clear the three inserters inside the untreatedZero branch) + Form C (one #-joined cell-interacted term per polynomial degree) + DEF-2 align-down (a dropped treated-cell dose SLOPE under untreatedZero is now fatal, reusing the existing bigcatastrophe path; `contreatcontrols` threaded into `wooldid_ccscreen` as `cc()`). **Verified pinned vs patched:** at poly1 the explicit engine moves 0.261604770 -> 30.227550127, landing bit-exactly on the default engine's value, which does not move; poly2 and poly3 agree across engines to 1e-10 and 2e-8 where pinned had them disagreeing wildly. DEF-2: poly3/poly4 went rc 0 + caution -> r(322) on BOTH engines, and poly2's caution disappeared entirely because Form C stopped creating the collinear term. safety(off) still returns rc 0 on all four |
| ~~R6~~ | `ELAS` | ship it | **MOOTED** 2026-09-20 by R2, and verified so. Every row in ELAS P1's table is a `semielasticity`/`eydx`/`eyex` arm: the defect is that reghdfe's `predict` omits the absorbed part, so an operator that DIVIDES BY THE PREDICTION differs from the explicit-dummy engine. The `dydx`/`dyex` controls showed error exactly 0 because they never divide by it, and all the dividing operators were removed. **Measured on the current build**, default reghdfe vs `unconditionalse robust`, worst relative gap over matched row NAMES (not indices - the documented verification trap): binary plain 1.6e-15, timetrends 4.1e-14, event study 2.1e-13, subgroups 1.6e-15, contreat dydx 6.7e-14, dyex 6.7e-14, dyex+timetrends 5.5e-14. `timetrends` carried a **13.12%** gap on the pinned build and is now noise. No code shipped; nothing left to repair |
| ~~R7~~ | `SMP-1` | ship it | **DONE** 2026-09-20. One insertion in `wooldid_bootrun`: rows outside `e(sample)` are removed from the resampling sample concept, which carries `wdb_C`, the nesting predicate, the re-id frame and `n_c` with it. Round 1's guard hardening adopted (a bare `gen byte = e(sample)` returns rc 0 with an ALL-ZERO marker when no sample is set, which would have died later with a misleading two-cluster message). **Verified with a seed-free oracle**, clean panel vs the same panel plus one wholly dead cohort: pinned moved `boot_df` 5->6, `boot_completed` 21->22 and the bootstrap SE by ~5%; patched is bit-identical on every cell. `N_clust` stayed 6 on both, so the defect was invisible in the reported cluster count |
| ~~R8~~ | `VCV-1` e(V) | rescale on lossy runs | **DONE** 2026-09-20, **shipped with VCV-2 as review required**. `wooldid_bootfin` runs two moment passes on DIFFERENT row sets — each standard error over the replications in which its own column was available, the covariance over the rows complete across every reported estimate — so `sqrt(diag(e(V)))` was not the printed standard error whenever any replication cell was lost, and the help's promise that a hand `test` or `lincom` reproduces the `bootdisplay(wald)` numbers was false there. The listwise diagonal is now rescaled to the reported column. Every mandated correction folded in: the **broadcast form** `V = (V :* vsc') :* vsc` rather than the outer product, so no second K x K temporary is allocated; the help says each **IDENTIFIED** diagonal, because a column with no listwise variance keeps the missing-variance sentinel and is deliberately not rescaled; and the regression check uses a **TOLERANCE, not `==`**. Rescaling by a positive diagonal is a congruence transform, so every listwise correlation and positive semi-definiteness are preserved exactly. VCV-2 ships alongside: **2a** the footer printed the literal text `e(boot_N_V)` and never its value (the line above it concatenates its macro outside the quotes, so this was a slip), and **2b** a caution plus a replay-surviving `e(boot_V_shortfall)` when the covariance rests on under half the completed replications — the package's own low-availability threshold, with the `e(df_r)` clause conditioned because `bootdf(normal)` blanks it. **Verified:** on a lossy run (196 columns, `e(boot_N_V)` 15 of 40 completed) `sqrt(diag(e(V)))` now matches the printed standard error to a worst relative difference of **2.9e-15** across 28 compared rows, and `e(boot_V_shortfall)` fires; on a clean run the count equals the completed count and no shortfall is posted |
| ~~R9~~ | bootstrap non-viable draws | **exclude AND caution**: drop draws at 1e-8, caution the user that X draws had SEs below 1e-5, reporting both counts. Never trim or shrink an SE to look reasonable | **DONE** 2026-09-20, layered on R4. Both thresholds come off ONE pass in `wooldid_bootviab`, which gained a `ctol` argument, and both are scale-free — read against the outcome's own standard deviation inside that replication's estimation sample. **The two cuts read different variances, on purpose.** Exclusion asks whether the replication produced any variance at all, so it reads the LARGEST diagonal. The disclosure is about the studentized tail, and a percentile-t statistic divides column by column, so it reads the SMALLEST positive diagonal: one collapsed column is enough to throw a t-ratio, and `BAND_STATUS.md` measured that as exactly how the band blows up (ordinary numerators, one denominator 33x to 120x too small). Columns margins dropped carry a zero or missing diagonal and are skipped rather than counted as collapsed. Nothing is trimmed, floored or rescaled; `ctol` only counts. New `e(boot_sesmall)` beside `e(boot_novar)`, and a caution that fires only when the looser cut catches something the stricter one did not, so the two counts are never reported twice. **Verified:** healthy design 0/0 and silent; near-degenerate design 5 excluded / 5 disclosed; and on an instrumented copy with the disclosure cut widened to 1e-1 (not shipped) the counts separate as designed — 30 vs 0 and 36 vs 5 — and the caution renders with both numbers. **The obvious objection to a min-based cut was raised in review and tested, twice.** A coefficient's standard error scales as sd(y)/(sd(x)*sqrt(N)), so a regressor measured in millions would carry a standard error at 1e-5 to 1e-8 of sd(y) on a perfectly sound replication and would trip a smallest-diagonal rule on every draw. Measured on the healthy design: adding `controls(big)` with `big = 1e6*rnormal()` leaves `e(boot_sesmall)` at 0 and the per-replication ratio BIT-IDENTICAL to the ordinary-scale run (relmin 3.4975e-02 in both) - controls never reach the verdict's `e(V)` at all. Scaling the DOSE, which does reach it, by 1, 1e3 and 1e6 likewise leaves the ratio identical to five figures (relmin 3.7991e-02). The metric sits three to four orders of magnitude above the cut on healthy data on both axes, so the disclosure has real headroom rather than resting on the fixtures' scale. **Worth the author's eye:** on every design measured the 1e-5 tier COINCIDES with the 1e-8 exclusion, because the continuous per-column near-degeneracy that actually drives the band sits orders of magnitude above 1e-5 (design H's decisive denominator was 3.05e-02). So the disclosure catches collapse, which is what was asked for, and not the near-degeneracy next door |
| ~~R10~~ | `GTE-1` | blank them | **DONE** 2026-09-20. Two SCALE-FREE floors in `wooldid_gtx` (`se < 1e-10*|est|`, and `v <= 1e-12*max(cell variances)`), a count in `e(gt_ndegen)`, a caution at the posting site and a note above the printed table (the note above the table was removed on 2026-09-21 at the author's decision, `1323749`; the caution and `e(gt_ndegen)` stand). Both rounds' text corrections folded in: the "robust covariance is rank-deficient" parenthetical is gone (false under an iid VCE and under an exact fit, so both mechanisms are now named), and `permute()` is out of the remedy list (it blanks these columns outright). **Verified pinned vs patched** on a 40-state/20-year/3-wave design: 10 of 39 cells went from `se ~ 1.5e-15` (printed as ordinary rows with p of order 1e-45 and zero-width intervals) to MISSING, min real se 1.57e-02, `e(gt_ndegen) = 10`. **Re-bless note:** under `boottstat` this also moves the p of every HEALTHY cell, because the studentized quantile changes once residue replications stop being counted |
| ~~R11~~ | `FIX-2` | hard `r(459)` for all four types; **tau fixed at ~1e-5**, widen the bracketing beyond one design | **DONE** 2026-09-20. Built from `FIX2_FINAL.md` — nothing existed. New `ccref` mode in `wooldid_amescreen`, generalising the shipped `bycoh` screen: every type inserts a baseline dose slope and differs only in what it lets that slope vary over (`basic` over nothing, `byTreated` over ever-treated status, `byTtre` over the onset value, `byTime` over the calendar period), and whatever the grouping the slope is identified only by dose variation among that group's `_t == 0` reference rows. Unconditional, with no `safety()` gate, matching the byCohort template the ruling cited. All four review corrections shipped: **R1** it fires under `contreatwithin` (verified 459 on two types), reversing the proposal — `contreatwithin` does not rescue the misspecification and byCohort's screen already fires there; **R2** the hard leg uses a RELATIVE tolerance of 1e-5 rather than `sd == 0`, and the `bycoh` call site moved with it so the two stay in step, with its near-tier text corrected to stop blaming `regtol()` (the estimate in that band is deterministic, not noisy); **R3** every requested type is screened, compound spellings included; **R4** the byTime message names the offending calendar periods and both measured remedies, and never says the dose lacks variation. **What the bracketing sweep did and did NOT establish — stated plainly, because the ruling asked for the bracketing to be widened and it was only half done.** On three independent designs `byTime` refuses at relative jitter 0, 1e-7 and 1e-6 and stops refusing at 1e-5. That is a fact about the PREDICATE and it is close to tautological, since 1e-5 is the threshold the screen codes; it is a useful consistency check that the relative yardstick behaves the same way across designs, and nothing more. What it does NOT do is reproduce `FIX2_FINAL` R2's substantive claim — that the reported estimate stays pinned to the exactly-collinear wrong answer out to a relative sd of about 1e-5 and becomes correct above it — which is the finding that would independently justify the constant. The fixtures used here do not recover their own truth at ANY jitter (the `basic` arms sit at -0.0127 against a truth of 1.0 below the cut and scatter to -54, -398 and -1080 above it), so they could not show a wrong-to-right transition and cannot speak to where the constant belongs. **1e-5 is the author's ruled value and ships on that authority**, with the single-design bracketing in `FIX2_FINAL` behind it; it is not corroborated here. **False positives: zero.** On an ordinary panel with a live comparison pool all seven specs return rc 0 (`basic`, `byTreated`, `byTtre`, `byTime`, `byCohort`, `untreatedZero`, and the compound `basic byTime`). On the degenerate panel `byTreated`, `byTtre`, `byTime` and `basic byTime` all refuse while `basic` and `untreatedZero` correctly proceed — `basic` is exactly right on these panels and `untreatedZero` is the remedy the message points at. Independent confirmation: the message named 5 affected periods on a design where `FIX2_FINAL`'s formula `T - max(onset) + 1` predicts exactly 5. **One cost note for the record, not a change:** like the byCohort screen it copies, `ccref` carries no `wdrepmode` gate, so it runs inside every `permute()` and `bootstrap()` replication rather than at rep 0 only. That is correct under the ruling — what is at stake is point identification, and a replication can be degenerate where the observed run is not — but on a large panel with `byTime` it is two `gegen` passes per requested type per replication |
| ~~R12~~ | `G3b` cohort id 0 | stop the classifier reading it | **DONE** 2026-09-20. One term deleted. Verified: pinned differs on 28 of 112 cells under cohort relabelling, patched is identical, and patched leg A lands on the pinned relabelled (correct) numbers |
| ~~R13~~ | `SMP-3` | **disclosure only** — rare, mostly user specification error | **DONE** 2026-09-20 as SMP-4, written from the proposal text (fixprop §0: "absent from all six builds", and "as written kills the run"). Sequenced after R7 for the stated reason: the straddle count is defined IN THE ESTIMATION SAMPLE, which is exactly what SMP-1 changed, so building it first would have pinned a count against a sample concept about to move. Three parts: `wdb_nstrad` computed at the existing nesting-predicate site in `wooldid_bootrun` (one `gegen tag` off the `wdbcmin`/`wdbcmax` already in memory — SMP-3 hunk (a) was going to supply it and was ruled down, so this is the proposal's own named fallback); a caution when it is positive; and `e(boot_reid)` + `e(boot_nstraddle)` posted from the SHELL after `wooldid_bootwire` (verified safe: bootwire and bootjnt both `ereturn repost`, never `post`, so existing entries survive). Round 2's correction adopted: the "`e(boot_nstraddle)` > 0 and `e(boot_reid)` == 0 is precisely the configuration FCM-1 describes" sentence is NOT shipped — round 1 showed it unreachable. The caution's 20-29% figure is sourced to FCW-01 / `mechanisms.md` ("one straddling unit out of 96 cuts the att_i bootstrap SE by 20-29%"), not to SMP-3's contested 10-14%. **Verified in both directions:** clean nesting panel → `e(boot_nstraddle)` 0, `e(boot_reid)` 1, silent; one straddler → 1, 0, caution naming 1; three straddlers → 3, 0, caution naming 3; and neither scalar is posted without `cluster()` or without a cell-weighted estimand. Four help edits: both statements of the two-regime dichotomy (sthlp:241 and the `bootstrap()` entry — fixing only one would have left the wrong rule documented at the other), `bootdf()`'s *G* now says "represented in the estimation sample", and the stored-results list gains both scalars. The `egen newi = group()` remedy carries `, missing` in both the caution and the help, in step with R16's two sites |
| ~~R14~~ | Roth / `PTP` | new work only; existing deviations stand | **DONE (no code)** 2026-09-20. The ruling resolves to nothing to build, and the check that it does is worth recording. wooldid departs from Roth in two SHIPPED places — `F(k, df_r)` where Roth uses chi-squared, and `invttail` where Roth uses a hard 1.96, the latter being the author's own earlier ruling BS-11 — and the gap is a small-G phenomenon: true size of a nominal-5% cutoff at k=3 is Roth 0.1301 / wooldid 0.0784 / the declined proposal 0.0500 at nu=15, but 0.0534 / 0.0513 / 0.0500 at nu=279. "New work only" leaves both in place. The three proposed `pretrendpower()` legs were already closed by the earlier "mimic Roth" ruling, and Leg 3 (the singleton-cell caution) was ruled DROP and sits in the take-no-action table below — note it would have newly cautioned a large share of ordinary runs, since `i` = unit is the default call. No site in this slate adds new Roth work, so nothing here changes |
| ~~R15~~ | `SG-1` phantom rows | cheap conjunction detector | **DONE** 2026-09-20. **Scope corrected first:** an earlier draft of this row tied R15 to the bootstrap's phantom CLUSTERS and resequenced it behind R7. That was a conflation — the ruled question was `fastnose` + `contreat()` + `subgroups()`, where fast mode PRINTS event-study rows the full run drops as non-estimable. Phantom ROWS in the displayed table; independent of R7. Spec: `author_followups.json` entry 7. Built at `ado:3456`: under `nose` a column is condemned iff BOTH its value is zero to 1e-12 AND no cohort in that (subgroup, horizon) cell with nonzero estimand weight has a live treated-cell dose slope in `e(b)`; the verdict is written as an `o.` into the colname, which `wooldid_sglbl` already keys the SI-7 filter on and slices `tab`/b/V by, so nothing downstream changes. Liveness is computed ONCE per `horizitem` (a done-flag local, not `sgcase == first of sglist` — the `cap confirm matrix` above can skip the first subgroup) via one Mata pass, `wooldid_sg1live`. **Verified pinned vs patched on the four designs that discriminate:** on D1 the pinned fast run published `pre2 pre1 post0 post1 post2` for sg1 where the full run has `post0 post1 post2`, and three columns for sg2 where the full run has one — patched fast matches the full run exactly on both. D6 (a subgroup mixing a live-slope cohort with an all-dropped one, where the naive any-dropped-term rule false-positives) and W0/W1 (`controls(c.dose#c.z)`, where the liveness half alone false-positives) are UNCHANGED on both builds — no estimable row suppressed. Full mode is bit-inert. A `permute()` + `contreat()` + `subgroups()` run returns rc 0 on both builds; that design does not trigger SKEP-02's rep-0 abort, so the permute leg shows no regression rather than a repair |
| ~~R16~~ | split units | refuse with `r(459)` | **DONE** 2026-09-20. SG-2 with round 2's domain, not round 1's: the screen reads the relative-period variable (`_t != 0`, i.e. rows carrying a treatment dummy), because the `t >= ttre - esprelength` window form was measured to re-open a silent wrong number under `esfixedbaseperiod`. `relt()` and `rm()` threaded into `wooldid_sgdom`; reuses `f13sgimin`, so no extra gegen. Residual admission stated in-code: a label missing only on `_t == 0` rows is admitted. **Verified:** a partially labelled unit went rc 0 (silently share-weighted) -> **r(459)**, while fully labelled, WHOLLY unlabelled, and missing-only-on-never-treated all still pass, and safety(off) still proceeds. Also fixed the pre-existing help/code remedy: bare `egen group()` blanks rows and drops them from the sample, so both sites now say `, missing` and note that splitting `i` redefines the estimand |
| ~~R17~~ | `level()` | add wooldid's own refusal | **DONE** 2026-09-20. L1+L2+L3 with both rounds' corrections. Verified: 8 legal `level()` values run with exact `e(level)`, all 8 session levels run (the 1768-level bug), 6 illegal values give r(198), a stale `WOOLDID_LEVEL0` no longer leaks |
| ~~R18~~ | `G1c` | tighten to the measured allocation | **DONE** 2026-09-20 — **as G1a + G1c**, because the ruling was given on an option whose text reads "Ships with the primary G1a fix", and the frame decides the answer: read through the pinned bound the overhead is `si8ilev + si8tlev + 200`, read through G1a's it is just 200 on the absorbing path, and the same break point means a different multiplier in each. The bisection the proposal asked for was run in full — `fixprop/G1C_MEASUREMENT.md`. Findings: the axis is treated cells, not the event window (300 cells in a 2000-wide window runs; 1700 in the same window dies); units are not allocated at all on the absorbing path (800 cells runs at 140, 600 and 1200 units); the break sits between si8cells 900 (runs) and 927 (r(907)) on reghdfe and ppmlhdfe alike, bracketing the multiplier in (1.85, 2.05] and the constant in (194, 248]; and `contreat()` breaks in the SAME place, so the factor is a **flat 2 replacing** the `(1 + contreat)` doubling, never multiplying it — charging 4 per cell would refuse a contreat run at 500 cells that completes. `+200` is kept and must be: at or below 194 it under-charges. The reg/poisson branch is left exactly as it was, because its arms are governed by two constraints rather than one product and there is no single measured allocation to tighten to; a residual false positive there is recorded rather than patched (800 cells + 1200 units: guard refuses at 2260, `safety(off)` completes). **On the landed build:** the G1 headline (5,200 units, maxvar 5000) returns rc 0 in 1.18 s where pinned refused with "at least 5458", `e(wdidmainresults)` bit-identical in `%21x` to the `safety(off)` run, and the `customvce` leg still refuses. Zero false positives across twenty absorbing arms; four late r(907)s become an early named refusal and two false refusals are removed |
| ~~R19~~ | perm cellthresh 2a | two string edits; **leave the statement qualitative**, no number | **DONE** 2026-09-20. Six sites, not two: 3 runtime cautions + 3 help passages. Scope ruled in-code; the continuous-treatment ame remedy is explicitly OUT (different statistic). Verified the caution renders with the qualifier |

## Ruled: take no action

| item | ruling |
|---|---|
| perm-inert | do nothing |
| structural zero | leave it |
| `jointtests` | leave it alone |
| `PTP` Leg 3 | drop it |
| clean margins operators (`dydx`/`dyex` linear, and their exponentiated-scale counterparts) | leave them entirely alone |
| help split | **DEFERRED** — the chunk sits below the "DO NOT WORK ON CHUNKS BELOW THIS LINE" marker at `wooldid_ado_review.md:8`. Notes kept; work deferred on the author's instruction |

## Landing requirement — DISCHARGED 2026-09-20

Run pristine TWICE (`test/run_batt_A`, `test/run_batt_B`), and the two passes compare
bit-identical to each other at `--rtol 0`, so nothing below is draw-stream noise.

Against the previous baseline there are 110 difference groups and they decompose exactly,
with **zero lines unaccounted for**. The classes, named before the compare output was opened:

| class | what | where |
|---|---|---|
| C1 | `macro wdid_semielas` REMOVED | 58 configs (R2). The 59th never reaches that macro |
| C2 | `scalar boot_novar` / `boot_sesmall` ADDED, both 0 | the 16 `boottstat` configs (R4, R9) |
| C3 | `scalar boot_reid` / `boot_nstraddle` ADDED | 10 clustered cell-weighted bootstrap configs (R13) |

**No number moved anywhere** — no matrix cell, no existing scalar, no `rc`, no `varlist`.
Five classes that were expected and did not appear are therefore also recorded, because their
absence is information: no battery fixture holds a degenerate `untreatedZero` reference pool
(R5), a `gteffects` residue cell (R10), a dead cohort or singleton in a bootstrap frame (R7),
a `fastnose` + `contreat()` + `subgroups()` arm (R15), or a bootstrap run lossy enough to move
`e(V)` (R8). The last is the load-bearing one: **VCV-1's rescale is exactly inert on clean
runs**, which the proposal asserted and this measures, since every clean bootstrap config's
`e(V)` is bit-identical rather than differing in the last places.

One internal cross-check falls out of C3 and is worth keeping: the three configs reporting
`boot_nstraddle = 3` are exactly the three reporting `boot_reid = 0`. The new disclosure
identifies which runs took the keep-`i` rule, on the project's own fixtures, without being
told.

Not rerun here, and not claimed: the cert suite, `oldnew`, `coverage` and `degeneracy`. The
user asked for the battery. `configs_ext.csv` carries `byTime` (x3), `byTreated` (x2),
`byTtre` (x1) and two compound `contreatcontrols` spellings whose `rc_new` may now be 459
under R11; those are `oldnew` rows and the next full run should adjudicate them.

---

# The non-number-moving batch (2026-09-20)

The nineteen rows above are what the author RULED on. The blue-sky review produced 56
proposals in all; only 13 moved a number and so only those needed rulings. This section is
the rest: items that adversarial review verified, that move no reported number, and that were
therefore never put to the author. They are built here.

## Built

| item | what it does | mechanism |
|---|---|---|
| **LAT L1** (4 edits) | the all-zero-covariance cautions state the cause that actually applied | caution text |
| **LAT L4** | the `lattice()` discussion gains the support requirement, distinct from the functional-form threat already described | help |
| **VCV-3** | `unconditionalse` is tested FIRST in the vce ladder, so `e(vce)`/`e(vcetype)` name the stage the reported SEs come from and `e(vce_stage1)` names the one underneath | `e()` |
| **GTE-4** | `e(gt_ncontained)` counts `e(gteffects)` cells that lie wholly inside one cluster, with a caution | `e()` + caution |
| **BSAVE-1** | the studentized path posts `e(boot_jointWobs)`/`e(boot_advWobs)`, and `bootsave()` reads them from `e()` rather than from the worker's locals | `e()` |
| **PTP Proposal 2** | singleton-cell detection behind the rank-deficiency caution, recorded in hidden `e(roth_cellmax_*)` | `e()` + caution |
| **G2** | the identifier-domain assert keys on `safetyinput`, not the `safety` that CS-98 rewrote to `off` under `contreatwithin` | screen |
| **SMP-2/BSAVE-3** | the small-cluster caution moves from `wooldid_bootctl` to `wooldid_bootrun`, where it can key on the realized `smp` | caution siting |
| **G1b** | r(907) is named as Stata's variable limit rather than left bare | error text |
| **G4** | `%7.0f` -> `%7.0g` for the relyr/numpre column | display format |
| **FIX-PERM-3** | the `permute()` nesting refusal names the smallest offending cluster and explains the design-vs-estimation sample distinction | error text |
| help | BAND-3, BAND-1b, DOC-1 (x3), IMP P2 (x3), GTE-5 (x3), plus stored-results entries for every new `e()` above | help |

## Closed without building

* **FIX-PERM-2a** — already shipped. Both `sthlp:729` and the runtime caution at `ado:1672`
  already carry the qualification ("...but only where that finer clustering leaves enough
  distinct clusters inside each treated cohort-by-period cell"). Moot.
* **FIX-PERM-2b** — the numeric clusters-per-cell predicate. Ruled out by the author's standing
  instruction to leave that statement qualitative.
* **ELAS P2**, **ELAS P3** — moot after R2 removed the options they concerned.
* **GTE-2** — shipped via R10.
* **BAND-1b** items 1, 3, 4 — already true of the shipped text.
* **LAT L5** — a "do not fix" with new evidence: the `regress` refit removed on 2026-09-08
  should stay removed, because s0.3 measures its answer at `se/sd ~ 0.03`. Its second half
  (stop steering readers to `customvce(robust)`) is what L1 edit 4 implements.

## Held for ruling

**LAT L3**, **FIX-PERM-1**, **GTE-3** — each adds a screen whose firing is a new event in a
user's log, and the evidence behind each is the proposer's fixtures rather than a measured
rate on ordinary data. Written up for a one-line ruling in
`HELD_FOR_RULING_2026-09-20.md`.

## LAT L1 — all five branches exercised

The proposal's fixture was not preserved, so it was rebuilt from its description (counties
inside states, dose near-constant within county, `i` = state, so a 3-bin lattice cuts treated
cells into 1-2 observation pieces).

| arm | expected | observed |
|---|---|---|
| contreat + lattice + `robust` | new unclustered text WITH the lattice clauses | as expected |
| contreat, no lattice, `robust` | new unclustered text, no lattice clauses | as expected |
| contreat + lattice + `cluster()` | R1 fires, vzero defers, R4 points at a caution that really printed | as expected |
| contreat + `cluster()`, no lattice | branch 1, the original clustered contreat text | **not reached** - R1 fired here too |
| `timetrends` + `cluster()`, no contreat | branch 2, the original clustered text verbatim | as expected |

Branch 1 is the shipped text, unchanged; the only edit to it is the added `clusvars != ""`
condition on its gate, which is the fix. No fixture reached it, and that is recorded rather
than papered over.

The bug edit 4 repairs is visible on the two bare-`robust` contreat arms: the old text sent
the reader to a cell-containment caution that was never printed, then offered "a
finer-than-cell clustering" to a user who had clustered nothing and "conditional robust
standard errors" to a user who had just asked for them. Both are gone.

**Deliberate change of substance, named so it can be reversed:** edit 4 also reorders R4's
remedy list, demoting conditional robust standard errors below `bootstrap()`/`permute()`.
This touches the AMESE R4 message ruled on 2026-08-24. The evidence is the frozen-design Monte
Carlo in `fixprop/proposals_LAT.md` s0.3, which puts the conditional robust standard error at
0.017 of the estimator's own conditional sampling sd on the design that prints this caution.
It is one local (`wdr4alt`) if the author prefers the old order.

## The R11 `configs_ext` follow-up — CLOSED, nothing to change

`1c09efa`'s message named a loose end: `configs_ext.csv` carries eight `contreatcontrols`
rows — `byTime` (x3), `byTreated` (x2), `byTtre` (x1) and two compound spellings — whose
`rc_new` might now be 459, because R11's `ccref` screen is unconditional, is **not** gated on
`safety()`, and is deliberately not carved out under `contreatwithin`. Six of the eight carry
`safety(off)` and two carry `contreatwithin`, so none of the usual exemptions would have
saved them.

Run on the current build through the suite's own runner (`mkrun.py` on a filtered
`configs_r11.csv`, both shards clean):

| config | recorded `rc_new` | observed |
|---|---|---|
| ext_cont_ext_controls_bytreated | 0 | **0** |
| ext_cont_ext_controls_byttre | 0 | **0** |
| ext_cont_ext_controls_cohort_time | 0 | **0** |
| ext_cont_ext_ccnoabsorb | 0 | **0** |
| ext_char_within_bytime_off | 0 | **0** |
| ext_cov10_cont_within_vary | 0 | **0** |
| ext_val_untreatedzero_combo | 198 | **198** |
| ext_val_within_bytime_safety | 198 | **198** |

Eight of eight match. R11 fires on none of them, because the reference pools in
`continuous_control_compact`, `continuous_control_dense` and `binary_balanced` all carry
usable dose variation — which is what the screen asks about. **No `configs_ext.csv` edit is
needed and no ledger row is owed.** The two 198s never reach R11 at all; they are refused
earlier, at validation.

## R5 / DEF-2 — SCOPE CORRECTION (2026-09-20), a regression `1c09efa` shipped

**Correction to the R5 row above.** DEF-2's "engines align down" was built as a single line in
`wooldid_worker`'s multicollinearity block:

```stata
if "`cc'" == "untreatedZero" local bigcatastrophe 1
```

That makes **every** dropped interacted term fatal under `untreatedZero`. The ruling, and the
comment directly above the line, say a dropped treated-cell **dose slope**. The gap inverted the
screen:

| run | dose in treated pre-periods | `1c09efa` |
|---|---|---|
| **obeys** `untreatedZero`'s premise | zero, as documented | **r(322), refused** |
| **breaks** that premise | nonzero | rc 0, allowed |

Honouring `untreatedZero` makes the pre-period dose columns structurally zero; they are dropped
as arithmetic, not as a rank failure, and the guard could not tell the difference. It fires on
any compliant `untreatedZero` run that requests pre-treatment horizons (`att`-only and
post-only runs are silent). It also contradicts `sthlp:889` ("expected zero doses in untreated
periods ... do not trigger this") and the R6 header below it, which claims `targetflag` scoping
exempts those cells — untrue under `esfixedbaseperiod`, where they carry `1o.targetflag`.

**Attribution.** `f8b4a64` rc 0; `1c09efa` rc 322; removing only that line from the current
build returns rc 0 with `ame b = -.03330222, se = .06590857`, bit-identical to `f8b4a64`.

**The fix** narrows to the ruling's own words: (a) the dropped token is a dose term, via the
`#`-terminated idiom the R6 loop uses, and (b) the flagged cell carries a nonzero dose.

**Verified in both directions, because half of it is the dangerous half.** Un-refusing the
compliant case proves nothing on its own — a screen narrowed until it fires on nothing is a
silent revert of the ruling, and no battery config is on this path.

* **Positives all hold.** 22 `untreatedZero` arms on the adversary's own `d8`/`d9` fixtures
  (recovered, and their `datasignature` re-verified against a rebuild from `mkdata2.do`/
  `mkdata3.do`) return the identical rc on the narrowed build, on `62b121c` and on `1c09efa`:
  r(322) everywhere, on `customvce(ols)`, `customvce(robust)`, `unconditionalse robust` and
  `contreatpoly(2)` — precisely the arms the ruling was written about. `safety(off)` reproduces
  the documented truths exactly (d8 `att` 2.254205175163; d9 `att` 1.314259241731,
  `ame` 0.729333333333, the documented 0.066 deficit).
* **Condition (b) is the whole discriminator**, measured per flagged cell rather than argued:
  the refusing cells carry 40 of 40 rows with nonzero dose; C2's 60 released cells carry 0. All
  60 of C2's tokens *are* dose terms, so condition (a) alone would have refused it too.
* **The unit conversion is load-bearing and points the right way.** `d9` hangs on one token; had
  the count written `+ `min_t'` rather than `+ `min_t' + 1`, it would have read 0 and released a
  case the ruling requires refused.
* **`test_amescreen.do`: AMESCREEN_ALL_PASS**, 16 of 16.

**A second class is released, named here rather than left to be discovered.** A dropped
treated-cell slope on an **interactive control** (not the dose) under `untreatedZero` now
returns rc 0 with the generic caution, where `1c09efa` refused. This is inside the ruling's
words, and measured on the pre-ruling build `c1d8979` those arms are rc 0 — so the correction
*restores* prior behaviour rather than inventing new behaviour, and `1c09efa`'s refusal there
was collateral. If the author wants interactive-control slopes fatal too, condition (a) is the
one line to widen.

## Two pre-existing defects found while verifying the above — NOT fixed here

1. **`lattice()` with `verbose` off disables the DEF-2 test entirely.** The outer gate reads
   `wordcountproblem <= 3 | "`lattice'" == "" | "`verbose'" != ""`, so with a lattice active and
   `verbose` off the `> 3` branch never runs. Measured identically on both builds: `d8`/`d9` +
   `untreatedZero customvce(robust) lattice(byCohort 2)` returns rc 0 with exactly the DEF-2 and
   DEF-3 outputs the guard exists to refuse. **R5's "every engine refuses where the default one
   does" therefore does not hold under `lattice()`.**
2. **`verbose` changes a return code.** The same command is rc 0 without `verbose` and r(322)
   with it. A diagnostic option should not decide whether a run errors.

Both predate the scope correction and are untouched by it.

## Process: what let this ship

`1c09efa` and `e2b5fad` were verified on the regression battery alone. The battery captures
`e()`, `rc` and varlists for its own 59 configs and **never executes a `test_*.do`**, so a red
certification suite survived two commits and one "does anything remain open" audit. `test/`
holds 43 targeted cert suites and `test/suite/run_full.sh cert` runs all of them in about 23
minutes. It should have been run at both checkpoints, and is run now.

## DEF-2 out of the anti-spam gate, and interactive controls back in (ruled 2026-09-20)

Two author rulings on the two items the scope correction left open.

**1. "You shouldn't put a hard gate in an unrelated anti-spam device. Let's fix that."** The
DEF-2 refusal sat inside the `wordcountproblem > 3` display branch, whose
`"`lattice'" == "" | "`verbose'" == "verbose"` gate exists only to stop a meaningless caution
repeating on lattice runs, where dropped interacted terms are routine. The refusal inherited
that gate. The refusal is now its own branch, gated on `untreatedZero` alone; the advisory
keeps its gate. The full-data `count` above is widened for `untreatedZero` only, so no other
run pays, and the comment justifying that skip is corrected, since it no longer holds.

A second gain from the move: at the new site `problem_t` is still in internal `_t` units, the
same units the count on the line above uses, so the `+ `min_t' + 1` undo of the display
branch's relative-horizon conversion **is deleted outright**. That was the step whose
direction, if reversed, silently released `d9`.

**2. "We can widen things for that untreatedZero + interactive controls case."** A dropped
treated-cell slope on an interactive control is fatal again. The match is against `interlist`,
the demeaned tempvars that actually enter the model — `interactivecontrols` holds the user's
raw names and would never appear in a dropped token, so matching on it would have made the
widening a silent no-op.

### Measured, targeted (no full battery for this change, on the author's instruction)

| arm | before | after | meaning |
|---|---|---|---|
| `d8` untreatedZero customvce(robust) **lattice(byCohort 2)** | rc 0 | **322** | the fix |
| `d9` same | rc 0 | **322** | the fix |
| `d8` same **+ verbose** | 322 | 322 | `verbose` no longer decides the rc |
| `d8`/`d9`, no lattice, customvce(robust)/customvce(ols)/unconditionalse/default | 322 | **322** | hoist broke nothing |
| `dz` untreatedZero + interactivecontrols(z) | rc 0 | **322** | ruling 2 |
| **`d6` healthy + lattice(byCohort 2)** | rc 0 | **rc 0** | no false positive |
| **`d6` healthy + lattice(all 3)** | rc 0 | **rc 0** | no false positive |
| `d6` healthy, no lattice | rc 0 | rc 0 | unchanged |
| `d6` compliant untreatedZero + esf pre-horizons (the C2 shape) | rc 0 | **rc 0** | the scope correction still holds |

The false-positive question raised when this was proposed — a dropped dose slope under a
lattice may be innocuous, since the remaining bins can still identify the slope — is answered
on the fixtures available: **zero false refusals** across three healthy lattice arms. That is
not a proof over all lattice designs, and it is stated as what it is.

### Battery: predicted IDENTICAL, 0 groups, before the compare

Both `untreatedZero` battery configs (`cont_uz`, `cont_uz_ss`) specify `safety(off)`, and the
screen's whole call site is gated `if "`safety'" != "off"`, so neither reaches it. No battery
config pairs `untreatedZero` with `lattice()`. The change therefore cannot move a battery cell.
