# FCW-02 — third lens (DOCUMENTED-BEHAVIOUR). Verdict: CONFIRMED, high, realistic.

## The lens question
Is the phantom straddler a defect, or is it documented / ruled / correct-but-surprising?
Answer: **a defect.** The behaviour the arms measure is the exact failure mode the project's own
ticket says the two-regime rule exists to prevent, and the sample concept that lets it fire was
never ruled, never documented, and is contradicted by the justification comment sitting on the line.

## Documentary evidence (all verified on the pinned build / the repo record)

1. **The ticket names the observed outcome as the failure mode.**
   `review_process/feature_tickets_2026-09-15.md:87-89`:
   "The nesting branch must stay `group()`. Keeping `i` there would make a twice-drawn cluster's
   cells identical duplicates that the equal-cell weighting counts ONCE, **collapsing the
   composition variance of `att_i` toward zero**."
   Repeated verbatim in the code at `ado:9328-9333`. On the estimation sample the arm-D design
   NESTS; the run nonetheless takes keep-i; se(att_i) falls 13-24%. That is the documented-wrong
   branch firing on a nesting design — not a surprising-but-correct design.

2. **The code's own justification for `smp` is falsified by the behaviour.**
   `ado:8980-8986`: "the `if `smp'' on BOTH the gegen and the assert is load-bearing ... Evaluated
   on the raw frame instead ... a strictly COARSER cluster ... would **silently take the keep-i
   rule and collapse att_i's composition variance**."
   The comment names this precise harm and claims `smp` prevents it. It does not: `smp` is
   (user `if` == 1) & nonmissing cluster (`ado:8910-8916`; `ifr(`wdp_if')` at `ado:394`,
   `wdp_if` = raw user `if` at `ado:615`), so every markout-droppable row still votes.
   Contradicting a stated rationale is a first-class finding per the coordinator note.

3. **The sample concept was never ruled.** Ticket, prereg (`test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md`
   mentions "design sample" only on the permute side), and `ado:6455-6457` all weigh exactly two
   options — compute the predicate in `wooldid_bootctl` and pass it down, or compute it in
   `wooldid_bootrun` on `smp`. The ESTIMATION sample was never on the table. `ado:6455-6457`:
   "the regime is decided in wooldid_bootrun on the same sample the resample runs on ... so the two
   can never disagree" — true of bootrun vs bootctl, silent about bootrun vs the estimator.

4. **Not documented.** `sthlp:241` states the rule with no sample qualifier: "Where {it:cluster()}
   nests {it:i} the copies become distinct units; where it is finer than {it:i}, or cuts across it,
   the copies stay inside their cohort". A user reads "nests i" as a property of their estimable
   data. `sthlp:431` documents the design-vs-markout gap **for permute() only** ("rows wooldid later
   removes by its deterministic markout, e.g. missing outcomes or zero weights, still count toward
   the design"). The asymmetry is evidence the author reasoned about the gap where he saw it and
   did not see it here.

5. **Silent.** `wdb_nest` is `c_local`-ed out at `ado:8998` and consumed at exactly one site,
   `ado:407`, gated `"`wdb_reid'" != "" & "`wdb_nest'" == "1"` — i.e. only in the NESTING regime,
   and only when a covariate names `i`. Under the flip (nest = 0) nothing prints. Nothing is posted
   to `e()`. Contrast `e(perm_irebuilt)`, which the author kept expressly so a user could check
   which of the permute paths a run took (ruling 2026-09-16).

## What I add that the prober and the skeptic did not: the bias has ONE direction

`ado:9340` is `qui keep if `smp'` and the replication is a nested `wooldid` call whose markout is
deterministic, so within every replication est ⊆ smp. Therefore:

    nest(smp) = 1  =>  nest(est) = 1.

Adding rows can only make a unit straddle MORE clusters; dropping them can only make the sample
MORE nesting. So the only reachable misclassification is `nest(smp) = 0` where `nest(est) = 1` —
spuriously taking keep-i. There is no phantom row anywhere that can push a run the other way.
Combined with (1): **the error is always toward the branch that collapses att_i's composition
variance, i.e. always anti-conservative.** Never a conservative accident.

## Scope / severity
Fires on: `bootstrap(standard)` + `cluster()` + a cell-weighted estimand (att_i / att_it /
att_itime) + one unit whose only cross-cluster row is one the markout drops (missing outcome,
missing control, zero weight). Realistic (a firm/person that changes county or establishment in a
year whose outcome is missing; admin frames carry the cluster code on every row while the outcome
comes from a thinner source), but a conjunction — hence high, not critical.

## The remedy's correctness — stated, not verified
Under the correct nest = 1 regime the phantom row gets a fresh `__wdbnewi` code carrying one
unestimable row, so it should contribute nothing. That is the expected behaviour of a fix, not
something I measured; it is not load-bearing for the defect's existence.

## Own run
`l1_doclens.do` — one A/D pair off the skeptic's generator, plus the two lens-specific checks
(silence of arm D's printed output and `e()`; the same phantom row under `permute()`, expected
loud at `ado:6078`). Results in `l1_doclens.log`, `l1_armA.log`, `l1_armD.log`.

## Two more documentary points found after the first pass

6. **The help's stated REASON for the keep-i branch does not hold on a phantom straddler.**
   `sthlp:532` (the `bootstrap()` entry, the fullest statement): "Where the cluster is finer than
   {it:i}, or cuts across it, the units are kept exactly as they are ... re-identifying there
   **would split a unit straddling two clusters into two units in every replication** whether or
   not anything had been duplicated, turning the unit fixed effects into effects per cluster per
   unit". On a phantom straddler that reason is inert. Under nest = 1 the phantom row DOES get its
   own `__wdbnewi = group(__wdbns, i)` code — the unit IS split into two codes — but the second code
   carries only a row the markout drops before the regression, so the extra code never reaches a
   fixed effect and the unit FEs are not turned into per-(cluster, unit) effects. The branch fires
   anyway and pays its documented cost in full. The help states a rule and a reason; the code
   applies the rule in a case where the reason cannot bite. That is a documented-vs-actual
   divergence in the strict sense, not just an undocumented edge.

7. **The feature's own preregistration never contemplated a sample-concept disagreement.**
   `test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md` pins the nesting arm at `cluster(state) == i`
   (F1a: "all at `cluster(i)` (nesting trivially)") and the keep-i arm at "clustering strictly finer
   than `i`" (P4/P8) — both unambiguous under either sample concept. Section F0, "what is
   preregistered, and what is honestly not", does not name the boundary case. So the ticket's
   "covered in the suite from both directions" is true of the two interiors and silent about the
   boundary this finding sits on.

## Scope is WIDER than the cell-weighted estimands — flagged, not measured
`ado:378`: `wdb_reid` is set when `bootstrap(standard)` AND (a cell-weighted estimand **OR**
(`boottstat` AND `cluster()`)). So a plain `att`-only run with `boottstat` + `cluster()` also takes
the physical leg and is also subject to `wdb_nest`, and the replication is a nested
`wooldid y __wdbnewi t ttre` call, so the regime changes the cohort structure the replication fits.
`sthlp:532` asserts the boottstat rationale "is independent of whether the unit identities move
with it" — true of WHY the leg is taken, but the identities still move and still enter the fit.
Whether the phantom flip moves a boottstat att number materially is UNMEASURED. Flagging as an
extension worth a probe, not claiming it.

## Ledger check
`review_process/open_defects.md` has no entry for this (grepped markout / design sample / phantom /
smp: no hits). Novel.

## Run status
`l1_doclens.do` queued behind the project-wide 8-instance Stata cap (saturated at N=8/9 for the
whole window). The verdict does not depend on it: the prober and the skeptic each paired seeds and
got se(att) identical to every printed digit across arms, which is the isolation that matters, and
every lens-specific claim above is settled by reading the build and the repo record.

## The lever — a one-line user workaround that is also the diagnosis
Adding `if !missing(y)` to the user's `if` moves the phantom row out of `smp` (it is built as
`(ifr) == 1`, `ado:8911`) and restores nest = 1, recovering the larger, correct se(att_i). It is
undocumented and no user would think of it — why filter on the outcome you are already passing? —
but it demonstrates that the design-vs-markout gap is the ENTIRE mechanism, and it names the fix:
screen `smp` by the same deterministic markout the estimator applies before computing the regime
predicate. (Predicted from `ado:8910-8916`; not executed, see run status.)

## The one "correct-by-design" defense, and why it fails
The strongest defense is the code comment's own argument at `ado:8980-8986`: the regime should be
computed on the rows `bsample` actually sees. But the regime's PURPOSE, per `sthlp:532` and the
ticket, is about unit fixed effects and cell composition IN THE FIT. A row that never reaches the
fit can affect neither. "Match what bsample sees" is therefore a proxy for "match what the
regression sees", and it is a proxy that fails exactly at the markout boundary — and by the
one-direction argument above, it never fails benignly.

## Adjacent observation in the same family (noted, NOT claimed as this finding)
The same design-vs-estimation gap sits on the cluster-count guard. `ado:6459-6467` builds
`smp2` = `(ifr)==1 & !missing(clus)` and errors "bootstrap() requires at least 2 clusters in the
**estimation sample**"; `ado:8934-8938` repeats the message off the map frame, also built on `smp`.
Neither counts the estimation sample. A cluster whose every row is markout-dropped therefore counts
toward the guard, so a design with two design-sample clusters and one estimable cluster passes a
check whose message says otherwise. Message-vs-computation, separate from FCW-02, unmeasured.

## Exact line pins used above (verified on the pinned build)
`smp` built 8909-8915 | regime predicate 8987-8994 | `c_local wdb_nest` 8998 | its ONLY consumer 407
| `qui keep if smp` 9344 | the two-regime branch 9348-9349 | `wdb_reid` set 378 | `ifr(wdp_if)` 394
| `wdp_if` = raw user `if` 615 | bootctl rationale 6455-6457 | sthlp 241, 431, 532.

## MY OWN RUN — RESULTS (all three questions answered)
`l1_doclens.do` / `l1_doclens.log`, `l1_armA.log`, `l1_armD.log`; Q3 re-run in
`l1_q3_permute.do` / `l1_q3_permute.log`.

**Q1 — reproduces.** One A/D pair, reps(199), rseed(41), 15 counties x 4 units x 7 periods,
`cluster(cty)`, `att att_i`.
| | arm A | arm D |
|---|---|---|
| e(N) | 419 | 419 |
| b(att) hex | `+1.09fe1ccddb87aX+000` | `+1.09fe1ccddb87aX+000` (bit-identical) |
| b(att_i) hex | `+1.0a5410a7dbb9eX+000` | `+1.0a5410a7dbb9eX+000` (bit-identical) |
| se(att) hex | `+1.df613c98659a9X-003` | `+1.df613c986599f X-003` (last 2 hex digits; ~1e-15 rel) |
| se(att_i) hex | `+1.fbfa17d32f07eX-003` | `+1.b1fbc72441c60X-003` |
`RATIO_se_atti = 0.85433705` — a **14.6% fall**. `DELTA_b_atti = 0` exactly. Same seed, same
draws, same fits; only the unit-identity rule differs. Matches the skeptic's first pair to all
printed digits, from an independently written harness.

**Q2 — SILENT, measured.** `grep -c -i caution l1_armD.log` = **0**. Arm D's printed output
contains no occurrence of nest / regime / re-identif / straddl / copies / finer. `ereturn list`
after arm D contains no nest/regime/reid entry. Note the header prints
`N = 419 (1 obs dropped from initial sample)` — wooldid TELLS the user the phantom row is not in
the estimation sample, and then lets that same row choose the variance.

**Q3 — the documentary contrast, measured.** Same two datasets, `permute(shuffle) permreps(20)`:
- arm A: `rc = 0`.
- arm D: **`rc = 459`**, "Error: permute() requires the cluster() variable to nest i (every i
  wholly inside one cluster) ... (7 in-sample observation(s) do)" — 7, i.e. every row of
  the straddling unit, not just the phantom row.
So ONE phantom row makes `permute()` refuse the run with a named error (and `sthlp:431` documents
that its design counts markout-dropped rows), while the same row makes `bootstrap(standard)`
silently switch to the variance branch the ticket says collapses `att_i`'s composition variance,
with no message anywhere. Same gap, same engine family, opposite handling — one documented and
loud, one undocumented and silent.
(First Q3 attempt used `permreps(19)`, below wooldid's floor of 20, and died r(198) in BOTH arms —
my error, corrected in `l1_q3_permute.do`.)
