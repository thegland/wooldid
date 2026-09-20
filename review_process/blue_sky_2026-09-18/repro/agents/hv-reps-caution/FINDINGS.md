# hv-reps-caution — FINDINGS

## HV-REPS-01  (VERDICT: CLAIM CONFIRMED, and it is a regression; consequence narrower than stated)

**Title.** The reps() p-value-resolution Caution is unreachable without `cluster()` — a
regression introduced 2026-09-13 by commit 813c414 (RF-2), and a doc-vs-code divergence.

**Claim.** The three-legged caution block at ado:6473-6478 in `wooldid_bootctl` is nested
inside `if "`clus'" != "" {` (opens ado:6459, closes ado:6501). No `cluster()` -> no caution,
whatever `reps()` is. `wooldid.sthlp:535` and `:634-636` promise it unconditionally.

**Static evidence (three independent routes).**
1. Printed ado:6455-6502 with line numbers. The only `}` between the abort block's close at
   6468 and 6501 are 6478 (the caution block's own), 6499 and 6500.
2. Mechanical brace/condition sweep (`sweep_cautions.py`, output `sweep_all.txt`) that
   strips comments and quoted-string contents and does NOT trust indentation:
   `ado:6475/6476/6477  under 6459:if "`clus'" != "" { , under 6473:if `reps' > 0 {`.
3. Call site ado:4941 passes `clus(`cluster')` verbatim from the user option; no default
   substitution, so the local is empty whenever `cluster()` is omitted.

**Provenance — this is a REGRESSION (brief category 2).**
`git log -S"Moved below every aborting" -- wooldidinstall/wooldid.ado` -> 813c414
(2026-09-13, "Fix three disclosure defects the assessment found"). The diff moves the block
from top level (parent-file line 6014, between the CS-20 comment and `* companion-value
validation`; parent saved at ./prev_wooldid.ado) into the cluster branch. The commit message
says "Verified both ways: ... a healthy run still prints it exactly once" — that healthy run
must have carried `cluster()`.
The block exists because of CS-20 (author ruling 2026-09-07 "caution, no floor"), whose own
in-code comment (parent:5996-6013) states the problem it was shipped to fix: "reps() carried
neither, so reps(20) printed percentile p-values on a .10 grid in silence", and notes that
the conditional footer `wooldid_pfloor` "prints only when a displayed p is exactly 0 ... so a
low-rep run whose p-values all happen to be non-zero was told nothing at all."

**Empirical evidence.** `repro_reps_caution.do` / `.log`, 100 units x 10 years, three treated
cohorts, 40 never-treated, `reps(10)`. Count of the fragment `Caution: with reps`:

| leg | command | caution |
|---|---|---|
| 1 | `bootstrap(standard) reps(10)` | **0** |
| 2 | ... `cluster(unit)` | 1 |
| 3 | `bootstrap(standard) reps(10) boottstat` | **0** |
| 4 | ... `cluster(unit)` | 1 |
| 5 | `bootstrap(standard) reps(10) bootdisplay(wald) jointtests` | **0** |
| 6 | ... `cluster(unit)` | 1 |
| 7 | `bootstrap(bayes) reps(10)` | **0** |
| 8 | `bootstrap(standard) reps(10)` on an outcome with a real effect | **0** |

3 for 3 with `cluster()`; 0 for 5 without. All nine legs returned rc 0.

**CONSEQUENCE SURFACE — narrower than the seeded claim, and one leg is worse.**
The seeded claim says "an unclustered low-reps bootstrap gets no warning at all." That is
true of the *Caution*, but only fully true on one of the three display families, because a
ROUTINE (unconditional) footer line in `wooldid_display` at ado:11784-11785 also names the
grid:
- **percentile/bc/bca (the default), unclustered:** footer prints "The rank-based coefficient
  p-values have resolution 2/10 = 0.200: a p printed as 0 means p below that ...". The bare
  number is disclosed. What is lost is the judgment the Caution adds: that 0.200 is **above
  .05**, that therefore *no* coefficient p-value can reject at 5% whatever the data say, and
  the "use at least 40 replications" guidance. PARTIAL mitigation.
- **boottstat, unclustered:** footer prints "...add-one ranks over 10 completed replications,
  so their resolution and their smallest attainable value are both 2/(n + 1) = 0.182".
  PARTIAL mitigation, same loss.
- **`bootdisplay(wald)` without `boottstat`, unclustered: NOTHING AT ALL.** ado:11784-11785 is
  an `if ... else if "`wdbdsp'" != "wald"` chain with **no wald branch**, and
  `wooldid_pfloor` (ado:10783) exits immediately under `bootdisplay(wald)` without boottstat.
  Leg 5's complete output is in the log: the whole "Bootstrap Inference" footer contains no
  occurrence of the words "resolution", "grid", "add-one" or any floor. So
  `bootstrap(standard) reps(10) bootdisplay(wald) jointtests` with no `cluster()` is the
  configuration in which the 1/(reps+1) = 0.091 floor on the joint / equality / adversarial
  p-values is stated **nowhere** — not by the caution, not by the footer, not by pfloor.

**Who this hurts / how realistic.** Realistic. An unclustered bootstrap is a documented,
first-class leg ("clusters ... or observations when no cluster() is given", sthlp:532), and
`reps()` has no minimum precisely so that small values are used to time a specification — the
scenario CS-20 was written for. A user timing at `reps(10)`, reading a printed p of .200 or
.400 in the default display, is told the grid but not that the grid has made rejection
impossible; under `bootdisplay(wald)` they are told nothing.

**Expected.** `sthlp:535`: "A rank-based p-value cannot be finer than its own grid, and
{cmd:wooldid} prints a Caution naming that grid whenever it is coarser than .05". `sthlp:634-636`:
"{cmd:reps()} and {cmd:permreps()} both print a caution when the floor they imply sits above
.05, naming the number." Neither sentence is conditioned on `cluster()`.

**Fix location.** RF-2's stated intent was "below every aborting validation, above the other
bootstrap cautions". On the unclustered leg there are no aborting validations inside the
cluster block at all, so the faithful landing spot is immediately after ado:6501 (the close of
the cluster branch), at program top level. The reason this survived review is that ado:6473-6478
is indented two spaces — the same indent as the enclosing `if` at 6459 — so it reads as
top-level.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-reps-caution/repro_reps_caution.do`
severity medium | confidence certain | realistic true | surface bootstrap()/reps()/bootdisplay(wald)

---

## Permute side: CLEAN (negative result)
`wooldid_permctl`'s counterpart caution at ado:5899-5912 is at program top level (mechanical
sweep: no enclosing condition but its own `if `permreps' > 0 & "`permexact'" == ""`). It could
not be stranded by a cluster conditional in any case, because `permute()` REFUSES a run
without a single-variable `cluster()` (ado:5946-5949), and `permreps()` has a hard floor of 20
(ado:5885-5892) that bounds the hazard the caution describes.

## Sweep of every other Caution in both engines: CLEAN (negative result)
`sweep_cautions.py` emits, for all 143 `Caution`/`Warning`/`Note` display lines in the ado, the
chain of enclosing block headers computed from brace depth after stripping comments and quoted
strings. Judged on subject match, the only mismatch is HV-REPS-01. In `wooldid_bootctl` /
`wooldid_permctl` specifically:
- ado:6442 bayes cell-weighted caution — under `if "`bootstrap'" == "bayes" & "`cellwtd'" != ""`. Correct.
- ado:6479 fewer-than-10-clusters — under `if "`clus'" != ""`. Correct (subject IS clusters).
- ado:6498 fewer-than-6-treated-clusters — under clus + ttrevar + `r(N) < 6`. Correct.
- ado:6056 / 6064 / 6070 permute cautions — top level. Correct.
- ado:6119 / 6122 subgroup cautions — under `if "`sg'" != ""`. Correct.
A false positive I chased and cleared: ado:4077 `if `evpostrc' == 0` nested inside ado:4053
`if `evpostrc' != 0` is NOT dead code — `evpostrc` is reassigned at ado:4076.

---

## Boundary control (negative result — the caution is CORRECT when it is reachable)
`probe2_boundary.do` / `.log`, same panel:

| leg | command | `Caution: with reps` |
|---|---|---|
| B1 | `reps(39) cluster(unit)` | 1  (2/39 = .0513 > .05 — fires, correct) |
| B2 | `reps(40) cluster(unit)` | 0  (2/40 = .05, not > .05 — silent, correct) |
| B3 | `reps(39)` no cluster | **0** — the defect at the exact setting sthlp:538 warns about |

`repro_reps_caution.do` leg 9 (`reps(400)`, unclustered, real effect) is the other control:
no caution (grid .005), and `wooldid_pfloor` correctly prints "p = 0 above means p < 0.005".

---

## HV-REPS-02  (companion; LOW)

**Title.** The permute engine still prints its resolution caution for a run that then aborts —
the exact "premature caution" disclosure defect that 813c414 fixed on the bootstrap side.

**Claim.** In `wooldid_permctl` the permreps() resolution caution (ado:5899-5912) executes
BEFORE the `wordcount("`clus'") != 1` refusal at ado:5946-5949. So `permute(shuffle)
permreps(20)` with no `cluster()` advises the user about the p-value resolution of a run that
never executes, then exits r(198).

**Observed** (probe2_boundary.log, leg `WDLEG_B4_PERMUTE_NOCLUS`, rc 198):
```
Caution: with permreps(20) the smallest attainable per-coefficient permutation p-value is
0.095 (2/(permreps + 1)), above .05, ...
Error: permute() requires cluster() with a single variable - the unit whose treatment-onset
assignment is randomized.
```
**Expected.** Commit 813c414's own statement of the rule: "a user whose run never executed was
advised about its p-value resolution. Moved below every aborting validation." The permute
counterpart was not moved. Note the irony: the bootstrap move is what created HV-REPS-01, and
the identical ordering defect it was fixing survives in the other engine.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-reps-caution/probe2_boundary.do`
severity low | confidence certain | realistic true (forgetting `cluster()` on a first permute
run is an ordinary user slip) | surface permute()/permreps()

---

## Also checked, CLEAN
- The other two fixes in commit 813c414 landed correctly on the pinned build: only one
  `error _rc` survives (ado:1471, the deliberate `qui cap` ppmlhdfe site the commit message
  names), and the `bootsave(name, replace)` comma refusal is live at ado:6282-6289.
- sthlp:532's "or two observations on the unclustered leg" refusal is NOT stranded even
  though ado:6465's copy is inside the cluster branch: `wooldid_bootrun` sets
  `local clusv `clus'` and falls back to an observation id when it is empty (ado:8918-8921),
  so the `wdb_C < 2` refusal at ado:8956 covers the unclustered leg. The only divergence is
  that it fires at run time rather than "when the options are screened" as the help says —
  cosmetic, not filed.

---

## HARDENING OF THE "worst leg" PARAGRAPH (probe3_waldjoint.do / .log)
The first `bootdisplay(wald)` leg (repro leg 5) had no event-study window, so its joint tests
were absent and the silence was partly "nothing to disclose about". Rerun WITH a window
(`esprelength(2) espostlength(2)`), where a joint p really does print:

| leg | command | joint p printed | `Caution: with reps` |
|---|---|---|---|
| C1 | `reps(18) bootdisplay(wald) jointtests espre/post(2)` no cluster | **0.6842** = 13/19 | **0** |
| C2 | ... `cluster(unit)` | 0.4211 = 8/19 | 1 |
| C3 | `reps(10)` ... no cluster | **0.8182** = 9/11 | **0** |

18 is the largest `reps()` at which the one-sided add-one floor 1/(n+1) = .0526 still exceeds
.05. On C1 the printed joint pre-trend p is a one-sided add-one share that cannot go below
.0526 whatever the data say, and the run's complete output contains no statement of that floor
anywhere: the caution is behind the cluster branch, ado:11784-11785's footer chain has no
`wald` branch (it is `if boottstat ... else if "`wdbdsp'" != "wald" ...`), and
`wooldid_pfloor` exits at ado:10783 under `bootdisplay(wald)` without `boottstat`. The
footer's only remark about the family is descriptive ("p = add-one share of replications at or
beyond the observed statistic") with no number. C2 is the positive control: identical run,
`cluster(unit)` added, caution present.

## Two observations (not filed as findings)
- **Fix-location note corrected.** RF-2 asked for "below every aborting validation, above the
  other bootstrap cautions". Once the block leaves the cluster branch those two cannot both
  hold: either before ado:6459 or after ado:6501 is correct, and the author picks. The
  substantive point is that on the unclustered leg there is no aborting validation in
  `wooldid_bootctl` for the block to sit below.
- **Off-by-one between the two engines' caution texts, for the identical formula.** The
  permute caution (ado:5911, seen in probe2 leg B4) says "Use at least 40 replications for a 5%
  test" for a 2/(permreps+1) floor; the boottstat caution (ado:6475) says "use at least 39" for
  the same 2/(reps+1) formula; sthlp:470 says "below 39 a 5% rejection is unreachable". At 39,
  2/40 = .05 exactly. Conservative in both directions, not wrong, but inconsistent.
- The caution is keyed on `reps()`, not on what the table ends up displaying: in repro leg 6
  (clustered, wald, no ES window) it fired about joint p-values the run did not print. That is
  its design, not a defect.
