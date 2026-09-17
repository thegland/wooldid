# Open defects ledger

Running list of defects found while working the ruled slates, held for a fix pass at the
end rather than fixed inline. Opened 2026-09-10 at the author's instruction ("keep the
defects found in a separate to do file which we will address at the end").

Rules for this file. A defect goes in OPEN when it is a genuine bug or a real deviation
that (a) is not what the item under work was about, or (b) would move a reported number
and so wants its own verified commit. Anything that is a judgment call for the author
goes to `held_for_ruling_2026-09-09.md` instead, not here. Anything fixed gets moved to
CLOSED with its commit, so this file is also the record of what the sweeps bought.

---

## RULINGS ON THE DOCUMENTATION REVIEW'S ITEMS (author, 2026-09-11)

All fifteen were put to the author as multiple choice and all fifteen were ruled. Build
order below is by risk, not by item number: the harness first (so everything after it is
judged by the stricter gates), then the changes that move no number, then the three that do.

| id | ruling |
|---|---|
| F033 cells | **Apply the add-one.** `2*(min+1)/(n+1)` in the studentized cell branch, matching `wooldid_boottp` and `wooldid_bootpwt`. One convention for every recentred rank-based p. Moves reported cell p-values; needs a re-bless. |
| D-1 | **Carve-out like `fastnose`.** Keep the 6-column cfx shape under `permute()` even with `cfxplotnose`, so the cells get real permutation inference. |
| D-2 | **Fix the wording.** "no permutation inference", naming which columns survive; the analytic SE stays, as on every other permutation surface. |
| D-8 | **Keep the refusal, say why in the caution.** The bootstrap critical value is a maximum over the family taken replication by replication; the caution must name the horizon that caused the refusal. |
| D-7 | **Caution naming the families** that got no band, replacing the all-or-nothing one. |
| D-10 | **Caution on both paths** in `wooldid_permsaveasg` (no assignment block written / N cluster columns omitted). |
| D-6 | **Split the caution** so `wdptsmiss` is spent: too few completed replications versus no usable observed analytic SE. |
| D-3 | **Rename to a `wooldid_` name** at all four `mata: M = ...` sites. |
| D-5 | **Delete all four inert sites.** All behaviour-neutral; the explanatory comments go with them. |
| D-4 | **Honour `dots()`** on the physical-expansion leg, so the option means one thing everywhere. |
| D-9 | **Reseed per family** in the band loop, so a family's band depends only on that family. Moves band endpoints; needs a re-bless. |
| D-11 | **Widen the regex and add both savers** to `mkrun.py`'s file pin. |
| D-12 | **Remove `compare.py`'s absolute floor at all rtols**, matching `compare_oldnew.py`. |
| D-13 | **Make UNRECOGNIZED fatal and add the `fw` bucket** in `coverage_check.py`. |
| D-14 | **Tighten `check_validate_sync.py`** to a real assignment. |

---

## OPEN — to address at the end

Nothing.

The fifteen items the documentation review of 2026-09-10 turned up were ruled on 2026-09-11
and all fifteen are built - see CLOSED below.

---

## HELD FOR THE AUTHOR (not defects — judgment calls)

Nothing. This list was stale: all five items it carried were ruled on 2026-09-10 and are
recorded as built or confirmed in the table at the head of
`held_for_ruling_2026-09-09.md` — PS-2 (name the exporters only) built, PS-62's sentinel
diagonal confirmed as built, BS-15's cut confirmed at 1 - 1e-2, PS-45's band built as a
percentile-t under boottstat, PS-46 built. The F033 cell p-value entry added on 2026-09-10
was ruled on 2026-09-11 and built the same day (see CLOSED).

---

## CLOSED

**RETRACTED the same day: `si8demand`'s under-prediction is the documented design, not a defect
(opened and closed 2026-09-16).** I opened this as a bug after watching a plain 200-unit x
20-period panel measure `cells=2004 ilev=200 tlev=20 demand=2424` against `c(maxvar) = 5000`,
stay silent, and then die on Stata's own generic `maxvar too small`. The measurement is correct
and reproduced on a build carrying one added `dis`. **The conclusion was wrong.** The block's own
comment states the intent outright (`wooldid.ado:1342-1345`):

> *"The preflight uses the EXACT treated interaction-cell count (distinct (_i,_t) pairs, x2 under
> contreat for the dose slopes), **a lower bound on the regression's variable demand**, so an
> early refusal can only fire when the late r(907) was certain anyway."*

`si8demand` is deliberately a LOWER bound. The design accepts false negatives - runs that will hit
r(907) and are not pre-empted - in order to guarantee zero false positives, because a false
positive would refuse a specification that would actually have estimated fine. Given that the
guard's whole value is a better message for a failure that happens either way (a case-1 guard,
see `guard_cost_2026-09-16.md`), trading recall for precision is the right way round: a spurious
refusal is worse than a generic error message.

So the observed miss is the design working as specified, not failing. **What the episode actually
shows is how wide the gap can be** - the true requirement exceeded the lower bound by more than
2x - which is worth knowing when judging what the guard is worth, and is recorded in the
guard-cost report for that purpose. No code change is wanted.

I should have read the eleven lines of comment above the block before filing. The same comment
also explains the `safetyinput` gate (P6, author ruling 2026-09-07) that this session's rep-gate
work extended.

**A10's caution now tested on both branches - opened and closed 2026-09-16.** Closed by three
message-assertion legs in `test_bootstrap_gates.do` (m7/m8/m9), not the two the opening entry
proposed. The gate is a conjunction - `wdb_reid' non-empty AND `wdb_nest' == 1 AND a covariate
naming `i' - so two legs would have left a vacuous pass available: a gate that never fires at all
also satisfies "fires here, silent there" if the silent leg is the only negative. The three legs
pin each conjunct separately, on the fixture's `s = ceil(unit/4)` (s nests unit):

| leg | regime | covariate | expected |
|---|---|---|---|
| m7 | `i=unit`, `cluster(s)` - nesting | `fe(unit)` | caution MUST appear |
| m8 | `i=s`, `cluster(unit)` - finer | `fe(s)` | MUST NOT appear (the T2/T3 suppression) |
| m9 | `i=unit`, `cluster(s)` - nesting | none | MUST NOT appear (the `wdb_ikeyed' conjunct) |

All three run at rc 0 and the file prints `BOOTGATES_ALL_PASS`.

**And the legs immediately found a defect in the harness they were written for - the more
important half of this entry.** `logh` asserts a needle appears in a log it opens, runs a command
into, then closes and reads. Stata echoes every command into EVERY open log, so the
`logh <marker> <needle>' line itself was written into the log a moment before `logh' closed it -
carrying the needle. **The needle always matched its own echo. All six pre-existing `logh`
assertions had been passing vacuously and could not have failed.** It surfaced only because the
new negative helper `lognh` has the opposite polarity, so the same echo made m8 and m9 fail
where it had been silently rescuing everything else.

Fixed by skipping SMCL command-echo lines (`{com}. ` and continuations `{com}> `); genuine output
is unprefixed or carries `{txt}`/`{res}`/`{err}`. A second latent bug fell out of the same fix:
`gettoken` leaves the separating blank on the remainder, so the needle carried a leading space and
could never have matched a message starting at column 0 - invisible while the echo was matching,
fatal the moment it stopped. Both helpers now `strtrim()` the needle.

Re-run after the fix: all six original assertions still pass, so their CLAIMS were true - only the
test was. That is this slate's fifth instance of one shape, and the first on the negative side:
a check credited with catching a failure it structurally could not catch. The general lesson stands
unchanged - calibrate an instrument on a case whose answer you already know, in BOTH directions.

Eleven other files under `test/` both open a log and read a file. All eleven were then checked and
**none shares this shape**, for two independent reasons: (1) structurally, each closes its log
BEFORE calling its assertion helper, so the helper's own invocation line cannot reach the log it
reads - `logh` was unique in closing the log inside the helper; and (2) empirically, no assertion
needle in any of the eleven occurs in that file's own command text, which is the only other way an
echo could satisfy a needle. (A literal scan flags `if flag == 1` in `test_edges.do`,
`joint_horiz` in `test_permute3.do` and two log filenames; all are commands, comments or e()
names, none an assertion needle.) The residual exposure is therefore zero rather than unmeasured.

**RNG bleed in `test_bootstrap_gates.do` — opened and closed 2026-09-16 (`REDACTED_SHA`).**
Cases without their own `rseed()` draw from whatever RNG state earlier cases left, so one can go
red after an edit that changes how many draws an EARLIER case consumes, reading as a regression in
whatever was just touched. Found because a strengthened G9 asserting 50 of 50 completed
replications passed in isolation and failed at 49/50 inside the file.

**The opening report overstated it and the corrected scope is the useful part.** "41 of 99" counted
`bootstrap(` occurrences, comments included, not calls. There are 29 unseeded CALLS, and they are
not one population:

- **19 assert a nonzero rc** - parse-time refusals (`bootstrap(wild)`, a missing `reps()`,
  `bootalpha(nonsense)`, `dots(-2)`). They error before any resampling and never touch the RNG, so
  seeding them is churn that changes nothing.
- **2 must stay unseeded**: G31 at lines 745 and 750 runs the same command twice with no seed,
  captures `e(boot_rngstate)` and asserts `mreldif(G31D1, G31D2) != 0`. The shared ambient state IS
  the test. A seed there would not tighten the case, it would delete what it checks.
- **8 genuinely run and assert an rc or a log message**, and those were the real exposure: the two
  `bayes` + `att_it` runs under `safety(off)`, the G7/G7b nesting controls, G10, G15b, and the three
  `logh` message assertions, whose printed text can shift with the draws.

Those 8 now carry distinct seeds (70101-70108). Verified by running the leg end to end:
`BOOTGATES_ALL_PASS`, no failures. Necessary rather than ceremonial, since changing a seed changes
the draws and a case pinned on a caution's wording can legitimately flip - which is precisely the
G9b failure mode.


### The help-file audit — 55 gaps written, 1 ado defect fixed, 2026-09-11
Second half of the documentation instruction. 14 agents worked FROM THE ADO (options off the
syntax line, `e()` names harvested from every `ereturn` site, behaviours read out of the
engines) and checked each against the help: 579 items, 59 claimed gaps, 55 confirmed after
adversarial re-verification, 3 downgraded, 1 refuted.

- **`14c30d3`** — the 55 confirmed gaps written up across the option entries, the inference
  sections and Stored Results.
- **`69641f0`** — the same prose rewritten in the file's own style. The gaps were documented
  correctly but in the register three days of code-comment writing had left behind; a
  coverage audit checks coverage, not voice. Measured against `2a2da17` (the sthlp stood
  unchanged 2024-07-31 to 2026-07, so that commit IS the "12 months ago" reference): all-caps
  emphasis at 9 sites where 2024 uses caps only for acronyms; ` - ` dash parentheticals at
  13.6 per 1000 words against 2024's 0.26; four authorial asides; and one
  `{help x:{cmd:y}}` nesting neither version of the file had ever used. All four now zero in
  the new text, the pre-existing campaign-era text untouched. Also in that commit: the
  stray-companion screen at `wooldid.ado:5965` refuses nine boot/perm companions with an
  error and only three said so, so `reps()`, `rseed()`, `boottstat`, `bootdf()`,
  `bootdisplay()` and `bootsave()` gained the clause.
- **the ado commit below** — audit item #57, the one finding that was not a help change.

**#57, the permutation off-centre caution (`wooldid.ado:7950`).** It cited, as its worked
example, "the ame rows under lattice() with contreatcontrols(basic), where the baseline dose
slope is not partialled out (see the standing issue on the reference cell)". SI-20 FIXED that
condition on 2026-09-05 (`1df8613`: `basic` now always absorbs `i.lattice##c.(ringer)`; the
measured reference margin went from 0.27..2.00 to exactly 0 on all 59 affected configs) and
the ticket moved to `wooldid_solved_issues.md`. So the caution was sending a user after a
cause that no longer exists and a section that no longer stands. The mechanism is kept, the
stale example and the dangling pointer are gone, and a comment above the line records why
and warns against re-attaching an unmeasured example. Straight text defect; no number moves.
It had been named only in `14c30d3`'s commit message, which made it a standing item nothing
in this file recorded — that is the reason it is written up here rather than only fixed.


### The documentation review's fifteen items — ALL BUILT 2026-09-11
Ruled by the author on 2026-09-11 (multiple choice, all fifteen answered; the ruling table is
above). Built in three risk-ordered commits, each gated.

**`ed95674` — the harness (D-11..D-14), first, so everything after it is judged more strictly.**
- D-12 `test/compare.py`: the 1e-12 absolute floor is gone at every tolerance, matching
  `compare_oldnew.py`. The selftest now pins the property with the two cases that used to
  pass: `close(0.0, 1e-300)` and `close(1e-13, 5e-13)`.
- D-13 `coverage_check.py`: an UNRECOGNIZED option token FAILS the run; `fweight` has its own
  bucket.
- D-14 `check_validate_sync.py`: the assignment regex requires the name to BE the thing
  assigned. Tightening it immediately surfaced `item`, which the worker sets fourteen times
  as a `foreach` loop variable and the old pattern had been suppressing for the wrong reason;
  `foreach`/`forvalues` are recognized now, so the script is OK on its own terms rather than
  by accident.
- D-11 `mkrun.py`: the file pin covers every file-writing option at every legal abbreviation,
  and `dumpres` takes a list of prefixes. `save(...)` used to slip it silently; `permsave()`
  and `bootsave()` were outside it altogether.

**`9bb24c6` — the seven that move no number (D-2..D-8, D-10).**
- D-7/D-8: `wooldid_bootband` reports WHY it refused, naming the horizon; `wooldid_esband`
  collects the declined families and names them. The refuse-the-whole-family policy stays, as
  ruled, and the caution now carries its reason.
- D-2: the drift caution says what is true — no permutation inference, p and interval blank,
  the se column analytic, the t blank unless `permtstat`. The replay line matches.
- D-6: the permutation blanking caution is split, so `wdptsmiss` is finally read: too few
  completed replications (raise `permreps()`) versus no usable observed analytic SE (a
  property of the specification).
- D-10: both silent exits in `wooldid_permsaveasg` now report.
- D-3: the four one-liner Mata temporaries are `wooldid_Mtmp`.
- D-4: `dots()` is honoured on the physical-expansion leg, using the documented cadence.
- D-5: the four inert sites are deleted, each leaving a note saying what went and why it
  could not fire.

**`82934b3` — the three that move numbers, each with before/after evidence.**
- **F033 on the studentized CELL p-value.** Now `2*(min+1)/(n+1)`, one convention for every
  recentred rank-based p. The effect was larger than the estimate in the finding: on the cert
  fixture at `reps(40)`, the OLD build produced **22 of 50 cell p-values at exactly 0** with a
  minimum of 0; the new build produces none at zero and none below `2/(completed+1)`.
- **D-1, `permute()` + `cfxplotnose`.** `cfxplotnose` now takes the same `permute()` carve-out
  `fastnose` has. Old build: the first cfx matrix had `colsof = 1`, no `e(permdist_*)`, and a
  missing cell p. New: `colsof = 6`, `permdist_` posted, a real permutation p.
- **D-9, the band's seed.** Derived per family from the family's NAME rather than set once
  before the loop. Same fixture and spec, `e(esband_crit_att_it)` with `att` banded before it
  versus `att_it` drawing first: before `2.474072413106145` vs `2.492870068598342`; after
  `2.47972450315949` in both.

**Pinned.** None of the 504 dual-build configs reaches these three paths, so the changes would
have been invisible to the suite. Three cert legs were added and each was verified to FAIL on
the pre-change build: `test_esband.do` B10 (the band's critical value is keyed to the family,
not the draw order), `test_permute3.do` Q33 (`cfxplotnose` no longer truncates under
`permute()`), `test_bootstrap_tstat.do` T12 (no studentized cell p is zero or below the
add-one floor).


### PS-3 — built 2026-09-10, after a wrong diagnosis
**Ruling.** *"Per-column formats. Keep the table structure, give each column an explicit
format."* **Author ruling 2026-09-10 on the four options put to him: the conventional
formats, as ruled** - accepting that a p of 1.84e-21 prints as 0.000, which is what every
comparator does. The magnitude is not lost from the data, only from the printed cell;
`e(wdidmainresults)` still carries the full value.

**Built.** `wooldid_mlfmt` builds matlist's `cspec()` and `rspec()` from each matrix's own
names, so no call site hard-codes a column count: statistics two decimals, p-values three,
`relyr` an integer, everything else a general format so no estimate loses precision. The
p-value rule is EXACT-match, because the summary-statistics table has `p25` and `p75`
columns that are percentiles. The rules are placed where the default put them - none above
the header, one under it, one at each equation-group boundary, none at the foot - so the
structure is unchanged and only the number formatting moves.

**The earlier failure was mine, not matlist's.** The first attempt was abandoned on a probe
that passed a `nooutput` option matlist does not have, so every row-spec length errored
regardless, and this file recorded "matlist's contract fights these tables" as a finding. It
does not. `rspec()` takes one character per horizontal rule and the count is rows + 2,
independent of equation groups - measured on all three tables. **A single failing probe is
not a property of the tool.**

### The eleven audit items closed 2026-09-10
From `review_process/research_2026-09-10/sweep_audit.md`. PS-3 is the twelfth and is above.

- **PS-20** the Ibragimov-Muller banner: the header now names IM inference and its reference
  distribution instead of the analytic cluster specification the run did not use.
- **PS-64** the joint pre-trend caution now fires under `pretrendpower()` as well as
  `jointtests`, and names whichever entry point asked for it - "one statistic, one rule, in
  both entry points", as the ruling says.
- **CS-98** resolved by NAMING rather than by moving: the two screens that run before the
  `contreatwithin` resolution now read `safetyuser`, snapshotted once at the top of the
  validator, so no reader has to know a line number to know which value a gate means. Hoisting
  the resolution instead would have extended the relaxation into a screen it has no
  relationship with, and those two screens terminate the run, so they were never in CS-98's
  hazard class - a message overturned by a later rewrite. No behaviour moves.
- **CS-60** the subgroup joint p is a scalar at all THREE posting sites (the audit found two;
  the permutation engine had the same pair as the bootstrap one).
- **CS-76/CS-115** the bootdisplay() entry no longer claims every family is computed on every
  run - a claim it contradicted five sentences later.
- **CS-118** `bootdisplay(bca)` added to the `e(boot_stat)` value list; `e(boot_display)`'s
  legal values enumerated.
- **CS-120** `e(wdid_fcncell)` and `e(wdid_fcskip)` documented, and the Replay paragraph no
  longer disowns the two cell-collapse names the `fastcollapse` entry tells users to read.
- **P8(ii)** the Standard Errors synopsis split into labelled sub-runs (analytic /
  shared-by-both-engines / permutation-only / bootstrap-only); `ties` states that it requires
  `bootstrap()`; `permside` and `permrawp` carry their abbreviation markers; `regtol` shows
  its argument.
- **P7 accident 3** both inline compatibility lists reshaped: reasons on every item, and the
  two refusals `permute()`'s list omitted (multi-way `cluster()`, `customvce()`) added.
- **CS-25/40** the control-group refusal added to the inventory's class rows, and the "gates
  fifty-eight sites" claim replaced with what the section actually delivers.
- **SI-A1** (author ruling 2026-09-10, a note and only a note): the `cluster()` entry states
  that clustering below the level of `i` is the intended practice, what goes wrong at or above
  it (pointwise standard errors understated to roughly 0.66-0.75, joint tests oversized), and
  why there is no plug-in correction.


### C-1..C-4. The four defects the batch-5 shared assessment surfaced (8280a2e, 2026-09-10)
1. `cleanmatrices` runs `clear mata` inside the worker, dropping wooldid's own Mata
   kernels; PS-62 had put a Mata call on the mainline analytic path, so a plain
   `att cleanmatrices` fit died `r(3499)`. Guarded. **No cert suite caught this** — only
   the 501-config dual-build comparison did. `test_stored` S9 now pins it.
2. `wooldid_bootwire` relabelled `e(wdid_vcvstate)` as `analytic-diagonal` on the
   repost-failure fallback, even when the matrix left standing was the spliced one.
3. The PS-62 splice ran inside every bootstrap/permute replication, where nothing reads
   `e(V)`'s off-diagonals. Gated on the replication marker.
4. `e(wdid_vblocks)`/`e(wdid_vblockn)` describe the hidden `e(V_analytic)` under
   `bootstrap()`, not the reported `e(V)`. The help now says so.
