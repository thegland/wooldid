# Documentation review of the code files — 2026-09-10

**Ask.** *"Review the codebase and check for sections that are low on explanatory/descriptive
comments... robustly commented enough to ensure that users freshly encountering the project
for the first time always had excellent context about the code. Ideally, the comments would
also caution about common or particularly tricky to spot pitfalls and subtleties of the code
(or of the behavior of called upon features) that normally might complicate either coding or
results assessment."*

**Method.** Three measurement passes (ado programs, Mata functions, Python harness) to find
the thin regions, then a 12-agent read-only *find* pass — one agent per region, each asked
two questions: where would a competent Stata programmer meeting this code for the first time
be unable to tell what a block does or why; and where does a trap sit uncommented. Agents
returned structured findings and wrote nothing. The comments themselves were written by
hand, in one voice, from those findings plus the traps this campaign has measured. 127
findings came back; all were acted on, folded into a neighbouring comment, or recorded in
`open_defects.md`.

## What the measurement said

| region | code lines | comment density |
|---|---|---|
| ado programs, all | 10,812 | 26.4% |
| Mata block, all | 2,476 | 24.1% |
| `wooldid_permjnt` | 320 | 13.0% |
| `wooldid_esplots` | 2,810 | 18.6% |
| `wooldid_display` | 408 | 19.0% |
| `wooldid_permwire` | 253 | 19.7% |
| `wooldid_bootjnt` | 529 | 20.3% |
| `wooldid_sgchain`, `wooldid_pfloor` | 25, 26 | 0% |
| Mata `bootcellkey`, `bootwrowx`, `bootadvr` | 15, 20, 53 | 0% |
| `test/suite/coverage_check.py` | 190 | 0.5% |
| `test/compare.py` | 122 | 0.8% |
| `test/suite/mkrun.py` | 215 | 2.3% |
| `test/suite/compare_oldnew.py` | 386 | 3.0% |

Density turned out to be the wrong ranking and was not used as one. `esplots` is the largest
thin region but a misreading there produces an ugly figure; `permjnt` and `bootjnt` compute
p-values from quadratic forms under three metric conventions, and a misreading there ships a
wrong test. The work was ordered by consequence. The 0%-density programs were mostly not
undocumented at all — their headers had been *separated from them* (below).

## Commits

| commit | scope |
|---|---|
| `45d1a0a` | the inference engines: permjnt/permcellsx, bootjnt, permrun/permwire, bootrun/bootwire, the engine screens and savers |
| `fdc35bf` | the Mata kernels: permutation conventions, bootstrap moments, the bands, the adversarial statistic |
| `3e0ea38` | presentation, plots, screens; five headers reunited with their programs |
| `a16b6b1` | the Python test harness |

All four are comment-only apart from three text corrections, each a case of documentation
contradicting the code: `wooldid_esplots`' replay gate said it required contract `'1'` while
testing for `'2'`; the help said `histogramcohorteffects` returns one row per cohort-by-period
cell under `att_it`/`att_itime` when both margins calls are `over(i)`; and the F038 block's
last sentence had been truncated mid-clause.

## The structural finding

Five documentation blocks had been separated from the programs they describe by the
compiled-size extractions out of `wooldid_worker`. The effect was worse than a missing
comment: several programs looked undocumented at their definitions while others carried two
paragraphs about their neighbours, and the blocks' own "below" references pointed across an
intervening program.

- the F6 permutation-engine block sat above `wooldid_bootctl`, which it never described →
  moved to `wooldid_permrun`, which had no header at all
- the F038 p-value-floor block was pushed off `wooldid_pfloor` by the PS-3 insertion of
  2026-09-10 and left above `wooldid_mlfmt` → restored, completed, and its three rulings named
- the NFADV adversarial block and the BS-8/SI-7 subgroup-labelling block were stacked above
  `wooldid_ccscreen` → each moved to its own program
- the E1 constraint-chain block sat above `wooldid_fcol` → moved to `wooldid_sgchain`, the
  file's only zero-comment program

Two headers had also drifted from their code and were rewritten rather than appended to: the
adversarial block still documented the free-intercept `alt==0` statistic with
`q = rank(V) - rank(M)` that BS-12 replaced on 2026-09-08, and `wooldid_bootctl`'s header plus
its "C-stage refusals" block described two refusals that no longer exist and recorded the
`pretrendpower` refusal as both permanent and lifted.

## The pitfalls now recorded in the code

Written at the site where each bites, not in a document:

- `ereturn matrix name = mat` and `ereturn repost V=` MOVE their source (measured: rc 0
  before, r(111) after) — at the three sites that manufacture it
- `ereturn repost V=` does not re-run `ereturn post`'s screen, so `e(V)` leaves the bootstrap
  wiring deliberately indefinite with a `-1e100` sentinel diagonal
- `_rc` is reset by the very next command — and, separately, plain `local` assignments do
  *not* touch it, which is why one `cap confirm` ... `if _rc` pair separated by two locals is
  correct and reads like a bug
- `st_matrix(name, J(0,0,.))` deletes the Stata matrix rather than creating an empty one
- Mata passes by address, so `wooldid_advW`'s symmetrization mutates its caller's `V`; it is
  safe only because `(V+V')/2` is bitwise idempotent
- a bare `mata: M = ...` leaves a global `M` in the user's session
- `matlist`'s `rspec()` is one character per rule, count rows + 2, and a mispaired format
  string aborts a *completed* estimation
- `twoway` joins points in dataset order, not x order
- the positional alignment contract shared by four e(b) walks, and which of them is guarded
- the four-way p-value convention integer, and the two places the expression is duplicated
- `do <path>` is the one load style under which wooldid cannot recover its Mata kernels —
  now stated on the harness side as well as the ado side

## What came out of it that is not documentation

Fourteen items are recorded in `review_process/open_defects.md` under OPEN, and one — the
F033 add-one missing from the studentized **cell** p-value — is in
`held_for_ruling_2026-09-09.md`, because fixing it moves reported numbers. Each is also
marked at its code site, so none of them is invisible to the next reader.

## Gate

Comments are free against Stata's per-program compiled-size ceiling, so the size headroom is
unchanged. The ado was reloaded after every group, and the full gate was run at the end:

- `verify.sh`: ado loads, 20-leg branch fingerprint identical, battery `TOTAL_DIFF_GROUPS=0`
- 504-config dual build: `unexplained=0 explained=5118 additions=16024 stale_ledger=0
  unscoped=0 prov_fail=0`, 504 of 504 compared
- cert 43/43, newfeat 7/7, battery IDENTICAL vs baseline, `COVERAGE_OK`
- every Python harness file compiles; `check_validate_sync.py` passes all three of its checks

## One landing error, found and fixed

The comment blocks were inserted by an anchored script whose line hints came from the finder
agents' PRE-EDIT line numbers. After ~400 lines of insertions, one anchor (`local ++k`) had
its intended duplicate outside the search window, so the block resolved to the only match
nearby — in `wooldid_icscreen` — and shipped there in `3e0ea38`. Found by an audit that took
every inserted block's first line, located it, and printed its ENCLOSING PROGRAM against the
intended one; 102 blocks checked, one wrong. Moved to `wooldid_cfxplots` in the follow-up
commit. **The lesson for the next pass of this shape: verify the landing PROGRAM, not the
landing line — a stale hint plus a duplicated anchor resolves silently to the wrong site.**


## Outcome, 2026-09-11: all fifteen ruled and built

The fifteen items above were put to the author as multiple choice and all fifteen were ruled.
They were built in three risk-ordered commits - `ed95674` (the harness, first, so the rest is
judged more strictly), `9bb24c6` (the seven that move no number), `82934b3` (the three that
do) - and the ruling table and the closure record live in `open_defects.md`, which now has
nothing under OPEN and nothing under HELD.

Three of the fifteen change reported values, and none of the 504 dual-build configs reaches
those paths - so they would have been invisible to the suite. Each was therefore verified with
a before/after probe and pinned with a cert leg that was checked to FAIL on the pre-change
build:

| item | before | after |
|---|---|---|
| F033 on the studentized CELL p | 22 of 50 cell p-values exactly 0, minimum 0 | none zero, none below 2/(completed+1) |
| D-1 `permute()` + `cfxplotnose` | first cfx matrix `colsof = 1`, no `e(permdist_*)`, cell p missing | `colsof = 6`, `permdist_` posted, a real permutation p |
| D-9 per-family band seed | `crit_att_it` = 2.474072413106145 or 2.492870068598342, depending on what drew first | 2.47972450315949 either way |

The F033 number is the one worth noticing: the finding estimated the gap at about 2/(n+1), and
on the cert fixture it was that plus a long tail of exact zeros - nearly half the cells were
printing a p of 0 that the coefficient surface could not print.

New cert legs: `test_esband.do` B10, `test_permute3.do` Q33, `test_bootstrap_tstat.do` T12.


---

# Phase 2 (2026-09-11): the help file

Second half of the same instruction: *"review it to confirm that every new feature and option
we've added to wooldid is documented there... and when writing documentation, please aim to
maintain the original style of the help file (ie, the style it used 12mo ago)."*

## What "12 months ago" is

`wooldid.sthlp` stood unchanged between 2024-07-31 and 2026-07, so `2a2da17` **is** the
12-months-ago file and is the style reference. Its measurable habits, against which the new
prose was checked:

| marker | 2024 (`2a2da17`) | the campaign-era file it joins | new text as first written | new text as shipped |
|---|---|---|---|---|
| all-caps words | acronyms only (ATT, AME, OLS, BJS, CFX, STATA) | some | 9 emphasis sites | 0 |
| ` - ` dash parentheticals per 1000 words | 0.26 | 6.81 | 13.6 | 0 |
| parentheses per 1000 words | 29.1 | - | 60.3 | 62.2 (16 pairs added; dashes went to commas, colons, semicolons and full stops more often than to parentheses, precisely because this rate was already double the 2024 one) |
| `{ul:}` for emphasis | 20 uses | - | - | the device the two surviving emphases use |
| a directive nested inside link text (`{help x:{cmd:y}}`) | never | never | 1 | 0 |
| authorial asides ("and deliberately so", "the case worth checking first") | never | - | 4 | 0 |

## The audit

14 agents, find then adversarially verify: 579 documentable items enumerated from the ado
(options, `e()` objects, cautions, refusal conditions), each candidate gap put to an
independent verifier that was asked to refute it. 55 confirmed gaps survived, written up in
15 edit batches as `14c30d3`. The style pass above is the follow-up commit.

## Two substantive additions made while re-reading

The stray-companion screen at `wooldid.ado:5965` refuses `reps()`, `ties`, `dots()`,
`rseed()`, `bootalpha()`, `bootdf()`, `bootdisplay()`, `boottstat` and `bootsave()` with an
error on a run without the master option. Only `ties`, `dots()` and `bootalpha()` said so;
`reps()`, `rseed()`, `boottstat`, `bootdf()`, `bootdisplay()` and `bootsave()` now do too.
`rseed()`'s clause says "a run with neither", matching the screen's
`& "\`permute'" == ""` guard rather than overstating it as bootstrap-only.

## Claims spot-checked against the ado after writing

A wrong sentence in user-facing help is worse than a missing one, so the load-bearing new
claims were re-derived from the source rather than taken from the verifier agents:

- `e(wdid_fcguard)` and its "see {opt fastcollapse}" pointer: the target text exists (sthlp
  line 384, the `fastcollapse` entry). Not dangling.
- `fe(state)` noticed but `fe(i.state)` not: `wooldid.ado:400-401` compares whitespace tokens
  with `==`, so a factor-prefixed token cannot match. Confirmed.
- `cfxplottypes()` strips a leading `+` and dedups: `wooldid.ado:5047` strips, `: list uniq`
  at the end of the loop dedups. Confirmed.
- `e(esband_ndraw)` on the simulated engine only: `wooldid.ado:10500` posts it under
  `"\`mode'" != "bootstrap"`, the same condition that makes `e(esband_engine)` read
  `analytic-simulated`. Confirmed.
- the 22-character cfx-name threshold: `wooldid.ado:6455` guards `strlen("boot_\`ftag'_\`nm'") <= 32`
  with `ftag` at most 4 characters (`pctl`, `wald`, `bc`), so the budget for `nm` is
  32 - len("boot_") - 4 - 1 = 22. Confirmed.
- "one digit tighter than reghdfe's own 1e-8": agrees with the claim the same entry has
  carried since 2024 ("By default, reghdfe sets a value of 8").

## What the render check actually proves

`help wooldid` returning rc 0 in batch mode proves the file was found and parsed, not that
every tag painted - batch `help` does not write the rendered page into the log. The real
check is the per-line brace balance (4174/4174, zero unbalanced lines, no tag spanning a
newline), and that is what is being asserted.

## Item #57: the one audit finding that was not a help change

The permutation off-centre caution in the ado cited, as its worked example, "the ame rows
under lattice() with contreatcontrols(basic), where the baseline dose slope is not partialled
out (see the standing issue on the reference cell)". **SI-20 fixed exactly that** on
2026-09-05 (`1df8613`: `basic` now always absorbs `i.lattice##c.(ringer)`, and the measured
reference margin went from 0.27..2.00 to exactly 0 on all 59 affected configs), and the
ticket moved to `wooldid_solved_issues.md`. So the caution was pointing a user at a cause
that no longer exists and at a "standing issue" section that no longer stands. Straight text
defect, no number moves; fixed rather than held.
