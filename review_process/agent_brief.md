# Standing brief for wooldid review/verification agents

Read this in full before doing anything. It is the canonical brief for sub-agents working
on wooldid. Process update 2026-08-18: fixed-format requirements (full battery for every
task, mandatory coverage grid on every report) are retired — scale your evidence to the
question you were asked. The HARD CONSTRAINTS below still bind absolutely.

## The codebase
- Repo: `/Users/hegland/Documents/Research/wooldid_dev` (branch `bugfixes`).
- Program under review: `wooldidinstall/wooldid.ado` — Wooldridge-style staggered DiD via
  reghdfe/ppmlhdfe + margins. Once a single ~3,500-line program, it is now a family of
  subprograms (thin `wooldid` shell, `wooldid_worker`, plus helpers — run
  `grep '^prog def' yourcopy.ado` for the current list).
- Help file: `wooldidinstall/wooldid.sthlp`. Worked example: `wooldid_simple_example.do`.
- **Line numbers in issue descriptions refer to an older commit and drift as fixes land.**
  Locate code on your pinned commit yourself; never trust quoted line numbers.

## HARD CONSTRAINTS
1. **Repo is READ-ONLY for you.** Never modify, create, or delete any file inside it. No
   state-changing git commands (no add/commit/checkout/stash/reset).
2. **Pin to a commit; never copy the working tree.** The working tree changes while you run
   (the coordinator implements other fixes in parallel). Obtain source with
   `git -C <repo> show <sha>:wooldidinstall/wooldid.ado > yourcopy.ado`. Your task prompt
   names the sha(s).
3. **All writing/testing happens in YOUR working directory** (given in your task prompt;
   create it if needed).
4. **At most ONE Stata instance at a time**, and check the shared cap before launching:
   `pgrep -x stata-mp | wc -l` — if 8 or more are running, wait and re-check. Prefer one
   consolidated do-file over many small runs.

## Comparing two builds
Prefer SEPARATE Stata processes, one per build: put each build's `wooldid.ado` in its own
directory, `adopath ++ <dir>` in each process, dump results to disk, compare offline
(`test/compare.py` parses the `%21x` dumps). This avoids rename hazards entirely.

Fallback — two builds in ONE session require renaming every program in the old copy or the
session silently shares whichever helper loaded last:
`sed 's/wooldid_worker/wooldid_old_worker/g; s/wooldid_anticgrid/wooldid_old_anticgrid/g; s/wooldid_anticesf/wooldid_old_anticesf/g; s/wooldid_ptpcheck/wooldid_old_ptpcheck/g; s/wooldid_ptp /wooldid_old_ptp /g; s/wooldid_display/wooldid_old_display/g; s/wooldid_esplots/wooldid_old_esplots/g; s/wooldid_validate/wooldid_old_validate/g; s/wooldid_sgchain/wooldid_old_sgchain/g; s/^cap prog drop wooldid_ptp$/cap prog drop wooldid_old_ptp/; s/^prog def wooldid_ptp,/prog def wooldid_old_ptp,/; s/wooldid_ptp,/wooldid_old_ptp,/g; s/wooldid_gtx/wooldid_old_gtx/g; s/wooldid_permctl/wooldid_old_permctl/g; s/^cap prog drop wooldid$/cap prog drop wooldid_old/; s/^prog def wooldid, eclass$/prog def wooldid_old, eclass/'`
(The helper substitutions match nothing on older builds, so the recipe is safe for
every era. Do NOT use a bare `s/wooldid/wooldid_old/g` — it would rewrite
`ereturn local cmd "wooldid"` and display strings, corrupting e()-level comparisons.)

## Liveness protocol (prevents wasted tokens and lost work)
- **Write findings to disk as you go**: keep a running `NOTES.md` in your working directory
  and update it after each meaningful result. If your session is interrupted (API errors
  happen), your resume — or a replacement agent — must be able to reconstruct state from
  disk alone.
- **Every Stata launch gets a deadline.** Decide up front how long the run should plausibly
  take; if it exceeds ~4x that, kill it, salvage what the log contains, mark the rest
  untested, and move on. Never idle indefinitely waiting for a completion event.
  (Known slow cells: poisson + subgroups + event-study margins can run >12 minutes.)
- **Always deliver a report.** A partial report with an honest untested list beats a perfect
  report that never arrives. Your final message IS the report (raw findings, not
  user-facing prose).

## Running Stata
- The `stata-mp` on PATH is a broken symlink. Use the full binary path:
  `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do yourfile.do`
- Batch mode writes `yourfile.log` to the current working directory and **aborts at the
  first error** — wrap expected failures in `capture noisily` and print sentinels
  (`dis "T1_RC=" _rc`). One command per line: Stata does NOT treat `;` as a separator.
  Batch logs hard-wrap at 80 columns regardless of `set linesize` — grep for short
  mid-line fragments.
- Stata 19 (StataNow) MP. `set maxvar 20000` if you hit r(907).
- **Use the tested harness template** at `review_process/harness_template.do` (copy it into
  your working directory) instead of writing dump/compare machinery from scratch — past
  agents lost hours to self-inflicted harness bugs.
- Installed packages: reghdfe, ppmlhdfe, ftools, gtools, boottest, honestdid.
  gtools is a HARD requirement (since `40414c7`): wooldid checks reghdfe/ftools/gtools at
  entry and exits r(198) if any is missing.
- Deterministic data: `set seed` AND `set sortseed` before generation.
- **Tolerance exception:** the plain reg/poisson estimation paths with a clustered VCE
  (`unconditionalse`, `customvce(cluster ...)`) are NOT bit-reproducible run-to-run —
  the same build run twice in one session can differ in the last 1-2 ulps of V/V_es_*/
  pwsgs/wdid*results (verified 2026-07-26, f11_v3). Two-build comparisons touching those
  cells need `--rtol 1e-12`; everything else stays at rtol 0. For old-vs-new comparisons,
  the project-wide indifference threshold is rtol 1e-10 (the author does not care about
  smaller differences unless they signal a deeper issue).

## Scaling your evidence (replaces the old mandatory coverage grid)
- **You do NOT need to run the full test battery or fill the full assessment grid for
  every task.** Run what probes your specific question. The testing-and-assessment suite
  has clearly labeled components — `battery` (bit-exact vs baseline), `cert` (the
  targeted test_*.do suites), `oldnew` (dual-build v1.3-vs-HEAD comparison, subsettable
  with ONLY=core/trig/covars/adm/plots/screens/patho/ext*), `compare` (re-judge existing
  dumps, no Stata), `newfeat` (HonestDiD / randomization-inference / identity
  cross-checks), and `coverage` — see `test/suite/README.md` (including its "Which
  subpart should I run?" table). Run the relevant subparts rather than the whole thing;
  the full suite (~20-30 min) is the campaign-checkpoint tool.
- Build targeted test cases for the specific issue at hand; do not limit yourself to
  datasets already in the suite.
- The coverage grid (`review_process/assessment_grid.md`) is a tool, not a toll: use it
  when your report makes breadth claims ("no regressions across the option surface");
  skip it for narrow tasks. Either way, **close your report with an honest list of
  relevant-but-untested areas** — untested cells are fine; hiding them is not.
- Calibration (from the author): failures on *realistic* inputs matter most — messy but
  plausible data (imbalanced panels, missing periods, tiny cohorts, administrative ID-style
  codes) and unusual option combinations. Deliberately absurd inputs (1e32 subgroup values)
  are low priority. Exception: Stata's 32-character name limit is a genuine real-world
  hazard — treat name-construction overflows as realistic.
- Pathological test datasets are pre-built: run `test/mkpathological.do` (from `test/`) or
  use `test/data/dta_path_*.dta` if present. Every suite dataset is classified in
  `test/suite/data_classes.csv` (wellbehaved/stress/pathological). Read results by class:
  a pathological-data case only asks "correct results or observable failure?" — it is not
  feature-quality evidence, and a defect found there should also be probed on well-behaved
  input before you diagnose it (see `test/suite/README.md`, "Role of pathological data").

## What wooldid does (orientation)
`wooldid y i t ttre, options` — y outcome, i cohort/unit id, t time, ttre treatment-onset
time (missing = never treated). Builds treated-cohort × relative-period dummies, estimates
via reghdfe (or ppmlhdfe with `poisson`, or plain reg/poisson when
`unconditionalse`/`customvce` given), then computes estimands (att, att_it, att_i,
att_itime) and event-study horizons via margins with specially constructed weights. Two
treatment modes matter for risk analysis: **binary** (default) and **continuous**
(`contreat()` + `contreatcontrols()` + optional `lattice()`); plus subgroups(), weights
[aw/pw], `esfixedbaseperiod`, `impvalues`, `summarystats`, plots. Four distinct sample
concepts exist — user-if, screened `ifvar`, realized `e(sample)`, margins subpops — and
"which sample does this code use" is a recurring bug axis.
Features added in the 2026 campaign (no old-build analogue): `fastnose` (point-estimates-
only fast mode), `permute()` randomization inference, `gteffects` per-cell effects,
`pretrendpower()` Roth-style pre-trend diagnostics, event-study VCV returns
(`e(b_es_*)`/`e(V_es_*)`), `e(escomp)` composition diagnostics, result replay on a bare
`wooldid` call, and `anticipation()`/onset-shift handling.
