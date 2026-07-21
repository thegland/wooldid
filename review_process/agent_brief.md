# Standing brief for wooldid review/verification agents

Read this in full before doing anything. It is the canonical version of the brief used in the
B1–B23 campaign, updated with the process rules adopted 2026-07-21.

## The codebase
- Repo: `/Users/hegland/Documents/Research/wooldid_dev` (branch `bugfixes`).
- Program under review: `wooldidinstall/wooldid.ado` — a single ~3,500-line eclass program
  implementing Wooldridge-style staggered DiD via reghdfe/ppmlhdfe + margins.
- Help file: `wooldidinstall/wooldid.sthlp`. Worked example: `wooldid_simple_example.do`.
- **Line numbers in issue descriptions refer to an older commit and drift as fixes land.**
  Locate code on your pinned commit yourself; never trust quoted line numbers.

## HARD CONSTRAINTS
1. **Repo is READ-ONLY for you.** Never modify, create, or delete any file inside it. No
   state-changing git commands (no add/commit/checkout/stash/reset).
2. **Pin to a commit; never copy the working tree.** The working tree changes while you run
   (the coordinator implements other fixes in parallel). Obtain source with
   `git -C <repo> show <sha>:wooldidinstall/wooldid.ado > yourcopy.ado`. Your task prompt
   names the sha(s). To run two versions in one Stata session, rename one program:
   `sed 's/^cap prog drop wooldid$/cap prog drop wooldid_old/; s/^prog def wooldid, eclass$/prog def wooldid_old, eclass/'`.
3. **All writing/testing happens in YOUR working directory** (given in your task prompt;
   create it if needed).
4. **At most ONE Stata instance at a time**, and check the shared cap before launching:
   `pgrep -x stata-mp | wc -l` — if 6 or more are running, wait and re-check. Prefer one
   consolidated do-file over many small runs.

## Liveness protocol (prevents wasted tokens and lost work)
- **Write findings to disk as you go**: keep a running `NOTES.md` in your working directory
  and update it after each meaningful result. If your session is interrupted (API errors
  happen), your resume — or a replacement agent — must be able to reconstruct state from
  disk alone.
- **Every Stata launch gets a deadline.** Decide up front how long the run should plausibly
  take; if it exceeds ~2x that, kill it, salvage what the log contains, mark the rest
  untested, and move on. Never idle indefinitely waiting for a completion event.
- **Always deliver a report.** A partial report with an honest untested list beats a perfect
  report that never arrives. Your final message IS the report (raw findings, not
  user-facing prose).

## Running Stata
- The `stata-mp` on PATH is a broken symlink. Use the full binary path:
  `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do yourfile.do`
- Batch mode writes `yourfile.log` to the current working directory and **aborts at the
  first error** — wrap expected failures in `capture noisily` and print sentinels
  (`dis "T1_RC=" _rc`). One command per line: Stata does NOT treat `;` as a separator.
- Stata 19 (StataNow) MP. `set maxvar 20000` if you hit r(907).
- **Use the tested harness template** at `review_process/harness_template.do` (copy it into
  your working directory) instead of writing dump/compare machinery from scratch — past
  agents lost hours to self-inflicted harness bugs.
- Installed packages: reghdfe, ppmlhdfe, ftools, gtools, boottest, honestdid, `_gxtile`
  (egenmore itself NOT installed). gtools presence switches wooldid code paths
  (`gtoolscheck`): to exercise the NON-gtools path, redirect detection in your COPY
  (`sed 's/cap which gtools/cap which gtools_neutralized/'`) — never in the repo.
  Note: the non-gtools path is run-to-run non-bit-reproducible (~1 ulp bysort tie-break
  noise); compare it with a relative tolerance, or fix `set sortseed` for determinism.
- Deterministic data: `set seed` AND `set sortseed` before generation.

## Coverage discipline
- Your report MUST end with the coverage grid from `review_process/assessment_grid.md`,
  marking each relevant cell tested / untested / n-a. Untested-but-relevant cells are not a
  failure — hiding them is.
- Calibration (from the author): failures on *realistic* inputs matter most — messy but
  plausible data (imbalanced panels, missing periods, tiny cohorts, administrative ID-style
  codes) and unusual option combinations. Deliberately absurd inputs (1e32 subgroup values)
  are low priority. Exception: Stata's 32-character name limit is a genuine real-world
  hazard — treat name-construction overflows as realistic.
- Pathological test datasets are pre-built: run `test/mkpathological.do` (from `test/`) or
  use `test/data/dta_path_*.dta` if present. See the grid file for what each covers.

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
"which sample does this code use" is a recurring bug axis; see S2 in the review file's
structural-suggestions section.
