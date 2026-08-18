# wooldid test assets

Two layers live here:
- **The regression battery** (this file, below): the per-fix bit-exactness tool.
- **The comprehensive testing & assessment suite** (`suite/` — see `suite/README.md`):
  the campaign-level harness comparing current wooldid against published v1.3 (git
  `main`) at rtol 1e-10 with a documented expected-divergence ledger, plus new-feature
  cross-checks (HonestDiD, external randomization inference, identities) and the
  option-surface coverage assertion. Run it once after a large battery of changes;
  run its labeled subparts for targeted work.
- **The targeted cert suites** (`test_*.do` + `check_validate_sync.py`): self-checking
  per-feature certification scripts (`<NAME>_ALL_PASS` sentinel); the suite's `cert`
  component runs them all.

# wooldid regression battery

A bit-exact regression harness built during the 2026-07 bugfix campaign (B1–B23), in which
it validated every change: it runs wooldid across a fixed set of configurations on frozen
deterministic data and dumps **every** e() scalar, macro, and matrix cell in Stata's exact
`%21x` hexadecimal format, plus the return code and post-run variable list. Comparing two
runs therefore detects *any* change in returned results, down to the last bit.

## Layout
- `runbattery.do` — runs the 29 configurations against the repo's `wooldidinstall/wooldid.ado`.
- `data/` — frozen input datasets (binary, continuous, tiny, simple; deterministic seeds).
- `mkdata.do` — regenerates `data/` from scratch (run from `test/`). Only needed if the
  datasets are ever deliberately changed — doing so invalidates `baseline/`.
- `compare.py` — compares two output directories (`python3 compare.py BASE NEW [--rtol 1e-9]`).
  Parses `%21x` (hex mantissa **and hex exponent**) and compares numerically with a relative
  tolerance; missing-value sentinels and unparseable strings compare exactly. Self-tests on
  startup. Use `--rtol 0` for strict bit-exactness.
- `baseline/` — the committed reference outputs (as of commit `3cb0c41` — the B21
  gen-double re-bless, 2026-07-21; gtools active, Stata 19.0 MP).
- `mkpathological.do` — builds `data/dta_path_*.dta`, the pathological data shapes used by
  the coverage grid (`review_process/assessment_grid.md`, dimension D12): no never-treated
  units, no controls, cohorts lacking pre-/post-periods, missing relative years, single-obs
  ttre-year cells, deeply imbalanced panels with an fe() singleton inducer, poisson
  separation. Each shape's baseline wooldid behavior is recorded in the grid. These are for
  targeted/wave testing; they are not yet battery configs (roadmap item 3).

## Usage
```bash
cd test && mkdir run1 && cd run1
"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do ../runbattery.do
python3 ../compare.py ../baseline .
```
Expected: `29/29 IDENTICAL, TOTAL_DIFF_GROUPS=0` for a change that should not alter results.
When a fix intentionally changes output, the diff must match the predicted signature exactly
(e.g. B12's only signature was one variable leaving the poisson config's varlist) — then the
new outputs become the baseline: `cp run1/*.csv baseline/` and commit with the fix.

Notes:
- Run one battery at a time; each takes a few minutes.
- The battery is bit-reproducible. (gtools is a hard requirement since `40414c7`; the
  non-gtools code path — which had inherent ~1-ulp bysort tie-break noise — no longer
  exists. Comparisons against pre-refactor builds' non-gtools legs need a tolerance.)

## Build-out roadmap (gaps the campaign identified)
1. **Plot-path configs.** The battery never ran `makeplots`/`makecfxplots`/`histogramcohorteffects`,
   yet the cfx matrices are e() returns — a cfx config would have caught B7's truncation
   directly. Add: a continuous `makecfxplots` + `makecfxplotsbysg` config and a binary
   `makeplots` + `histogramcohorteffects` config (`set graphics off` for speed).
2. **A dependency-check negative leg.** (Replaced the old gtools-neutralized leg: since
   `40414c7` gtools is mandatory and the non-gtools path is gone.) Add a variant run that
   loads a copy of the ado with a package renamed in the entry check's foreach list and
   asserts the informative r(198) missing-package error fires.
3. **A drop-inducing config with documented expected values.** The current datasets produce no
   estimator drops, so B10-class weight-semantics changes are invisible here. Add a config with
   an `fe()` singleton cascade whose correct (post-B10) values are pinned in the baseline.
4. **An `if`-condition config** (exercises B23's restriction path and the `(if)==1` semantics).
5. **Author-priority messy-data configs**: deeply imbalanced panels (units spanning 2–12
   periods), very few observations per cohort, scattered missing y/controls.
6. **Provenance stamps.** Have `dumpres` write the repo commit hash and date into each dump so
   stray baselines are self-identifying.
7. **rc expectations.** Encode each config's expected rc in the runner and flag deviations
   loudly, independent of the baseline comparison.
8. **Targeted-suite preservation.** The campaign's per-fix targeted do-files (extended-missing
   scenarios, unconditionalse combos, cfx column checks, graph-survival checks, etc.) currently
   live in a session scratchpad; port the durable ones here as `targeted_*.do` cert scripts.
