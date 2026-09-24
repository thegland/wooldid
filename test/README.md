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
- `runbattery.do` — runs the 59 configurations against the repo's `wooldidinstall/wooldid.ado`.
- `data/` — frozen input datasets (binary, continuous, tiny, simple; deterministic seeds).
- `mkdata.do` — regenerates `data/` from scratch (run from `test/`). Only needed if the
  datasets are ever deliberately changed — doing so invalidates `baseline/`.
- `compare.py` — compares two output directories (`python3 compare.py BASE NEW [--rtol 1e-9]`).
  Parses `%21x` (hex mantissa **and hex exponent**) and compares numerically with a relative
  tolerance; missing-value sentinels and unparseable strings compare exactly (a Stata missing
  is dumped as `+1.0000000000000X+3ff`, `.a` as `+1.0010000000000X+3ff`, ...; since 2026-09-13
  those decode to the `.`/`.a` sentinels, as in `suite/compare_oldnew.py`). A dump's `files|`
  line - the files a config saved, pinned by name - compares as a string; dumps written before
  2026-09-13 have no such line and read as "no files". Self-tests on startup. Use `--rtol 0`
  for strict bit-exactness.
- `baseline/` — the committed reference outputs for the 59 configs (gtools active, Stata
  19.0 MP): the 37 configs that predate 2026-09-13, the 11 permute()/bootstrap() configs of
  the 2026-09-13 suite audit (perm_sgjt_adv .. level_analytic) and 11 more added 2026-09-16
  (`3ce24d4`, `07d7c6b`, `373e764`), the two later blocks blessed additively. Last re-blessed
  in full at `1c09efa` (2026-09-20, run pristine twice): 58 of the 59 files changed, by
  stored-result additions and removals only (no cell moved), and every file now carries a
  `files|` line; `perm_atti` was byte-unchanged. Later re-blesses touch named configs only:
  `e2b5fad` (aw_uncond, boot_ts_gt, boot_ts_sgadv, gt_basic), `aae69dd` (perm_atti) and
  `929ae6c` (band_analytic, cont_lat, perm_cfxnose, perm_tstat, ptp_two); the full record is
  `git log -- test/baseline`. Earlier history: `3cb0c41` (the B21 gen-double re-bless,
  2026-07-21), then a re-bless on 2026-09-01 when the cell collapse became the default
  (`fastcollapse` / `fullprecisionslow`): every file gained the hidden macro `wdid_fcol`;
  the 17 collapse-eligible configs were bit-identical to the previous baseline; `cont_lat`,
  `cont_uz`, `cont_uz_ss` moved at 1e-16..1e-10 in the shared session but were
  bit-identical to the pristine build when each was run alone (session-order contamination
  from the earlier configs' collapse, not a change to the contreat path, which does not
  collapse).
- `mkpathological.do` — builds `data/dta_path_*.dta`, the pathological data shapes used by
  the coverage grid (`review_process/assessment_grid.md`, dimension D12): no never-treated
  units, no controls, cohorts lacking pre-/post-periods, missing relative years, single-obs
  ttre-year cells, deeply imbalanced panels with an fe() singleton inducer, poisson
  separation. Each shape's baseline wooldid behavior is recorded in the grid. These are for
  targeted/wave testing; they are not yet battery configs (roadmap item 3). Every dataset
  here and in `suite/data/` is classified wellbehaved/stress/pathological in
  `suite/data_classes.csv`; a pathological-data test asks only "correct results or
  observable failure?", never whether a feature performs well — see `suite/README.md`,
  "Role of pathological data".

## Usage
```bash
cd test && mkdir run1 && cd run1
"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do ../runbattery.do
python3 ../compare.py ../baseline .
```
Expected: `59/59 IDENTICAL, TOTAL_DIFF_GROUPS=0` for a change that should not alter results.
When a fix intentionally changes output, the diff must match the predicted signature exactly
(e.g. B12's only signature was one variable leaving the poisson config's varlist) — then the
new outputs become the baseline: `cp run1/*.csv baseline/` and commit with the fix.
When configs are ADDED (appended at the END of `runbattery.do`, so no earlier dump can move),
bless additively: require exactly `N IDENTICAL` for every pre-existing config plus the new names
as `NEW config (no baseline)` - compare.py exits 1 until they are blessed, which is the expected
pre-bless state - then copy ONLY the new files (`cp run1/<new>.csv baseline/`). `cp run1/*.csv`
re-blesses every config and hides an unintended move. Run the battery a second time and compare
the new files run to run at `--rtol 0` before committing them, in one commit with `runbattery.do`.

Notes:
- Run one battery at a time. The timing below was measured on 2026-09-13, when the battery held
  48 configs; it has since grown to 59, so treat these as a floor rather than a current figure
  (measured on two fresh processes: 275-292 s for the first 37 plus 90-108 s for the
  permute()/bootstrap() block, 366-400 s in all; a loaded machine took 450 s for the first 37
  alone).
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
