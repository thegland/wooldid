# wooldid comprehensive testing & assessment suite

The campaign-level assessment harness. It answers two questions:

1. **Old-vs-new fidelity** — does current wooldid (repo HEAD) reproduce published
   wooldid v1.3 (`git main`) point estimates and SEs at relative tolerance **1e-10**
   across the full feature surface and a ladder of datasets, except where a documented,
   intentional bugfix explains the difference?
2. **New-feature validity** — do the features with no old-build analogue (fastnose,
   permute, gteffects, pretrendpower, F2 VCV returns, escomp, replay, anticipation)
   satisfy external benchmarks (HonestDiD), external re-implementations (randomization
   inference), and internal identities?

**Cost calibration:** the full run (`./run_full.sh all`) takes hours and is meant to run
ONCE after a large battery of changes — not per fix. Every component (and `ONLY=` subsets
of the dual-build component) runs standalone; per-fix work should run only the relevant
subparts. Sub-agents: see "Which subpart should I run?" below.

## Components

| component | what it does | builds | cost |
|---|---|---|---|
| `data` | builds all suite datasets into `suite/data/`: our generators (`mkdata_suite.do`) + the colleagues' portable fixture generators (`extfix/`, 12 datasets) | — | ~2 min |
| `oldnew` | THE dual-build comparison: stages old (git `main`) and new (HEAD) builds, runs every config in `configs.csv` (+ `configs_ext.csv` if present) on both in separate Stata processes, dumps full e() surfaces in %21x, compares with `compare_oldnew.py` | both | the heavy one (scales with config count; shards run in a job pool) |
| `cert` | the 10 targeted `test/test_*.do` suites + `check_validate_sync.py` — HEAD certification (each self-checking, `*_ALL_PASS` sentinel) | new | ~20 min pooled |
| `newfeat` | cross-checks in `suite/newfeat/`: randomization-inference external re-implementation, HonestDiD integration, feature identities (anticipation shift, escomp recompute, F2 diag, F7 formulas, fastnose/gteffects spot identities) | new | ~10 min |
| `battery` | the standing 29-config bit-exact regression vs `test/baseline` (`--rtol 0`) | new | ~5 min |
| `coverage` | static assertion: every option in the wooldid syntax line is exercised by ≥1 config (`coverage_check.py`) | — | seconds |

`oldnew` subcomponents (subset with `ONLY=core,trig,...`):
`core` (the 29 battery configs), `trig` (divergence-ledger trigger configs), `covars`
(controls/fe/interactivecontrols/ccc/cwm — options with no coverage elsewhere), `plots`
(plot machinery, graphics off), `screens` (the rc-contract deck: inputs where old and new
return different error codes by design), and — mined from the colleagues' bundle —
`extgrid`/`extlattice`/`extcase`/`extval` (see `mine_external.py`).

## Usage

```bash
cd test/suite
./run_full.sh all                 # the full campaign run (hours)
./run_full.sh data oldnew         # dual-build comparison only
ONLY=screens ./run_full.sh oldnew # one oldnew subcomponent
./run_full.sh cert                # HEAD certification suites only
STATA_JOBS=4 ./run_full.sh all    # throttle the pool (hard project cap: 8 Stata instances)
```

Outputs land in `suite/out/` (gitignored): per-side dumps (`run_old/`, `run_new/`),
`oldnew_compare.txt`, `cert.log`, `newfeat.log`, `coverage.txt`, `pool.log`
(per-shard OK/TIMEOUT/INCOMPLETE), and a cumulative `SUITE_REPORT.md`.

## Comparison semantics (`compare_oldnew.py`)

- **rc gate:** if the builds' rcs differ, or either is nonzero, the config is judged on
  the rc pair ALONE, against the `rc_old`/`rc_new` expectations in `configs.csv` (an
  erroring old build leaves residual internal `margins` e() state — content would be
  noise). Blank expectation = 0.
- **Content (both rc 0):** shared objects compared at `--rtol 1e-10` with **no absolute
  floor** (deliberately stricter than `test/compare.py`'s rtol>0 mode, which carries a
  hidden 1e-12 absolute floor). New-only objects (the F2/F8/F11 additive surface, hidden
  `wdid_*` macros, feature returns) are informational ADDITIONS. Old-only objects are
  REMOVALS = failures unless ledgered (the old e() surface is a strict subset of the new
  one, verified by audit).
- **Ledger (`ledger.csv`):** the expected-divergence register. Every row cites the fix
  that explains the difference. Rows with a numeric `rtol` only *bound* a difference
  (storage-noise classes); rows with empty `rtol` explain any difference in the matched
  object (reserved for documented LARGE fixes). **A ledger row that never fires on a
  matching config is STALE and fails the run** — the ledger cannot rot silently.
- The verdict per config is `MATCH`, `EXPLAINED`, `RC_EXPLAINED`, or `DIFFERS`; the run
  passes iff nothing is UNEXPLAINED and nothing is STALE. **Every UNEXPLAINED difference
  is a finding**: either a real regression, an undocumented behavior change (document it
  or fix it), or a suite bug. Never widen a ledger bound to make a run pass — that
  inverts the tool's purpose.

## Tolerance policy

- Project indifference threshold: **rtol 1e-10** (author directive: differences below
  this don't matter unless they signal something deeper).
- Documented carve-outs above 1e-10 (bounded ledger rows): svmat/predict float→double
  (~4e-8 rel on IM estimates, tail-amplified p), B21 gen-double (~1.2e-8 rel),
  centered-timetrends (~1e-8 rel on 4-digit years).
- `unconditionalse`/`customvce(cluster)` are 1–2 ulp non-reproducible run-to-run on BOTH
  builds — comfortably inside 1e-10, but never compare those cells at `--rtol 0`.

## Which subpart should I run? (sub-agents)

- Touched estimation/weights/margins code → `battery` (bit-exact vs baseline) + the
  `cert` suites nearest the change; add `ONLY=core ./run_full.sh oldnew` if the change
  claims old-parity.
- Touched validation/screens/error text → `ONLY=screens ./run_full.sh oldnew` + the
  matching cert suite (`test_advguard`, `test_sgsplit`, ...).
- Touched a new feature (fastnose/permute/gteffects/roth/replay) → its cert suite +
  `newfeat`.
- Added an option → `coverage` (and add a config for it).
- Claiming "no old-vs-new regressions anywhere" → the full `oldnew` run; nothing less
  supports that claim.

## Maintenance

- **Adding a config:** one row in `configs.csv` (inline data prep goes in `preps.json`).
  Set `rc_old`/`rc_new` if either build is expected to error; set `slow=1` for
  poisson-family or wide-ES-subgroup configs (isolated shards, longer timeout).
- **Adding a ledger row:** only with a citation (commit sha / solved-issues entry) and
  the tightest bound that fires. Config-glob narrowly.
- **`configs_ext.csv`:** generated by `mine_external.py` from the colleagues' validation
  bundle (`wooldid_external/wooldid_test_comparison_code`, untracked); the generated CSV
  is committed so the suite never depends on the bundle at run time.
- The frozen battery assets (`test/data/`, `test/baseline/`, `test/runbattery.do`) are
  unchanged by this suite; the battery remains the per-fix bit-exactness tool.

## Known limits / future work

- The colleagues' independent/metamorphic **oracle library** (their `tests/cases/oracles.do`)
  is not ported — a candidate future component.
- Console output is not diffed (message texts intentionally changed program-wide: C12
  Error:/Caution: prefixes); e()/rc/varlist comparison is the contract.
- HonestDiD anchors F2's VCV returns (rc 0 + Original-CI identity). F7's power/MDT
  numbers have **no external package analogue** — they are checked against
  first-principles formulas, not an external implementation.
- No large-panel (10^5+ obs) performance/precision configs yet.
