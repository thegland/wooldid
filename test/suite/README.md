# wooldid comprehensive testing & assessment suite

The campaign-level assessment harness. It answers two questions:

1. **Old-vs-new fidelity** — does current wooldid (git HEAD) reproduce published
   wooldid v1.3 (`git main`) point estimates and SEs at relative tolerance **1e-10**
   across the full feature surface and a ladder of datasets (minimal → rich → realistic
   administrative-scale → pathological), except where a documented, intentional bugfix
   explains the difference?
2. **New-feature validity** — do the features with no old-build analogue (fastnose,
   permute, gteffects, pretrendpower, F2 VCV returns, escomp, replay, anticipation)
   satisfy external benchmarks (HonestDiD), external re-implementations (randomization
   inference), and internal identities?

**Cost calibration (measured on this machine):** the full `./run_full.sh all` run is
**~20–30 minutes** at pool=6 — run it at every campaign checkpoint. Individual
components and `ONLY=` subsets run in seconds-to-minutes and are the per-fix tools.
The orchestrator **exits nonzero if any component fails** and appends every verdict to
`out/SUITE_REPORT.md`.

## Components

| component | what it does | builds | cost (measured) |
|---|---|---|---|
| `data` | builds the frozen datasets in `suite/data/` if absent: ours (`mkdata_suite.do`), the realistic admin panel (`mkdata_admin.do`), the colleagues' portable fixtures (`extfix/`) | — | ~3 min (skipped when present) |
| `oldnew` | THE dual-build comparison: stages old (git `main`) and new (**git `HEAD`** — never the working tree), runs every config on both in separate Stata processes, dumps e() + e(sample) digests in %21x, compares with `compare_oldnew.py` | both | ~15 min pooled |
| `compare` | re-judges EXISTING dumps against the current ledger/configs (no Stata) — the triage loop tool | — | seconds |
| `cert` | the targeted `test/test_*.do` suites + `check_validate_sync.py` — HEAD certification (fresh logs enforced; column-1 `*_ALL_PASS` + no `*_FAILS=` grading) | new | ~2 min pooled |
| `newfeat` | cross-checks in `suite/newfeat/` (7 scripts): external randomization-inference re-implementation (exact, sampled, unequal-count, t-stat, probabilities), HonestDiD integration + the 11-case honestdid grid, the Roth `pretrends`-package crosscheck, the reghdfe common-onset oracle (b_es/V_es at machine precision), the 22-config fastnose parity grid, feature identities incl. the V_es total-variance check | new | ~8 min |
| `battery` | the standing bit-exact regression vs `test/baseline` (`--rtol 0`) | new | ~3 min |
| `coverage` | static assertion: every syntax-line option exercised by ≥1 config, with **clean-rc gating** (an option covered only by erroring configs is flagged ERR-ONLY) and dead allow-zero detection; allow-zero = `update,permprobs` only. Also joins configs to `data_classes.csv` (warns PATHO-ONLY, fails on table drift — see "Role of pathological data") | — | seconds |

`oldnew` subcomponents (subset with `ONLY=core,trig,...`): `core` (the 29 battery
configs), `trig` (divergence-trigger configs), `covars`, `adm` (the realistic
administrative panel: FIPS-style counties, calendar years, unequal cohorts, unbalanced
entry/exit, lognormal population weights — carries the only coverage anywhere of
`customvce(cluster)`, poisson+pweight, clean factor-variable `fe()`, subgroups×summarystats,
unconditionalse+pweight, and the saveplots×subgroups filename pin), `plots`, `screens`
(rc-contract deck), `patho` (the pathological shapes), `boot` (the bootstrap-SE deck:
both engines, the re-id leg, weighted, subgroups, joint tests, observation-level), and the mined
`extgrid`/`extlattice`/`extcase`/`extval`.

A separate DEEP-VALIDATION component lives in `mc/` (own orchestrator `mc/run_mc.sh`,
own README): preregistered Monte-Carlo size/power calibration of the statistical
claims — RI p-value exactness, adversarial line-test size, pretrendpower power, and
bootstrap standard-error calibration (`mc/mc_bootse.do`).
Deliberately NOT part of `all` (~1.5 h); run it when inference machinery changes.

## Usage

```bash
cd test/suite
./run_full.sh all                 # the full campaign run (~20-30 min)
./run_full.sh data oldnew         # dual-build comparison only
ONLY=screens ./run_full.sh oldnew # one oldnew subcomponent (judges only those configs)
./run_full.sh compare             # re-judge existing dumps after a ledger/config edit
./run_full.sh cert newfeat        # HEAD certification only
STATA_JOBS=4 ./run_full.sh all    # throttle the pool (hard project cap: 8 instances)
```

Outputs land in `suite/out/` (gitignored). `pool.log` records per-shard OK/TIMEOUT/
INCOMPLETE with durations. On full (non-`ONLY`) runs the dump dirs are cleared first —
mixed-vintage dumps cannot be compared.

## Comparison semantics (`compare_oldnew.py`)

- **Universe accounting:** the configs to judge come from the config CSVs; a configured
  case missing a dump, a stray dump, a truncated dump (missing its `end|` sentinel), or
  mixed provenance stamps within a side are all failures. The summary prints
  `COMPARED N of M` and the exact `prov` (build label + git sha) pair.
- **rc gate:** if the builds' rcs differ, or either is nonzero, the config is judged on
  the rc pair alone against `rc_old`/`rc_new` (an erroring old build leaves residual
  internal `margins` e() state). NOTE: both-error parity pins assert the **return code
  only** — error messages are out of contract.
- **Content (both rc 0):** shared objects at rtol 1e-10 with **no absolute floor**.
  Stata missing codes decode to sentinel strings (`.`/`.a`/…), so a missing-code flip
  is a string difference that no bounded ledger row can absorb. The **e(sample)
  digest** (`esamp|` lines: N + per-variable sums over the realized sample) catches
  sample-composition changes even when b/V coincidentally agree. `files|` lines pin
  saveplots output names. New-only objects are informational ADDITIONS; old-only
  objects are REMOVALS = failures unless ledgered.
- **Ledger (`ledger.csv`):** every row cites the fix that explains the difference.
  LARGE rows (empty bound) come first and are scoped to named configs or narrow globs;
  bounded rows are set at **~10× the measured ceiling of their class** and only absorb
  finite numeric diffs within the bound — they are noise-vs-wrong-answer separators,
  NOT a forward tolerance policy (forward regressions are the rtol-0 battery's job).
  Rows that stop firing go STALE and fail the run; rows whose glob matches nothing are
  UNSCOPED (fatal on full runs). `--verbose` prints each bounded row's worst fired
  reldiff against its bound.
- **Accepted blind spots (documented):** the total-divergence trigger configs
  (`trig_tt_big`, `trig_allest_drop`, `trig_extmiss`, `trig_mintsent_off`,
  `trig_esfsup_ref`, `cov_kitchen`, `cov_icbad2`) carry unbounded `*` rows by design — they exist to prove a documented
  LARGE fix fires, not to constrain its magnitude. Everything else is bounded or named.

## Tolerance policy

- Project indifference threshold: **rtol 1e-10** (author directive).
- Documented carve-out classes above 1e-10, all bounded in the ledger with measured
  ceilings recorded in the reason column: svmat/predict float→double (~4e-8 rel
  estimates, tail-amplified p), B21 gen-double (~1e-8 rel; interactivecontrols
  regressor copies up to ~1.3e-4), timetrends-poisson IRLS (run-varying: within 1e-10
  one run, ~1.2e-5 the next — bounded 1e-4). Degenerate near-zero-SE cells amplify ULP
  noise run-varyingly (pwsgs observed at 2.2e-7 and 1.1e-4 across runs), so those
  bounds are calibrated across runs, with the evidence in each row's reason. Timetrends
  on the reghdfe path at calendar-year scale agrees **within 1e-10** (no carve-out).
- `unconditionalse`/`customvce(cluster)` are 1–2 ulp non-reproducible run-to-run on
  BOTH builds — inside 1e-10, but never compare those cells at `--rtol 0`.

## Role of pathological data (author doctrine, 2026-08-29)

Every dataset (and inline-prep config) the comparison configs run on is classified in
`data_classes.csv` — `wellbehaved` / `stress` / `pathological`, one row per dataset
plus `config:<name>` rows for inline-prep cases and deliberate per-config overrides
(an override wins over the dataset lookup). The `coverage` component consumes the
table and fails when it goes out of sync (unclassified config or stale row), so the
marking stays live as configs are added.

How to read a result by class:

- **`pathological`** data is engineered so a trustworthy answer is impossible (no
  controls, no never-treated units, single-obs cells, separation, sentinel codings,
  knife-edge estimability). A config on such data answers ONE question: *does wooldid
  handle this input correctly — true correct results, or an observable failure?* The
  catastrophic outcome it guards against is plausible-looking wrong numbers. It is NOT
  evidence that a feature functions well; do not cite a passing patho config as
  feature-quality evidence.
- **`wellbehaved`** (and, with judgment, `stress`) data is where feature-quality
  evidence lives: a diff there speaks to whether the feature works as intended on
  useful input.
- The boundary is admittedly imperfect, and edge-case failures can reveal deeper
  errors — which is why both kinds stay in the suite for the features that matter.

Two working rules follow. **Coverage:** features should tend to be tested on BOTH
kinds of data. The `coverage` component warns `PATHO-ONLY` when an option's clean-rc
coverage exists only on pathological-class data (never shown working where a good
answer is expected — fix by adding a well-behaved config). The inverse — no
pathological exposure — is reported informationally but is NOT a defect: per the
effort-allocation principle (`wooldid_standing_issues.md` preamble), not every feature
needs a pathological probe. **Diagnosis:** when an error surfaces on pathological
data, probe the same feature on well-behaved input before diagnosing; combining the
two signal streams (how it fails on bad data, how it behaves on good data) is not
guaranteed to add value but is a good first path to investigate.

Out of scope for the table: the `mc/` and `newfeat/` components generate purpose-built
calibration data inside their own scripts, and the `screens` deck's subject is
option-pathology (bad option combinations, usually on fine data) — the config
`component` column already carries that purpose.

## Which subpart should I run? (sub-agents)

- Touched estimation/weights/margins code → `battery` + the nearest `cert` suites; add
  `ONLY=core ./run_full.sh oldnew` if the change claims old-parity; `ONLY=adm` for a
  realistic-data pass.
- Touched validation/screens/error text → `ONLY=screens ./run_full.sh oldnew` + the
  matching cert suite.
- Touched a new feature → its cert suite + `newfeat`.
- Touched the bootstrap feature → `ONLY=boot ./run_full.sh oldnew` + the four
  `test/test_bootstrap_*.do` cert suites (and `mc/mc_bootse.do` if the inference math
  changed).
- Added an option → `coverage` (and add a config for it).
- Edited the ledger or configs → `./run_full.sh compare` (seconds, no Stata).
- Claiming "no old-vs-new regressions anywhere" → the full `oldnew` run.

## Maintenance

- **Adding a config:** one row in `configs.csv` (inline prep → `preps.json`). Set
  `rc_old`/`rc_new` when either build is expected to error; `slow=1` for poisson-family
  or wide-ES-subgroup configs.
- **Adding a ledger row:** only with a citation and the tightest bound that fires,
  ordered into the correct section (LARGE before bounded).
- **Regenerating mined configs:** `mine_external.py` writes `configs_ext.gen.csv` /
  `preps_ext.gen.json` and NEVER touches the curated `configs_ext.csv`/`preps.json`
  (which carry hand-verified rc pins) — diff and merge deliberately. The miner needs
  the untracked `wooldid_external/` bundle; the suite itself never does at run time.
- The `msg_old` column in `configs_ext.csv` is triage-only metadata (the colleagues'
  Stata 17 old-build error messages); no tool consumes it.
- Frozen battery assets (`test/data/`, `test/baseline/`, `test/runbattery.do`) are
  untouched by this suite.

## Known limits / future work


- Both-error parity pins compare rc only (messages out of contract; console output is
  never diffed — message texts changed program-wide by design).
- The colleagues' oracle/metamorphic library is not ported (future component).
- `.gph` files are compared by NAME only (`files|`), never by content (embedded
  timestamps).
- `ccc_absorb()` is exercised but errs identically on both builds (defect ticket in the
  certification report); `cov_icbad` under-triggers the B21(b) collinearity class
  (levels ~1e6 are float-representable — a ~1e8 variant would pin it old-side).
- Performance is not compared beyond shard durations in `pool.log`.
