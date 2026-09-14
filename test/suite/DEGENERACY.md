# The `degeneracy` component

A dumps-driven assertion that runs after `oldnew` (no Stata, ~3 s): for every clean-rc config on
well-behaved or stress data it checks that the surface the config was written to test is still
**live** - not b-only, not `-1e100` SEs, not a blanked `pwsgs`, not an empty event-study window -
and compares the set of degeneracies it finds against the curated `degeneracy_expected.csv`.

Why it exists: `oldnew` compares HEAD against v1.3 cell by cell, so a surface that is dead on
**both** builds compares equal and passes, and `coverage_check.py` (static) still credits the
option. That is the "looks covered, is degenerate" class the Stage A audit (2026-09-13) set out to
find, and this component makes it mechanical. As shipped, `degeneracy_expected.csv` holds 175
rows: 131 BY_CONSTRUCTION, 13 PART_RULED, 11 HIDDEN, 11 DATASET_DESIGN, 9 RULED.
0 of those rows carry an `UNRULED` citation and are waiting on the author. A further
12 are `PINNED_DEGENERATE`: ruled 2026-09-13 to be KEPT as pins of their degenerate
surface rather than re-hosted or silently accepted. Those still WARN on every run, on purpose -
the warning is what keeps them visible rather than forgotten - but they are no longer pending.
Count the UNRULED citations, not the HIDDEN class, when you want to know what is still open.

Files: `degeneracy_scan.py` (the scanner), `degeneracy_expected.csv` (the contract), this note.
The `run_full.sh` hook (`do_degeneracy`, component name `degeneracy`, after `oldnew` in `all`)
writes `out/degeneracy.txt` and `out/degeneracy/`.

## What it reads, per config

- the `%21x` dump in `out/run_new/<config>.csv`: every result matrix row (`wdidmainresults`,
  `wdidsgresults_*`, `pwsgs`, `gteffects`, `cfx*`) is LIVE / INFDEG (estimate finite, an inference
  cell missing or `-1e100`) / VOID (estimate missing); `e(V)` and `e(V_es_*)` diagonals; the ado's
  own flags (`perm_*`, `boot_*`, `es_ndegen_*`, `wdid_sgdegen*`, `wdid_vcvstate`); the realized
  `b_es_*` window against the requested `esprelength/espostlength`; exact-zero `b` cells; `N_clust`;
  joint / adversarial / permutation / bootstrap objects present vs requested;
- `configs.csv` + `configs_ext.csv` (options canonicalised with `coverage_check.py`'s abbreviation
  rule), `data_classes.csv` (class), `ledger.csv` (LARGE rows that already name a dead object),
- the shard logs `out/wd_new/*/*.log`: column-0 `Caution:` / `Note:` / `Error:` lines after each
  config's `. * ---- <cfg> (<comp>) ----` marker, 80-column wrap re-joined (a raw grep is useless:
  every shard log opens with the echoed ado source).

## Classes

| class | meaning |
|---|---|
| CLEAN | no flag |
| BY_CONSTRUCTION | only `ADV_IMPOSSIBLE_WINDOW` / `ADV_POOLED_NOP` / `JT_WINDOW_SHORT` / `SG_ONE_LEVEL`: the config's own options make the surface impossible |
| INFO | only `NOSE_BY_DESIGN` / `NCLUST_TINY` / `BES_ABSENT(impvalues\|permute)` / `SGDEGEN(:report)`: values present, or no variance by design; never a stem |
| RULED | every substantive flag's e() object is named by a LARGE ledger row THIS FILE CITES for the config (see the citation gate below) |
| PART_RULED | the ledger names some of the dead objects |
| HIDDEN | nothing rules it |
| DATASET_DESIGN, DATA_EDGE | curated only (the CSV): the data was built to lack the cells; the scanner sees HIDDEN or PART_RULED |

Flags that name no e() matrix of their own (`WINDOW_SHORT`, `B_ZERO`, `ES_NDEGEN`, `VES_ABSENT`,
`BES_ABSENT`) count as a degeneracy of `mat b`, so only a LARGE ledger row on `b` (or `*`) rules
them (critique A4).

## The contract: `degeneracy_expected.csv`

Columns `config,class,stems,citation[,ncaution]`; `#` rows are comments. One row per clean-rc
good-data config that is **allowed** to show the listed flag stems (flag names without their
parentheses). The discipline mirrors the ledger's:

- a stem that fires and is not listed -> `DEGEN_FAIL <cfg>: NEW degeneracy [...]`
- a listed stem that no longer fires -> `DEGEN_STALE <cfg>: STALE expectation [...]` (retire it)
- a row whose config is gone from the config files, or is pathological-class / expected to error
  -> `DEGEN_STALE` (it can never be asserted)
- a config with **no dump** -> `DEGEN_FAIL` for every config regardless of class (a timed-out
  shard, or an `ONLY=` subset `oldnew` run; never silently skipped - critique A2)
- a good-data config whose rc differs from `rc_new`, or a config with no data class -> `DEGEN_FAIL`
- malformed rows (duplicate config, unknown class, empty stems) -> `DEGEN_BADROW`
- never failing: `DEGEN_WARN` (class moved vs the CSV; a cited ledger id that no longer exists;
  a pinned `ncaution` that moved; degen-only options and option PAIRS over the current dumps; the
  count of `UNRULED` rows) and `DEGEN_NOTE` (old-vs-new caution diff; skipped digest).

**The citation gate (author ruling 2026-09-14).** For RULED/PART_RULED the `citation` column is
not documentation - it is the INPUT to the classifier. `degeneracy_scan.py`'s `ledger_known`
credits coverage only from the ids cited here (`candidates`). Before this it credited ANY ledger
row matching the same e() MATRIX, which promoted `sg_degpre` to RULED on the strength of
`F1SG-degpre-m-b`, a row about a positional insertion in `b` that says nothing about the
surface in question, and `ext_cov10_rcs_pweight_ddd_c1` on the same shape. Nothing in `kind`,
`pattern` or the matrix name separates those from the honest cases - the difference lives only
in the prose - so a human has to assert the tie, and this column is where. It was NOT done by
qualifying the ledger's `pattern` with the dead stem: `pattern` is also `compare_oldnew.py`'s
cell-match key (`compare_oldnew.py:216`), so narrowing it there would have stopped those rows
explaining the cells they exist to explain. Consequence to expect: adding a flag to a config no
longer silently self-rules because some unrelated row touches the same matrix - you must cite
the row, and a `DEGEN_NOTE` lists configs where the whole ledger WOULD have credited more.

`citation` is free text by convention: ledger ids separated by `;` for RULED/PART_RULED,
`BY_CONSTRUCTION`, `DATASET_DESIGN: <why>`, or `UNRULED[: <what>]` for a degeneracy that still
needs an author ruling (re-host, twin config, ledger row, or accept). The optional `ncaution`
column pins HEAD's per-config console-caution count; blank = not pinned. `expectations.csv` in
the output directory is the scanner's GENERATED version of this file (with `ncaution` filled) -
copy rows from it when curating, never assert against it.

Pathological-data configs are never asserted for stems (README, "Role of pathological data");
they appear in `summary.txt` as information only.

## Running it

```
python3 test/suite/degeneracy_scan.py --dumps out/run_new --configs configs.csv --configs configs_ext.csv \
    --data-classes data_classes.csv --ado ../../wooldidinstall/wooldid.ado --out out/degeneracy \
    --expect degeneracy_expected.csv [--ledger ledger.csv] [--logs-new out/wd_new] [--logs-old out/wd_old] \
    [--dumps-old out/run_old]
```

`--ledger` defaults to `ledger.csv` beside the first `--configs`; `--logs` is a synonym of
`--logs-new`; any other flag is a usage error (exit 2). Without `--expect` it scans only and exits
3. The last stdout line is `DEGEN_CHECK_OK` (exit 0) or `DEGEN_CHECK_FAIL: ...` (exit 1), preceded
by one `DEGEN_SUMMARY` line of counts. Outputs in `--out`: `degeneracy_table.csv` (one row per
config, every signal, CLASS, FLAGS), `option_inverted.txt` (per-option live-coverage column and
the DEGEN-ONLY option / option-PAIR lists), `cautions.csv` (+ `cautions_old.csv`,
`caution_diff_old_new.csv` with `--logs-old`), `summary.txt`, `expectations.csv`.

It needs the dumps of a FULL `oldnew` run: after `run_full.sh oldnew ONLY=<subset>` every other
config is NO_DUMP and fails by design. Run it on `out/run_new` only - the v1.3 dumps trip era
flags (v1.3 posts no `b_es_*`), which is exactly what proves the assertion is not vacuous (376
FAIL / 41 STALE against `out/run_old`).

## Maintenance

- A new config on good data that shows any non-INFO flag needs a row (the run fails NEW until it
  has one); a config whose degeneracy is fixed fails STALE until its row is retired.
- A rename in the ado of an e() name the surface contract keys on (`surface()` in the scanner)
  needs a scanner update - the same class as `coverage_check.py`'s syntax-line parser.
- `WINDOW_SHORT` cannot tell "data edge" from "the package dropped a horizon"; that is what the
  citation column records.
- The old-vs-new caution diff is a permanent block while `oldnew` stages `main` = v1.3 (237 configs
  up on 2026-09-11); it is reported as one NOTE line and a CSV, not per-config WARNs. The
  regression signal for cautions is the optional `ncaution` pin.
- The census sentence above ("holds N rows: ...", and the UNRULED count) is MACHINE-READ by
  `degeneracy_scan.py` and a disagreement with `degeneracy_expected.csv` is a `DEGEN_FAIL`.
  Keep its shape: `holds <N> rows: <n> <CLASS>, ...` on one sentence, and
  `<n> of those rows carry an ...UNRULED...`. Rewording freely is fine; dropping the sentence
  downgrades the check to a WARN.
