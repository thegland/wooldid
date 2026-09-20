# Predicted full-suite result at ff7d978 — written BEFORE any component was read

Components: data, oldnew, cert, newfeat, battery, coverage, degeneracy.

| component | prediction | why |
|---|---|---|
| **data** | skipped | 25+ frozen `.dta` present, so `do_data` returns without regenerating. Regeneration would invalidate `baseline/`. |
| **cert** | **44/44 PASS** | It was 41/44, failing only `test_s6`, `test_stored`, `test_si89`, all on `wooldid_vzanchor() not found` after `cleanmatrices` cleared Mata. That helper is gone. All four suites re-run green individually. |
| **battery** | **IDENTICAL, 0 groups** | `baseline/` was re-blessed at `e2b5fad`; the three commits since move no `e()`. Confirmed twice already, at `65466b0` and at `90c8d43`. |
| **oldnew** | unknown — first run against any of this week's work | The dual-build comparison against published v1.3 with its divergence ledger. The R-batch changed real behaviour (R1, R11, R13, R15 and others), so divergences are expected; the question is whether each is ALREADY LEDGERED. An unledgered divergence is a finding, not a pass. |
| **newfeat** | unknown | HonestDiD / external randomization-inference cross-checks. Nothing this week touched those paths, so no movement expected. |
| **coverage** | unknown | The option-surface coverage assertion. The help audit found all 83 options documented, which is a different check but points the same way. |
| **degeneracy** | unknown | Scans the degeneracy expectations table. R11's escalations and the DEF-2 narrowing both changed which designs refuse, so `degeneracy_expected.csv` may name rows that moved. |

## The falsifier

Any cert suite still red is a defect. Any battery group at all is a defect. For
`oldnew`, only an UNLEDGERED divergence is a defect. For the rest, any failure is
adjudicated on its own terms before being called expected — the rule this week has
repeatedly proved right is that a red suite is red until measured otherwise, and that
"it was probably already like that" must be demonstrated by running the older build,
not assumed.
