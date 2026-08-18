# Pre-registered signature for the first full suite run (2026-08-18, HEAD 8e25cfa)

Locked BEFORE the run. Any deviation not listed here is stop-and-investigate.

## oldnew (431 configs: 73 ours + 358 mined)
- core/trig/covars/plots/screens: the 23 data-independent configs repeat their validated
  verdicts (0 unexplained). The data-dependent trigger configs now fire their ledger rows:
  trig_extmiss (A6), trig_sgadv_deeppre (rc 2000->0), trig_tt_year/_pois (bounded 1e-6),
  trig_tt_big (any — old trends inert), trig_bigcodes (A10/A11 any), trig_imdegen (A16),
  trig_mintsent_off (A15), trig_mintsent_def (rc 0->459), trig_singleton (A7),
  cov_* (A7/A8b), scr_strclus (rc 109 parity), scr_permexact_t / scr_permstrata (198->0).
- Mined ext configs — EXPECTED UNEXPLAINED CLUSTERS to triage (first run = triage run;
  ledger rows will be added with citations, precisely scoped, after inspection):
  1. adversarial p/line diffs of the svmat storage class (~1e-8..1e-6 rel) on ext grid
     cases (all 216 carry adversarial; only esf references produce the fits).
  2. B16-class LARGE adversarial changes in below-window regime cells (esrel(-2) with
     esprelength < 2, i.e. f2-reference h11-class rows) and gap-adjacent shapes.
  3. e(summarystats) B1/B2 diffs wherever summary stats run (covered by config-* rows).
  4. The 9 blank-rc_new extval configs get their new-side rcs pinned from this run.
  5. ext_contract10_summary_control: guaranteed VALUE divergence at rc parity (A1).
  6. Possible ext parity surprises from `set emptycells keep` preps and oldsyntax cases.
- No config should show UNEXPLAINED matnames diffs after tempvar canonicalization except
  documented label fixes (C-SGlab class) and dimension changes from documented fixes.
- pois-family configs show the _ppmlhdfe_d varlist removal (E-B12) where poisson runs
  clean on both builds.

## cert
11/11 PASS (10 test_*.do sentinels + check_validate_sync).

## newfeat
3/3 ALL_PASS (verified standalone by the builder; this re-runs them via the harness).

## battery
29/29 IDENTICAL vs committed baseline at --rtol 0 (ado untouched since e5445c4).

## coverage
COVERAGE_OK (verified pre-run).

## Run mechanics
- Expect pool.log all OK; poisson+subgroup+wide-ES shards may TIMEOUT at 90 min —
  if so, mark those configs untested and report (documented slow class), do not retry
  blindly.
