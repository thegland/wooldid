# SKEPTIC notes — refuting NFC-01 (nofastcollapse vs default e(V) divergence)

## Plan
- r1.do : own micro panel (30 states x 25 individuals x 8 years, i=state), arms
  NOISE / NOCTRL / Z / AGE / INC / CELLCONST / ONEROW / SHUFOFF / SHUFON.
  Metric split: max rel on ALL of e(V) vs DIAGONAL only vs e(wdidmainresults).
- r2.do : mechanism OUTSIDE wooldid — plain reghdfe + margins on full data vs a
  hand-collapsed one-row-per-cell dataset, plus an EXACT analytic oracle c'Vc.
  Decides whether this is wooldid code or a property of Stata's margins.
- r3.do : (if needed) severity — does any printed digit / decision move.

## Status
- r1 launched.

## Results so far (2026-09-17)
- **r1 (my micro panel, i=state, 30 states x 25 ppl x 8 yrs)**: NFC-01 does NOT
  reproduce.  default vs nofastcollapse max rel on e(V): no controls 0.00,
  controls(z) 0.00, controls(age) 6.84e-16, controls(inc, mean 4e4) 6.83e-16,
  cell-constant control 0.00, aggregated one-row-per-cell 0.00.  Row shuffle
  (collapse off) 2.45e-15, (collapse on) 2.45e-15.  SE[1] hex-identical.
- **r2 (plain reghdfe + margins, hand-built collapse, no wooldid)**: margins'
  delta-method V equals the EXACT oracle c'Vc to 0.00e+00, on full data AND on
  the hand-collapsed one-row-per-cell data with the FIRST row's covariate, at
  mean|xb| up to 3504.  So this is NOT a generic property of Stata's margins.
- **r3 (my own data, the prober's SHAPE: 150 units, 3 cohorts + never, i=cohort)**:
  IT REPRODUCES.  controls(xpct=50+8*x1): 6.59e-06 on e(V) AND on its DIAGONAL;
  no controls 5.65e-16; jointtests irrelevant (6.59e-06 either way).
  LADDER (control = c + 8*x1, only |xb| moves):
     c=0    3.27e-07   c=5  3.27e-07   c=50  6.59e-06   c=500 3.28e-05  c=5000 4.73e-05
  So even a MEAN-ZERO unit-scale control gives 3.3e-07 here, 3 orders above rtol 1e-10.
- Prober's own log p13.log: arm F (row shuffle, collapse OFF) 2.92e-07, arm G
  (shuffle, collapse ON) 2.28e-07 — their repro HEADER says 7e-16 / 9e-8, which
  contradicts their own logged numbers.  On their data the collapse-vs-nocollapse
  gap (5.5e-6) is ~20x the pure row-order gap (2.9e-7), but the row-order gap is
  itself 3 orders above rtol.
- Open: r4 (which cell moves, and self-consistency c=0 vs c=50 within an arm),
  r5 (is the trigger i granularity? same micro data with i=state vs i=cohort).

## r4b (cohort panel, 1500 rows) — the decisive run
PART 1, controls(xc = 50 + 8*x1), default vs nofastcollapse, max rel on e(V):
   att                                   1.92e-16  CLEAN
   att esprelength(1) espostlength(1)    6.47e-16  CLEAN
   att esprelength(2) espostlength(2)    5.35e-07  FAIL
   att esprelength(3) espostlength(3)    1.51e-06  FAIL
   att att_i att_it (no ES)              5.41e-16  CLEAN
   att att_i att_it ES 3/3               6.59e-06  FAIL
   att esprelength(5) espostlength(5)    4.65e-07  FAIL
 => a within-cell-varying control is NOT sufficient: an EVENT-STUDY WINDOW with
    >= 2 pre/post horizons is also required on this panel.
PART 3, which cells move (ON vs OFF, c=50): the PRETREATMENT rows.
   pretreatment att_it se .11036045039018183 vs .11036008692391645  rel 3.29e-06
   pretreatment att    se .11064425394958512 vs .11064417035467353  rel 7.56e-07
   main att_i          se .12688804712927898 vs .12688805058582747  rel 2.72e-08
   pre2                se .14038386747878376 vs .14038385846673826  rel 6.42e-08
PART 4 (MY OWN test, the prober never ran it): adding a pure CONSTANT to the
   control cannot change any reported quantity in exact arithmetic.
   collapse OFF (nofastcollapse): NO cell moves above 1e-9.        <- exact
   collapse ON  (the default):    se/t/ci move up to 3.46e-06.     <- inexact
 => the DEFAULT path is the inaccurate side; the prober's "not solely
    attributable to the collapse / shared instability" framing is wrong here.

## r6 (loose ends, failing panel)
   A1 row shuffle, collapse OFF in BOTH   1.32e-15   <- uncollapsed path is order-INDEPENDENT
   A2 row shuffle, collapse ON  in BOTH   6.59e-06   <- collapsed path is order-DEPENDENT
   B  cell-CONSTANT control, ON vs OFF    2.81e-15   <- within-cell variation necessary
   C  micro panel, 3 onsets, ES 3/3       6.72e-16   <- still clean on the other shape
   D  PRINTED table differs: pretreatment att_it se .11036045 vs .11036009,
      pretreatment att se .11064408 vs .11064417, ci -.24568926 vs -.24568854,
      ub .18727346 vs .18727275, att_i ub .72469403 vs .72469402.

## r8 (one-line patched COPY of the ado: hand controls() to the dose-averaging branch)
   P0 patch INERT, ON vs OFF              6.59e-06
   P1 patch ARMED, ON vs OFF              4.84e-16   <- REMEDY WORKS
   P2 patch ARMED, shuffle, ON in both    1.38e-15   <- order dependence gone
   P3 patch INERT, shuffle, ON in both    6.59e-06
   P4 patched-ON vs OFF, e(b)             1.31e-15   <- point estimates unaffected
 => mechanism: the retained row carries the FIRST row's RAW control value; margins'
    variance depends on it at ~1e-6.  Carrying the cell's weighted MEAN (what the
    code already does for the contreat dose) closes the gap to machine precision.

## Final state
- FINDINGS.md written.  Verdict on NFC-01: PARTIAL (real, reproduced independently;
  trigger and attribution corrected; severity medium).
- SK-01 (sort-order dependence of the printed SEs under the default collapse) is a
  SEEDED sharpening of the prober's arm G, not an unseeded find.
- UNRESOLVED: the third ingredient that makes the cohort panel fire and the micro
  panel immune.  Scale, i granularity, cohort count and ES window were all varied.
