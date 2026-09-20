# fp-estimands — first-principles oracle for att / att_i / att_it / att_itime + ES horizons

Agent: FP-ESTIMANDS. Method: INDEPENDENT REIMPLEMENTATION (numpy OLS + from-scratch weights).

## Files
- `oracle.py` — the reimplementation. Rebuilds `_t`/`_i`/`targetflag` in BOTH reference
  regimes, builds the saturated-cell design matrix (unit FE + time FE + (g,h) cell dummies),
  solves by `lstsq`, rebuilds the four estimand weights from raw data, and does the margins
  aggregations. **Reads none of wooldid's e() output.**
- `wddump.do` — dumps `e(b)`, `e(b_es_*)`, `e(gteffects)`, `e(escomp)` to csv in %21x + %.17e.
- `genA.do` / `run1.do` / `cmpA.py` — Design A (unit-level, unbalanced, 3 cohorts, 8 periods).
- `genB.do` / `genC.do` / `run2.do` / `cmpBC.py` — Designs B (multi-row-per-cell admin shape)
  and C (cohort entering the panel at/after its own onset).

## Build assertion
`which wooldid` = the pinned build, `*! version 1.4 2026-08-19` (run1.log:114-119). Good.

## STATUS

### Oracle validated on Design A (run1.log, cmpA.py)
30 units, t=1..8, cohorts {4,6,never}, 8 rows dropped to unbalance. All four estimands,
both pretreatment rows, every ES horizon, six option regimes. Worst relative difference:

| run | spec | worst reldiff |
|---|---|---|
| A1 | unweighted | 2.03e-14 |
| A2 | `[aw=w]` | 1.23e-14 |
| A3 | `[pw=w]` | 1.23e-14 |
| A4 | `[aw=w] customweightmultiplier(cwm)` | 2.75e-13 |
| A5 | `esfixedbaseperiod esrelativeto(-1)` | 1.45e-13 |
| A6 | `customweightmultiplier(cwm)`, no sample weights | 1.79e-13 |

=> the estimand machinery on a well-behaved unbalanced panel is CLEAN to 3e-13.
A2 vs A3: aw and pw give the identical point estimates (both match the same oracle) — the
invariant holds.

### Seeded leads settled
- **`[fw=]` is unreachable for a user.** ado:620-631 refuses it, `error 198`. run1.log:196
  `FWRC=198`. The internal `wdrepmode` path demotes fw to aweight for everything except the
  regression. Nothing to probe. NOT a finding.
- **cwm enters the numerator of all four estimand weights and the denominator of none** —
  CONFIRMED by code and by the oracle matching at 2.7e-13 with that formula. But it is
  *documented*: sthlp:1107 ("the weights normally used to produce the requested type of
  treatment effect (ATT, ATT_i, etc.) will be multiplied by this variable"). The only gap is
  that the `att_it`/`att_i`/`att_itime` phangs (sthlp:220-236) spell the weight out as
  `w/W_it`, `w/W_ipre|post`, `w/(W_it*T_iprepost)` with no cwm factor, and keep calling the
  result "the simple, unweighted mean". Low-severity doc note, filed as FP-CWM-DOC.
- **ESF window**: `wooldid_validate` ado:5152-5154 bumps `esprelength` by one under
  `esfixedbaseperiod` when `abs(esrelativeto) <= esprelength`, so `esprelength(3)` reports
  pre4/pre3/pre2 (pre1 being the reference). Correct behaviour; my first oracle missed it.

### FINDINGS so far
See FINDINGS.md.

---

## ORACLE VALIDATION TABLE (final) — wooldid vs from-scratch reimplementation

Pass threshold 1e-9 (reghdfe converges at `regtol(9)`); everything came in far tighter.

| run | design | spec | worst reldiff |
|---|---|---|---|
| A1 | 30 units, 8 periods, unbalanced, 3 cohorts | unweighted | 2.03e-14 |
| A2 | " | `[aw=w]` (w varies within unit) | 1.23e-14 |
| A3 | " | `[pw=w]` | 1.23e-14 |
| A4 | " | `[aw=w] customweightmultiplier(cwm)` | 2.75e-13 |
| A5 | " | `esfixedbaseperiod esrelativeto(-1)` | 1.45e-13 |
| A6 | " | `customweightmultiplier(cwm)` unweighted | 1.79e-13 |
| B1-B5 | 14 states x 7 yrs, **2-5 rows per (i,t)**, w varies in-cell | unweighted / aw / pw / cwm / esf | <= 7.26e-15 |
| C1-C4 | 36 units incl. 6 entering at their own onset | unweighted / aw / defaults / esf | <= 1.68e-13 |
| F1 | 40 units, 9 periods, **mid-window gaps + onset-period gap + out-of-window cohort** | unweighted | 1.85e-12 |
| F2 | " | `[aw=w]` | 1.42e-12 |
| F3 | " | `anticipation(1)` | 9.99e-14 |
| F4 | " | `oostreatedcontrols(3)` | 5.03e-13 |
| F5 | " | `esfixedbaseperiod esrelativeto(-2)` | 3.73e-13 |

Quantities compared each time: `main_att`, `main_att_it`, `main_att_i`, `main_att_itime`,
the four `pretreatment_*` rows, and every realized `ES_att_*` / `ES_att_it_*` horizon.
Also checked each time against a SECOND, "documented-semantics" oracle built straight from
sthlp:216-236 (simple mean over cells; simple mean over cohorts of each cohort's
weighted effect; simple mean over cohorts of each cohort's simple mean of cells) — that
agreed too everywhere except under `customweightmultiplier()`, which is FP-CWM-DOC.

### e(escomp) — CLEAN
`cmp_escomp.py`: `ncohorts`, `nunits`, `nclusters`, `nobs`, `sumw_att`, `wshare_att`,
`wshare_att_it`, `relyr`, on 8 configurations across designs A/B/C. Worst reldiff 8.74e-16,
and the help's own claim ("each column sums to 1 across the reported horizons") holds to
1e-15. Nobody had checked this matrix.

### The cell collapse (`nofastcollapse`) — CLEAN
`run6.do` / `genI.do`: 16 units x 7 periods, 2-5 rows per cell, `[aw=w]`, with a control
that VARIES INSIDE THE CELL, a cluster that CROSSES `i`, an `fe()` that varies inside the
cell, `impvalues`, `subgroups()` and `jointtests`. Collapsed vs `nofastcollapse`, comparing
every number in `e(b)`, `e(b_es_*)`, `e(gteffects)`, `e(wdidmainresults)` and `e(escomp)`
(~1000 numbers per pair): worst 1.32e-14. The help's "no answer changes" holds.

### Duplication invariance — CLEAN
`expand 2` on one cohort with the sample weight halved on the copies leaves every estimand
and horizon unchanged at 9.25e-15. Without the weight adjustment only the `att` family
moves (by up to 15%), which is correct: `att` is observation-weighted, so doubling a
cohort's rows doubles its weight.

### Poisson — CLEAN, and it is the control that proves FP-SEMIELAS-FE
`ppmlhdfe` IS called with `d` (ado:1462,1468), so its prediction carries the fixed effects.
Result (`run5.log:2023-2247`): `poisson poisexpresults` agrees between the ppmlhdfe engine
and the explicit-dummy `poisson` engine to ~1e-9, with and without `semielasticity`; and
`poisson` without `poisexpresults` returns exactly the `poisexpresults semielasticity`
number, which is the BS-6 remap being self-consistent. Where `d` is passed the engines
agree; on the linear `reghdfe` path, where it is not, they do not.

### [fw=] — unreachable
`ado:620-631` refuses a user-typed fweight with `error 198` (run1.log:196). The seeded lead
is moot.

## UNTESTED
- `contreat()` / `lattice()` / the ame family (a different oracle: average marginal effects
  of a dose, not a binary onset). FP-SEMIELAS-FE almost certainly extends to
  `contreatelasticitytype(eydx)` and `(eyex)` on the reghdfe path by the same mechanism —
  both divide by the prediction — but I did not run it.
- `timetrends` against the oracle (the trend terms would have to be rebuilt in the
  reimplementation); it was exercised only through the reltime code path.
- `interactivecontrols()`, where the per-row effect is no longer cell-constant and
  `att_it`'s "simple mean of cell effects" reading needs re-deriving.
- `subgroups()` against the oracle (only the engine-split reach test was run).
- The inference surface (SEs, p-values, intervals). Point estimates only, throughout.

## LATE CORRECTIONS (post-review)
- FP-NOPRE: the header DOES print "N = 240 (36 obs dropped from initial sample)" and
  sthlp:73 DOES say "cohorts treated before the sample begins are dropped". The finding was
  narrowed to FP-NOPRE-LATEENTRY (the screen is per-unit against the unit's OWN first
  in-sample period, so a late-entering unit with an in-window onset is dropped too) and
  downgraded to LOW. Verified run11.log:230-246.
- FP-CTZERO: `verbose` on the explicit arm names the dropped term --
  6o._i#3o._t#1o.targetflag#co.d2, one cell's dose slope (run11.log:125-127). The default
  arm prints nothing even under `verbose` (run11.log:172-222). Identical R2 proves the same
  FIT, not the same dropped term; that distinction is now stated in FINDINGS.md.
