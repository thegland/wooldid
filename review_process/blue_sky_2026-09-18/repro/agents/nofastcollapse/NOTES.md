# nofastcollapse differential-testing probe — NOTES

Agent: NFC (nofastcollapse / bit-identity of fast vs slow cell collapse)
Method: DIFFERENTIAL TESTING, pushed into cells existing evidence does not cover.
Build asserted: .../scratchpad/bluesky/build/wooldid.ado (confirmed via `which wooldid` in mkdata.log:52)

## Infrastructure
- `lib_dump.do`  : `wddump` (dumps every e(scalar)/e(macro)/e(matrix) in %21x hex to a text file)
                   `wdcell` (runs one cell twice — arm A = collapse default, arm B = nofastcollapse —
                             from a FRESH `use` of the same .dta each time, and records
                             datasignature/_N/c(k) before and after BOTH arms)
- `cmp.py`       : parses two dumps, compares scalars/macros/matrices elementwise at a
                   relative tolerance; flags presence/dimension/rowname differences separately
                   (a horizon present in one arm only will not show as a numeric diff).
- `mkdata.do`    : builds data/dC.dta (i=cohort, fat cells), data/dU.dta (i=unit, 4 rows per
                   unit-period, RCS), data/dCm.dta (messy: entry/exit/gaps/missing y and x).
                   Uses the recon generator gen_core.do. NO repo test data, NO repo suite.

## Pass criterion
The ado's own comment (wooldid.ado:5605) and sthlp:330 say the collapse is NOT bit-identical —
(S+c1)-(S+c0) leaves ULP residue. Project indifference threshold is rtol 1e-10. The sthlp makes
two harder, MEASURED claims I am testing against:
  - binary: bounded by ~1 ulp of the linear prediction; se/V at ~1e-16 relative.
  - continuous: "no cell moved by more than 1e-9 in relative terms ... largest relative move on a
    quantity of ordinary size was 1.9e-10".
Anything above those is a contradiction of a documented measurement = first-class finding.

## Structural reading done before running (what is in the cell key)
key = `_i' `_t' `targetflag' `ifvar' `t' `ttre' `subgroups' [+ `cluster' iff impvalues]
      [+ `i' and the dose/ringer means iff contreat]
Every margins `over()` in the worker stripes on _t, cluster, subgroups, _i or cfxlatticeplot;
every subpop() keys on ifstatement(=`ifvar'==1), targetflag, t, ttre, _t, subgroups.
cfxlatticeplot is the only one NOT in the key — and makecfxplots is refused without contreat()
(ado:5267) while contreat()+makecfxplots is on the ineligible list (ado:5078), so it is closed.
`coarsecohortcontrols()` expands to i._i#i._t#(x) — NOT interacted with targetflag — so it cancels.
margins df is always passed explicitly as df(`df_r'), so it cannot pick up the collapsed row count.
`local fcnoes` is cleared at every one of the 8 `restore, preserve` sites, so the fcdeg fallback
cannot leave margins running noesample on full data.

## Status
- p01 (noise floor + base sweep, i=cohort: plain/robust/cluster/aw/fw/pw/impvalues x2/controls/subgroups) — QUEUED
- p02 (contreat x contreatcontrols modes, poisson, cwm with exact zeros, subgroups) — QUEUED
- p03 (timetrends, ccc, fe(unit), esfixedbaseperiod+adversarial, anticipation+oos, fastnose,
       histogram/summarystats/gteffects, if-restriction, esrelativeto, messy panel, RCS) — QUEUED
- p04 (state: retired option names, hidden results, datasignature/sortorder leak, user preserve,
       error raised INSIDE the collapse window, reuse after an error) — QUEUED
- p05 (rank-deficient ES covariance in a NON-subgroup run — the one place the help says the
       collapse cannot change WHAT is reported; permute/bootstrap forwarding; ineligible list;
       large-mean outcome; the "row order moves it 19x more" calibration claim) — QUEUED

HARNESS NOTE: the machine is at the 8-Stata cap held by other agents; runs are queued behind it.

## RESULTS (running log)

### p01 base sweep (dC, i=cohort)  — harness noise floor EXACTLY 0 in both arms
| cell | worst rel | verdict |
|---|---|---|
| c01 plain (4 estimands + jointtests + ES) | 2.89e-14 | CLEAN |
| c02 robust | 3.62e-14 | CLEAN |
| c03 cluster(region) | 2.29e-14 | CLEAN |
| c04 [aw=w] | 1.45e-14 | CLEAN |
| c05 [fw=fwt] | r(198) both arms (fweights refused) | n/a |
| c06 [pw=w] | 7.78e-15 | CLEAN |
| c07 impvalues cluster(region) — cluster CROSSING i | 2.98e-15 | CLEAN (confirms sthlp's 1e-12 claim) |
| c08 impvalues cluster(clnest) — cluster NESTING i | 4.28e-16 | CLEAN (ditto) |
| **c09 controls(x1 x2)** | **5.34e-09 on e(V)** | **FAIL — NFC-01** |
| c10 subgroups(sgc) | 1.44e-14 | CLEAN |

### p07 controls() dissection (dCc) — see FINDINGS.md NFC-01
controls(x1) 1.8e-7 | controls(xpct) 5.5e-6 | controls(xdollar) 2.0e-6 | +robust 2.0e-6 |
+cluster(region) 2.6e-6 (joint p 2.5e-8) | contreat+controls 8.8e-7
MECHANISM ARMS, all CLEAN: cell-constant control 3.6e-14, fe(region) 3.6e-14,
controls(i.region) 5.0e-14, no controls 2.9e-14.
PRINTED TABLE DIFFERS (p07.log, DISPLAY_A vs DISPLAY_B): displayed se .16194921 vs .1619492,
displayed CI -.30483918 vs -.30483917, displayed ES estimate .16215286 vs .16215285.

### p02 continuous treatment / poisson / weights-with-zeros / subgroups (dC)  — ALL CLEAN
k01 contreatcontrols(basic) 6.3e-15 | k02 byTime 5.6e-13 | k04 byTtre 5.4e-15 |
k05 ccnoabsorb 8.8e-15 | k06 contreatelasticitytype(dyex) 5.9e-15 | k07 contreat+[aw] 2.6e-14 |
k08 poisson 2.4e-14 | k09 poisson+contreat 8.9e-15 | k10 poisson+dyex 8.9e-15 |
k11 poisson+semielasticity 1.2e-14 (collapse DOES engage — the BS-6 remap verified) |
k12 customweightmultiplier with exact zeros 1.5e-14 | k13 same + contreat 8.2e-14 |
k14 subgroups+contreat 1.6e-14 | k15 subgroups+cluster(cohort) 4.6e-16 with
e(wdid_fcguard)==fallback — the fcdeg screen FIRED and the fallback reproduced nofastcollapse.
(k03 contreatcontrols(byCohort) with i=cohort is r(459) in both arms — unestimable by design.)

### p03 eligible-list spot checks / messy panel / repeated cross-sections — ALL CLEAN
m01 timetrends 1.5e-14 | m03 fe(unit) 7.2e-14 | m04 esfixedbaseperiod+adversarial 2.2e-14 |
m05 anticipation(1)+oostreatedcontrols(2) 5.5e-15 | m06 fastnose 3.5e-14 |
m07 histogramcohorteffects+summarystats+gteffects 2.9e-14 | m08 `if` restriction 1.4e-14 |
m09 esrelativeto(-2) 7.2e-15 | m10 messy panel + controls(x1 x2) 2.85e-13 |
m11 messy + contreat + controls 2.2e-14 | m12 messy + cwm zeros 1.7e-13 |
m13 repeated cross-sections (4 rows/cell) EXACTLY 0 | m15 RCS + subgroups + cluster EXACTLY 0.
NOTE m10: the SAME option (controls) that fails at 5e-6 on dCc is at 2.85e-13 on dCm.  The
magnitude is dataset-dependent, which is what cancellation residue looks like; the failure is
reproducible per dataset (noise floor 0), not run-to-run noise.

### p11 inference-path amplification (does the 1e-6 residue survive into a decision?)
| cell | result |
|---|---|
| v01 `permute(shuffle) permexact cluster(cohort) controls(xdollar)` | **4.26e-14 CLEAN**, `e(permdist_b)` included. A rank-based permutation p does NOT move. |
| v03 `bootstrap(standard) reps(200) rseed(1) cluster(region) controls(xdollar)` | bootstrap b/se/ci CLEAN; the RETAINED ANALYTIC `e(V_analytic)` moves **2.32e-07** — NFC-01 again. |
| v04 `impvalues cluster(region) controls(xdollar)` | **1.85e-15 CLEAN**. impvalues sets `nose`, so there is no delta-method VCE to carry the residue — a strong corroboration of the mechanism. |
| v05 `impvalues` on a NOISELESS design | b 6.0e-10, a far-tail p 2.5e-8 (pathological input; the amplification the help documents) |
| v06 noiseless design + controls | b_es 1.98e-09 relative on a structurally-zero pre-coefficient of -2.3e-4; the ABSOLUTE move is 4.6e-13 against a prediction of ~4010, i.e. ~0.5 ulp — this one IS inside the help's stated bound and is NOT a finding. |
(v02 r(459) both arms: cluster(region) crosses the cohort onsets, so permute refuses. By design.)

## Static verification (no Stata needed), all CLEAN
- The 2026-09-16 consolidation: `grep -n "\`fastcollapse'\|\`fullprecisionslow'\|FASTCOLLAPSE\b\|
  FULLPRECISIONSLOW" wooldid.ado | grep -v NOFASTCOLLAPSE` returns NOTHING. Neither retired name
  survives as an option, a macro reference or a branch; both `syntax` lines (ado:555, ado:4572)
  carry only `NOFASTCOLLAPSE`. wooldid.sthlp mentions neither.
- All EIGHT `restore, preserve` sites in the fcdeg fallback (ado:2822, 2878, 2915, 3031, 3049,
  3093, 3148, 3170) are followed on the next two lines by `local fcnoes` and `local fcsg 0`, so
  the fallback cannot leave a later `margins` running `noesample` against restored full data, and
  cannot re-arm the screen. The single plain `restore` is at ado:3496 under the same `fcon==1`
  guard as the `preserve` at ado:1876.
- `e(wdid_fcncell)` is posted as the cell count when the collapse ran and as `.` otherwise
  (ado:4452-4453), matching sthlp's description.
- ONE loose documentation statement found: sthlp says "The full data is restored before any
  plotting or display work", but `wooldid_hist` (histogramcohorteffects) is invoked at ado:1944,
  INSIDE the preserve window, and the cfx `marginsplot` calls sit at ado:2492-2625, also inside it.
  No measured consequence: the histogram is built in its own frame from a margins matrix, and
  m07_hist matched at 2.9e-14 including e(histogramestimates); the cfx sites are unreachable under
  the collapse because contreat()+makecfxplots is on the ineligible list. Recorded, not filed.

## Remaining queue (runq.sh): pmin (noise floor WITH controls + row-shuffle framing +
## coarsecohortcontrols/ccc_absorb/interactivecontrols), then pALL (repro + p12 + p04 state +
## p10 rows-per-cell ladder), p05, p06, p08.
## The machine has been held at 8-9 Stata instances by other agents for most of this session,
## so these are queued behind the cap.

## OPEN / latent, being measured by dx02.do
`wooldid_fcol` picks the cell's surviving row with
    qui gen long `ord' = _n
    qui gcollapse (sum) `wv' (min) `pick' = `ord', by(`key') merge replace
    qui keep if `ord' == `pick'
If gcollapse types `pick' as FLOAT (Stata's default `c(type)` is float), then once the estimation
sample exceeds 2^24 = 16,777,216 rows the comparison is long-vs-float and a cell can keep two
rows or none - a silent weight double-count or a vanished cell.  dx02.do measures the realized
type and the comparison behaviour at ord values just above 2^24 without needing a 17M-row panel.
I can already rule out the twin risk on the WEIGHT side: if `(sum) `wv'' came back float the
no-control cells would have disagreed at ~1e-7 (the B21 signature), and they disagree at 1e-14.

## STATE AS OF THE LAST WRITE (a replacement can resume from here)
DONE and analysed: p01 (10 cells), p02 (15), p03 (15), p07 (11 controls cells), p11 (6) -
57 paired arms, plus the 4 noise-floor runs. All results are in FINDINGS.md.
ONE FINDING: NFC-01 (continuous `controls()` -> e(V)/se/CI/p disagree up to 5.5e-6, printed digits
change).  Mechanism isolated by three clean arms.  Confidence `certain` (see the four-line
argument in FINDINGS.md - the decisive one is that F/N/df_r/ll/rmse/rss are BIT-IDENTICAL between
the arms in every failing cell, so the regression stage is not the source).
ONE MINOR: NFC-02, sthlp's "full data is restored before any plotting or display work" vs
wooldid_hist at ado:1944 inside the preserve window.  Cosmetic.

STILL QUEUED (blocked all session behind other agents holding 8-9 Stata instances):
  p09   = the repro file itself (validates it runs; its section 0 is the noise floor WITH the
          control present, the one measurement that would have been needed had the bit-identical
          regression scalars not already settled determinism)
  pcrit = pmin (noise floor with controls + ROW-SHUFFLE FRAMING CELL + coarsecohortcontrols/
          ccc_absorb/interactivecontrols) + the adversarial 1e-12 rank-truncation flip test
          (cluster(cohort), 4 clusters vs 9 horizons, + controls(xs4)) + the ROWS-PER-CELL LADDER
          (4/8/32/128/512 rows per cell) which is the pre-registered severity-upgrade criterion
  pdup  = does the collapse key contain a DUPLICATE variable when the user clusters on the period
          or onset variable under impvalues, or on `i' under impvalues+contreat?  (ado:1891-1896
          appends `cluster' and `i' to a key that already holds `t' and `ttre'.)  UNTESTED.
  dx02  = the storage type gcollapse gives `(min) pick = ord' - the >2^24-row latent risk
  p05   = permute/bootstrap forwarding, ineligible-list identity, large-mean outcome
  p06   = 1-row-per-cell shape, where the collapse is a pure re-sort

## Build assertion
`which wooldid` at the top of every probe resolves to
  .../scratchpad/bluesky/build/wooldid.ado   (*! version 1.4 2026-08-19)
verified in mkdata.log, p01.log, p02.log, p03.log, p07.log, p11.log.  The PLUS-path v1.3 was
never in play: `nofastcollapse` parsed without r(198) in all 51 paired cells, and
`e(wdid_fcol)`/`e(wdid_fcncell)` (which v1.3 does not post) came back on every collapse arm.

## LATE RESULTS (pcrit + p09) — see FINDINGS.md "LATE RESULTS" for the full write-up
- noise floor on dCc WITH controls: EXACTLY 0, three ways. Determinism settled.
- FRAMING: row permutation with collapse OFF = 2.08e-13; collapse vs nofastcollapse = 1.40e-06;
  row permutation with collapse ON = 8.70e-08.  sthlp's "a row permutation moves it ~19x MORE
  than the collapse does" is inverted by ~1.3e5.
- rows-per-cell ladder 4/8/32/128/512: 2.4e-15 ... 2.85e-13 with the control, 0 without.
  NO growth with cell size -> pre-registered severity upgrade did NOT fire, severity stays medium.
- adversarial 1e-12 rank truncation with a rank-deficient clustered ES covariance: numbers move
  1.69e-06 but NO rank flip, no horizon dropped, no fcguard asymmetry.
- subgroups + controls: 4.81e-06 on e(V), subgroup joint tests move 1e-8, no structural change.
- coarsecohortcontrols(c.xdollar) 1.38e-06 FAILS; ccc_absorb(xdollar) 5.76e-15 clean;
  interactivecontrols(xdollar) EXACTLY 0 (ineligible in both arms).
- The FIRST self-contained repro (hand-built 120-firm panel) did NOT reproduce: all arms ~6e-16.
  The effect is DESIGN-dependent. repro v2 (repro/repro_NFC01_controls_V.do, also queued as
  p13.do) regenerates the recon panel with gen_core.do so that it does fire, and adds the
  row-permutation arms.

## Memory/state angle — coverage summary
The `wdcell` harness ran `datasignature set, reset` before and `datasignature` after BOTH arms of
every paired cell.  Across 74 paired arms: **zero** signature changes, and `_N` and `c(k)`
unchanged everywhere.  So on the success path the collapse's `preserve`/`restore` leaves the
user's data byte-identical, including its observation order (datasignature is order-sensitive),
and the `recast double` of the contreat variable and the `sort key` inside the window do not
survive.  What remains untested empirically is the ERROR path (an error raised between the
preserve at ado:1876 and the restore at ado:3496), nesting inside a user's own `preserve`, and
the r(198) message a user gets for the retired `fastcollapse`/`fullprecisionslow` spellings -
all three are in p04.do, which stayed queued behind the 8-instance cap.  Statically those paths
look sound: Stata auto-restores a program-scoped preserve on error exit, the single plain
`restore` at ado:3496 sits under the same `fcon=="1"` guard as the `preserve`, and all eight
`restore, preserve` fallback sites reset `fcnoes` and `fcsg` on the next two lines.

## CORRECTION (final) — retracted an over-reach
I initially claimed sthlp:357's "a row permutation moves the same near-zero coefficient about 19
times more than the collapse does" was inverted by 4e10.  That compared a coefficient CLAIM
against a VARIANCE measurement.  Recomputed on the actual near-zero coefficient (b_es_att[pre2]):
row permutation 2.06e-14 vs collapse 1.33e-15 -> ratio 15.5 against the stated ~19.  **The
documented calibration is CORRECT.**  The finding is a uniform SCOPE GAP: coefficients obey every
documented bound (1e-16..1e-15); variances do not (1e-7..5e-6).  Plus the genuinely new arm-G
consequence: the collapse makes e(V) row-order dependent (6.63e-16 -> 8.70e-08) where it is not
otherwise, while leaving coefficients' row-order behaviour untouched (2.06e-14 -> 2.13e-14).

## Candidate design properties that make it fire (NOT isolated — the obvious next bisect)
Firing: dC / dCc (gen_core).  Not firing: my hand-built 120-firm panel, the cs4..cs512 ladder.
Two differences stand out and each is one cell to test:
  (a) gen_core's truetau varies by cohort AND horizon; the non-firing DGPs have a constant effect.
  (b) dC's y contains -0.25*x2, which is OMITTED from controls(xpct)/controls(x1); the non-firing
      DGPs have no omitted covariate.
Also candidates: gen_core's common period shock _gt, and unequal cohort sizes.
