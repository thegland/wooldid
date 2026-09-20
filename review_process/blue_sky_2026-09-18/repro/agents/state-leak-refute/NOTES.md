# SKEPTIC notes -- refuting the "state-leak" findings

Rule followed: every claim rebuilt on MY OWN data and MY OWN do-files. I never ran the
prober's scripts. DGP used throughout: 40 "states" x 12 years (2005-2016), 3 treated
cohorts of 10 (onset 2009/2011/2013), 10 never-treated states, unit + year effects,
true ATT 0.40, sd(eps)=0.5.

Files:
- rA_missctrl.do / .log      -- finding 01, eight legs
- rB_verbose_and_remedy.do / .log -- finding 01 exact-hex + window remedy; finding 02 full
- rC_state_catchall.do / .log -- NEW: settings/globals/RNG/sort/frames catch-all
- rD_mata_wooldid.do / .log  -- finding 03, wooldid alone (own process)
- rE_mata_reghdfe.do / .log  -- finding 03, reghdfe alone (own process)
- rF_mata_gtools.do / .log   -- finding 03, gtools alone (own process)

## state-leak-01 -- results (rA, rB)
- A1 intact: rc 0, N=480, N_dropped=0.
- A2 y missing on all 120 never-treated rows: r(459), message as claimed, PRECEDED by
  "Caution: estimability checks found time periods exist in sample without a control unit,
   causing them to be dropped."
- A4 never-treated rows PHYSICALLY DELETED, default safety: **identical r(459) and identical
  message**. So the message is not a mis-diagnosis: it is the same message the genuinely
  control-less design gets, and by that point the screened sample genuinely has no control.
- A3 (missing y) + safety(off): rc 0, header prints
  "N clusters = 30; N = 240 (240 obs dropped from initial sample)" -- i.e. HALF the dataset,
  stated in the header. The prober's "NO caution of any kind" is false for the
  dropped-observation disclosure the help promises (sthlp:898/927).
- A5 (rows deleted) + safety(off): rc 0, N=240, N_dropped=120.
- rB: e(wdidmainresults) from A3-equivalent vs A5-equivalent -> MAXABSDIFF = +0.0X-3ff,
  BITIDENTICAL=1. safety(off) on the missing-y data returns EXACTLY the hand-cleaned answer.
  Not a wrong number.
- rB B3: the message's SECOND named remedy -- "modify the sample window manually" --
  WORKS on this data: `wooldid y st yr onset if yr<=2012` -> rc 0, N=240, N_dropped=80,
  att .21855514. The 2013 cohort becomes never-treated. B4 (hand-restricted dataset)
  is identical. So "both remedies the message offers are wrong for this data" is FALSE.
- A7 (zero weight on never-treated) and A8 (i missing on never-treated) -> same r(459):
  the behaviour is the general markout->screen ordering, not outcome-specific.
- A6 partial missingness (30% of control rows): rc 0, N=436, N_dropped=44. Handled.
VERDICT: PARTIAL, severity low. Only the message-wording half survives.

## state-leak-02 -- results (rB part 2)
- PLAIN: c(k) 8 -> 8, no _est_ variable, estimates dir empty.
- VERBOSE: c(k) 8 -> 9, _est____wooldidfullmodel present, estimates dir = ___wooldidfullmodel.
  Observation CONFIRMED.
- `est drop ___wooldidfullmodel` (the documented handle on the retained estimate) removes
  the variable too: c(k) back to 8.
- A subsequent PLAIN run removes it: POST2C_ESTVAR=0, c(k)=8. It does not accumulate.
- Three back-to-back verbose runs: c(k) = 9, 9, 9. No growth.
- **The prober's proposed remedy is destructive.** Drop _est____wooldidfullmodel and then
  `est restore ___wooldidfullmodel`: rc 0, but `count if e(sample)` returns **N = 0**
  (it was 480). Silently emptying e(sample) of the model the help tells users to restore.
  "the variable can be re-created" is false -- Stata's estimates system owns that marker.
VERDICT: PARTIAL, severity low -> at most a one-line help nit. Not a leak.

## state-leak-03 -- pending (rD/rE/rF running)

## state-leak-03 -- SETTLED: gtools owns it (rD/rE/rF, three separate processes)
Each arm starts `mata: mata clear` and plants `mata: k = 12345`.
- wooldid alone: k -> 3; mata globals __00000L __00000M k; scalars __gtools_gfile_ghdfeabs,
  __gtools_gfile_gregvcov.
- reghdfe alone: k -> 12345 (INTACT); no mata globals; no gtools scalars.
- gtools alone (gegen by(st) + gcollapse by(st yr), NO wooldid, NO reghdfe): k -> 2 (CLOBBERED);
  mata globals __000000 __000001 k; BOTH __gtools_gfile_* scalars.
The value left in k equals the by-variable count (2 with by(st yr), 3 inside wooldid). gtools'
by-key machinery owns the name. REFUTED as a wooldid defect; the prober set exactly this
strike condition themselves.

## NEW probes (my own half)
rC_state_catchall: plain run against a session at `set emptycells keep` / `set level 80` /
`set varabbrev off` / `set cformat %9.4f`. Only change in the whole session: the two gtools
scalars. emptycells program-scoping HOLDS (contradiction check on B18 passed). level() restored
on success, on r(459), and on parse-stage r(198); no WOOLDID_LEVEL0 survives. Plain run does NOT
move c(rngstate). Error path leaves only the documented WOOLDID_LASTRC.

rG_microdata_state: 12,000-row microdata (40 individuals per state-year, so the fast collapse
fires). G1PRE vs G1POST and G2PRE vs G2POST: **NO DIFFERENCES AT ALL** in level, emptycells,
c(k), c(N), sortedby, frame, datasignature, estimates dir, globals, frames dir, varlist; row
order preserved (`assert myrow==_n`). Fast collapse vs nofastcollapse: max|diff| over
e(wdidmainresults) = exactly 0 (bit-identical).

rH_micro_state2: xtset survives a wooldid call (panelvar/timevar and _dta[iis]/_dta[tis]
unchanged). A caller's own `preserve` straddling a wooldid call works and restores byte-exactly.
verbose on microdata: the retained marker marks 12000 rows = c(N) = e(N) (micro scale, correct -
`est sto` fires before the collapse preserve), and `est restore` + `count if e(sample)` = 12000.
NOTE: H2/H4/H5 "ROWORDER=0" lines are a HARNESS ARTIFACT of mine (myrow was generated before a
later `bys`/`egen` re-sorted the saved file, and those arms did not re-`sort myrow`); the clean
row-order evidence is rG's G1/G2 and rI's I8/I9, which pin the order first.
`permute(onset)` is correctly refused r(198) with a message naming permute(shuffle)/permute(probabilities).

## rI_esample -- e(sample) is sound on every path I tried (NEGATIVE RESULT)
sthlp:1426 promises e(sample) characterizes the estimation sample. It does. Under the fast
collapse the `ereturn post ... esample()` happens while the working data is collapsed, so this
was the natural place for a scale mismatch. There is none:
  unit-level default: e(N)=300, count if e(sample)=300
  microdata default (collapse fires): e(N)=12000, count if e(sample)=12000
  microdata nofastcollapse: 12000 / 12000
  microdata with 40 rows made missing: e(N)=11960, count=11960, e(N_dropped)=40
  unit-level with if(): e(N)=290, count=290, e(N_dropped)=0
  microdata with if(): e(N)=11600, count=11600, e(N_dropped)=0
  reghdfe reference on the same data: 12000 / 12000
(e(N_dropped)=0 under if() matches the documented definition at sthlp:1411 -- rows excluded by
the if were never part of the "initially specified sample".)
I8 permute(shuffle) and I9 bootstrap(standard): rc 0, caller's datasignature UNCHANGED, row
order pinned (`assert myrow==_n` passes), sortedby preserved, `frames dir` = default only,
no wooldid globals left.
NOTE: rC's log ends at its ARM 7 because of a typo of mine (`qui runiform()` is not a command);
arms 6-9 of that file were rebuilt and run in rG/rH/rI.

## rJ_roworder -- RUNNING: is the answer invariant to the caller's physical row order?

## rJ_roworder -- RESULT: same rows, different physical order (rtol threshold 1e-10)
  J1 unit, shuffled, att+att_i+att_it+att_itime      3.563e-15   CLEAN
  J1 unit, reverse-sorted, same estimands            8.212e-15   CLEAN
  J2 unit, [aw=wt] + es(3,3) + jointtests            3.732e-08   ** ABOVE THRESHOLD **
  J3 microdata, default fast collapse                9.910e-13   CLEAN
  J4 microdata, nofastcollapse                       9.910e-13   CLEAN
  J5 bootstrap(standard) reps(50) rseed(4242)        5.264e-15   CLEAN
  J6 permute(shuffle) permreps(50), same set seed    7.512e-15   CLEAN
  J7 same order, different `set sortseed`            0 (BIT-IDENTICAL)
maxabs in J2 is 1.23e-08 on a cell of magnitude ~0.33, i.e. ~0.3 float-ulp -- the signature of a
single-precision store, the same class as the B21 float-to-double fix. rK_isolate_aw.do bisects
it: att-only vs +es vs +jointtests, aw vs pw vs fw vs unweighted, double vs float weight
variable, and exactly-representable integer weights.

## HARNESS CAVEATS a verifier should know
1. rC_state_catchall.log ends at its ARM 7 because of MY typo (`qui runiform()` is not a
   command, r(199)). Arms 6-9 of that file were rebuilt and run in rG/rH/rI. The log on disk is
   from a second launch of the same (already-fixed) do-file; the sentinels are identical.
2. rA_missctrl.log's `%21x e(wdid_att)` / `e(wdid_att_se)` lines are VOID -- those e() scalar
   names do not exist, so every one prints `+1.0000000000000X+3ff` (Stata's hex for missing).
   The exact-hex comparison is in rB, off `e(wdidmainresults)`, and the displayed tables in rA
   carry the same numbers to 8 significant figures.
3. rH's H2/H4/H5 `ROWORDER=0` lines are MY artifact: `myrow` was generated before a later
   `bys`/`egen` re-sorted the file that was then saved, and those arms did not re-`sort myrow`
   first. The clean row-order evidence is rG (G1/G2) and rI (I8/I9), which pin the order before
   the call; both pass.
4. rH's H8 arm is void: it issues `clear all`, which drops the matrices it then tries to compare.

## rK_isolate_aw -- bisection of the J2 row-order sensitivity
  K1  att only, unweighted                       9.100e-16
  K2  att only, [aw=wt]                          6.915e-16
  K3  att + es(3,3), unweighted                  1.448e-14
  K4  att + es(3,3), [aw=wt]                     3.732e-08   <-- fires
  K5  att + es(3,3) + jointtests, unweighted     1.448e-14
  K6  att + es(3,3) + jointtests, [aw=wt]        3.732e-08   <-- same magnitude
  K7  att only, [pw=wt]                          6.915e-16
  (K8 [fw=wt] aborted the batch: non-integer frequency weights are refused, as they should be;
   K9-K12 did not run. runpair is capture-wrapped in rL.)
It needs BOTH an event-study window AND weights. The moving cells are exactly the `pre3` row's
se / t / p / lb_95ci / ub_95ci -- one horizon's STANDARD ERROR and what is derived from it. The
point estimates do not move above 1e-10.
rL_attribute.do tests whether this is reghdfe's documented iterative-absorption tolerance
(wooldid default regtol(9) = tol(1e-9), sthlp:1117 -- one digit tighter than reghdfe's own 1e-8):
tighten to regtol(12)/regtol(14), loosen to regtol(6), and bypass reghdfe entirely with
unconditionalse / customvce() (which run regress on explicit dummies and take no tolerance).

## rL_attribute -- ATTRIBUTION OF THE ROW-ORDER SENSITIVITY (final)
  L1  [aw=wt] regtol(9) default        3.732e-08   fires
  L2  [aw=wt] regtol(12)               3.732e-08   UNCHANGED -> not the absorption tolerance
  L3  [aw=wt] regtol(14)               3.732e-08   UNCHANGED
  L4  [aw=wt] regtol(6)                4.793e-09
  L5  [aw=wt] unconditionalse          1.407e-13   CLEAN (regress, no HDFE)
  L6  [aw=wt] customvce(cluster st)    1.790e-14   CLEAN (regress, no HDFE)
  L7  [pw=wt]                          3.732e-08   fires
  L8  [aw=iwt] integer weights         1.019e-13   CLEAN
  L9  [aw=one] constant weight         4.870e-15   CLEAN
  L10 unweighted                       1.448e-14   CLEAN
  L11 [aw=wt] robust                   2.496e-08   fires
  L12 [aw=wt] fe(st)                   1.662e-08   fires
Conclusion: HDFE path + non-representable weight values. NOT regtol. ~0.3 float ulp. wooldid's own
source has no `gen float` and builds weights `gen double` (ado:1644), so the store is upstream in
reghdfe on this evidence. Filed as skept-roworder-01, severity LOW.

## rM_floatcheck -- MECHANISM SETTLED
M0  untyped `gegen x = total(y), by(g)` type = DOUBLE (gtools); explicit `gegen double` bit-identical;
    plain `egen` would be float. So ado:1619-1621 are NOT a surviving B21 float store. NEGATIVE.
M1  att+es+[aw=wt], c(type)=float    6.270e-13     (different draw of the same DGP shape as rL!)
M2  same under `set type double`     6.270e-13     -> `set type` is irrelevant
M3  att_it/att_i/att_itime+es+aw     1.681e-13
M4  same under `set type double`     1.681e-13
M5  cell-weighted, unweighted        3.252e-12
M6  cell-weighted, integer aw        1.304e-14
M7  reghdfe ALONE, [aw=wt]   e(V) 8.079e-10 ; e(b) 2.234e-13
M8  reghdfe ALONE, [aw=iwt]  e(V) 0 (EXACT)  ; e(b) 1.241e-13
M9  reghdfe ALONE, unweighted e(V) 0 (EXACT) ; e(b) 2.738e-13
=> reghdfe's WEIGHTED CLUSTERED VCE is row-order dependent at ~1e-9 and exactly order-invariant
   with integer or no weights. wooldid inherits it and, on an ill-conditioned pre-horizon,
   amplifies it to 3.7e-8 (rL) -- but on another draw of the same DGP shape only to 6.3e-13 (M1).
   So the magnitude is dataset-specific; the root cause is upstream and now isolated.
