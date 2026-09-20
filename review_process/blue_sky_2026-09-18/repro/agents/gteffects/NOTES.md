# gteffects probe — NOTES (agent: GTE, method: aggregation identity)

## Harness
- Build asserted: `which wooldid` -> .../bluesky/build/wooldid.ado  (p01_smoke.log:230)
- Helpers: `gtlib.do` (mkpanel generator + gtcheck/gtdefcheck wrappers),
  `gtmata.do` (compiles gtident/gtidentpre/gtidentes to .mo in the workdir;
  workdir is on adopath so the .mo files resolve).
- e(b) estimand names are `main_att`, `main_att_it`, `main_att_i`, `main_att_itime`,
  `pretreatment_*`, `ES_<est>_post<k>` / `ES_<est>_pre<k>`.
- `i` must be POSITIVE nonzero integers (r(411) otherwise) -> generator emits `cohp = coh+1`.

## Code facts established by reading (pinned build)
- gteffects builder: ado:1749-1773 (runs BEFORE wooldid_fcol at ado:1902, i.e. on the
  UNCOLLAPSED data, reading the stage-1 e(b)/e(V) directly).
- Mata worker: `wooldid_gtx`, ado:13466-13537. regex
  `^([0-9]+)\.<_i>#([0-9]+)\.<_t>#1\.<targetflag>$`; any factor marker (o./b./bn.) fails it.
- cell weight cols = `gegen total(<w_*>) if e(sample)&targetflag==1, by(_i _t)`.
- `_i` IS the user's `i` value (ado:928), so the table's `i` column is the user value.
- `w_att = w*cwm` (ado:1650); `w_att_it = w*cwm/w_per_it` with
  `w_per_it = total(esamp*w), by(i t)` — note the DIVISOR carries NO cwm (ado:1619,1645).
- main att margins subpop: ado:2030 `ifstatement & targetflag==1 & t>=ttre`.
- onset non-constant within `i` is REFUSED at ado:782 (so the "arbitrary tagged row"
  hazard for the onset/caltime columns is unreachable) — NEGATIVE RESULT.
- default-safety dropped-coefficient screen: `wooldid_ccscreen` ado:11818+; fires only when
  the coefficient NAME carries `o.`; errors r(322).

## Runs
- p01_smoke.do  — build assertion + colname discovery. OK.
- p02_identity.do — 12 configs x {4 estimands post, 4 pre, 12 ES horizons}.
  **280 identity checks, 0 FAIL, worst reldiff 6.4e-13** (that one under nofastcollapse,
  ES_att_it_pre3, 3.3e-14; the 6.4e-13 is a near-zero pre-trend row).
  Configs: balanced i=unit / i=cohort; 12% interior gaps; entry+exit attrition;
  aweights at unit level (i=unit and i=cohort); unbalanced+aw+3 rows per cell;
  nofastcollapse (collapse oracle — agrees with the collapsed run bitwise on att);
  timetrends; esfixedbaseperiod; anticipation(1); level(90); biennial grid.
  Definitional checks all clean: w_att_it == 1 per post cell; sum w_att_i == 1 per cohort;
  sum w_att_itime == 1 per cohort.
  => The core identity is SOLID. See FINDINGS.md for what is left.

## Next
p03: (a) `verbose` + `estimates restore ___wooldidfullmodel` as an INDEPENDENT oracle for
the est/se/p/ci columns under iid/robust/cluster/level(90)/fastnose;
(b) option identity sweep (subgroups, fe, controls, interactivecontrols, poisson,
impvalues, unconditionalse, cluster, pw, esrelativeto, oostreatedcontrols);
(c) customweightmultiplier varying WITHIN a cell (predicted: att_it cell weights != 1);
(d) out-of-window post rows (espostlength(1) with 6 realized post periods);
(e) safety(off) + an omitted treated cell;
(f) fastnose e(b_es_att)/e(numpre_att) (known_bydesign S7 item 11) and permute blanking.

## Results so far (2026-09-17)

| probe | what it tests | verdict |
|---|---|---|
| p02_identity | 12 core configs x 4 estimands x pre/post/ES = 280 identity checks | **0 FAIL**, worst 6.4e-13 |
| p04_options  | 20 option arms | fails ONLY under interactivecontrols (time-varying), semielasticity, impvalues — all documented |
| p09_options2 | 12 more arms incl. poisson, ccc, customvce, subgroups(i=unit), anticipation+oos | clean; poisexpresults fails (documented) |
| p05_truth    | EXACT-truth oracle: zero-noise DGP, cell estimate vs the cell's true mean effect | **MAXABSDEV 1.6e-14 / 5.1e-15 / 4.0e-15 / 1.0e-14**; NOBSMISMATCH 0; caltime-onset==relyr on every row |
| p03_oracle   | se/p/lb/ub recomputed independently from the restored stage-1 e(b)/e(V)/e(df_r) | **reldif EXACTLY 0** on b, se, p, lb, ub in all 10 runs (iid / robust / cluster(unit) G=80 / cluster(st) G=8 / level(90) / level(97.5) / fastnose / impvalues / i=unit 401 cells) |
| p08_size     | mkmat row cap -> "found no coefficients" caution | mkmat OK at 20,000 rows; 492-cell table built and exact. CLEAN |
| p06/p10      | control sample ends early; cohort treated in period 1 | wooldid's estimability screen drops the control-less periods first; NO omitted cell arises; identity holds |
| p07 (aborted) | subgroup rebuild | my bug: `e(GTEFFECTS)` uppercase -> re-run as p12 |

### Decoded df_r from p03 (all match the stage-1 regression exactly)
- iid: 766 ; cluster(unit) G=80: 79 ; cluster(st) G=8: 7 ; impvalues+cluster(st): 7
  (so the cells DO use the analytic stage-1 df under impvalues, as sthlp says)

### Seeded leads settled
- "fastnose still carries analytic SEs even though fastnose promises none" — STALE. The
  help explicitly says gteffects is "Fully available under fastnose (the regression's own
  VCE is intact...)". p03 O6 confirms the SEs are exactly the regression's own.
- "gteffects under permute() is blanked — check the user is told" — display note exists at
  ado:11726 and is existence-gated on e(perm_gt_blank); to be confirmed in p13 M6.
- "the regex silently drops cells carrying a factor marker" — MECHANISM CONFIRMED by
  reading (ado:13479), but I could NOT reach it on realistic input: level 0 exists for
  both _i and _t by construction so no treated cell is ever the base, and an omitted cell
  always carries `o.` which wooldid_ccscreen catches with r(322) under default safety.
  p13 forces one to measure what happens under safety(off).

### §7 item 11 (known_bydesign) — ANSWERED
`fastnose` DOES post e(b_es_att) and e(numpre_att) and does NOT post e(V_es_att)
(p06_dropped.log, G3_HAS_*). So sthlp:1325 listing fastnose under "the covariance is not
returned but the coefficient vector still is" is ACCURATE, and the 2026-09-13 ruling
item 11 ("remove fastnose from that list") was not built — or rested on a wrong premise.

## p11_window (inwindow column) — CLEAN
`gtwincheck`'s own e(b) parse was broken (Stata expands `=regexs(1)` before the `if`
runs), so read the manual comparison instead: for every arm the set of relyr values with
inwindow==1 equals exactly the set of ES_att_pre*/post* horizons in e(b).
  W1 default espre2/espost2      -> {-2,-1,0,1}                 matches pre2,pre1,post0,post1
  W2 esfixedbaseperiod           -> {-3,-2,0,1}                 (base at -1 is excluded)
  W3 esf + esrelativeto(-2)      -> {-4,-3,-1,0,1,2}            (base at -2 excluded)
  W4 espostlength(0)             -> {-2,-1}, no post row in window  (help says exactly this)
  W5 esprelength(0)              -> {0,1,2}
  W6 espre6/espost8 on a 12-year panel -> {-4..-1, 0..7}; the missing -5/-6 are announced by
     "Caution: reference period for some units forced (for want of panel length)..."
  W7/W8 anticipation(2)          -> matches
  W9 biennial (tstep 2), espre4/espost4 -> {-4,-2,0,2}; the odd slots are empty, which is the
     documented "contiguous integer slots" convention (BS-36)

## repro_cwm.log — GTE-1 confirmed on a clean self-contained case
  (A) no multiplier : per-cell w_att_it min=max=1 ; per-cohort w_att_i min=max=1
  (B) with customweightmultiplier(pop), pop CONSTANT within the cohort-by-year cell:
      per-cell w_att_it 1.43607 .. 4.92387 (3.43x spread) ; per-cohort w_att_i 2.62936 .. 4.75544
      and NO caution is printed.

## p12_sg2 — subgroup rebuild WORKS (positive result worth telling the author)
With subgroups(sg) and i=unit, every subgroup estimand in e(b) is reproduced from the
SINGLE full-sample e(gteffects) table by restricting to that subgroup's cells:
  main_SG0_att      reported +1.1fe4366047e22X+000   rebuilt +1.1fe4366047e21X+000 (1 ulp)
  main_SG0_att_i    reported +1.21798d7b376b8X+000   rebuilt +1.21798d7b376b9X+000 (1 ulp)
  main_SG1_att      reported +1.75a3730595ff4X+000   rebuilt +1.75a3730595ff5X+000 (1 ulp)
  main_SG1_att_i    reported +1.75b2e1978398eX+000   rebuilt +1.75b2e1978398fX+000 (1 ulp)
(att_it == att and att_itime == att_i here because i is the unit.)
So "with subgroups() a single full-sample table is returned" is not a limitation: since
the subgroup variable is REQUIRED to be constant within the (i,t) cell, the cells
partition cleanly and the subgroup aggregates rebuild exactly. Worth documenting.

## p12 S2 — weight columns for estimands that were NOT requested
  request att_i      -> w_att_it POPULATED (372/372), w_att_i populated, others missing
  request att_itime  -> w_att_it POPULATED (372/372), w_att_itime populated, others missing
  request att_it     -> only w_att_it
  (ado:1643 generates w_att_it whenever any of att_it/att_itime/att_i is asked for.)
sthlp says the weight columns are "missing for estimands not requested". Benign (the
values are the correct att_it weights) but the sentence is not accurate.

## THE LATE-ENTRY DROP (ado:12018-12059) — mechanism read, measurement in p15
  tempvar mint ; gegen `mint' = min(`t') if `ifstatement', by(`i')
  gegen byte `nopre' = max(`ttre' <= `mint') if `ifstatement', by(`i')
  replace `ifvar' = 0 if `nopre' == 1        <- NO MESSAGE
`mint' is per-`i', so this drops any treated `i' whose OWN first observed period is at or
after its onset - not only "cohorts treated before the sample begins" (sthlp:73). The
anticipation() version of the identical drop twelve lines above (ado:12045) DOES print a
caution naming the unit and observation counts. Under timetrends the rule is stricter
(needs 2 pre-periods) and IS documented ("drops cohorts without sufficient pre-treatment
data", sthlp:1055) though still silent at runtime.

## p13_omit — I could NOT reach an omitted treated cell on a completing run
Three attempts, all blocked by a DIFFERENT wooldid mechanism:
  M1/M2  fe(i.st#i.year) with the treated cohort entirely inside one state,
         default safety AND safety(off): both r(459) with the named BS-33 error
         "the fitted model carries no estimable treatment-effect coefficient at any of
          the requested event-study relative periods ... the cohort-by-period treatment
          terms were absorbed or dropped ... check absorbed fixed effects".
  M3/M4  controls(z) where z is the EXACT indicator of one treated cohort-by-period
         cell: the run completes and the omitted term is `o.z' - Stata drops the
         CONTROL, not the cell dummy, because the cell interactions come earlier in the
         varlist (ado:1476). The table keeps all 7 cells and the identity holds to
         1.4e-16; M5 (no z at all) returns the identical numbers.
  p06/p10 realistic collinearity shapes: the estimability screen removes the offending
         periods/cohorts before the regression, so no `o.' cell is ever produced.
=> the "silently dropped cell" the seeded lead describes is REAL AS A MECHANISM
   (ado:13479) but I found no route to it that ends in a completed run.

## p13 M6 — permute(): the user IS told (seeded lead CLEAN)
e(perm_gt_blank)=1, and the run prints
  "The gteffects cells are keyed by cohort onset, which the replications reassign, so
   they have no permutation distribution: their inference columns are blank."
The POINT estimates survive and still rebuild every aggregate exactly (250 rows, all
identity checks ok, worst 6.3e-15). Nothing hidden.

## p14_messy — all clean, plus one help-sentence miss
N1 (`if' whose expression is missing on 200 rows), N3 (a singleton treated cohort),
N4 (a cohort treated in the LAST observed period), N5/N5b (cohorts treated after the
sample end, with and without oostreatedcontrols(4)), N6 (45% of rows missing at random),
N7 (i=unit + aweights + entry/exit attrition + cluster(state)): every identity check ok,
JOINMISS 0, CALTIME check 0, inwindow consistent. N2: fweights are refused by design
("wooldid does not accept fweights. Specify aweights or pweights.").

**N8**: sthlp (gteffects paragraph) says "under {it:timetrends} the coefficient-level
variances are typically not identified and the {it:se} column is missing".
MEASURED on a 90-unit panel with i = cohort and the default iid VCE:
  SE|N8_timetrends|rows=14|se_missing=0|se_exactly_zero=0|p_missing=0
i.e. every cell has a standard error. (N9 plain and N10 fastnose likewise 0 missing.)
The "typically" hedge partly covers it, and the known all-zero-covariance case is
timetrends + cluster() AT i — but on the cohort-level i the package itself recommends
(SI-8), the sentence as written does not describe what happens.

## Verifier caveats (read these before judging any log)
1. `p07_subgroup.log` ABORTS with r(111) at `matrix GT = e(GTEFFECTS)` — my own
   case bug, not a wooldid defect. The corrected re-run is `p12_sg2.do`.
2. `gtwincheck` (in gtlib.do) prints a garbage `ebhoriz=[. . 1]` list: Stata expands
   `` `=regexs(1)' `` before the `if regexm(...)` on the same line runs. The inwindow
   verdicts in NOTES were taken from a MANUAL comparison against the ES_att_* names in
   e(b), not from that line.
3. `p05_truth.log` shows PRE_* "FAIL" lines with reldiff ~1e-5. Both sides there are
   numerical zeros (~1e-10 on a zero-noise DGP where the true pre-effect is 0), so the
   relative metric is meaningless; the absolute difference is ~1e-15. The exact-truth
   verdict is the TR| line (MAXABSDEV), which is 1e-14 or better in all four arms.
4. p04's F5/F6 and p03's O9 poisson arms failed r(411) (my generator's y can be
   negative); the poisson arms that matter were re-run as p09 H2/H3 on exp(y/4).

## p15_lateentry — GTE-3 CONFIRMED, and the magnitude is large
One 70-firm register 2000-2010, 30 firms treated in 2004, only the first-observed
year of ten of them varied (i = unit):
  entrants seen from   e(N)  N_dropped  treated firms in e(gteffects)   att
      2000              770      0            30                      2.5372
      2003              740      0            30                      2.5775
      2004              660     70            20                      1.7223
      2005              660     60            20                      1.7401
  same data with i = COHORT, entrants from 2004:  730, 0 dropped,      2.4337
No caution anywhere in the 2004/2005 blocks (grep 'Caution' over the block: the only
hit in the 2003 block is the unrelated reference-period-forcing one, so cautions are
not suppressed in this run).
=> 33% swing from the silent drop; 41% swing between two legitimate choices of `i`
   on IDENTICAL data, with nothing in the output to attribute it.

## UNTESTED / not reached
- `contreat()` mode: gteffects is refused there by design (ado:4695) — not probed.
- A cell count above ~11,000 pushed through wooldid itself. Standalone `mkmat`
  succeeded at 20,000 rows, so the `cap mkmat` failure path at ado:1765 (which would
  print the FALSE caution "gteffects found no cohort-by-relative-period main-effect
  coefficients to report") is almost certainly unreachable, but I never ran a
  12,000-cell wooldid to prove it.
- `boottstat`, `bootdisplay(bc)` and `bootdisplay(bca)` cell families (only the
  default percentile family was exercised, under bootstrap(standard) reps(25)).
- Under `poisson poisexpresults` the identity fails as documented; I did not
  separately verify that the CELLS themselves are on the xb scale.
- `permstepdown` / `simultaneousbands` companions of the cell table.
- The esfixedbaseperiod reference-period screen at ado:855-860 drops units that do
  not OBSERVE the reference period, also with no message except when the anticipation
  attribution at ado:860 fires. Same class as GTE-3; not separately measured.
