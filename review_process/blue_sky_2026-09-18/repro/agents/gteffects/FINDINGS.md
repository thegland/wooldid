# FINDINGS — gteffects (agent GTE; method: AGGREGATION IDENTITY)

Pinned build asserted with `which wooldid` on every run
(p01_smoke.log:230 -> .../bluesky/build/wooldid.ado).
Workdir: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects

---------------------------------------------------------------------------
## THE HEADLINE IS A NEGATIVE RESULT, AND IT IS A STRONG ONE
---------------------------------------------------------------------------

I attacked the documented weighted-mean identity
    aggregate == sum(estimate * w_<est>) / sum(w_<est>)
five different ways and could not break it anywhere the help does not already
say it breaks.

1. **Internal identity, 280 checks, 0 failures** (p02_identity.do). 12 core
   configurations x {4 estimands over relyr>=0, the same 4 over relyr<0, 12
   event-study horizons}. Worst relative difference 6.4e-13, on a pre-trend row
   whose value is ~1e-7 (pure cancellation); typical 1e-16.
   Configurations: balanced i=unit and i=cohort; 12% interior missing rows;
   entry+exit attrition; unit-level aweights under both i choices; unbalanced +
   aweights + 3 rows per cohort-period cell; `nofastcollapse`; `timetrends`;
   `esfixedbaseperiod`; `anticipation(1)`; `level(90)`; a biennial time grid.
   WHAT THIS WOULD HAVE CAUGHT: a cell dropped by the name regex (ado:13479),
   a bad raw-matrix join, a wrong relyr encoding, a mismatch between the cells'
   subpopulation and the margins subpopulation, and a weight column read off the
   wrong variable.

2. **The collapse oracle.** gteffects extracts the cells at ado:1753, BEFORE
   `wooldid_fcol` collapses the data at ado:1902; the aggregates are computed by
   `margins` AFTER it. So the identity under `nofastcollapse` vs the default is
   an independent read of the same regression against the collapse. C6 vs C7 on
   unbalanced + aweighted + 3-rows-per-cell data: att agrees BITWISE
   (+1.4b7188f3dabdaX+000 both ways) and every identity check passes both ways.

3. **Exact-truth oracle on the cell VALUES** (p05_truth.do). Zero-noise DGP
   whose effects lie exactly in the span of (unit FE, period FE, cell dummies),
   so each cell coefficient must equal that cell's true (weight-)mean effect.
   MAXABSDEV 1.6e-14 (i=unit, 144 cells), 5.1e-15 (i=cohort), 4.0e-15
   (i=cohort + aweights + 10% gaps), 1.0e-14 (short window, 9 out-of-window post
   periods). `nobs` matched a hand count on every row; `caltime - onset ==
   relyr` on every row of every run in every probe (0 violations, thousands of
   rows).

4. **Independent oracle on the inference columns** (p03_oracle.do). `verbose`
   stores the stage-1 model as `___wooldidfullmodel`; I restored it and
   recomputed every cell's se/p/lb/ub from its own e(b)/e(V)/e(df_r).
   **Relative difference EXACTLY 0 on b, se, p, lb and ub in all 10 runs**:
   iid, robust, cluster(unit) G=80, cluster(state) G=8, level(90), level(97.5),
   fastnose, impvalues+cluster, i=unit with 401 cells. The df actually used is
   the stage-1 one (766 iid; 79 at G=80; 7 at G=8; 7 under impvalues), which is
   what the help promises ("the e(gteffects) cells ... use e(df_r_analytic)").

5. **Option sweep, 32 further arms** (p04_options.do, p09_options2.do). The
   identity fails ONLY where the help says it will: time-varying
   `interactivecontrols()`, `semielasticity`, `poisexpresults`, `impvalues`.
   It HOLDS under fe(), controls(), controls(i.t##c.x), coarsecohortcontrols(),
   cluster(), customvce(), pweights, subgroups(), poisson (xb scale),
   esrelativeto(), oostreatedcontrols(), anticipation(), espostlength(0), a
   short window with out-of-window post rows, jointtests+adversarial, and
   customweightmultiplier().
   A nice confirmation of a documented corner: with `interactivecontrols(x)`
   where x is UNIT-level and i=unit (so the cell slopes are degenerate and the
   help says the run reduces to `controls(i.t#c.x)`), the identity is EXACT
   (max reldiff 7.9e-15, H6) — the approximation really is confined to the
   time-varying case.

---------------------------------------------------------------------------
## FINDINGS
---------------------------------------------------------------------------

### GTE-1 (found by me) — customweightmultiplier() silently turns att_it, att_i
### and att_itime from the documented SIMPLE means into weighted means
severity low-medium | confidence certain | realistic true | surface
customweightmultiplier() x {att_it, att_i, att_itime}
repro: repro_cwm.do

sthlp `att_it`: "The final ATT_it estimate is the simple, unweighted mean of all
of the cohort-period-specific estimates."  sthlp `att_i`: "the simple,
unweighted mean of the cohort-specific treatment effects."
With `customweightmultiplier(v)` neither sentence holds, because the numerator
carries the multiplier (`w * v`, ado:1631/1637/1650) while the normaliser does
not (`W_it`/`W_ipost` = `total(esamp*w)`, ado:1619-1620). Each cell therefore
enters att_it with weight = its own weighted mean of v, not 1.
MEASURED (p04_options.log): per-cell w_att_it totals 1.2420 to 3.2259 and
per-cohort w_att_i totals 1.5440 to 2.4191 with a multiplier that is CONSTANT
WITHIN the cohort-by-period cell (arm F16) — so within-cell variation is not the
cause. Without the multiplier the same design gives exactly 1 and 1.
Nothing at runtime says so, and the cwm help paragraph (sthlp:1107) describes
only the numerator, so the two help paragraphs disagree with each other.
HONEST CAVEAT: cross-cell reweighting is the only thing a multiplier CAN do to
att_it (the cell derivative is constant within the cell), so this may well be
the intended "bespoke treatment effect". Read as a missing qualifier on the
att_it/att_i/att_itime sentences rather than as a wrong number.
IN THE CODE'S FAVOUR: `e(gteffects)` is the only surface on which a user can see
this at all — the identity still holds, and the weight columns show the
distortion. That is the table doing its job.

### GTE-2 (seeded, known_bydesign S7 item 11) — fastnose DOES post
### e(b_es_att)/e(numpre_att); sthlp:1325 is right and the ruling was not built
severity low | confidence certain | realistic true | surface fastnose x e(b_es_*)
repro: repro_fastnose_esatt.do  (measured in p06_dropped.log, G3_HAS_*)

Measured on the pinned build: under `fastnose`, e(b_es_att) present, e(numpre_att)
present, e(V_es_att) ABSENT (and likewise for att_it). That is exactly the
behaviour sthlp:1325 describes when it opens "the covariance is not returned but
the coefficient vector still is" list with {cmd:fastnose}. The 2026-09-13
ruling item 11, which was to REMOVE fastnose from that list, is therefore either
unimplemented or was based on a wrong premise. This resolves the open question
in known_bydesign S7 in the direction OPPOSITE to the one it anticipated ("if it
does not, this is a live doc-vs-code divergence" — it does, so the doc is fine
and the ruling is the loose end).

---------------------------------------------------------------------------
## SEEDED LEADS THAT DID NOT SURVIVE CONTACT
---------------------------------------------------------------------------

- "gteffects under fastnose still carries ANALYTIC standard errors even though
  fastnose promises no standard errors at all" — STALE. The help explicitly
  says gteffects is "Fully available under fastnose (the regression's own VCE is
  intact even though the margins stage skips variances), making gteffects the
  one inference surface that survives fast mode." p03 O6 confirms the SEs are
  bit-exactly the regression's own. Not a divergence.

- "the Mata regex SILENTLY DROPS any cell whose name carries a factor marker" —
  the MECHANISM is real (ado:13479: `^([0-9]+)\._i#([0-9]+)\._t#1\.targetflag$`
  matches nothing with `o.`/`b.`/`bn.`), but I could not reach it on realistic
  input:
    * a treated cell can never be the BASE level: `_i` is 0 for every
      non-target row and `_t` is 0 for the reference period (ado:928-929,
      ado:1009), so level 0 always exists for both factors and is always the
      base, and no treated cell name carries `b.`;
    * an OMITTED cell always carries `o.`, and `wooldid_ccscreen` (ado:11818+)
      scans colnames e(b) for exactly that and raises r(322) under default
      safety, before the table is built (the screen at ado:1713 runs before the
      builder at ado:1753);
    * two realistic collinearity designs — a control sample that stops being
      collected four years before the treated one (p06/p10 L), and a cohort
      treated in the very first observed period, i.e. with no pre-treatment rows
      at all (p10 K) — produce NO omitted cell: wooldid's estimability screen
      drops the control-less periods first with a caution, and the identity then
      holds to 1e-14 (p10, K3_OMITTED_TERMS=0).
  So the silent-drop path requires `safety(off)`; p13_omit.do forces one to
  measure the consequence. (Result recorded in NOTES.md when it lands.)

---------------------------------------------------------------------------
## CLEAN (probed hard, nothing found)
---------------------------------------------------------------------------
- mkmat row cap at ado:1765 (`cap mkmat` failing would leave `gtn`=0 and print
  the false caution "gteffects found no cohort-by-relative-period main-effect
  coefficients to report"): mkmat succeeds at 20,000 rows on this Stata even
  though c(matsize) reads 400, and a 492-cell table builds and reconciles
  exactly (p08_size.do). Dead concern.
- onset varying within `i`, which would make the table's `onset`/`caltime`
  columns come from an arbitrary `gegen tag()` row: REFUSED at ado:782.
- `nobs` column: matched a hand `count if e(sample) & i==.. & year==..` on every
  row of every truth run (NOBSMISMATCH 0).
- Weight-definition checks the identity cannot see: per post cell w_att_it is
  exactly 1; per cohort sum(w_att_i) is exactly 1; per cohort sum(w_att_itime)
  is exactly 1 and its cells are equal within the cohort (1/n_i) — under
  unbalanced panels, aweights, pweights, multi-row cells and gaps.

---------------------------------------------------------------------------
## MORE FINDINGS (measured after the first draft)
---------------------------------------------------------------------------

### GTE-3 (found by me) — a treated `i` that enters an unbalanced panel at or
### after its own onset is dropped SILENTLY, and the help describes the rule at
### the wrong level
severity medium (magnitude argues for high) | confidence certain |
realistic true | surface sample construction / which cohorts reach e(gteffects)
repro: repro_lateentry_silentdrop.do   (measurement: p15_lateentry.do)

ado:12018-12020 computes `mint' = `gegen min(t), by(i)' — the first period each
`i' is OBSERVED, not the first period of the sample — and ado:12055-12059 then
does
      gegen byte nopre = max(ttre <= mint) if ifstatement, by(i)
      replace ifvar = 0 if nopre == 1
with no display line of any kind. sthlp:73 documents this as "cohorts treated
before the sample begins are dropped". In an unbalanced panel that is not what
happens: a firm treated in 2004 that only enters the register in 2004 or 2005,
in a 2000-2010 panel, is dropped although the sample begins in 2000 and the firm
has seven post-treatment years. It never appears in e(gteffects), so the
cohort-by-period table gives no hint either, and the reported att is then the
att over the self-selected subset of treated units that happen to have
pre-treatment coverage.
MEASURED (p15_lateentry.log). One 70-firm register, 2000-2010; 30 firms treated
in 2004; ONLY the first-observed year of ten of them is varied:

  entrants first seen  e(N)  N_dropped  treated firms in e(gteffects)   att
       2000             770       0                30                 2.5372
       2003             740       0                30                 2.5775
       2004             660      70                20                 1.7223   <-
       2005             660      60                20                 1.7401   <-
  same data, i = COHORT instead of i = unit, entrants first seen 2004:
                        730       0                (1 treated cohort)  2.4337

So the headline att moves 33% (2.5775 -> 1.7223) purely because ten treated firms
were silently removed, and 41% (1.7223 -> 2.4337) between two equally legitimate
choices of `i` on IDENTICAL data. The only trace in the printed output is the
generic "(70 obs dropped from initial sample)" count, which is shared with
`if'-drops and markout drops; grepping the whole 2004 block for "Caution" finds
none. (The 2003 arm, where nothing is dropped, DOES print the unrelated
"reference period for some units forced" caution, so cautions are not suppressed
in this run generally.)
THE INCONSISTENCY IS INTERNAL, NOT JUST DOC-VS-CODE: the identical drop
performed twelve lines earlier for `anticipation()' (ado:12045) prints
  "Caution: anticipation(#) places the effective onset of N unit(s) at or before
   their first usable pre-period; those unit(s) and their M observations are
   dropped from the sample."
so the package already knows this drop is worth naming. Under `timetrends' the
rule is stricter (two pre-periods needed) and IS documented at sthlp:1055
("drops cohorts without sufficient pre-treatment data to estimate a time
trend"), though still silent at runtime.

### GTE-4 (found by me) — weight columns are populated for estimands that were
### not requested
severity low | confidence certain | realistic true | surface e(gteffects) columns
repro: p12_sg2.do section S2

sthlp for gteffects: the weight columns are "missing for estimands not
requested". Measured (p12_sg2.log):
   request att_i      -> w_att_it populated on all 372 rows (att_it not requested)
   request att_itime  -> w_att_it populated on all 372 rows (att_it not requested)
   request att_it     -> only w_att_it
   request att        -> only w_att
Cause: ado:1643 generates `w_att_it' whenever ANY of att_it/att_itime/att_i is
requested, and the builder at ado:1758-1762 populates a column for every weight
variable that exists. The values are correct att_it weights, so this is benign —
a user who takes the sentence literally would just be surprised, and one who
uses the column gets the right answer. Fix is a one-line help edit or a gate.

---------------------------------------------------------------------------
## A POSITIVE RESULT WORTH DOCUMENTING
---------------------------------------------------------------------------
The single full-sample table returned under `subgroups()` IS enough to rebuild
every subgroup estimand exactly, because the subgroup variable is required to be
constant within the (i,t) cell, so the cells partition cleanly (p12_sg2.log):
  main_SG0_att    reported +1.1fe4366047e22X+000  rebuilt +1.1fe4366047e21X+000
  main_SG0_att_i  reported +1.21798d7b376b8X+000  rebuilt +1.21798d7b376b9X+000
  main_SG1_att    reported +1.75a3730595ff4X+000  rebuilt +1.75a3730595ff5X+000
  main_SG1_att_i  reported +1.75b2e1978398eX+000  rebuilt +1.75b2e1978398fX+000
(one ulp each). The help currently presents "a single full-sample table is
returned" as a limitation; it could say instead that restricting the table to a
subgroup's cells reproduces that subgroup's aggregates.

### GTE-5 (found by me) — "under timetrends ... the se column is missing" is not
### what happens on a cohort-level i
severity low | confidence certain | realistic true | surface gteffects x timetrends
repro: p14_messy.do section N8 (self-contained apart from the generator include;
the numbers are printed as SE|N8_timetrends|...)

The gteffects help paragraph says "under {it:timetrends} the coefficient-level
variances are typically not identified and the {it:se} column is missing
(degenerate zero variances are always reported as missing, never zero)".
MEASURED on a 90-unit, 10-year panel with i = cohort (3 cohorts) and the default
iid VCE: all 14 cells carry a standard error; 0 missing, 0 exactly zero, 0
missing p-values. Values run .475 to 3.04 and grow with the horizon, which is
what a cohort-trend model should produce.
The "typically" hedge partly covers this, and the known degenerate case is
`timetrends` with `cluster()` AT the level of i (documented at sthlp:1058, and in
the project's own record). But on the cohort-level i that the package itself
recommends (the >100-distinct-i caution steers users there), and with the default
VCE, the sentence does not describe the behaviour. No wrong number is produced —
p03's independent oracle shows the ses are exactly the regression's own — so this
is a help-text scope issue, not an inference defect.

---------------------------------------------------------------------------
## MORE CLEAN
---------------------------------------------------------------------------
- p14_messy.do: an `if` whose expression is MISSING on 200 rows; a singleton
  treated cohort (one unit); a cohort treated in the LAST observed period;
  cohorts treated after the sample end, with and without oostreatedcontrols(4);
  a panel with 45% of rows deleted at random; i=unit + aweights + entry and exit
  attrition + cluster(state). Every identity check passes, JOINMISS 0, the
  caltime-onset-relyr identity holds, and inwindow matches e(b) in each.
  (fweights are refused by design.)
- p12 S3: under `bootstrap(standard) reps(25) cluster(unit)` the point estimates
  still rebuild every aggregate exactly (372 cells, worst 3.8e-15) — the
  bootstrap rewiring touches only se/p/lb/ub, as documented.
- p12 S4: `unconditionalse robust` — identity exact on att.
- p13 M6: under `permute(shuffle)` e(perm_gt_blank)=1, the run PRINTS "The
  gteffects cells are keyed by cohort onset, which the replications reassign, so
  they have no permutation distribution: their inference columns are blank", and
  the surviving point estimates still rebuild every aggregate (250 cells, worst
  6.3e-15). The seeded "check the user is told" lead is clean.

### GTE-3 — the self-contained repro's own output (repro_lateentry_silentdrop.log)
  entrants first observed 2000: e(N)=770  N_dropped=0   firms in table 30  att 2.598970
  entrants first observed 2003: e(N)=740  N_dropped=0   firms in table 30  att 2.532414
  entrants first observed 2004: e(N)=660  N_dropped=70  firms in table 20  att 1.818684
  entrants first observed 2005: e(N)=660  N_dropped=60  firms in table 20  att 1.715990
Nothing but the observation count distinguishes the last two runs from the first two.

---------------------------------------------------------------------------
## REPRO INVENTORY (all four have been RUN; their .log sits beside them)
---------------------------------------------------------------------------
  repro_cwm.do                  -> repro_cwm.log                  (GTE-1)
  repro_fastnose_esatt.do       -> repro_fastnose_esatt.log       (GTE-2)
  repro_lateentry_silentdrop.do -> repro_lateentry_silentdrop.log (GTE-3)
  repro_weightcols.do           -> repro_weightcols.log           (GTE-4)
      requested [att]     : w_att=180  w_att_it=0    w_att_i=0    w_att_itime=0
      requested [att_it]  : w_att=0    w_att_it=180  w_att_i=0    w_att_itime=0
      requested [att_i]   : w_att=0    w_att_it=180  w_att_i=180  w_att_itime=0
  p14_messy.do section N8       -> p14_messy.log                  (GTE-5)
Each of the four repro_*.do files is standalone: it sets its own adopath to the
pinned build, asserts it with `which wooldid`, generates its own data and
includes nothing.
