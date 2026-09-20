# FINDINGS -- agent FP-ESTIMANDS (first-principles oracle / independent reimplementation)

## FP-SEMIELAS-FE -- `semielasticity` divides by a prediction that excludes the absorbed fixed effects  [HIGH, certain, NOT seeded]

Surface: `semielasticity` on the default (reghdfe) engine. By extension
`contreatelasticitytype(eydx)`/`(eyex)`, and every output reusing the same margins call
(event-study horizons, histogramcohorteffects, makecfxplots, subgroup rows).

CLAIM. wooldid's default engine is `reghdfe` with `absorb(i.i i.t ...)`, so the unit and
time fixed effects are not in e(b) and `predict, xb` returns `_cons + Xb` with the fixed
effects omitted. `margins , eydx(targetflag)` averages ln(yhat|tf=1) - ln(yhat|tf=0) over
the subpopulation, so on that engine the semi-elasticity measures the effect against the
sample GRAND MEAN, not against each observation's own predicted outcome. The identical
model with explicit `i.i i.t` dummies -- what customvce()/unconditionalse do (ado:1484-5)
-- divides by the true fitted value and returns a different number. No caution. The help
(sthlp:243-245) promises "the average percentage change ... between the two sets of
predicted outcomes".

REPRO: repro_semielas_fe.do  (in this directory; self-contained)
  "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_semielas_fe.do

OBSERVED (repro_semielas_fe.log:267-268, 413-414, 472-505):
  ARM1 (30 units, y ~ 10..230, constant level effect 5, no calendar trend)
    LEVEL   reghdfe= 5.198696458  explicit= 5.198696458  ratio= 1.000000
    SEMIEL  reghdfe= 0.048760539  explicit= 0.080343215  ratio= 1.647710   <-- 65% apart
  ARM2 (CONTROL: every unit at the same level; nothing material absorbed)
    LEVEL   reghdfe= 4.876573741  explicit= 4.876573741  ratio= 1.000000
    SEMIEL  reghdfe= 0.039676977  explicit= 0.039776647  ratio= 1.002512   <-- agree
  mechanism, by hand on the same fit:
    CORR(reghdfe xb , true fitted value) = -0.001479
    SD(reghdfe xb | post)      =    1.77908
    SD(true fitted value|post) =   56.88270
    HAND mean[ln(xb_hdfe)-ln(xb_hdfe-b)] =  0.048760539  <- reproduces the DEFAULT answer
    HAND mean[ln(yhat)   -ln(yhat   -b)] =  0.080343215  <- reproduces the EXPLICIT answer

A second independent design (dataE.dta, 24 units, esprelength(3) espostlength(4)) splits the
same way: main_att = 0.036668142715654 on the default engine vs 0.051863777168292 under
customvce(robust), factor 1.414, and every event-study horizon splits too (resE1.csv vs
resE2.csv).  run4.log:144-149,199,227 shows margins' own eydx return equals the hand-computed
log-difference on each engine to 1e-13, so the operator is correct; its input prediction is not.

EXPECTED. The level (dydx) answers agree bit-for-bit across engines (ratio 1.000000) because
the fixed effects cancel in a difference. eydx divides by the prediction, so they do not
cancel -- but the two engines estimate the same model and must still agree. The
explicit-dummy arm is correct: its denominator is the fitted outcome. The author already
knew margins needs the whole prediction: ado:1459 says "ppmlhdfe requires d for predict,
which margins requires" and passes `d` to ppmlhdfe (ado:1462,1468) but NOT to reghdfe
(ado:1476). Under plain dydx that omission is harmless, which is why it survived.

REALISTIC: yes. A semi-elasticity is the natural summary exactly for a strictly positive,
widely dispersed outcome (revenue, spending, caseload, earnings) -- the shape on which the
error is largest -- and the default engine is what almost everyone runs.
SEVERITY: high.  Not in known_bydesign.md; sthlp:243-249 discusses only the poisson remap
and the cell-collapse ineligibility.

---

## FP-CWM-DOC -- the att_it/att_i/att_itime weight formulas in the help omit customweightmultiplier()  [LOW, certain, SEEDED -- and the seeded worry is NOT confirmed]

The seeded lead was that cwm enters the numerator of all four estimand weights and the
denominator of none, so a cell-varying cwm would silently stop att_it/att_i/att_itime being
"the simple, unweighted mean" of sthlp:220-221.  SETTLED: the code is coherent and the
option's own help entry authorises it.  sthlp:1107: "the weights normally used to produce
the requested type of treatment effect (ATT, ATT_i, etc.) will be multiplied by this
variable" -- exactly w_item*cwm.  In a one-row-per-cell panel w_att_it = cwm, so att_it
becomes the cwm-weighted mean of cell effects, which is the only placement under which cwm
does anything at all for that estimand.  My oracle reproduces all four estimands and every
horizon under cwm to 2.7e-13 on both a unit-level and a multi-row-per-cell design (cmpA.py
A4/A6, cmpBC.py B4).

What remains is a doc gap: the three phangs at sthlp:220-236 spell the weights out as
w/W_it, w/W_ipre|post and w/(W_it*T_iprepost) -- no cwm factor -- and keep calling the
result "the simple, unweighted mean", which under cwm it is not.
REPRO: run1.do (A4/A6) + cmpA.py.  SEVERITY low.

---

## FP-NOPRE-LATEENTRY -- the "no pre-period" drop is per-unit, not per-sample, and the help says otherwise  [LOW, certain, NOT seeded]

Surface: pooled (default) reference regime, wooldid_reltime ado:12047-12055.

CORRECTED after checking: wooldid DOES report a count -- the header prints
"N = 240 (36 obs dropped from initial sample)" -- and sthlp:73 DOES document a drop:
"cohorts treated before the sample begins are dropped".  What survives is narrower.

CLAIM. The screen is `ttre <= mint` where `mint` is THAT UNIT's own first in-sample
period, not the first period of the sample.  A unit that enters the panel late and whose
onset falls comfortably INSIDE the sample window is therefore dropped too, which sthlp:73's
"treated before the sample begins" does not describe.  No caution names the cause -- while
the anticipation-driven version of this same screen both prints one and is documented as
printing one (sthlp:824: "units whose effective onset falls at or before their first usable
period are dropped (with a caution reporting the count)").

REPRO: repro_nopre_silent.do; confirmation in run11.log:230-246.
OBSERVED. Design C spans t = 1..8.  Six units have ttre = 3 -- inside the window -- and are
first observed at t = 3.  They are dropped (e(N) 240 of 276 rows, e(N_dropped) = 36,
run3.log:1755-1761).  The header gives the count with no cause; the only printed caution is
about a different rule ("reference period for some units forced ...").  The same run under
anticipation(3) DOES print "Caution: anticipation(3) places the effective onset of 9 unit(s)
at or before their first usable pre-period; those unit(s) and their N observations are
dropped from the sample" (repro_nopre_silent.log:513).
EXPECTED. Either sthlp:73 should say "at or before their first observed period" rather than
"before the sample begins", or the unconditional screen should print the caution its
anticipation twin already prints.
SEVERITY low (doc-vs-behaviour; the numbers are right for the sample used).  REALISTIC yes:
late panel entry is ordinary in administrative data, and such a unit is silently excluded
from an estimand the user believes covers it.

---

## NEGATIVE RESULTS
The core estimand machinery is CLEAN.  An independent reimplementation (numpy OLS on an
explicitly built design matrix, weights rebuilt from raw data, margins aggregations redone by
hand) reproduces att, att_it, att_i, att_itime, both pretreatment rows and every event-study
horizon to <= 2.7e-13 across 10 option regimes on 3 designs, including a multi-row-per-cell
administrative shape where all four estimands genuinely differ.  See NOTES.md.

---

## FP-SEMIELAS-FE — additional evidence (reach and a third engine arm)

THIRD ENGINE ARM. `unconditionalse robust` (also explicit dummies, ado:1484-5) returns
0.080343215 on ARM1 and 0.039776647 on ARM2 -- identical to `customvce(robust)` to every
printed digit (repro_semielas_fe.log:309,495). Two independent explicit-dummy engines agree
with each other and with the hand computation; only the default reghdfe engine differs.

REACH. Every reported row inherits it, not just `main_att` (run5.log, resH1.csv vs resH2.csv,
`att semielasticity subgroups(sg) esprelength(2) espostlength(3)` on dataE):

```
row                     reghdfe(default)      explicit dummies      ratio
main_att                +3.674580812270e-02   +5.170423388990e-02   1.4071
main_SG0_att            +3.374902343623e-02   +5.332147792872e-02   1.5799
main_SG1_att            +3.974259280916e-02   +5.008698985108e-02   1.2603
ES_att_post0            +3.888884323085e-02   +5.483201515028e-02   1.4100
ES_att_post1            +3.340497179932e-02   +4.637639897929e-02   1.3883
ES_att_post2            +3.918408359892e-02   +4.960809063008e-02   1.2660
ES_att_pre1             -3.204115488237e-03   -9.268265152963e-03   2.8926
ES_att_pre2             -5.304291784924e-03   +1.683047578474e-03  -0.3173   <-- SIGN FLIP
ES_att_SG0_pre2         -5.345136705720e-03   +8.249688177697e-04  -0.1543   <-- SIGN FLIP
ES_att_SG1_pre2         -5.263446864128e-03   +2.541126339177e-03  -0.4828   <-- SIGN FLIP
pretreatment_att        -4.254203636581e-03   -3.792608787244e-03   0.8915
```
Three PRE-TREATMENT event-study coefficients change sign with the engine. Under
`semielasticity` the reported pre-trend is therefore engine-dependent in direction, not
only in magnitude.

NO WARNING. The only message a default `semielasticity` run prints is the cell-collapse
caution, whose text is "the run is slower, and **no answer changes**" (run7.log L1 block).

A CORRECT PATH EXISTS but is not the default: `customvce(robust)` works for the
cell-weighted estimands too (`att_it semielasticity customvce(robust)` returns
0.05494737436613 against the default's 0.03980378341365, run7.log:213-276), while
`unconditionalse` is refused for them by design.

---

## FP-CTZERO-ENGINE — `contreatcontrols(untreatedZero)`: the two engines return opposite signs, and only one warns  [MEDIUM, UNCERTAIN, NOT seeded — reported as a lead]

**Surface:** `contreat()` + `contreatcontrols(untreatedZero)`. This is the AME family, which
is not my assigned target; I found it while testing whether FP-SEMIELAS-FE extends to the
continuous-treatment elasticity forms. I could NOT determine which side is right.

**Claim.** On a low-noise, well-identified design the default (reghdfe) engine and the
explicit-dummy engine return `main_att` and `main_ame` of OPPOSITE SIGN from the same data
with the same R2, and wooldid's own multicollinearity caution fires on the explicit-dummy
arm only.

**Repro:** `repro_contreat_engine_split.do` (self-contained).

**Observed** (repro_contreat_engine_split.log:143-284):
```
untreatedZero DEFAULT(reghdfe)     att=  4.3577330  ame=  2.6678204  R2=0.310242  N=882
untreatedZero EXPLICIT(customvce)  att= -5.0726170  ame= -3.7078336  R2=0.310242  N=882
basic         DEFAULT(reghdfe)     att=  4.3781039  ame=  2.3722126  R2=0.310788  N=882
basic         EXPLICIT(customvce)  att=  4.3781039  ame=  2.3722126  R2=0.310788  N=882
```
Both `untreatedZero` rows report p < 0.001 (t = 5.64/5.94 and -3.51/-4.29). The true dose
effect in the DGP is +3.0. Only the EXPLICIT arm prints "Caution: a treatment period
coefficient interacted with something ... was dropped for multicollinearity" plus "Note: the
multicollinearity caution above was triggered 1 times total". The DEFAULT arm prints nothing, under `verbose` either (run11.log:172-222).

WHICH TERM. Re-run with `verbose`, as the caution instructs: the explicit arm names
"6o.__000005#3o.__000006#1o.__000007#co.d2" -- i.e. the DOSE SLOPE of the single cell
(_i = 6, _t = 3), omitted as collinear (run11.log:125-127). Identical R2 shows the two
engines reach the same FIT; it does NOT show they dropped the same term, and I did not
establish what reghdfe did with that slope. A contreat-focused reviewer should start there.

**Expected.** With `contreatcontrols(basic)` on the SAME data the two engines agree to the
last bit (ratio 1.00000 on both `att` and `ame`), so this is specific to `untreatedZero`
rather than to the engines. Either the screen should fire on both paths or the estimate
should not be engine-dependent; a user on the default path currently gets a sign-flipped
"significant" estimate with no warning of any kind.

**Confidence:** uncertain — I do not know the `untreatedZero` model formulation (ado:1313)
well enough to say which sign is correct, and a contreat-focused reviewer should adjudicate.
The ASYMMETRY OF THE WARNING is the part I am confident about.
**Severity:** medium. **Realistic:** yes as a specification (untreatedZero is a documented,
substantively-motivated choice); the particular collinearity may or may not be.
