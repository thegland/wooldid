# fp-estimands-refute -- SKEPTIC verdicts on the FIRST-PRINCIPLES ORACLE findings

Build asserted in every log: `.../scratchpad/bluesky/build/wooldid.ado`.
Every number below is from MY OWN data and MY OWN do-files. I read the seeded prober's
scripts but never ran them; each verdict rests on a repro I built from the claim's own terms.

Self-contained repros in this directory:
- `repro_semielas_own.do` / `.log`  -- FP-SEMIELAS-FE + two extensions
- `repro_uz_ringer.do`     / `.log`  -- FP-CTZERO-ENGINE corrected + the new ringer defect
- `repro_nopre_own.do`     / `.log`  -- FP-NOPRE-LATEENTRY
- supporting: `pA_semielas.do`, `pB_contreat.do`, `pC_mixed.do`, `pD_ringer.do` (+ logs)

---

## 1. FP-SEMIELAS-FE -- CONFIRMED (severity kept HIGH, certain)

Own DGP: 40 units, 10 periods, unit-level intercepts alpha_i = 60 + 2*(i-20), calendar trend,
staggered onset t=4 / t=7, constant level effect 8, sigma 1.

| arm | default reghdfe | customvce(robust) | ratio |
|---|---|---|---|
| level (dydx) | 8.219294310  `+1.070475805f158X+003` | 8.219294310  `+1.070475805f170X+003` | 1.000000 (2 ulp) |
| semielasticity (eydx) | 0.116514051 | 0.162023562 | **0.7191** |
| eydx, homogeneous levels (control) | 0.125533755 | 0.125907468 | 0.99703 |

(second, independent seed in `pA_semielas.log`: 7.885943102 vs 7.885943102 on the level;
0.111828326 vs 0.155101018 on eydx, ratio 0.7210; homogeneous control 0.99777.)

Mechanism reproduced with NO wooldid at all -- `reghdfe y post, absorb(i.i i.t)` + `margins,
eydx(post)` against `reg y i.i i.t post` + `margins, eydx(post)`, same data:
- dydx **1.000000000** exactly (8.216934249 both sides); eydx 0.118799295 vs 0.133531105
- `corr(reghdfe prediction, true fitted value) = -0.354195138`
- `SD(reghdfe prediction) = 3.94` against `SD(true fitted value) = 21.64`
- hand oracle `mean(b/xb_hdfe) = 0.118799295` == the reghdfe eydx to 9 dp;
  `mean(b/xb_full) = 0.133531105` == the explicit eydx to 9 dp.

Arithmetic check of the claimed mechanism on the `pA` design: mean(y)=69.36, mean(post)=0.3575,
so reghdfe's `_cons` = 66.5 and ln(74.5/66.5) = 0.1136 -- which is the 0.1118 the default arm
reports. The reported semi-elasticity is the effect measured against the SAMPLE GRAND MEAN.

The only message a default `semielasticity` run prints is the cell-collapse caution, whose
text ends "the run is slower, and **no answer changes**".

Two notes of my own (neither is a correction of the seeded write-up, which was right on
both points): my wooldid-free demo uses a CONTINUOUS `post`, where `margins, eydx` is
`mean(b/yhat)`; wooldid's own call is on the FACTOR `targetflag`, where it is the mean log
difference the seeded finding describes -- same denominator, same defect. And
`unconditionalse` ALONE is refused r(198) ("requires specifying either robust or clustered
standard errors"), so the reachable spellings are `unconditionalse robust`,
`unconditionalse cluster()` and `customvce()`.

### Extension A (measured, not inferred): contreatelasticitytype(eydx)/(eyex)
The seeded finding said "by extension". Measured, on a 30-unit continuous-treatment panel:

| operator | ame default | ame customvce(robust) | ratio |
|---|---|---|---|
| dydx | 1.972616603 | 1.972616603 | 1.000000 |
| eydx | 0.020997076 | 0.027770557 | 0.7561 |
| eyex | 0.051142250 | 0.067165897 | 0.7614 |

`att` (a dydx quantity here) is 5.005706073 on every arm. So exactly the ey* operators split.

### Extension B (NEW, mine): the default arm CRASHES where the explicit arm is correct
Design: treated units high-baseline (alpha 26.5..34), controls low-baseline (~6), effect -20,
sigma 0.3. **Every y > 0 (min 5.47). Every true fitted value > 0 (min 5.98).** But reghdfe's
`_cons + Xb` is NEGATIVE on 80 of 400 rows (min -3.57), because
`_cons` = mean(y) - b*mean(post) = 16.4 while b = -20.

- `wooldid y i t ttre, att semielasticity` -> **r(3598)** with a raw Mata traceback:
  `variable __00001O not found / st_data(): 3598 Stata returned error /
   __marg_cr_eval(): - function returned error / deriv__call2user_v(): ...`
  An internal tempvar name and a Mata stack; no wooldid concept is named, and nothing
  tells the user that the prediction the option divides by has gone negative.
- `wooldid y i t ttre, att semielasticity customvce(robust)` -> **-1.107274232**,
  which is right: mean ln((alpha-20)/alpha) over the treated-post cells is about -1.1.

**Provenance of the r(3598), checked (`pE_checks.log` E2).** The traceback is Stata's OWN
`margins` behaviour, not something wooldid generates: with no wooldid at all, on the same
data, `reghdfe y i.post, absorb(i.i i.t)` then `margins, eydx(post)` prints the identical
stack (`variable __000013 not found / st_data(): 3598 ...`), while the continuous spelling
`reghdfe y post` + `margins, eydx(post)` gives the informative r(459)
"could not calculate numerical derivatives -- discontinuous region with missing values".
So what wooldid owns here is (a) driving `margins` onto a prediction that has no business
being negative -- which is FP-SEMIELAS-FE -- and (b) not trapping the resulting r(3598)
into a wooldid-level diagnostic. That downgrades the crash itself to LOW; the root cause
stays HIGH.

---

## 2. FP-CTZERO-ENGINE -- PARTIAL: real split, WRONG CULPRIT, severity raised to HIGH

### (a) The seeded harness overstates its own design
Its `ui = 100 + 20*runiform()` is generated AFTER `expand`, so it is per-ROW noise with SD
about 5.8, not a unit intercept. "Low-noise, well-identified" is false; the effective residual
SD is ~19x the stated 0.3. I rebuilt it with a genuinely unit-level intercept.

### (b) On clean data the split is real -- and the DEFAULT arm is the CORRECT one
Truth is computable: dose ~ U(0.2,2.7), true dose effect 3.0, no pure onset shift, so
ame_true = 3.000 and att_true = 3.0 * E[dose | treated-post].

`repro_uz_ringer.log` (truth att 4.2397, ame 3.000):

| arm | att | ame |
|---|---|---|
| untreatedZero, **default reghdfe** | 4.2713432 | **2.9672819** |
| untreatedZero, customvce(robust) | 0.1577321 | **0.0565105** |
| untreatedZero, unconditionalse robust | 0.1577321 | **0.0565105** |
| untreatedZero, **poisson** (default engine, xb scale; truth 0.5653 / 0.400) | 0.0771772 | **0.1003708** |
| basic, default (control) | 4.2712002 | 2.9819023 |
| basic, customvce(robust) (control) | 4.2712002 | 2.9819023 (bit-identical) |

`pB_contreat.log`, independent seed (truth 4.2521 / 3.000): untreatedZero default
4.1976525 / 3.0022080; untreatedZero customvce -0.1027635 / -0.0318788.
`pC_mixed.log` C1, third seed (truth 4.3919 / 3.000): default 4.4161636 / 2.9908976;
customvce -0.1335864 / -0.1169391; unconditionalse robust identical to customvce.

So the seeded framing -- "a user on the DEFAULT path gets a sign-flipped significant estimate
with no warning" -- is INVERTED. The default path is correct on every seed. The explicit-dummy
path reports a true effect of 3.0 as an insignificant approximately 0 (p = 0.605 / 0.089 / 0.282).
"Opposite sign" is seed-specific; the invariant is "collapsed to ~0".

### (c) The mechanism is not rank deficiency -- it is the DF-1 "ringer"
`wooldid ... verbose` prints the three regressions (see `pC_mixed.log`, `pD_ringer.log`):

```
DEFAULT : reghdfe y  i._i#i._t#i.TF#c.d2  i._i#i._t#i.TF , absorb(i.i i.t) tol(1e-9)
EXPLICIT: reg y i.i i.t  i.__000009##c.(__000008)  i._i#i._t#i.TF#c.d2  i._i#i._t#i.TF , vce(robust)
POISSON : ppmlhdfe ycount i._i#i._t#i.TF#c.d2 i._i#i._t#i.TF ,
                 absorb(i.i i.t  i.__000009##c.(__000008)) tol(1e-9) d
```

`i.__000009##c.(__000008)` is the DF-1 "ringer" -- a POOLED dose control.
- ado:1233 installs it whenever
  `strpos(contreatcontrols,"basic") | unconditionalse | customvce() | poisson`.
  The test does NOT exclude `untreatedZero`.
- ado:1300 puts it in `contreatcontrolinserterFE`.
- ado:1313 -- the `untreatedZero` branch -- rewrites `contreatinj` ONLY. It never clears
  `contreatcontrolinserterFE`, which ado:1485 (reg/poisson) reads directly and ado:1462
  (ppmlhdfe) reads through `contreatcontrolinserterABS`.

The DF-1 comment justifies the ringer as "a numerical twin of columns the margins-visible
block already carries through its ## chain (c.`contreat`, ...). It adds no rank." That premise
holds under `basic`. Under `untreatedZero` there IS no margins-visible baseline dose term --
sthlp:1033: "By not including a g() function at all, identification of marginal effects of
Contreat_jt will occur using only within treated cohort-period variation." So the ringer is
not redundant: it is a real pooled dose slope, it absorbs the entire dose effect, and
`margins` switches off only `targetflag`, leaving the pooled slope in BOTH predictions, where
it cancels out of the reported contrast.

Demonstrated with NO wooldid, on the same data, by rebuilding the design by hand
(`repro_uz_ringer.do`, MECHANISM block):

```
reg y i.i i.t        i.cell i.cell#c.d2  ->  mean effect on treated-post =  4.308945
reg y i.i i.t c.d2   i.cell i.cell#c.d2  ->  mean effect on treated-post =  0.195330
                                    pooled coefficient on d2 =  2.910771
```
(`pD_ringer.log`, other seed: 4.307445 vs -0.054503, pooled coefficient 3.081415.)
The pooled term takes the whole effect; what margins reports is the remainder.

The multicollinearity caution the seeded finding leans on is a DOWNSTREAM SYMPTOM -- once the
pooled slope is in, one cell dose slope becomes redundant -- not the cause. It fires on the
arm that is wrong, but its text ("This could imply your results are not identified, depending
on exactly what was dropped") does not tell the user the estimate has been zeroed out, and a
verbose re-run names one cell, not the pooled control that did the damage.

### (d) NEW extension: the DEFAULT engine is hit too, via `poisson`
`poisson` is in the ado:1233 trigger list and its ringer goes INSIDE `absorb()`. With
`contreatcontrols(untreatedZero) poisson` on the default engine the true xb-scale ame of 0.400
comes back as 0.1004 (and 0.0737 on the other seed). So this is not "customvce vs reghdfe";
the trigger set is `untreatedZero x {unconditionalse | customvce() | poisson}`.
(The poisson arm also takes the rc-3200 fallback, which re-inserts the ringer as an explicit
regressor -- either placement leaves it margins-invisible.)

POISSON CONTROL ARM (`pE_checks.log` E1), same data, same xb-scale truth of 0.400 / 0.5653:

| poisson arm | att | ame |
|---|---|---|
| contreatcontrols(basic), default engine | 0.5442619 | **0.4377304** |
| contreatcontrols(basic), customvce(robust) | 0.5442619 | 0.4377304 (bit-identical) |
| contreatcontrols(untreatedZero), default engine | 0.0771772 | **0.1003708** |

So the xb-scale truth IS recoverable under poisson and the two engines agree there; the
untreatedZero poisson arm is the deviant, not the scale or the DGP.

### (e) Not already ruled
`review_process/cell_collapse_contreat_2026-09-02.md:118-120` records that
"`contreatcontrols(untreatedZero)` **alone** puts the ringer nowhere - it only rewrites
`contreatinj`". That is correct, and is precisely why the COMBINATION was never examined.
Nothing in `known_bydesign.md` or the review record covers
untreatedZero + unconditionalse / customvce / poisson.

---

## 3. FP-NOPRE-LATEENTRY -- CONFIRMED (low/certain)

`repro_nopre_own.do`. Panel t=1..8, 30 units; units 11-16 have `ttre = 4` (comfortably inside
the window) but are first observed at t=4. ado:12053 computes `mint = min(t) by i` and drops on
`ttre <= mint`, so those six units go: `e(N) = 192` of 222 rows, `e(N_dropped) = 30`. The
header prints "N = 192 (30 obs dropped from initial sample)" and NO caution of any kind fires
(grep for "Caution" over the whole arm returns nothing). Rebuilt with the same units observed
from t=3 -- one pre period -- nothing is dropped and `e(N) = 228`. Because those six carry a
true effect of 6 against the others' 2, the reported att moves 2.1609926 -> 3.9429456, i.e.
the run silently estimates over a different population than the user specified.

sthlp:73 says only "cohorts treated before the sample begins are dropped". A help-wide grep
for "first observed", "no pre-treatment", "without a pre" and "first usable" finds only
sthlp:824 -- the `anticipation()` twin of the same screen, which DOES print a caution
(ado:12036-12046, "places the effective onset of N unit(s) at or before their first usable
pre-period; those unit(s) and their N observations are dropped"). So the per-unit-versus-
per-sample gap is genuinely undocumented and the unconditional screen is genuinely silent
about the cause. Dropping the unit is CORRECT behaviour (no within-unit pre-period exists),
so this stays a documentation/UX low.

---

## 4. FP-CWM-DOC -- CONFIRMED (low/certain), with one nuance

sthlp:222-236 spells `att_it`, `att_i` and `att_itime` as `w/W_it`, `w/W_ipre|post` and
`w/(W_it*T_iprepost)` with no `customweightmultiplier()` factor and calls each result "the
simple, unweighted mean". sthlp:1107 authorises the multiplication generally. The nuance the
seeded write-up does not state: the **att** phang at sthlp:220 ALREADY names cwm explicitly
("multiplied by the customweightmultiplier() variable if one is specified - that multiplier
enters the margins weights only and is deliberately not folded into the regression weight").
So it is exactly the three cell-weighted phangs that omit it, not the block as a whole.

Numeric confirmation on my own one-row-per-cell panel (`pD_ringer.log` D4): `att_it` =
4.552181151 without cwm and 4.548217942 with it, on identical data -- so the reported quantity
is demonstrably not the simple unweighted mean of the cell effects those three entries
describe. (I did not re-derive the exact weight formula; the seeded oracle's finding that the
CODE is coherent and only the printed formulas are short is consistent with everything I saw.)

---

## 5. WHAT I ATTACKED AND DID NOT BREAK (honest negatives)

- **wooldid_ccscreen blind on reghdfe?** I hypothesised that the screen (ado:11818), which
  scans `: colnames e(b)` for "o.", might be structurally dead on the default engine because
  reghdfe strips omitted terms. Measured directly (`pB_contreat.log` B3): reghdfe KEEPS them --
  `1b.g#co.z 2o.g#co.z 3o.g#co.z 4o.g#co.z 5.g#c.z ...` for a structurally-zero interaction,
  and `z o.zdup _cons` for a duplicated regressor. The screen can see omissions on both
  engines. REFUTED (my own hypothesis).
- **`contreatcontrols(basic)` engine parity**: bit-identical across default and customvce on
  three independent seeds. The untreatedZero defect does not generalise to the other control
  types.
- **Is the untreatedZero collapse just "opposite signs"?** No. Across three seeds the
  explicit arm returned +0.158, -0.103 and -0.134 for att and +0.057, -0.032 and -0.117 for
  ame against a truth of ~4.25 / 3.00. The invariant a user will see is **collapsed to ~0**
  (an insignificant null), not a sign flip -- the sign is whatever the residual noise leaves.
- **`dydx` / level estimands across engines**: agree to 1-2 ulp on every design I built,
  including the continuous-treatment `contreatelasticitytype(dydx)` arm (5.005706073 /
  1.972616603 identical). The engine split is confined to the ey* operators and to
  untreatedZero.
