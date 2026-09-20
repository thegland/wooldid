# fp-contreat-A — third lens (DOCUMENTED-BEHAVIOUR)

## Documentary record (complete; searched pinned build sthlp + ado + repo review_process)

1. **The help defines the quantity unambiguously.** sthlp:905-907 — eyex = "% change in y/%
   change in contreat"; eydx = "% change in y/unit change in continuous treatment".
   sthlp:243 — semielasticity = "% change in y / unit change in x", "implemented using the
   eydx option of margins".
2. **The ONLY documented default-path-vs-regress-path divergence is a STANDARD ERROR one.**
   sthlp:264: "the default estimation path and the regress-based path used by customvce()
   and unconditionalse apply different rules for fixed effects nested within the clustering
   variable, so the same specification can return slightly different STANDARD ERRORS
   depending on which path it takes." Point estimates are not mentioned anywhere.
   Corroborating: sthlp:532 (bootstrap) "point estimates are unchanged"; sthlp:1052 tells
   the user to write fe(i.state) "on any specification you may later run both ways" — the
   help's model is that the two paths are the same specification.
3. **No statement anywhere that the linear-path prediction excludes the absorbed effects.**
   grep over wooldid.sthlp for partial/prediction/absorbed: sthlp:249 and :371 say only that
   eydx "divides by the prediction" (in the context of the cell-collapse ban); sthlp:911 warns
   only that these elasticities need larger samples and non-sparse y. Nothing about what the
   prediction contains on the reghdfe path.
4. **The author's own written standard is the opposite.** ado:1459 comment, on the poisson
   twin: "ppmlhdfe requires d for predict, which margins requires" — and ppmlhdfe IS called
   with `d` at ado:1462/1468 while reghdfe at ado:1476 is not. ado:5024-5025: "poisexpresults
   is deliberately untouched: its scale is predict(n) = E[y] > 0, where eydx/eyex are the
   genuinely correct operators." The denominator the author treats as correct is E[y|x].
5. **Not ruled.** BS-6 (rulings_2026-09-07.md) is the POISSON xb-scale remap only; F032
   (research_2026-09-06) is the dyex analytic-Jacobian change; engine_parity.md is
   bootstrap-vs-permute and is irrelevant. Nothing in known_bydesign.md §1-6 or §7.
   review_process/assessment_grid.md D4 lists reghdfe / ppmlhdfe / plain reg (customvce) /
   plain poisson / unconditionalse as parallel legs of the SAME dimension to be cross-checked.

=> "documented or ruled" is answered NO. What remains is "correct-but-surprising design".

## My own discriminating test (probe L1_docbehav.do, section B)
Same engine (reghdfe) both times, same column space: fe(i.zone) [absorbed] vs
controls(i.zone) [main equation]. dydx and att must match to machine precision. If eydx/eyex
MOVE, the reported number depends on where the user wrote a term, not on the DGP — so it is
not an estimand at all, and "surprising but correct" is unavailable as a defence.

## Status
- do-file written: L1_docbehav.do
- waiting on the 8-process cap before launching

## MEASUREMENTS (mine, independent)

Fixtures: 20 states x 5 counties x 10 years (i = state, so treated cells hold 5 rows and
the within-cell dose slope is identified). State FEs carry most of y's level. N=1000.
Logs: L1c_docbehav.log, L1d_place.log, L1e_place2.log in this directory.

### (1) Engine parity, default(reghdfe) vs customvce(robust)(plain regress) — L1c
Same model, same data; only the variance estimator is supposed to differ.
    dydx  AME   DEF 0.898106454827   CVE 0.898106454827   RELGAP  0.0000000000
    dyex  AME   DEF 1.726089073045   CVE 1.726089073045   RELGAP  0.0000000000
    eydx  AME   DEF 0.018256146345   CVE 0.016007056918   RELGAP +0.1405
    eyex  AME   DEF 0.033719596671   CVE 0.030100558255   RELGAP +0.1202
    att         DEF -5.195339239851  CVE -5.195339239851  RELGAP  0.0000000000
BINARY side, no contreat at all:
    semielasticity att  DEF -0.125530252771  CVE -0.111536386352  RELGAP +0.1255
    plain att           DEF -5.193389348631  CVE -5.193389348631  RELGAP  0.0000000
Second seed (L1d): eydx +2.7%, eyex +2.5%, semielasticity +12.6%.
Third fixture (L1e): eydx +8.1%, eyex +9.6%, semielasticity -0.87%.
=> magnitude is design-dependent, measured 0.5% to 14% here; prober 13.7%, skeptic 2-4%.
   dydx / dyex / att are bit-identical across engines EVERY time.

### (2) THE DECISIVE TEST — placement invariance, reghdfe BOTH times (L1e)
`wave` = a randomly assigned 4-level county-year program indicator with a real effect on y,
orthogonal to treatment by construction. Written as fe(wave) [-> absorb()] or
controls(i.wave) [-> main equation]. Same column space, SAME FIT: R2 = 0.95987187894692 to
all 14 printed digits in both legs.
    att   FE 2.662412725378  CT 2.662412725378   RELGAP -0.0000000000
    dydx  FE 1.431779348211  CT 1.431779348211   RELGAP -0.0000000000
    dyex  FE 2.903627086455  CT 2.903627086455   RELGAP -0.0000000000
    eydx  FE 0.026051474007  CT 0.026098908955   RELGAP -0.0018100000
    eyex  FE 0.052408039570  CT 0.052433400290   RELGAP -0.0004800000
    semielasticity (binary) FE 0.048772828536  CT 0.048814629448  RELGAP -0.0008563
Identical data, identical fit, identical att and dydx to machine precision — and the three
divide-by-the-prediction quantities MOVE by 5e-4 to 2e-3 relative, 7-10 orders of magnitude
above the project's 1e-10 indifference threshold. The number therefore depends on where the
user wrote a control, not on the data-generating process. That forecloses "correct but
surprising": no correct estimand has that property. Magnitude scales with the share of y
absorbed; the unit FEs (always absorbed) are the big one, which is the 8-14% engine gap.

### (3) Mechanism, measured directly (L1c section D)
    full regress prediction  mean 50.80472236  sd 8.92652749  min 24.64040010
    reghdfe predict, xb      mean 50.80472236  sd 1.19252961  min 49.99234985
    y                        mean 50.80472236  sd 9.12496822  min 22.56026840
reghdfe's xb carries a _cons, so its MEAN equals mean(y) but its dispersion is gone: the
denominator is approximately the constant ybar. So the default path reports roughly
dydx/ybar — a semi-elasticity evaluated at the grand mean — where the help promises the
average of dydx/yhat_i. (The prober's "percentage of a partialled-out residual index" is
therefore wrong in detail; the skeptic's correction is right and I confirm it.)
Consequence the skeptic named and I endorse: the error is unbounded only where the partial
index nears zero, which is realistic for outcomes centred near zero (net flows, growth
rates, differenced outcomes) rather than for levels.

## VERDICT: CONFIRMED, high, realistic.
Not documented (sthlp:264 documents the only default-vs-regress divergence as STANDARD
ERRORS; sthlp:905-907/243 define the quantity as a % change in y; nothing says the
denominator omits the absorbed effects). Not ruled (BS-6 and F032 are the poisson-scale
rulings only; nothing in known_bydesign). Not correct-but-surprising (placement test).
Contradicts the author's own written standard at ado:1459 and ado:5024-5025.

SCOPE CORRECTION worth flagging to the author: the wider surface is `semielasticity` on the
plain BINARY default path — no contreat(), no continuous treatment, an ordinary option on an
ordinary staggered DiD run. It is affected by the same mechanism and I measured it three
ways.

## CAVEATS I am flagging myself (read before quoting my magnitudes)

1. **L1c/L1d carry a fixture confound.** In those two do-files `region = mod(_n,4)+1` and
   `coh = mod(_n,4)` — region IS the cohort partition, so `zfe` contains a cohort-by-post
   component collinear with treatment. That is why att came out -5.19 on a DGP with a +1.3
   dose effect, and why `fe(zone)` died r(322) ("a treated cohort-period dummy coefficient
   was dropped for multicollinearity"). The ENGINE-PARITY comparison is still valid there —
   both engines see identical data and fit identical models, and dydx/dyex/att are
   bit-identical while eydx/eyex/semi are not — but the 12-14% headline is measured on a
   misspecified fixture. **My clean magnitudes are L1e's**: eydx +8.1%, eyex +9.6%, binary
   semielasticity -0.87%, placement moves 5e-4 to 2e-3. (L1c is not unrealistic — real
   panels do carry unmodelled confounds — but the number should be quoted from L1e.)
2. **The sign of the gap is design-dependent and that is the same single mechanism.**
   Binary semielasticity was +12.5% on L1c and -0.87% on L1e. The gap has two parts: a
   Jensen term from the lost dispersion of the denominator, and a covariance term between
   the cell-level effect and the cell's level of y. When cell effects are heterogeneous and
   correlated with y's level the covariance term dominates and can flip the sign. One
   mechanism, two terms, either sign. (This is the objection the skeptic correctly raised
   against the prober's binary evidence; it does not undo the finding, and the placement
   test in L1e is immune to it because there the data and the fit are held fixed.)
3. **"Inherited upstream convention" is not a defence.** Stata's own `margins` after `areg`
   or `xtreg, fe` has the identical trap: `predict, xb` omits the absorbed effects, so any
   ey* operator divides by a partial index. One could call wooldid's behaviour an inherited
   convention. ado:1459 answers that: wooldid KNOWINGLY works around exactly this trap on
   the ppmlhdfe twin, passing `d` with the comment "ppmlhdfe requires d for predict, which
   margins requires". The linear path is the one place the author's own workaround is
   missing, not a place where the convention was adopted on purpose.
