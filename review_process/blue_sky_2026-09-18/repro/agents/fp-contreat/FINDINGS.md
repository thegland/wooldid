# fp-contreat — findings (INDEPENDENT REIMPLEMENTATION of the continuous-treatment AME family)

Method: EXACT-FIT ORACLE. Every probe generates a panel whose outcome is produced with zero
residual by a model nested in wooldid's continuous-treatment specification, so E[y|x] == y,
the margins-visible dose derivative is the planted beta_it, and every reported quantity has a
closed form I compute myself. No reconstruction of wooldid's tempvar regression is needed.

---

## FINDING A — contreatelasticitytype(eydx)/(eyex) divide by a prediction that omits the
## absorbed fixed effects on the DEFAULT (reghdfe) path

- repro: `repro_A_elasticity_absorbedFE.do`  (also probe04/05/06)
- surface: `contreat()` + `contreatelasticitytype(eydx|eyex)`, default reghdfe estimation
- severity: HIGH   confidence: CERTAIN   realistic: YES   seeded: NO (found by this method)

CLAIM. On wooldid's default estimation path the reported elasticity / semi-elasticity of the
continuous treatment is `dydx / xb`, where `xb` is reghdfe's linear prediction EXCLUDING every
absorbed fixed effect (i.i, i.t, fe(), and the absorbed dose control). The help defines these
as "% change in y" (sthlp:905-907). They are not: they are a percentage of a partialled-out
residual index. The same estimand on the plain-regress paths (`customvce()`, `unconditionalse`),
where the fixed effects are explicit regressors, is exact.

OBSERVED (repro_A, 20 states x 5 counties x 10 periods, staggered onsets, exact fit, R2 = 1):

| quantity | cluster(unit)  [reghdfe]| customvce(cluster unit) [reg] | closed-form truth |
|----------|-------------------------|-------------------------------|-------------------|
| dydx     | 1.320000000000          | 1.320000000000                | 1.320000000000    |
| dyex     | 2.647372379686          | 2.647372379686                | 2.647372379686    |
| eydx     | **0.016412789649**      | 0.019020203082                | 0.019020203082    |
| eyex     | **0.032813176169**      | 0.038004324224                | 0.038004324224    |

The reghdfe path is 13.7% too small on both elasticity rows. dydx and dyex are exact on both
paths (they do not divide by the prediction), which isolates the denominator as the cause.

MECHANISM, established by a graded experiment (`probe06_eydx_scale.do`): scale the absorbed
part of y by s and re-measure the relative error of the reghdfe path against the closed-form
truth.
    s = 0    (no fixed-effect variation at all) -> eydx rel.err 1.8e-09   (i.e. exact)
    s = 0.1                                      -> 1.72e-03
    s = 1                                        -> 1.30e-02
    s = 3                                        -> 2.48e-02
The error is zero exactly when there is nothing for reghdfe to absorb, and grows with the
share of y's level that is absorbed. In real panels the unit fixed effects usually carry MOST
of the level of y, so the error is not small: it is bounded only by how large the absorbed
level is relative to the residual index, and the denominator can in principle cross zero.

EXPECTED. Under an exact fit the conditional mean IS y, so the semi-elasticity is
mean(beta_it/y) and the elasticity mean(beta_it*dose/y) over the treated post rows. Both
estimation paths fit the identical model to the identical data; only the VCE differs. A point
estimate that moves when the user switches from `cluster(unit)` to `customvce(cluster unit)`
cannot be right on both paths, and only the plain-regress one matches the closed form.

NOT covered by the known-by-design record: `known_bydesign.md` has nothing on this, and the
sthlp discussion of eydx (sthlp:245-249) concerns only (a) numerical failure on sparse y and
(b) the poisson xb-scale remap, which is a different scale question. The remap comment at
ado:5010 is explicitly about poisson; nothing anywhere says that off the poisson path the
elasticity denominator is the partialled-out index rather than the fitted mean.

---

## FINDING B — `lattice()` + continuous treatment: the all-zero VCV is a property of the
## ABSORBING ENGINE, not of the design. Plain `regress` on the identical model and the
## identical data returns an ordinary robust VCE.

- repro: `repro_B_lattice_zeroVCV.do`   (built from probe10_p1like.do / probe11_lat_diag.do)
- surface: `contreat()` + `lattice()` on the default reghdfe path, `robust` or `cluster()`
- severity: HIGH   confidence: CERTAIN   realistic: YES
- seeded: YES (the diagnosis run that held_for_ruling_2026-09-13 item 5 asked for and never got)

CLAIM. On a realistic 12-state / 5-12-counties-per-state-year-cell panel with a unit-level
dose, every `lattice()` run on wooldid's default (reghdfe) estimation path returns an all-zero
covariance matrix and loses all 12 reported standard errors, while the SAME model fitted by
plain `regress` on the SAME data - identical N, identical df_r, identical rss, identical r2,
IDENTICAL point estimates to 8 decimals - returns finite standard errors. The zero is
therefore not the documented cell-containment geometry (which would afflict `regress` too),
and not a lattice-construction defect (a user-supplied lattice variable does it too).

OBSERVED (probe11_lat_diag.log; N = 2200, esprelength(2) espostlength(2),
contreatcontrols(basic), every leg rss = 3.317e+03, r2 = 0.642998, df_r = 1672 where the
VCE is not clustered, point estimates att = 0.78718487 and att_pre = 0.15042830 on EVERY leg):

| leg                                        | engine       | se on the att_pre row |
|--------------------------------------------|--------------|-----------------------|
| `lattice(all 3) robust`                    | reghdfe      | **missing** (V == 0)  |
| `lattice(all 3) robust regtol(14)`         | reghdfe      | **missing**           |
| `lattice(all 3) robust regtol(6)`          | reghdfe      | **missing**           |
| `lattice(all 3) ccnoabsorb robust`         | reghdfe      | **missing**           |
| `lattice(latu) robust` (user lattice var)  | reghdfe      | **missing**           |
| `lattice(all 3) cluster(st)`               | reghdfe      | **missing**           |
| `lattice(all 3) customvce(robust)`         | plain regress| **0.16205346**        |
| `lattice(all 3) unconditionalse robust`    | plain regress| **0.20168098**        |

The zero is invariant to regtol (6, 9 and 14 all give it), so it is not a tolerance artifact.
It is invariant to `ccnoabsorb`, so it is not caused by what is inside `absorb()` - reghdfe is
still the engine there. It fires under plain `robust` with no `cluster()` at all, which is the
case `known_bydesign.md` section 5 explicitly says is NOT documented and IS a finding.

It is also not confined to `lattice(all 3)`: `lattice(all 2)` and `lattice(all 5)` on the same
fixture do it too (probe10_p1like.log), with rss 3.637e+03 / 2.938e+03 and df_r 1795 / 1511.

SECOND HALF OF THE SAME OBSERVATION (worth its own weight): on this fixture the lattice runs
also produce WILD ame points with no standard error to warn the reader. The DGP's dose effect
is exactly +0.4 everywhere. Without a lattice wooldid returns ame 0.595 (se 0.286) and
pre-ame 0.142. With a lattice it returns, with every inference column blank:
    lattice(all 2): ame  7.013 , pre-ame -3.488
    lattice(all 3): ame  0.949 , pre-ame  3.702
    lattice(all 5): ame  5.789 , pre-ame -54.10
A user who ran only the default path would see a point estimate 17x the truth and no way to
tell it is noise. The help's lattice caveat (sthlp:1010-1018) is about effects loaded onto the
BIN INTERCEPTS; here the true dose response is linear with no bin structure at all, so that
caveat does not describe this.

EXPECTED. An all-zero covariance is documented at sthlp:1076 only for the cell-containment
geometry, i.e. continuous treatment with `cluster()` at or above `i`. Under `robust` there is
no containment argument, df_r is 1672 and the residual sum of squares is 3.3e+03, so the
variance is not algebraically zero - and the plain-regress leg proves it by computing it.

---

## FINDING C — the all-zero-VCV cautions blame a cluster VCE that was never requested, and
## point the reader at a caution that did not print

- repro: `repro_B_lattice_zeroVCV.do` (leg `lattice(all 3) robust`); text captured in
  `probe10_p1like.log`
- code: `ado:12983` (wooldid_amescreen, mode vzero) and `ado:1567`
- severity: MEDIUM   confidence: CERTAIN   realistic: YES   seeded: NO

CLAIM. When the estimation command returns an all-zero covariance under `contreat()`, wooldid
prints two cautions whose content is false whenever no `cluster()` was given:

1. `ado:12983` is gated ONLY on `"\`contreat'" != ""` and on the R1 cell-containment screen not
   having fired. It says: "The issue is likely related to continuous treatment effect estimation
   coupled with use of **analytic clustered SEs, with clusters occuring at a high level of
   aggregation**. Consider ... if finer clustering levels are non-viable." On the repro the VCE
   is `robust`; there are no clusters at all, so the named cause cannot be the cause and the
   named remedy ("cluster finer") is not available.
2. `ado:1567` then says "This is a property of the specification, **not a numerical accident**:
   **see the cell-containment caution above** for the mechanism and for the alternatives". No
   cell-containment caution printed - the R1 detect (`mode(detect)`) only runs over `clusvars`,
   which is empty under `robust` - so the pointer dangles; and the plain-regress leg of
   FINDING B shows the zero IS an artifact of the engine rather than a property of the design.

EXPECTED. Either the two cautions should be conditioned on a cluster-type VCE being in force,
or the unclustered case needs its own text. A user following the printed advice would spend the
afternoon looking for a finer clustering level for a run that has no clustering.

---

## FINDING D (low) — `fe()` with `customvce()`/`unconditionalse`: wooldid cautions, then silently
## enters the fe() variables as CONTINUOUS regressors

- repro: `probe08_customvce_dbg.do`, leg A
- severity: LOW   confidence: CERTAIN   realistic: YES   seeded: NO

`fe(unit)` under `customvce(cluster unit)` prints "Caution: unconditional SEs and custom VCE
specifications ... are not compatible with absorbed fixed effects. You must specify anything in
fe as controls, manually specifying any needed i. prefixes." and then runs
`reg y_dose i.i i.year unit  ...` - `unit` WITHOUT an `i.` prefix, i.e. a linear trend in the
identifier. The run completes rc 0 and reports R2 = 0.332 and ame = -0.178 where the correctly
absorbed twin reports R2 = 1.000 and ame = 1.370 (the truth). The caution is printed, so this is
not silent; but "warn, then fit a regression on the numeric value of a unit id" is a predictable
mis-read, and an error would cost nothing since wooldid already knows the combination is invalid.

### FINDING B — mechanism, reproduced OUTSIDE wooldid (probe14_reghdfe_repl.do)

The same three regressions, hand-built on the same fixture (cell = state x relative-year for
treated rows in the window, `lat` = dose tertiles, `ringer` = a double copy of the dose):

| command | max abs e(V) | df_r | rss | e(rank) | df_a |
|---|---|---|---|---|---|
| `reghdfe y_dose ib0.cell##i.lat##c.dose_it, absorb(i.st i.year i.lat##c.(ringer)) vce(robust)` | **0** | 2004 | 4100.2092 | 160 | 36 |
| `reghdfe y_dose i.lat##c.(ringer) ib0.cell##i.lat##c.dose_it, absorb(i.st i.year) vce(robust)` | **0** | 2004 | 4100.2092 | 165 | 31 |
| `reg y_dose i.st i.year i.lat##c.(ringer) ib0.cell##i.lat##c.dose_it, vce(robust)`            | **2923.571** | 2004 | 4100.2092 | 154 | - |

Identical residual sum of squares and identical residual degrees of freedom; only the engine
differs. So the zero is not a property of the design and not wooldid's own arithmetic - it comes
back that way from reghdfe. PROBABLE specific cause (not certain): reghdfe retains 6-11 more
columns than `regress` drops as collinear (rank 160/165 vs 154) on a design where the dose is
nearly constant inside a 1-4-county (cell x bin) group, so its X'X is singular where regress's
is not. Reporting the engine split as certain and the column-retention cause as probable.

A wooldid-side remedy exists and is cheap: the vzero screen already detects the all-zero V; on a
NON-clustered VCE it could refit with `regress` (which wooldid already knows how to drive - it is
the `customvce()` path) instead of discarding every standard error in the run.

### FINDING B — what it contradicts in the project record

`ado:1546-1565` records that the AMESE **R4 att-rescue refit** ("a re-estimation of the same
specification through regress with the fixed effects written out as explicit dummies") was
REMOVED on 2026-09-08 by author ruling, on the reasoning that

> What it bought was the subset of standard errors that survive a specification whose
> cell-level containment has already made the ame inference invalid

That premise holds for the clustered containment geometry the ruling was written about. It does
NOT hold for the case measured here: the VCE is `robust`, no cluster exists, the R1 containment
predicate cannot fire (it only runs over `clusvars`), df_r is 1672-2004 and the residual sum of
squares is 3.3e+03 - and `regress` computes an entirely ordinary robust variance for the same
fit. So on this design the removed rescue was buying real standard errors for a specification
that is NOT broken, and the run now discards every standard error it could have had.

The other stated reason for removing it - that the explicit-dummies solve is numerically noisier
(1.6e-8 vs 1e-13 under row reordering, research record research_2026-09-08) - is a cost worth
weighing against "no inference at all", and it does not apply to the POINT estimates, which
probe11 shows are identical across the two engines to 8 decimals.

### FINDING A — the precise code line, and the remedy that is already in the file

`ado:1476` (the linear path):
```
qui reghdfe `y' ... if `ifstatement' , absorb(i.`i' i.`t' `fe' ...) `vcespecifier' tol(1e-`regtol')
```
`ado:1462` / `ado:1468` (the poisson path):
```
qui cap ppmlhdfe `y' ... , absorb(i.`i' i.`t' `fe' ...) `vcespecifier' tol(1e-`regtol') d
```
The poisson engine is given `d`, so it carries the sum of absorbed effects and `predict`
returns the full prediction. reghdfe is NOT, so `predict` after it returns `xb` alone, and every
margins read that DIVIDES by the prediction - `contreatelasticitytype(eydx)`, `(eyex)`, and the
binary-side `semielasticity` - divides by an index that is missing i.i, i.t, fe() and the
absorbed dose control. The operators that do not divide by the prediction (`dydx`, `dyex`) are
unaffected, which is exactly the pattern measured.

---

# NEGATIVE RESULTS (probed hard, came back CLEAN)

These are reported as contributions. In each case the probe HAD power: a deliberately
mismatched arm moved the answer by 10-90%, so a defect of the shape being hunted would have
shown up.

1. **`contreatcontrols()` partials out exactly what it claims** - 54 configurations,
   `probe02_invariance.do`. Six planted nuisance dose functions g(dose) - none, a pooled
   delta*dose, delta_t*dose, delta_i*dose, delta*1{ever treated}*dose, delta_{ttre}*dose -
   crossed with nine contreatcontrols() spellings (basic, byTreated, byTtre, byCohort, byTime
   and the four "basic + byX" pairs). Every MATCHED pair returns mean(beta) to <7e-14 with
   R2 exactly 1; every MISMATCHED pair moves the ame by 0.68% to 91%. A drop-order inversion
   between the margins-visible `c.dose` and the invisible `ringer` copy - the DF-1 hazard the
   code comments describe at ado:1226-1300 - would have moved the ame by exactly delta = 0.9
   on at least one of the 54; none did.
2. **Engine-path invariance of the POINT estimates** - `probe03_engines.do`. The same battery
   under DEFAULT (reghdfe absorb), `ccnoabsorb` (reghdfe, controls explicit) and
   `customvce(robust)` (plain regress, fixed effects explicit): all 15 legs exact. In
   particular the base-level omission in `i.evertreated#c.(ringer)` that a pure `#` interaction
   implies on the explicit-regressor path does NOT cost the run its dose control.
   (`unconditionalse` alone is refused r(198) with "requires specifying either robust or
   clustered standard errors" - informative, correct.)
3. **All four estimand weightings and both weighting states** - `probe13_estimands_es.do`.
   att, att_it, att_i and att_itime, unweighted and `[aw=w]`, against weights I recomputed
   from the documented definitions (w/w_per_it, w/w_per_i, w/(w_per_it*t_per_i)): all eight
   main ame rows exact to 1e-15, R2 = 1 on every one. A weight-construction defect - a wrong
   `by()` grouping, a pre/post pooling error, a float store - would have separated at least one
   of the four from its recomputed truth.
4. **Event-study ame horizons** - same probe, with a planted PRE-TREATMENT dose slope
   (+0.50 at relyr -1, -0.40 at relyr -2) so the pre reads are non-trivial. Every horizon
   (-2, -1, 0, +1, +2) exact to 1e-15 weighted and unweighted.
5. **The four contreatelasticitytype() values are mutually consistent, and dydx/dyex are exact
   on every path** - `probe05_eydx_mech.do`. dydx = mean(beta) and dyex = mean(beta*dose) agree
   with the closed form to 1e-15 on reghdfe, plain regress and unconditionalse alike. Only the
   two operators that divide by the prediction disagree, which is what isolated FINDING A.

6. **The constant-dose seeded lead is CLEAN, in both its absurd and its realistic form** -
   `probe15_constdose.do`. Four dose shapes on the same noisy design:
     `vary`   (dose varies within every cell)                      -> ame .36585649 (se .17), identified
     `global` (dose = 2.5 on EVERY observation)                    -> **ame MISSING**, two cautions
     `bycell` (dose constant within each (i,t) cell - the realistic
               state-level policy-intensity case, sd .241 overall) -> **ame MISSING**, two cautions
     `byunit` (dose constant within unit over time, varying across
               units inside a cell)                                -> ame 1.10299961, identified (correctly:
                                                                      within-cell dose variation exists)
   The cautions name the mechanism exactly: "every treated cohort-period dose-slope term for the
   continuous treatment variable dose was dropped for multicollinearity ... The within-cell
   average marginal effects (ame) that wooldid reports for continuous treatments are NOT
   IDENTIFIED in this design - the regression's ame entries are structural zeros, not estimates"
   (wooldid_ccscreen R6, ado:12xxx). The att side still reports and says so. This is an
   informative refusal, not a silent number, and the realistic within-cell-constant case is
   handled identically to the degenerate global one. NOT a defect.

8. **contreatpoly(2), a user lattice with bin-specific slopes, and contreatwithin are exact**
   - `probe16_poly_lat_uz.do` / `probe17_within_collapse.do` inside `probe19_merged.log`.
     contreatpoly(2): ame = mean(b1 + 2*b2*dose) to 1e-15 on reghdfe and on plain regress.
     lattice(uservar) with bin-specific TRUE slopes: ame = mean(beta_{cell,bin}) to 4e-16.
     contreatwithin on a design with NO never-treated units at all: exact to 1e-15 both paths.
9. **The default cell collapse is exact** - `probe17`, noisy data, `[aw=w]`, esprelength(1)
   espostlength(2), all four estimands x {dydx, dyex}: seven of eight legs are BIT-IDENTICAL
   in hex to the `nofastcollapse` twin. The eighth (att_i + dyex) differs by 2 ulps
   (+1.03cf1b5e65078X+001 vs +1.03cf1b5e65076X+001, ~2e-16 relative) - far below the project's
   1e-10 line, reported only because the code comment claims bit-exactness.
10. **The pre-treatment ame row, and an unbalanced panel** - `probe20_pre_es.do` and
   `probe18_unbalanced.do`. The pretreatment ame (roweq `pretreatment`) matches the closed form
   for all four estimands weighted and unweighted, to <5e-14, on a design carrying a planted
   pre-trend dose slope. On a messy panel with 22% of unit-periods missing at random, all four
   estimands weighted and unweighted are exact to 1e-15 with R2 = 1 - the case where the cell
   and cohort weight denominators (w_per_it, w_per_i, t_per_i) actually differ from each other.

11. **`if`, `subgroups()`, missing dose, and the missing-contreatcontrols refusal** -
   `probe22_if_sg_miss.do`, exact-fit design, R2 = 1 on every leg.
     `if t <= 6`     -> ame 1.1300000000 = truth (7.9e-16)
     `if i != 2`     -> ame 1.1680000000 = truth (1.3e-15)   (drops a whole treated cohort)
     `if dose > 1.5` -> ame 1.1320000000 = truth (1.6e-15)   (an `if` on the dose itself)
     `subgroups(sg)` -> e(wdidsgresults_0) ame 1.3150000000 and e(wdidsgresults_1) ame
                        1.0262500000, each equal to its own subgroup's mean(beta) to 2e-16
     dose missing on 10% of rows -> ame 1.1414851485 = the truth recomputed over the
                        non-missing rows (1.6e-15)
     `contreat()` with no `contreatcontrols()` -> r(198), as sthlp:983 promises.
   A sample-concept defect - the ame reading e(sample) where it should read the subpop, or the
   margins weights not being renormalised on a restricted sample - would have shown up on the
   `if` legs or on the missing-dose leg; neither moved.

---

# UNTESTED / NOT REACHED

- `contreat()` under `permute()` and `bootstrap()`: the ame rows' permutation and bootstrap
  inference, and whether FINDING A's denominator problem propagates into replication-based
  intervals (it would, since every replication re-runs the same margins call).
- `impvalues` + `contreat()`: the IM ame path takes a different margins branch (`nose over(cluster)`)
  and was not measured against a closed form.
- `makecfxplots` / `e(cfxame_*)` per-bin point estimates against a per-bin closed form.
- `esfixedbaseperiod` / `esrelativeto()` interaction with the ame pre row.
- `coarsecohortcontrols()`, `interactivecontrols()` and `timetrends` combined with `contreat()`.
- `poisexpresults` (the one remaining path where margins' ey* operators still reach the engine;
  there the engine IS ppmlhdfe with `d`, so FINDING A should not apply, but it was not measured).
- The FINDING B zero-VCV boundary: I did not map where it starts and stops as a function of
  cell size, lattice count and dose within-cell variation. I established the engine split, not
  the frontier.
- UNVERIFIED HYPOTHESIS worth one probe by whoever picks this up: `ado:2652`
  `matrix ___C..._horiz = (f9tab[1..6, 1..colsof(r(b))]' , relyrcolumn)` concatenates the
  margins `over(_t)` output with `relyrcolumn`, which is built by regex over e(b)'s column
  names and SKIPS omitted (`o.`) coefficients. If a treated relative period keeps rows in the
  subpop while EVERY one of its coefficients is omitted, margins returns one more column than
  relyrcolumn has rows and the concatenation is non-conformable (r(503)), outside any `cap`.
  I could not engineer that state in the time available: absorbing a cell's intercept via fe()
  leaves its `#c.dose` slope, whose name still matches the regex, so the horizon stays in
  relyrcolumn and the counts agree. Not filed as a finding - no repro.

### FINDING A — second surface (same mechanism, outside my assigned target)

The binary-side `semielasticity` option takes `margins, eydx(targetflag)` through the same
prediction, so the same denominator problem applies there. probe06_eydx_scale.do measured only
the path DISAGREEMENT, not a closed form (margins' eydx on a factor variable is a discrete
difference of logs, and I did not derive its exact form):
    fixed-effect scale s = 0 : robust 0.322104771542 == customvce 0.322104771542  (identical)
    s = 1                    : robust 0.191049620678 vs customvce 0.191768105033  (3.7e-3 rel)
    s = 3                    : robust 0.104915119469 vs customvce 0.104942900858  (2.6e-4 rel)
Same signature as the contreat side - identical with nothing absorbed, divergent once there is.
Reported as a pointer, not as a measured error.

### FINDING B — the self-contained repro's own numbers (repro_B_lattice_zeroVCV.log)

Every leg: N = 2200, df_r = 1663, rss = 3249.87, r2 = 0.653964, and the point estimates agree
across engines to eight significant digits ON EVERY ROW, the ame rows included.

| leg | att se | att_pre se | **ame b** | **ame se** | ame_pre se |
|---|---|---|---|---|---|
| `lattice(all 3) robust`                 | . | . | 18.09100052 | **.** | . |
| `lattice(all 3) robust regtol(14)`      | . | . | 18.09100052 | **.** | . |
| `lattice(all 3) ccnoabsorb robust`      | . | . | 18.09100052 | **.** | . |
| `lattice(latu) robust` (user lattice)   | . | . | 18.09100052 | **.** | . |
| `lattice(all 3) customvce(robust)`      | .13937285 | .15399566 | 18.09100039 | **.43558016** | .74495755 |
| `lattice(all 3) unconditionalse robust` | .14966591 | .19849945 | 18.09100039 | **10.040394** | 3.49933 |
| no lattice, `robust`                    | .13513692 | .17420902 |  0.65200913 | .12109561 | .22023536 |

TRUE ame on this DGP = 0.4 exactly.

Two things follow, and the second is the one I would put in front of the author:

(a) The engine split is confirmed on a second, independently generated fixture: reghdfe
    returns an all-zero V, `regress` on the identical fit does not.
(b) The lattice ame point is **18.09 against a truth of 0.4 - 45x** - and NO path warns the
    reader properly. The default path prints cautions but no standard error at all. The
    `customvce(robust)` path reports se 0.436, i.e. "18.09 (0.44)", a 40-sigma miss presented
    as a precise estimate. Only `unconditionalse` (se 10.04) makes the number look like the
    noise it is. The cause is that with a unit-level dose plus 10% time noise, a dose-quantile
    lattice leaves almost no within-(cell x bin) dose variation, so each cell-by-bin slope is a
    near-zero-denominator ratio; the help's lattice caveat (sthlp:1010-1018) attributes
    discordant ames to effects being loaded on the BIN INTERCEPTS of the true f(), which is a
    different cause and does not warn about this one. I mark (b) probable rather than certain
    because a 45x miss with a small conditional SE could in principle be a legitimate
    conditional-variance statement about an ill-conditioned design rather than an SE defect;
    what is certain is that the shipped default path gives the user neither the SE nor a
    caution that names within-bin dose variation.
7. **All twelve documented `latticetype` values are exact** - `probe21_lattypes_poisson.do`.
   On an exact-fit design whose true dose response is cell-constant with NO bin structure,
   `lattice(<type> 3)` for all of `all byTime byTreated byTtre byCohort byTreatedByTime
   byTtreByTime byCohortByTime baseTreated baseTreatedRef baseTreatedPost baseTreatedPre`,
   each run unweighted and again `[aw=w] latticeignoreweights`, returns the closed-form ame
   (0.9325000000 unweighted, 0.9236953192 weighted) to <6e-15 on all 24 legs. A lattice
   CONSTRUCTION defect - wrong quantile group, wrong by(), a bin that swallows the reference
   period, a weight leaking into `latticeignoreweights` - would have moved at least one.
   (This is the same fixture family on which FINDING B's zero-VCV does NOT fire, which is
   itself informative: the zero needs the near-degenerate within-bin dose geometry, not a
   lattice as such.)

---

## FINDING E — `contreatcontrols(untreatedZero)` silently gains a baseline dose control on the
## plain-regress and poisson paths, and the reported ame is wrong (here: SIGN-FLIPPED)

- repro: `repro_C_untreatedZero.do` (first observed in probe16 inside `probe19_merged.do`)
- surface: `contreat()` + `contreatcontrols(untreatedZero)` + any of `unconditionalse`,
  `customvce()`, `poisson`
- severity: HIGH   confidence: CERTAIN (linear path measured; poisson twin pending)
- realistic: YES   seeded: NO

CLAIM. `untreatedZero` is defined as the control type that adds NO dose control at all -
ado:1311-1313 "if always zero when not treated, the only thing to do is not control for
anything! ... it is implemented via modifying the contreatinj statement and just not having
contreatcontrols of any kind" - and sthlp:995 says it "cannot be paired with any other type".
But the block that builds `basic`'s pooled control is gated at ado:1234 on

    if strpos("`contreatcontrols'", "basic") != 0 | "`unconditionalse'" == "unconditionalse" ///
       | "`customvce'" != "" | "`poisson'" == "poisson"

with no untreatedZero exclusion, so `ccringerhead = i.`null'##c.(`ringer')` is built and flows
into `contreatcontrolinserterFE` / `...ABS`, which the untreatedZero branch below never clears.
On an untreatedZero design the dose is exactly 0 outside treated cohort-periods, so that pooled
slope is NOT redundant with anything: it is identified off the treated post rows alone and eats
part of the dose response.

OBSERVED (probe19_merged.log, exact-fit design, dose == 0 on every untreated observation,
truth = mean(beta_it) over treated post rows = 0.9325000000):

| leg | reported ame | truth | rel. err |
|---|---|---|---|
| `contreatcontrols(untreatedZero) robust` (default reghdfe) | 0.9325000000 | 0.9325 | 2.9e-15 |
| `contreatcontrols(untreatedZero) customvce(robust)`        | **-0.1075000000** | 0.9325 | **1.115** |

A sign flip on a documented option pair, with no caution and rc 0. The default path is right,
so the two paths - which differ only in how the variance is computed - disagree by 111%.

EXPECTED. `untreatedZero` must take the same control set on every engine: none. The one-line
fix is to exclude it from the ado:1234 gate (or to clear
`contreatcontrolinserterFE`/`...ABS` inside the untreatedZero branch, which already rewrites
`contreatinj`).

CONFIRMED ON ALL THREE GATED PATHS (`repro_C_untreatedZero.log`; TRUE ame = 0.9325000000
on the linear design, TRUE dose semi-elasticity = 0.20 on the poisson design):

| leg | reported ame | truth | rel. err |
|---|---|---|---|
| `untreatedZero robust` (default reghdfe)         | 0.9325000000 | 0.9325 | 5.4e-16 |
| `untreatedZero ccnoabsorb robust`                | 0.9325000000 | 0.9325 | ~0      |
| `untreatedZero customvce(robust)`                | **-0.1075000000** | 0.9325 | **1.115** |
| `untreatedZero customvce(cluster unit)`          | **-0.1075000000** | 0.9325 | **1.115** |
| `untreatedZero unconditionalse robust`           | **-0.1075000000** | 0.9325 | **1.115** |
| `untreatedZero poisson robust`                   | **0.0187049681**  | 0.20   | **0.906** |
| `untreatedZero poisson ccnoabsorb robust`        | **0.0187049681**  | 0.20   | **0.906** |

The poisson case is the one I would worry about most in practice: `contreat()` with a dose that
is zero outside treatment and a count outcome is exactly what `untreatedZero` was written for
(sthlp:995 and the worked example at sthlp:1551), and the reported semi-elasticity comes back
an order of magnitude too small with no caution and rc 0. Note that `ccnoabsorb` does NOT
rescue the poisson path: under poisson the head is built and then, because `ccringertail` is
empty, it goes into `contreatcontrolinserterABS` too, so both spellings carry it.

### FINDING A — mechanism reproduced digit-for-digit outside wooldid (probe23_reghdfe_d.log)

On repro_A's own fixture, a bare reghdfe + margins with no wooldid involved:

    reghdfe y ib0.cell##c.dose , absorb(i.i i.t i.nullv##c.ringer) vce(robust) tol(1e-14)
    predict double xbh , xb
    summarize (y - xbh)        ->  mean  0.000000 ,  **sd 17.336133**
    margins , eydx(dose) subpop(if treated-post) noestimcheck nose
                               ->  **0.016412789649**

0.016412789649 is EXACTLY what wooldid reports for
`contreatelasticitytype(eydx)` on the default path in repro_A, against the closed-form truth
0.019020203082. The prediction margins divides by differs from the fitted mean by a quantity
with standard deviation 17.3 - the absorbed fixed effects - and the mechanism is therefore
certain, not inferred.

Remedy note: a bare `d` (which wooldid passes to ppmlhdfe) is REJECTED by reghdfe 6.12.3 with
r(198); reghdfe needs `savefe` / `d(newvar)` and a prediction built from it, so the fix is not a
one-token copy of the poisson line. Versions used throughout: wooldid 1.4 (pinned build,
sha c1d8979), reghdfe 6.12.3 08aug2023, ppmlhdfe 2.3.0 25feb2021, gtools 1.10.1, Stata 19.0 MP.

### FINDING E — the exact arithmetic, and what the user actually sees (probe24_uz_cautions.log)

The mechanism is not "some bias": it is an identity, and it is solver-order dependent.
The ringer column equals the dose everywhere, and under untreatedZero the dose is 0 off
treated-post, so the ringer lies exactly in the span of the cell#dose columns. On the
explicit-regressor path `regress` keeps the ringer (it is written first) and drops ONE cell
slope for collinearity - the last in the expansion. Writing beta_d for the pooled slope it
fits and beta_k for the true cell slopes, the fitted visible slopes are beta_k - beta_d and the
dropped cell reads 0, so

        reported ame  =  mean(beta)  -  beta_{the cell the collinearity screen dropped}

Measured, exactly: mean(beta) = 0.9325000000, and the dropped cell is (i = 2, relyr = 2) whose
true slope is 1.0400000000. 0.9325 - 1.04 = **-0.1075**, which is the reported number to every
digit. Which cell gets dropped is the engine's drop order, so the answer is not merely biased -
it is arbitrary.

On the poisson/absorb path the arithmetic is different and worse: the pooled slope is partialled
out FIRST, so the visible cell slopes estimate deviations from the pooled dose response. Under a
HOMOGENEOUS true effect (the poisson DGP has beta = 0.20 on every treated post row) those
deviations are zero, so the predicted reported ame is ZERO regardless of the true effect. The
measurement is 0.0187 with se 0.0752 and p = 0.803 - i.e. a real semi-elasticity of 0.20 comes
back as an insignificant zero.

WHAT THE USER SEES (each leg run noisily):
- `customvce(robust)`: a GENERIC caution - "a treatment period coefficient interacted with
  something (likely, a term implementing a continuous treatment variable or an interactive
  control variable) was dropped for multicollinearity ... This could imply your results are not
  identified, depending on exactly what was dropped" - and then `ame = -.1075` with
  **se = 3.216e-14, t = -3.3e+12, p = 0.000**. The caution does not name untreatedZero, does not
  say a baseline dose control was inserted, and the row it sits above looks infinitely precise.
- `poisson robust`: the ppmlhdfe rc-3200 fallback caution, then the same generic
  multicollinearity caution, then `ame = .01870497 (se .0752, p .803)`.
- default reghdfe `robust`: `ame = .9325` exactly, and NO caution at all.
So the failure is observable but mis-attributed, and on the linear path it is dressed in a
standard error of 3e-14.
