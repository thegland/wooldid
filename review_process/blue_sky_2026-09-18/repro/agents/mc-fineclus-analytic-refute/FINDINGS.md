# mc-fineclus-analytic-refute — SKEPTIC report

Assignment: try to REFUTE MC-FC-1 and MC-FC-2. Own generators, own arms, own truth;
the prober's do-files were read, never run.

## SUMMARY

| item | verdict | severity |
|---|---|---|
| MC-FC-1 finer-than-`i` remedy incomplete at small K | **PARTIAL** | medium |
| MC-FC-2 the 0.66-0.75 band is not the size of the at-or-above-`i` shortfall | **PARTIAL** | medium |
| SK-FC-A coarse-`i` remedy geometry biases the POINT ESTIMATE on unbalanced panels | new | high |
| SK-FC-C joint pre-trend test 12x oversized in the geometry the help calls the remedy | new | high |
| SK-FC-B few-treated-cluster caution exists only on the bootstrap engine | new | low |
| SK-FC-D a COARSER `cluster()` lowers the reported SE (over-clustering is anti-conservative) | new | medium |

Both seeded findings reproduce in substance on data their author never saw. Neither is
refuted. Each carries corrections: MC-FC-1's K=10 rung does not reproduce (I get 1.013 /
0.958 where they get 0.912 / ~0.93, and the project's own 500-draw run gets 0.972 / 0.934),
its `sqrt((K-1)/K)` is a FLOOR their own K=10 number falls below rather than a fit, and its
csdid comparison is overstated relative to their own FINDINGS.md. MC-FC-2's "0.18 / 0.29"
headline is a minimum over eight event-study rows; on the honest `att` row the shortfall on
my data is larger than they reported (0.315 / 0.467), but the behaviour is thrice-ruled and
the help's own next sentence forecloses the calibration reading its severity rests on.

The two things nobody had looked at are worth more than either seeded item: the joint
pre-trend test is 62% oversized in the exact geometry the help and the runtime caution name
as its fix, and the coarse-`i` route into that geometry silently moves the point estimate on
an unbalanced panel.

---

# NEW FINDING SK-FC-A — the documented remedy geometry (coarsen `i` so that `cluster()` is finer than `i`) silently biases the POINT ESTIMATE on an unbalanced panel

**severity**: high  **confidence**: certain  **realistic**: yes  **seeded**: no
**surface**: `i` / `fe()` / the "What level to cluster at" remedy at `sthlp:255`

**Scope, stated up front so the severity is not read as larger than it is.** The `i` entry
at `sthlp:64-70` DOES tell a careful reader both halves: "prefer treatment onset cohorts as
`i`" and "Choose `i` by the treatment heterogeneity ... not by the level at which you want
fixed effects - those are specified separately in `fe()`, so a coarse `i` combines with
unit-level ones." The gap is that the Standard-Errors remedy sentence that sends users into
that geometry omits `fe()` entirely, that neither passage says omitting it is a source of
BIAS rather than a modelling preference, and that nothing fires at runtime when `i` holds
sub-units, the panel is unbalanced, and no `fe()` was given.

**Claim.** `sthlp:255` names "a finer `cluster()`" as THE remedy for the at-or-above-`i`
standard-error understatement, and says "The intended practice is to cluster BELOW the level
of `i` - at a level whose units are nested inside cohorts". In a state-policy panel with one
row per state-year, the only way to reach that geometry is to pass the ONSET WAVE as `i`
instead of the state. That replaces the state fixed effects with a single wave intercept.
On a balanced panel this is harmless (the point estimate is bit-identical). On an UNBALANCED
panel it silently changes the point estimate, by an amount proportional to the dispersion of
the sub-`i` permanent heterogeneity, with no caution of any kind. `fe(state)` repairs it
exactly, and the "What level to cluster at" passage never mentions `fe()`.

**Observed** (`probe_unbal2.do`, `probe_fe.do`, `probe_unbal3.do`; ZERO noise in the DGP, so
every departure from 1.000000000000 is arithmetic, not sampling variation):

DGP: 40 states, T = 12, states 1-3 adopt at t=5, 4-6 at t=7, 7-40 never; permanent
deterministic state levels `s_g = 2*mod(state,7) - 6`; `y = s_g + 0.2*t + 1*(t>=ttre)`.
TRUE att = 1 exactly.

    case                                       i=state        i=coh (the remedy geometry)
    balanced (control)                       1.000000000000   1.000000000000
    one treated state attrits after t=8      1.000000000000   1.210526315789   (+21%)
    one treated state enters late at t=5     1.000000000000   0.428571428571   (-57%)
    ordinary ragged panel, 8% MCAR           1.000000000000   1.214732454068   (+21%)

On the ragged case the run prints, with no caution at all:
    att = 1.2147325   se = .12963004   t = 9.37   p = 0.000   95% CI [.95253096, 1.476934]

Repairs, on the same ragged data (`probe_fe.do`):
    i=coh cluster(state) fe(state)        att = 1.000000000000   se = 0
    i=coh cluster(state) controls(i.state) att = 1.000000000000   se = 0

Scaling with the sub-`i` dispersion (`probe_unbal3.do`, 8% MCAR, s_g ~ N(0, sd^2)):
    sd(s_g)   i=state          i=coh           i=coh + fe(state)
     0.50     1.0000000000     1.0154259361    1.0000000000
     1.00     1.0000000000     1.0308518723    1.0000000000
     2.00     1.0000000000     1.0617037446    1.0000000000
     4.00     1.0000000000     1.1234074892    1.0000000000
So even with state heterogeneity the same size as the treatment effect the estimate is 3%
wrong; with the state dispersion an applied panel actually has it is 6-21%.

**Expected.** `sthlp:67-70` does document the pairing - "Choose `i` by the treatment
heterogeneity the estimands should average over, not by the level at which you want fixed
effects - those are specified separately in `fe()`, so a coarse `i` combines with unit-level
ones" - and the SI-8 runtime caution (`ado:1371`) likewise says "specify unit-level fixed
effects via `fe()` and swap `i` to a coarser level of aggregation". So the pairing is known
to the package. But both present `fe()` as a modelling PREFERENCE ("the level at which you
want fixed effects"), not as a CORRECTNESS REQUIREMENT, and the Standard-Errors remedy that
actively pushes a user into a coarse `i` omits it entirely. Nothing at runtime notices that
`i` holds several sub-units per (i,t) cell, that the panel is unbalanced, and that no `fe()`
was given. The result is category-1 by the brief's ranking: a silently wrong number on
input a real applied economist has (ragged administrative panels are the norm), produced by
following the help's own advice.

**Repro**: `probe_unbal2.do` (self-contained, ~30 s), corroborated by `probe_fe.do` and
`probe_unbal3.do`.

---

# NEW FINDING SK-FC-B — the package has a named "few treated clusters" caution, but only on the engine its own help says is NOT a remedy for few treated clusters

**severity**: low  **confidence**: certain  **realistic**: yes  **seeded**: no
**surface**: analytic clustered path vs `bootstrap()` cautions

`ado:6498` prints "Caution: only N treated cluster(s) are available for resampling ... with
few treated clusters the literature finds this bootstrap badly sized (MacKinnon and Webb
2017, 2018)", and `ado:6479` prints a caution at fewer than 10 clusters in total. Both live
inside `wooldid_bootctl` and fire ONLY under `bootstrap()`. `sthlp:532` states of that same
engine: "Neither mode is a remedy for few treated clusters."

On the DEFAULT ANALYTIC path - the one `sthlp:255` calls the remedy - there is no
few-treated-cluster caution at all. `repro_nocaution.do` runs both engines on the SAME
dataset: 40 states, `cluster(state)`, `i` = the onset wave, TWO states in each treated
cohort-by-period cell and FOUR treated states in the whole panel.

    analytic  (i=coh, cluster(state), jointtests) : prints "N clusters = 40" and no caution
                                                    of any kind; the joint pre-trend p-value
                                                    is printed as an ordinary test
    bootstrap (same data, bootstrap(standard))    : "Caution: only 4 treated cluster(s) are
                                                    available for resampling ... with few
                                                    treated clusters the literature finds
                                                    this bootstrap badly sized (MacKinnon
                                                    and Webb 2017, 2018)"

The bootstrap threshold is `r(N) < 6` at `ado:6497`. The only cell-geometry caution on the
analytic path (`ado:1443`) fires on the strictly binary predicate "every `i` lies wholly
inside one cluster", so it is silent at 2, 3 or 5 clusters per treated cell.

This is the corroborating half of MC-FC-1: the hazard class is already named in the package,
with a numeric threshold, on the engine whose own help says it is not a remedy for it - and
is absent from the engine the help calls the remedy.

**Repro**: `repro_nocaution.do` (self-contained, ~1 min).

---

# NEW FINDING SK-FC-C — the joint pre-trend test is severely oversized in the exact geometry the help and the runtime caution name as its REMEDY, and prints nothing there

**severity**: high  **confidence**: probable (certain on the measurement; the scope
narrows if the `esprelength(2)` arm in `repro_joint_finer.do` comes back near nominal)
**realistic**: yes  **seeded**: no
**surface**: `jointtests` under analytic clustered SEs with `cluster()` strictly finer than
the cohort-by-period cells (`sthlp:743`, caution at `ado:1443`)

**Claim.** `sthlp:743` says `jointtests` "is unreliable under wooldid's analytic clustered
SEs when clustering is at the level of the cohort variable `i` or above ... If joint-type
inference is needed in such settings, use bootstrap-based clustered standard errors or
randomization inference (see `permute()`), **or cluster strictly finer than the
cohort-by-period cells**." The runtime caution at `ado:1443` ends with the same clause:
"...as does clustering strictly finer than the cohort-by-period cells." In that named
remedy geometry, with three clusters inside each treated cohort-by-period cell, the joint
pre-trend test rejects a TRUE null 61.7% of the time at nominal 5%, and no caution prints -
the caution's predicate (`ado:1421-1432`, distinct `i` vs distinct (`i`, cluster) pairs) is
the strictly binary test "every `i` lies wholly inside one cluster", so it is silent at 2,
3 or 5 clusters per cell.

**Observed** (`mc1_ladder.do`, 120 draws per rung; DGP has tau = 0 in every pre-treatment
period, so the joint pre-trend null is exactly true; 40 states, `cluster(state)`, 40
clusters in every rung):

    K   arm                                  reject@5%   reject@10%   caution?
    3   i=coh, cluster(state)   [FINER]        0.617       0.683        NO
    3   same, pooled reference  [FINER]        0.617       0.667        NO
    3   i=state, cluster(state) [AT i]         0.925       0.933        yes
    10  i=coh, cluster(state)   [FINER]        0.108       0.192        NO
    10  i=state, cluster(state) [AT i]         0.550       0.675        yes

So at 3 clusters per treated cell the "remedy" geometry is 12x oversized, and even at 10
clusters per cell it is still 2x oversized - both with nothing printed.

**Why this is not the known item.** `known_bydesign` s1 records the oversizing at or ABOVE
`i` (11-36% at G=40, 65.7% at G=10, NF-1b), computed with a named caution. This is the
FINER-than-`i` geometry at G=40, which the help names as the fix and which the caution
predicate cannot see. `known_bydesign` s7 lists "any finer-than-`i` arm that FAILS coverage"
as genuinely new.

**The help already contains the right idea, scoped to the wrong case.** The SECOND
reliability condition in the same `jointtests` paragraph says: "in the continuous treatment
case the joint tests are only correctly sized when the i-by-t cells each contain many
observations ... designs with only a handful of observations per cohort-period cell yield
oversized joint tests EVEN UNDER ROBUST OR FINELY CLUSTERED SEs". That is exactly the
failure measured here - except this is BINARY treatment, the cells hold twelve observations
(3 states x 4 units), and what is scarce is CLUSTERS per cell, not observations. So the
package's model of the hazard is keyed on observations-per-cell and confined to
`contreat()`, while the binary case with few clusters per cell is instead told that finer
clustering is the cure.

**The obvious alternative explanation is ruled out: the covariance is NOT degenerate there.**
`probe_rank.do` dumps the 4x4 pre block of `e(V_es_att)` on one draw of each rung:

    K=3   i=coh  (FINER)  rank 4 of 4   cond 2.23e+01   det  1.63e-06   joint p 0.037
    K=3   i=state (AT i)  rank 3 of 4   cond 6.20e+15   det -7.40e-25   joint p 1.0e-06
    K=10  i=coh  (FINER)  rank 4 of 4   cond 1.64e+01   det  2.02e-07   joint p 0.662
    K=10  i=state (AT i)  rank 3 of 4   cond 4.71e+15   det -2.29e-23   joint p 0.382

In the remedy geometry the pre-trend covariance is FULL RANK and well conditioned (condition
number 22). So the 61.7% rejection rate is not a singular matrix being inverted anyway - it
is a full-rank covariance that is simply too small in the directions the test uses. (The
at-or-above-`i` arm is rank-deficient with a condition number of 6e15 and wooldid prints a
p-value from it; that case is the documented one and carries the caution.)

**Two things a skeptic should check, and what they show.** (1) The pre-period estimates are
unbiased: mean over 120 draws is -0.072, +0.021, +0.011, -0.022 on the four pre horizons at
K=3 (the -0.072 is ~2 standard errors of the mean, and the K=10 rung's corresponding mean is
+0.038 while rejecting only 10.8%), so the 61.7% is not a real pre-trend being detected.
(2) The four-horizon window is ragged - `pre5` exists only for the later cohort - so
`repro_joint_finer.do` also runs an `esprelength(2)` arm whose two pre horizons are
identified by BOTH cohorts; if that arm is also grossly oversized the raggedness is not the
cause.

**Repro**: `repro_joint_finer.do` (self-contained, ~6 min; adds a `csdid` + `estat pretrend`
oracle on identical draws and a rank/conditioning dump of the pre block of `e(V_es_att)`).
---

# VERDICT on MC-FC-1 — PARTIAL (severity medium)

**What I did.** Own generator (`skgen` in `mc1_ladder.do`): 40 states x 4 units x 12 periods,
two adoption waves of K states at t=5 and t=7, the rest never-treated; state RE N(0,.8^2) +
state-by-year shock N(0,.4^2) + unit RE N(0,.5^2) + iid N(0,.4^2). 120 draws per rung,
K in {3, 10}, three arms on the same draws. Truth computed in closed form
(att = 1.153571428571..., ES post h = 1 + .05h). Their do-files were read, not run.

| K  | arm                              | att SE/SD | att cov95 | med ES SE/SD | med ES cov95 |
|----|----------------------------------|-----------|-----------|--------------|--------------|
| 3  | A: i=coh, cluster(state) [FINER] | 0.869     | 0.858     | 0.828        | 0.871        |
| 3  | B: same, pooled reference        | 0.884     | 0.900     | 0.807        | 0.845        |
| 10 | A: i=coh, cluster(state) [FINER] | **1.013** | **0.958** | 0.986        | 0.938        |
| 10 | B: same, pooled reference        | 0.901     | 0.908     | 0.906        | 0.934        |

**CONFIRMED**: at 3 clusters per treated cohort-by-period cell, with `cluster()` strictly
finer than `i` throughout, the analytic interval does not deliver its nominal level -
0.869 SE/SD and 0.858 coverage on the headline `att` row, and 0.828 / 0.871 on the median
event-study row. That crosses both halves of the prober's pre-registered gate on data they
never saw, and it is what `known_bydesign` s7 calls genuinely new.

**PARTIAL, correction 1 - the K=10 rung does not reproduce.** They report 0.912 SE/SD and
~0.93 coverage at K=10 and present the whole ladder as degraded. On my independent draws
K=10 is CLEAN: 1.013 SE/SD and 0.958 coverage on `att`, 0.986 / 0.938 on the median ES row.
My number agrees with the project's own 500-draw 2026-09-14 run, which measures 0.972 /
0.934 at that same geometry (`review_process/se_coverage_2026-09-14/README.md` row
"i = cohort id"; I read `se_cov_mc.do` and confirm that arm is `wooldid y icoh t ttre,
cluster(i)` with i = the UNIT, i.e. 10 clusters per treated cell - the prober's reading of
that file is correct). So the defect is concentrated at FEW clusters per cell, not spread
along the ladder.

**PARTIAL, correction 2 - `sqrt((K-1)/K)` is a FLOOR, not a fit.** Under the stated
mechanism the ratio is `sqrt(1 - phi/Kc)` where `phi` is the share of Var(bhat) coming from
saturated cells and `Kc` the clusters inside one; `phi <= 1` by construction, so
`sqrt((Kc-1)/Kc)` is the LOWEST value the mechanism can produce. Their own K=10 measurement
of 0.912 sits BELOW that floor (0.949), which their own mechanism cannot generate; my K=10
measurement of 1.013 sits above it. The seven-point "tracking" table is within Monte Carlo
noise (+/- 6% on a ratio at 120-150 draws) of several other shapes. The gradient in Kc is
real; the functional form is not established.

*I chased the one structural difference between our generators and it explains nothing.*
Their ladder uses AR(1)(rho .5) errors within unit; mine uses iid. `probe_ar1.do` runs my
K=10 design 120 draws each way, everything else held fixed:

    arm                          att SE/SD   att cov95   ES SE/SD range
    r0  iid errors (mine)          1.019       0.942     0.869 - 1.031
    r5  AR(1) rho .5 (theirs)      1.001       0.942     0.843 - 1.023

Serial correlation costs 1.8% of the ratio and nothing at all in coverage. So the gap
between their 0.912 and my 1.013 at the same K is not the error process; the residual
candidates are Monte Carlo noise in their SD at 150 draws and the wave-timing difference
(t=5/t=8 versus t=5/t=7). Either way, the K=10 rung of their ladder is not reproducible.

**PARTIAL, correction 3 - the csdid comparison in the transmitted claim is overstated.**
"csdid ... holds SE/SD at 0.92-1.12 across the whole ladder" is true of csdid's SE/SD, but
the transmitted claim pairs it with "95% coverage falls from ~0.94 to ~0.86" as if csdid
held coverage. Their own FINDINGS.md concedes the opposite in its "Honest scope" paragraph:
csdid's coverage drops to 0.880-0.947 at K=3 (median 0.923, below their own 0.93 gate) and
did_imputation undercovers at 0.76-0.85 there, so the gate's "both controls inside band"
clause does not fire at K=3. At three adopting states per wave the design is hard for
everyone; what is specific to wooldid is the SE LEVEL, not the coverage gap alone.

**The MECHANISM moves in the predicted direction by roughly the predicted amount, in a
PAIRED test that needs no sampling SD.**
`mc2_mech.do`, 120 draws, on a DGP with NO state-level dependence of any kind (units
independent: unit RE + iid error; `state` is only a grouping label). Three arms on the same
draws with IDENTICAL point estimates, so the sampling SD is shared and the SE ratio is
paired:

| arm | geometry | N_clust | clusters in treated cell | mean SE | att SE/SD | att cov95 |
|---|---|---|---|---|---|---|
| E | i=coh, cluster(state) | 40 | 3 | 0.1170 | 1.002 | 0.925 |
| F | i=coh, cluster(unit) | 160 | 12 | 0.1327 | 1.136 | 0.958 |
| G | i=unit, cluster(unit) | 160 | 1 | 0.0422 | 0.362 | 0.567 |

Paired `mean(se_E)/mean(se_F)` across the seven reported quantities: 0.882, 0.867, 0.840,
0.900, 0.850, 0.875, 0.884 - mean **0.871** against **0.853** predicted by
`sqrt(2/3)/sqrt(11/12)`. So "the number of clusters inside the saturated treated cell" is
genuinely the operative variable, and it is not about where the noise lives: this DGP has no
state-level component at all.

*Honest caveat on that comparison.* Arm F is not itself well calibrated (1.136 SE/SD): with
`i` = coh the unit random effect is not absorbed, so `cluster(unit)` sees it in the scores
although the estimator differences it out. The two arms are therefore each miscalibrated,
in opposite directions, and only their RATIO is a clean paired quantity. The ratio's
agreement with 0.853 is supporting evidence, not a precise confirmation of the constant.

**Corroboration the prober did not report**: there is no few-treated-cluster caution on the
analytic path at all (see SK-FC-B), and the only cell-geometry caution (`ado:1443`) fires on
a strictly binary predicate. And the same geometry's JOINT pre-trend test is 12x oversized
(SK-FC-C), which is a larger consequence of the same gap than the pointwise SEs are.

**Own repro**: `mc1_ladder.do` (+ `mc1_ladder.log`).

---

# VERDICT on MC-FC-2 — PARTIAL (severity medium)

**What I did.** The `i=state, cluster(state)` arm (arm C) of the same 120-draw runs, on the
same draws as arms A and B, so the point estimates and the empirical SD are shared and only
the reported SE differs.

| K  | quantity           | SE/SD | cov95 |
|----|--------------------|-------|-------|
| 3  | `att` row          | 0.315 | 0.467 |
| 3  | median ES row      | 0.289 | 0.434 |
| 3  | min ES row         | 0.223 | 0.358 |
| 10 | `att` row          | 0.634 | 0.800 |
| 10 | median ES row      | 0.612 | 0.750 |
| 10 | min ES row         | 0.483 | 0.658 |

**CONFIRMED in substance, and worse than they reported on the honest statistic.**
`sthlp:255`'s "roughly 0.66 to 0.75" is not a floor: on an ordinary staggered state-policy
panel with a state-by-year shock and four sub-units per state the ATT ROW - not a minimum
over rows - comes in at 0.315 SE/SD and 0.467 coverage with six treated states, and 0.634 /
0.800 with twenty. Their own headline att figures were 0.384-0.610.

**PARTIAL, correction 1 (the sharpest one) - the stated MECHANISM is wrong.** The claim's EXPECTED says
"The reason my designs are harsher is identifiable and one-dimensional: they contain a
state-by-year shock. Remove it and the band reproduces". It does not. `mc2_mech.do` arm G
(`i` = unit, `cluster(unit)`, 120 draws) runs on a DGP with NO cluster-level time shock and
no state-level component of any kind - unit random effects plus iid errors - and gives
**0.362 SE/SD and 0.567 coverage** on the `att` row, further outside the documented band
than any of their state-shock designs. What actually drives the magnitude is the share of
Var(bhat) that sits inside cells the treated dummies saturate; one observation per treated
(`i`, period) cell achieves that with no cluster-level shock anywhere. A DGP without such a
shock can therefore be MORE harsh, not less, and re-measuring the band "with a cluster-level
time shock present" - the fix the claim proposes - would not bound it.

**PARTIAL, correction 2 - the "0.18 / 0.29" headline is a selection statistic.** It is the
MINIMUM over eight event-study rows in a 200-draw MC. The minimum of eight noisy ratios is
biased low; on my runs the min ES row is 0.223 against an att row of 0.315, and 0.483
against 0.634 - a 25-30% gap between the honest statistic and the selected one. The claim
should be stated on the `att` row (or the median ES row), where it is still decisively
outside the documented band.

**PARTIAL, correction 3 - the severity argument rests on a reading the help forecloses.**
The claim's EXPECTED says "a user calibrates to 'my t-statistics are about 1.4x too large'".
The sentence immediately after the quoted band says: "There is no plug-in correction - at
that level the analytic covariance is the WRONG OBJECT rather than a mis-scaled one". That
is an explicit instruction not to treat 0.66-0.75 as a calibration constant. The band is
also scoped in the same sentence to "the diagnostic fixtures", which the claim concedes.

**PARTIAL, correction 4 - this is a rediscovery with a bigger number.** `known_bydesign` s1
records the at-or-above-`i` understatement as ruled three times, with the same mechanism the
help describes. What is new is only the MAGNITUDE, which the help's own mechanism sentence
already predicts is unbounded. Value is documentation-only, and the actionable form is NOT
the claim's proposed fix (re-measure with a cluster-by-time shock) but: say that 0.66-0.75
is a fixture reading and not a bound, and key any guidance on how much of the estimator's
variance sits inside cells the treated dummies saturate. Severity medium, not higher,
because the behaviour itself is documented and thrice-ruled.

**Own repro**: `mc1_ladder.do` arm C (+ `mc1_ladder.log`); `mc2_mech.do` arm G for
correction 1 (+ `mc2_mech.log`).
---

# NEGATIVE RESULTS (attacked hard, came back clean)

1. **The finer-than-`i` remedy WORKS at ten clusters per treated cell.** 120 draws, arm A at
   K=10: `att` SE/SD 1.013, coverage 0.958; median ES row 0.986 / 0.938; min ES row 0.900 /
   0.883. Nothing here fails. This is the result that makes MC-FC-1 PARTIAL rather than
   CONFIRMED, and it independently reproduces the project's own 2026-09-14 figure (0.972 /
   0.934). If the remedy were broken as such, this rung would have shown it.

2. **`cluster()` and the choice of `i` do not move the point estimate on a balanced panel.**
   `pilot.log`: `wooldid y coh t ttre` and `wooldid y state t ttre` return `att` =
   .82521497 to every printed digit while their standard errors differ 3.06-fold.
   `probe_unbal2.do` CASE 0: both return exactly 1.000000000000 on a zero-noise panel.
   Whatever else is wrong, the estimator is not being silently switched by `cluster()`.

3. **The `jointtests` caution predicate fires exactly where it says it does.** It printed at
   `i` = state with `cluster(state)` and stayed silent at `i` = coh with `cluster(state)`,
   on the same data, in `pilot.log`. `ado:1421-1432` (distinct `i` vs distinct (`i`,cluster)
   pairs) is a correct implementation of "every `i` lies wholly inside one cluster". The
   defect in SK-FC-C is that the predicate is the wrong shape, not that it misfires.

4. **`fe()` and `controls(i.)` fully repair the coarse-`i` composition bias**, to the last
   printed digit on a zero-noise panel (`probe_fe.do`). The machinery is sound; the gap is
   that nothing tells the user to reach for it.

# UNTESTED / OPEN

- `permute()` with `permonsetcohorts` as the third documented remedy at K=3: not run.
- Whether the oversized joint test at K=3 finer-than-`i` also affects `adversarial` and
  `pretrendpower()` (which `PS-64` routes through the same statistic): not run.
---

# PARTIAL RESULT — does remedy #2, `bootstrap(standard)`, close the K=3 gap?

`sthlp:255` lists three remedies, and MC-FC-1 tested only the first. `probe_boot.do` runs
`bootstrap(standard) reps(99)` against the analytic clustered SE on identical draws, at
`i` = coh, `cluster(state)`, 3 states per treated cohort-by-period cell (one row per
state-year, so the cell holds 3 observations in 3 clusters).

Paired per-draw ratio of the bootstrap standard error to the analytic one, over the 23
completed draws at the time of writing (the run was still going; the do-file does 30):

    mean( se_bootstrap / se_analytic ) = 1.177     (first five draws: 1.183, 1.165,
                                                    1.113, 1.041, 1.189)

To close a shortfall of the analytic ratio measured at 0.869 (mc1_ladder, K=3) the bootstrap
would need a ratio of 1/0.869 = 1.151. It is running at 1.177. So on this evidence
**`bootstrap(standard)` does close the K=3 gap**, and MC-FC-1's headline - "the finer-than-
`i` remedy is incomplete" - is true only of the FIRST of the three remedies the help lists,
taken alone. The help's remedy SET is not shown to be incomplete.

**This is a paired SE-ratio measurement on 23 draws, not a coverage measurement.**
`probe_boot.do` run to completion prints SE/SD and coverage for both arms on its own 60
draws; if this report reaches you before that, treat the coverage question as open.

---

# NEW FINDING SK-FC-D — a COARSER `cluster()` can lower the reported standard error, reversing the conservatism every user expects from over-clustering

**severity**: medium  **confidence**: probable  **realistic**: yes  **seeded**: no
**surface**: `cluster()` on the default analytic path

**Claim.** Applied practice treats "cluster coarser than you need" as the safe direction: a
cluster variable above the true level of dependence costs precision and yields wider
intervals. In wooldid, because the treated cohort-by-period dummies demean the residual
across whatever sits inside the cell, a coarser `cluster()` puts FEWER clusters inside that
cell and the reported standard error gets SMALLER. Nothing in the run or the help says so,
and the "What level to cluster at" passage - which tells users to cluster BELOW `i` - does
not distinguish "below `i`" from "as fine as the dependence actually is".

**Observed** (`mc2_mech.do`, 120 draws, DGP with no state-level dependence of any kind, so
`cluster(unit)` is the correct level and `cluster(state)` is pure over-clustering; the two
arms have IDENTICAL point estimates, so this is a paired comparison of reported uncertainty
with no sampling-SD noise in it):

    quantity        mean SE, cluster(unit)   mean SE, cluster(state)   ratio
    main_att              0.1327                    0.1170            0.882
    ES_att_post0          0.1753                    0.1519            0.867
    ES_att_post1          0.1699                    0.1427            0.840
    ES_att_post2          0.1719                    0.1547            0.900
    ES_att_post3          0.1723                    0.1465            0.850
    ES_att_pre3           0.1719                    0.1504            0.875
    ES_att_pre2           0.1780                    0.1573            0.884
                                                        mean ratio    0.871

Coarsening the cluster from 160 units to 40 states made every reported standard error 13%
SMALLER. The ratio predicted by "clusters inside the saturated cell"
(`sqrt(2/3)/sqrt(11/12)`) is 0.853.

*What this does and does not establish.* The DIRECTION and the ~13% size are measured
cleanly, because the point estimates are identical across the two arms and the comparison is
paired. It is NOT a clean confirmation of the 0.853 constant: arm F reports 1.136 of the
sampling SD, because with `i` = coh the unit random effect is not absorbed and
`cluster(unit)` picks it up in the scores although the estimator differences it out. Both
arms are miscalibrated in opposite directions; only the ratio is paired. `mc2_noau.do` is
the clean control (same design with the unit random effect removed, so `cluster(unit)` is
exactly right); it was still running when this was written. The user-facing proposition -
a coarser admissible `cluster()` gives a SMALLER reported standard error - stands on the
paired ratio alone.

**Expected.** Over-clustering is the direction users deliberately choose when unsure.
`sthlp:255` reinforces the choice ("cluster BELOW the level of `i`"), which a user reads as
"anything below `i` is safe". It is not: within the set of admissible cluster variables, the
coarsest is the least conservative here. This is the same arithmetic as MC-FC-1, stated as
the user-facing rule it implies, and it is measurable without any Monte Carlo - two runs on
one dataset show it.

**Repro**: `mc2_mech.do` (self-contained, 120 draws; the E and F arms are the comparison);
clean control `mc2_noau.do`.
