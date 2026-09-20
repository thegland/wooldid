# SKEPTIC verdicts on the fp-contreat findings

Every number below comes from a probe I wrote, on data I generated, with oracles I built.
The prober's do-files were read once for orientation and never executed.
All probes live in this directory with their `.log` files.

---

## fp-contreat-E — CONFIRMED (critical), and WORSE than claimed
repro: `R2_E_untreatedZero.do`, follow-up `R8_E_att_and_poisson.do`

Noisy panel (6 cohort groups x 12 units x 10 periods, 3 never-treated groups, dose == 0 off
treated cohort-periods, heterogeneous beta_it, N(0,0.5) noise), not the prober's exact-fit design.

| | ame |
|---|---|
| ORACLE A (`reghdfe y ib0.cell ib0.cell#c.dose, absorb(i.i i.t)`) — the intended model | 1.0162389082 |
| ORACLE B (`reg y i.i i.t c.dose ib0.cell ib0.cell#c.dose`) — pooled slope in front | -0.6495711425 |
| wooldid `robust`, `ccnoabsorb robust`, `cluster(unit)` | 1.0162389082 |
| wooldid `customvce(robust)`, `customvce(cluster unit)`, `unconditionalse robust` | -0.6495711425 |

The pooled dose slope in ORACLE B is 1.6658100507, which is exactly ORACLE A's slope for the
cell the collinearity screen drops (cell 15), and
`ORACLE_A_ame - b_dropped = 1.0162389082 - 1.6658100507 = -0.6495711425` to every digit.
DGP truth 0.916; the customvce se is 0.29-0.32, so the wrong answer is ~5 sigma off with the
sign flipped, and on noisy data the se looks perfectly ordinary (the prober's se of 3e-14 was
an artifact of their exact fit, which makes the defect LESS visible on real data, not more).

Scope correction, upward: the prober measured only the `ame` row. The att row moves too
(`R8`): DGP true ATT 2.5494633489; default path 2.6743762849 (se 0.180); `customvce(robust)` and
`unconditionalse robust` both -0.5331814282 (se 0.633) — 4.9 sigma and sign-flipped. So the
headline treatment effect, not just the marginal-effect row, is corrupted.

Mechanism correction on the poisson leg. The prober says the pooled slope "is partialled out
first, so the visible cell slopes estimate deviations from it and the reported semi-elasticity is
driven to zero regardless of the truth." That is the wrong mechanism. What actually happens
(`R2`, `R9`) is that ppmlhdfe fails r(3200) on wooldid's own one-level absorbed term and the
rc-3200 fallback moves the control into the main equation, where it is subject to the SAME
collinearity-drop mechanism as the reg path. On a heterogeneous-effect poisson DGP:
wooldid = -0.5506611369 == my ORACLE-POIS B (-0.5506611369); ORACLE-POIS A = 0.0454169223;
DGP truth 0.1587. Not zero. "Driven to zero" is only the special case of homogeneous effects.

Nothing catches it: `sthlp:906` states that the structural-zero ame screen is deliberately not
armed under untreatedZero ("expected zero doses in untreated periods ... do not trigger this").

---

## fp-contreat-A — CONFIRMED (high)
repro: `R4_spec_diag.do` (mechanism), `R3_A_elasticity.do` (wooldid legs), `R5_A_magnitude.do`,
`R5b_A_level.do` (magnitude)

`R4` pins the exact spellings outside wooldid:
`reghdfe y ib0.cell##c.dose, absorb(i.i i.t i.nullv##c.ringer)` reproduces wooldid's default path
to every digit, and `reg y i.i i.t i.nullv##c.ringer ib0.cell##c.dose` reproduces the customvce path.

* `dydx` is identical across engines (1.022234373986 on both, and on wooldid both ways) — the
  two engines fit the same model to the same data, so any elasticity gap is purely the denominator.
* `eydx`: 0.034876727614 (reghdfe/default) vs 0.034066172921 (regress/customvce) vs
  0.034803487530 (`ccnoabsorb robust`, a THIRD answer: ringer moves into the main equation so
  it enters xb, but i.i and i.t still do not). `eyex`: 0.073136451868 / 0.071262964505 /
  0.072724096817.
* The closed form that the help defines — mean over treated post rows of b_cell / the FULL
  prediction — is 0.034066172921 for eydx and 0.071262964505 for eyex, i.e. **exactly** the
  customvce/regress numbers and not the default-path ones. (`R3`, "CLOSED FORM" lines.)
* reghdfe's `predict, xb` has the same mean as y (27.395972 both) but a smaller sd
  (1.262981 vs regress 1.943938 vs y 2.626373). The absorbed effects are missing from the
  denominator, so the error is second-order in the dispersion of the absorbed part, not the
  first-order "percentage of a partialled-out residual index" the prober's wording implies.

Magnitude correction. The prober's "13.7% too small" is design-specific, not the rate.
`R5`/`R5b`: 2-4% on ordinary panels, rising as E[y] falls toward the dispersion of the absorbed
effects (eydx reghdfe 0.0215677240 at ymean 47 vs 0.0503237247 at ymean 13, against the regress
twin each time). It is unbounded only when the partial index approaches zero on some rows.

The binary `semielasticity` extension is CONFIRMED, though the prober's evidence for it was
self-contradictory (their gap shrank with their scale knob while the contreat gap grew, which
the shared-mechanism story cannot produce). On my data the gap is there and has the right sign:
`semielasticity` att = 0.087622452771 (robust) vs 0.085833510361 (customvce(robust)), 2.1%,
the same size and direction as the contreat-side eydx gap on the identical fixture.

Nothing in `sthlp:245-249` (semielasticity), `sthlp:905-907` (contreatelasticitytype) or
`known_bydesign.md` says the linear-path denominator omits the absorbed effects.

---

## fp-contreat-B1 — CONFIRMED (high), and more deterministic than claimed
repro: `R7_B_same_data_mc.do` part (a); fixture also in `R6_B_lattice.do`

My own fixture: 12 states x 5-12 counties x 10 years, county-level dose plus a
within-county time-variation knob, `i` = state, no `cluster()` anywhere, bare `robust`.

Six replications, each engine run on the very same saved dataset:

| rep | wooldid default ame / se | wooldid customvce ame / se | reghdfe max\|e(V)\| / rank / rss | regress max\|e(V)\| / rank / rss |
|---|---|---|---|---|
| 1 | 0.223331 / . | 0.223331 / finite | 0 / 230 / 1345.1673 | 645.59 / 177 / 1345.1673 |
| 2 | 1.365626 / . | 1.365626 / finite | 0 / 224 / 1349.8185 | 225.76 / 174 / 1349.8185 |
| 3 | 1.778866 / . | 1.778866 / finite | 0 / 229 / 1580.5488 | 108.98 / 177 / 1580.5488 |
| 4 | 1.535433 / . | 1.535433 / finite | 0 / 221 / 1490.2072 | 532.34 / 182 / 1490.2072 |
| 5 | 0.498829 / . | 0.498829 / finite | 0 / 226 / 1476.1284 | 222.51 / 190 / 1476.1284 |
| 6 | 15.635699 / . | 15.635699 / finite | 0 / 227 / 1481.6840 | 395.42 / 177 / 1481.6840 |

6/6, not "on some designs". `margins, dydx(dose)` on my own outside-wooldid fits returns the
same ame under BOTH parameterizations to six decimals (e.g. 0.223331 and 0.223331), which is
positive evidence that the ame is estimable — invariant to which collinear columns the
screen happens to drop — so the regress covariance is a legitimate object for it and the
reghdfe zero is a recoverable engine artifact, exactly as claimed.

Mechanism refinement (the claim stands, the framing is loose). The design IS rank-deficient:
reghdfe retains 40-55 more columns than regress drops (rank 221-230 vs 174-190) on an identical
fit. So the proximate cause is engine-side, but "not of the design" overstates it — the design
supplies the near-singularity and the engines differ in what they do with it. Also, the claim's
"identical N, df_r, rss" is off by one on df_r in my runs (854 vs 855 in `R6`); N, rss and the
point estimates are identical.

---

## fp-contreat-B2 — PARTIAL (high): the point-miss story is wrong, the inference story is worse
repro: `R7_B_same_data_mc.do` part (b)

12 independent draws from the same DGP (true dose effect +0.4 on every treated post row),
`lattice(all 3)`, `customvce(robust)`:

* ame: mean -0.342, SD 2.519, range [-6.091, +2.906]
* reported robust se: mean 0.341, range [0.305, 0.434]
* no-lattice control arm: ame mean 0.462, SD 0.146, reported se mean 0.110

Correction: there is no detectable BIAS at n=12. The mean is 1.0 standard errors of the mean
from the truth; the "45x the truth" observation is one draw from a very heavy-tailed sampling
distribution, and my own draws produced 0.919, 3.200, 0.148, 15.636, -6.091 and +2.906 on the
same DGP. So "the estimator is consistent ... so 18.09 against 0.4 is a miss" mis-diagnoses a
variance problem as a point-estimate problem.

But the finding is worse where it matters. The realized sampling SD is 2.519 while the
reported conditional robust SE averages 0.341 — a 7.4x understatement, so nominal 95%
coverage on this design is close to zero. The prober flagged exactly this ("se 0.44 on a 45x
miss") and marked it uncertain; the Monte Carlo settles it. The no-lattice arm understates by
only 1.3x, so the lattice is the aggravating factor, and the help's lattice caveat
(`sthlp:1010-1018`) is about the effect loading onto bin INTERCEPTS, which is not this.
The 1.3x no-lattice arm is essentially `known_bydesign.md` section 7's NF-3 territory
(ame-side inference oversized under `robust`, remedy R5 not implemented); the 7.4x figure is a
lattice-specific sharpening of that known-open item, not a rediscovery of it.

---

## fp-contreat-C — CONFIRMED (medium)
repro: `R6_B_lattice.do`, the "NOISY leg" section

Under `lattice(all 3) robust` with no `cluster()` token anywhere in the command, verbatim:

> Caution: estimation returned an all-zero variance-covariance matrix. The issue is likely
> related to continuous treatment effect estimation coupled with use of analytic clustered
> SEs, with clusters occuring at a high level of aggregation. Consider switching to use of
> bootstrap SEs, IM p-values, or randomization inference if finer clustering levels are
> non-viable.

immediately followed by

> ... This is a property of the specification, not a numerical accident: see the
> cell-containment caution above for the mechanism and for the alternatives (impvalues
> dispersion inference, a finer-than-cell clustering, or conditional robust standard errors).

No cell-containment caution appears anywhere above it (the only prior caution is the
fastcollapse one). One point the prober missed: the second caution's third remedy is
"conditional robust standard errors", which is precisely what the user already asked for with
`robust`. All three pieces of advice are unreachable on this run.

---

## fp-contreat-D — REFUTED (severity none) — documented behaviour
repro: `R8_E_att_and_poisson.do` part (iii); source at `wooldid.sthlp` lines 1051-1053

The help documents this exactly, with the same example:

> This holds on the default estimation path. Under `customvce()` or `unconditionalse`,
> estimation runs through `regress` ... and a BARE token there is read as a CONTINUOUS
> variable: `fe(state)` fits a linear term in the state code rather than a set of state
> indicators. Write the prefix explicitly - `fe(i.state)` - whenever those options are in use,
> or on any specification you may later run both ways.

And the runtime caution (ado:4637-4638) ends "...manually specifying any needed i. prefixes".
The documented remedy works: `fe(i.unit) customvce(cluster unit)` gives R2 = 0.960850 while
`fe(unit) customvce(cluster unit)` gives R2 = 0.725469 on the same data. The finding's stated
expectation — "wooldid already knows at parse time that this combination is invalid — it says
so" — is the opposite of what the help says: the combination is valid, the prefix is the user's
responsibility, and erroring would remove a documented capability. Residual: the caution's
phrase "specify anything in fe as controls" is looser than the help's "write fe(i.state)".
A wording nit, not a defect.

---

## NEW (mine): poisson + contreat() always takes the ppmlhdfe rc-3200 fallback
repro: `R9_poisson_fallback.do`

`contreatcontrols(basic)` and `contreatcontrols(untreatedZero)` both leave
`i.<a variable identically equal to 1> ## c.(ringer)` inside `absorb()` (ado:1234-1243, and
ado:1301-1303 clears that head only when a byX made `ccringertail` non-empty). ppmlhdfe fails
r(3200) on the one-level factor every time, and the user is told:

> Caution: initial estimation of ppmlhdfe failed, possibly due to inclusion of an absorbed
> categorical fixed effect with only 1 level ... If estimation still fails, the issue could be
> either with the fe() specified or a result of having few cohorts/time-periods.

The user specified no `fe()` and no one-level fixed effect: the one-level fixed effect is
wooldid's own. Two controls confirm the cause: `contreatcontrols(basic byTime)` (which clears
the one-level head) does NOT fire it, and a plain binary `poisson` run (no `contreat()`) does
NOT fire it. Consequences: a caution that sends the user to debug their own specification for
wooldid's construction; the documented performance penalty ("computationally less efficient")
on every poisson continuous-treatment run; and — the reason it matters here — it routes every
such run onto the explicit-regressor path whose collinearity drop order is the mechanism of
finding E.
