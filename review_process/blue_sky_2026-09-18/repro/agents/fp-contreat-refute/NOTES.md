# SKEPTIC notes — refuting fp-contreat findings A, B1, B2, C, E

All probes are my own; the prober's do-files were read once for orientation and never run.

## Verdict summary (as of the last probe)

| id | verdict | severity revised | evidence |
|----|---------|------------------|----------|
| fp-contreat-E  | CONFIRMED (mechanism PARTIAL on the poisson leg) | critical | R2 |
| fp-contreat-A  | CONFIRMED | high (magnitude design-dependent) | R3, R4, R5, R5b |
| fp-contreat-B1 | CONFIRMED, stronger than claimed | high | R6, R7 |
| fp-contreat-B2 | PARTIAL — variance not bias; the SE understatement is the real defect and is worse than claimed | high (was high/uncertain) | R7 |
| fp-contreat-C  | CONFIRMED | medium | R6 |
| fp-contreat-D  | REFUTED — documented at sthlp:1051-1053 | none | code read |

## R1 (discarded oracle, kept for the record)
`ib0.cell##c.dose` accidentally carries the pooled `c.dose` main effect, i.e. it IS the
buggy model. Useful by accident: reghdfe with the pooled dose gave exactly the customvce
number (-0.6495711425), which is itself independent evidence for E's mechanism.

## R2 — finding E, my own data + my own oracles
noisy panel, 6 cohort groups x 12 units x 10 periods, 3 never-treated groups, dose == 0
off treated cells, heterogeneous beta_it, N(0,0.5) noise.
- ORACLE A (intended model, no pooled dose slope): ame = 1.0162389082
- ORACLE B (pooled dose slope pasted in front, so one cell slope is dropped): -0.6495711425
- pooled dose slope = 1.6658100507 == ORACLE A slope of the dropped cell (cell 15)
- identity ORACLE_A_ame - b_dropped = -0.6495711425 EXACT
- wooldid `robust` / `ccnoabsorb robust` / `cluster(unit)` = 1.0162389082 (== ORACLE A)
- wooldid `customvce(robust)` / `customvce(cluster unit)` / `unconditionalse robust`
  = -0.6495711425 (== ORACLE B), se 0.29-0.32, DGP truth 0.916 -> ~5 sigma, sign flipped
- poisson: the ppmlhdfe rc-3200 fallback FIRES, so the control becomes an explicit
  regressor and the SAME drop mechanism applies. Heterogeneous-effect poisson leg:
  wooldid -0.5506611369 == my ORACLE-POIS B -0.5506611369; ORACLE-POIS A 0.0454169223;
  DGP truth 0.1587. So the prober's "driven to zero under poisson" is the right
  PREDICTION only under homogeneous effects and the WRONG MECHANISM.

## R3/R4 — finding A  (R3 was CORRECTED and rerun: its oracles now use R4's spellings)
R4 pinned the exact spellings: `reghdfe y ib0.cell##c.dose, absorb(i.i i.t i.nullv##c.ringer)`
reproduces wooldid's default path to every digit, and
`reg y i.i i.t i.nullv##c.ringer ib0.cell##c.dose` reproduces the customvce path.
- dydx identical across engines: 1.022234373986 both -> SAME model, SAME fit
- eydx: 0.034876727614 (reghdfe) vs 0.034066172921 (regress)
- reghdfe `predict, xb` has the SAME MEAN as y (27.395972) but a smaller sd
  (1.262981 vs regress 1.943938 vs y 2.626373): the absorbed effects are missing from
  the denominator, so the error is second-order in the dispersion, not first-order.
- binary-side `semielasticity`: 0.087622452771 (robust) vs 0.085833510361 (customvce),
  a 2.1% gap -> the SAME mechanism DOES reach the binary option (the prober's own
  binary numbers moved the wrong way with their scale knob; mine do not).
- R5/R5b: the gap is 2-4% on ordinary panels and grows as E[y] falls toward the
  dispersion of the absorbed effects (eydx reghdfe/regress: 0.02157/… at ymean 47,
  0.05032/… at ymean 13). The prober's 13.7% is a design where the absorbed part
  dominates; it is not a universal magnitude.

## R6/R7 — findings B1, B2, C
Fixture: 12 states x 5-12 counties x 10 years, county-level dose (+ optional
within-county time variation `tvar'), i = state, truth = +0.4 on every treated post row.
- B1: SAME dataset, 6/6 reps at tvar=1: reghdfe max|e(V)| = 0 EVERY time, regress
  max|e(V)| = 109-646, IDENTICAL rss to 4 dp, IDENTICAL margins-ame to 6 dp, ranks
  221-230 (reghdfe) vs 174-190 (regress). wooldid default path loses every SE;
  customvce(robust) keeps them on the same point estimates. No cluster() anywhere.
  The point estimate being invariant across two different collinearity drops is
  evidence the ame IS estimable, so the regress SE is a legitimate object.
- B2: 12-draw MC at tvar=1, lattice(all 3), customvce(robust):
    ame mean -0.342, SD 2.519, range [-6.091, 2.906]; TRUTH 0.400
    reported se mean 0.341, range [0.305, 0.434]
  => no detectable BIAS (mean is 1.0 se from truth) but the reported conditional
  robust SE is 7.4x too small. The "45x miss" is a heavy tail, not a systematic miss.
  No-lattice control arm: ame mean 0.462, SD 0.146, reported se mean 0.110 (1.3x).
- C: under bare `robust` with no cluster() anywhere, both cautions print verbatim as
  claimed, the "finer clustering levels" remedy does not exist, the referenced
  cell-containment caution never printed, and the third remedy the second caution
  offers ("conditional robust standard errors") is what the user already asked for.

## D — refuted from the source
sthlp:1051-1053 documents this exactly: "Under customvce() or unconditionalse ... a BARE
token there is read as a CONTINUOUS variable: fe(state) fits a linear term in the state
code rather than a set of state indicators. Write the prefix explicitly - fe(i.state)."
The runtime caution (ado:4637-4638) also says "manually specifying any needed i. prefixes".
The prober's expectation ("wooldid already knows the combination is invalid") contradicts
the help, which says the combination is valid with the prefix written out.


## Not run (boundary for a verifier)
- HC3 / leverage diagnosis of the lattice SE collapse (why the conditional robust SE is 7.4x
  too small): not attempted; `customvce(hc3)` and `predict, hat` would be the next step.
- finding A under `poisexpresults`, and under `poisson` (where eydx is remapped to dydx, so
  the defect should not reach it) — not measured.
- B1's caution text under `cluster(cty)` (a cluster FINER than i): the SE was lost there too
  in R6 at tvar=1, but I did not capture the caution wording on that leg.
- R6's quiet legs and its noisy leg draw DIFFERENT datasets (mkdata redraws). R7 fixes this by
  saving each draw to a tempfile; the R7 numbers are the same-data ones.
