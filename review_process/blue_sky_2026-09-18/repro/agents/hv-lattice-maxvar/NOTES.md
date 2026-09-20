# hv-lattice-maxvar -- verifying "the maxvar preflight omits lattice()"

## The claim under test
`ado:1371` builds `si8demand = si8cells*(1 + (contreat!="")) + si8ilev + si8tlev + 200`
and refuses with a named error if that exceeds `c(maxvar)`. It never reads `lattice()`.
Claim: a realistic lattice run therefore sails past the guard and dies at a bare r(907).
The guard's false-negative DESIGN is ruled (retracted defect, 2026-09-16) -- the open
question is the SIZE of the miss and what the user sees.

## Mechanism derived from the source (read, then to be measured)
- `ado:1476` regressor block is `i.`_i'#i.`_t'#i.`targetflag'`contreatinj'`
- `ado:1117` `contreatinj = `latticepart2hash'##`contreatpolypart'`
- `ado:1087` `latticepart2hash = ##i.`latticevar''`
=> the estimated term is a FULL FACTORIAL `cell ## lattice ## c.dose`.
   With L lattice levels and contreatpoly(p) the column count is ~ cells*(1+L)*(1+p)*...,
   i.e. the miss is MULTIPLICATIVE (~1+L), not additive. The preflight's `*2` covers
   only the p=1, L=1 case.
- `_i` (ado:928) = i on POST rows, 0 elsewhere; `_t` = relative-period index, 0 on the
  reference. So si8cells = # treated (cohort, post-relative-period) pairs.

## Harness facts measured
- Stata 19.0 MP, 8-core license. `c(maxvar)` default **5000**, `c(max_k_theory)` 5000.
  Stata/BE's ceiling is a FIXED 2048 and 2048 is also MP's `set maxvar` minimum, so
  `set maxvar 2048` is an exact emulation of a BE user's hard limit.
- Build asserted: `which wooldid` -> scratchpad/bluesky/build/wooldid.ado.
- gen_core's `dose` is 0 outside treatment; `contreatcontrols(basic)` + `lattice()` on it
  is degenerate (r(322), 51 collinearity cautions). Wrote my OWN generator instead:
  county panel, i = state, time-invariant county EXPOSURE defined on every row
  (Bartik/shift-share shape). No degeneracy.

## Probes
- `p0_orient.do` / `.log` -- orientation; flavor/maxvar; verbose regression string.
- `p1_bisect.do` / `.log` -- 25 states x 20 counties x 20 years, cells=213.
  fvexpand column counts; bisection of `set maxvar`; user-visible output at 2048;
  a 120-state design that TRIPS the preflight so wooldid prints its own si8demand
  (validates my reconstruction of the formula).
- `p2_default.do` -- (planned) 50-state county panel at the DEFAULT maxvar 5000.

## Status
- p1 launched. Results below when in.

## p1 RESULTS (2026-09-17) -- CLAIM CONFIRMED, UX half is TRUE-BUT-DIFFERENT
Design p1A: 25 states (i), 18 treated, 20 counties each, T=2000..2019, 10000 rows.
si8cells = 213, si8ilev = 25, si8tlev = 20  => **si8demand = 671**.

fvexpand column counts on the EXACT regressor block wooldid builds:
| regressor block                                | columns |
|------------------------------------------------|---------|
| `i._i#i._t#i.tf`                     (cells)   |   214   |
| `... ##c.exposure`        (what si8 models)    |   429   |
| `... ##i.lat5##c.exposure`  (lattice(all 5))   |  2579   |
| `... ##i.lat5##c.exp##c.exp` (+contreatpoly 2) |  3869   |

Measured `set maxvar` bisection, identical command, no `verbose`:
| maxvar | no lattice | lattice(all 5) |
|--------|-----------|----------------|
| 2048   | rc 0      | **rc 907**     |
| 2500   | -         | rc 907         |
| 3000   | -         | rc 907         |
| 12000  | -         | rc 0, 162 s, att = 1.0478 (design truth 0.8*E[exposure] ~ 1.02) |

So: preflight predicts 671; the run actually needs >3000. Miss factor **>4.5x**, and
structurally ~(1+L): 429 -> 2579 is x6.01 for L=5.

### What the user actually sees (verbatim, maxvar 2048, no verbose)
wooldid prints its fastcollapse caution, then STATA's own error, with no wooldid line:
```
maxvar too small
    You have attempted to use an interaction with too many levels or
    attempted to fit a model with too many variables.  You need to increase
    maxvar; it is currently 2048.  Use set maxvar; see help maxvar.
    ...try set emptycells drop ... use the c. operator ...
```
=> the claim's "no wooldid message" is literally true, BUT "uninformative" is WRONG:
Stata's r(907) text names the exact remedy (`set maxvar`). wooldid's OWN preflight
message (ado:1372-1373) never mentions `set maxvar` at all -- it says "Consider
specifying unit-level fixed effects via fe() and swapping i to a coarser level".
The early named error is therefore LESS actionable than the late raw one. That
inversion is the real UX finding.

### Note
Section 5 of p1 died r(199) -- `clear all` drops the user's own `mkpanel` program.
Rerun as p1b.

## p1b RESULTS -- reconstruction VALIDATED, contreatpoly shares the omission, P6 clean
- **Reconstruction exact.** Design B (120 states, 110 treated, 4 counties, T=20):
  my recon said cells=1271, si8demand=2882. wooldid printed, verbatim:
  "Error: this specification requires at least 2882 variables (about 1271 treated
   cohort-by-period interaction terms) but maxvar is currently 2048..."
  Both numbers match. The formula reconstruction is measured, not assumed.
- **P6 holds (negative result).** Same design + `contreatwithin` at 2048 still produced
  wooldid's own preflight error and rc 907. The preflight is NOT disabled by
  contreatwithin's internal safety rewrite. No contradiction of the ruling.
- **contreatpoly() is omitted by the same formula.** Design A, maxvar 2048,
  si8demand = 671 for BOTH rows:
  | spec | fvexpand cols | rc |
  |------|---------------|----|
  | `lattice(all 3)`                  | 1719 | **0**   |
  | `lattice(all 3) contreatpoly(3)`  | 3439 | **907** |
- **Death location (verbose, maxvar 3000):** the last thing wooldid says is
  "Wooldid is now estimating its regression." and then the `reghdfe` line; the r(907)
  comes out of `reghdfe` itself. So the failure is genuinely MID-RUN, which is exactly
  the outcome `sthlp:1167` says the preflight prevents.

## p2 aborted and was rerun as p2b
`fvexpand i._ii#i._hh#i._tf##i._lat##c.exposure` at the DEFAULT maxvar 5000 itself
returned r(907) -- direct evidence that the 50-state lattice block exceeds 5000
columns while si8demand predicts 1320. p2b reorders so the wooldid runs happen at the
untouched default and the counting happens afterwards at maxvar 20000.

## p3 (SELF-FOUND) -- is lattice the only omission?
`ado:1138` `coarsecohortcontrols` -> `i._i#i._t#(vars)` and `ado:1171`
`interactivecontrols` -> `i._i#i._t#1.targetflag#c.(vars)` are ALSO cells-multiplying
and ALSO absent from the ado:1371 formula. p3_family.do measures both with no lattice
at all.

## p4b RESULTS (self-found) -- where the blind spot lives
Same panel, only `i` changes. Identity `cols = 2*(n_cell+n_cell_lat)+2L+1` reproduces
fvexpand's 2579 exactly, so the counting is calibrated.
| i | si8cells | si8demand | cols nolat | cols lat5 | ratio |
|---|----------|-----------|------------|-----------|-------|
| state (20 rows/cell) | 213 | 671 | 429 | 2579 | 3.844 |
| county (1 row/cell)  | 4260| 9240| 8523| 17063 | 1.847 |
At i=unit the bound (9240) already exceeds maxvar 5000 so the guard fires; at i=state
the bound (671) is a seventh of maxvar and the guard is silent while the run needs
3000-3500. The blind spot is COARSE i + THICK cells -- the shape ado:1371 recommends.

## p1b bracket: BRACKET2_MAXVAR3500_RC=0
So design A's true demand is in (3000, 3500]. Preflight 671 -> miss factor 4.5x-5.2x.

## p3 RESULTS (self-found) -- CONFIRMED: the omission is a family, no lattice needed
Design A, si8demand = 671 on every row, maxvar 2048:
| spec (contreat, NO lattice)                   | cols | rc  |
|-----------------------------------------------|------|-----|
| baseline                                      | 429  | 0   |
| `interactivecontrols(x1..x8)`                 | 2293 | 907 |
| `coarsecohortcontrols(i.reg)` (12 regions)    | 2997 | 907 |
So patching si8demand for lattice() alone would leave three more multipliers uncovered.

## p1b maxvar-4000 arm killed as redundant (3500 already rc 0).

## FINAL STATE (everything needed is on disk; see FINDINGS.md for the writeups)
- `repro_lattice_maxvar.do` was REWRITTEN and RUN END-TO-END. The first version had the
  same bug that aborted p2: a bare `fvexpand` of the lattice block hits r(907) itself at
  the default maxvar, so batch Stata aborted before the wooldid runs. The rewrite counts
  columns with the nonempty-cell identity instead. Verified output:
    si8demand = 1320 | cols nolat = 1053 | cols lat5 = 6115 | RUN1_RC = 0 | RUN2_RC = 907
- `c(max_k_theory)` tracks `c(maxvar)` (5000 -> 20000 after `set maxvar 20000`), so
  maxvar is the single operative ceiling.
- p2b's optional maxvar-20000 estimability confirmation was left to run under a deadline;
  it is NOT load-bearing -- design A at maxvar 12000 (rc 0, 162 s, att = 1.0478 against a
  design truth of ~1.02) already shows the lattice spec class is estimable and correct.

## Findings shipped
F1 lattice+contreatpoly omitted from si8demand (seeded, CONFIRMED, medium/certain)
F2 the UX half: "no wooldid message" true, "uninformative" overstated (low/certain)
F4 the omission is a FAMILY: interactivecontrols + coarsecohortcontrols too (self-found,
   medium/certain)
F5 the granular-i caution's framing (self-found, low/probable, scope-limited)
F3 P6 holds under contreatwithin -> reported as a CLEAN negative, not a finding.
