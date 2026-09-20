# fp-contreat — INDEPENDENT REIMPLEMENTATION probes of the continuous-treatment AME family

## Method
EXACT-FIT ORACLE. Generate a panel whose outcome is produced with zero residual by a model
nested in wooldid's continuous-treatment specification
    y = a_i + g_t + g(dose) + 1{t>=ttre}*(kappa_it + beta_it*dose)
Then the fitted conditional mean equals y, every margins-visible dose derivative equals
beta_it, and the truth for each reported quantity is a closed-form weighted mean I compute
myself. No reconstruction of wooldid's tempvar-laden regression is needed.
Generator: gen_ct.do (arg = which nuisance dose function g() is planted).

## Log so far
- probe01_smoke: exact-fit design works; ame = .9325 = mean(betatrue) on treated post rows. R2 = 1.
- probe02_invariance: 54 configs, 6 nuisance shapes x 9 contreatcontrols() spellings.
  CLEAN: every matched pair reproduces the truth to <7e-14; mismatched pairs move 10-90%,
  so the test has power. contreatcontrols() basic/byTime/byCohort/byTreated/byTtre all
  partial out exactly what they claim on the default (reghdfe absorb) path.
- probe03_engines: same battery x {DEFAULT, ccnoabsorb, customvce(robust)} — ALL CLEAN.
  (unconditionalse alone is refused r(198): it requires robust or cluster. Informative.)
- probe04/05_elast: **HIT.** contreatelasticitytype(eydx) and (eyex) DISAGREE between the
  reghdfe path and the plain-reg path ON IDENTICAL DATA, and only the plain-reg path matches
  the exact-fit truth. dydx and dyex agree everywhere.
     eydx: robust 0.099345285547  vs truth/customvce 0.101147453290   (1.78% off)
     eyex: robust 0.196842161516  vs truth/customvce 0.196440383258   (0.20% off)
  Hypothesis: margins' prediction after reghdfe is xb EXCLUDING the absorbed fixed effects,
  so the "% change in y" divides by the wrong denominator.

- probe09/10/11 (lattice): **HIT.** all-zero VCV reproduced on a realistic 12-state fixture for
  lattice(all 2/3/5) AND a user lattice variable, under robust and cluster, invariant to regtol
  and ccnoabsorb; the SAME model on plain regress (customvce/unconditionalse) returns finite SEs
  with identical b, rss, r2, df_r. => the zero is an ABSORBING-ENGINE artifact. Also: the lattice
  ame points are wild (0.95 / 7.01 / 5.79 against a true 0.4) with no SE to flag them.
- probe13 (estimands x weights x ES horizons, exact fit): CLEAN to 1e-15.
  att/att_it/att_i/att_itime main ame all exact, weighted and unweighted; every ES_ame horizon
  exact including the two pre-treatment horizons with a planted pre-trend slope.
  LAYOUT NOTE: e(wdidmainresults) roweq is main / pretreatment / ES_att / ES_ame; the
  pre-treatment ame row sits under roweq "pretreatment", not "main".

- probe14 (reghdfe replication OUTSIDE wooldid): reproduces the zero V in three lines.
  reghdfe maxabsV=0 (rank 160, df_a 36) / reghdfe-with-explicit-controls maxabsV=0 (rank 165) /
  reg maxabsV=2923.571 (rank 154). SAME rss 4100.2092, SAME df_r 2004. reghdfe keeps 6-11 more
  columns than reg drops as collinear -> probable mechanism.
- FINDING C: the two all-zero-VCV cautions (ado:12983, ado:1567) blame "analytic clustered SEs"
  and point at a cell-containment caution that never printed, on a run whose VCE is `robust`.
- FINDING D (low): fe() + customvce()/unconditionalse -> wooldid cautions, then enters the fe()
  variable as a CONTINUOUS regressor (probe08 leg A: `reg y i.i i.year unit ...`).

## QUEUE STATE (runmaster.sh, background id bnyqdyxbl)
order: repro_B, probe15 (const dose), probe21 (lattice types + poisson),
       probe19 (=16 poly/latuser/untreatedZero + 17 within/collapse + 18 unbalanced),
       probe20 (pretreatment ame + non-att ES horizons), probe22 (if/subgroups/missing dose)
Stata contention from other agents has been at 8-9 all session; each launch waits up to 5 min
for a slot then proceeds.
