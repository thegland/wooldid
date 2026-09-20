# FUZZ agent — FINDINGS

## Method
Randomised, compatibility-aware option-combination fuzzing with **automated** assertions,
run in five arms. Nothing was eyeballed; every arm asserts a machine-checkable proposition.

| arm | what it does | runs |
|---|---|---|
| A / B | random sensible command lines over 19 generated datasets; assert on the whole e() surface | 250 |
| X (oracle) | zero-noise DGP where all four estimands have a closed-form design truth; assert equality | 160 |
| O (observe) | the same, with identification-CHANGING options added, to learn which move the number | 30 |
| P (metamorphic) | the same random command on a transformed twin (shuffle / tshift / idperm / yscale / yshift / ynegate / wscale / if-vs-keep / collinear control / rowdup / reverse sort) with the exact algebraic relation asserted cell-by-cell | 150 |
| R (refusal) | 80 deliberately invalid option combinations; assert the guard fires with a wooldid `Error:` | 80 |
| targeted | probes written to settle a specific documented claim | 9 |

Automated assertions on every rc-0 run: SE >= 0 and not missing where a number prints;
p in [0,1]; CI ordered and containing its point estimate; t == b/se on the analytic path;
CI half-width in standard errors sane; e(b)/e(V) conformable, symmetric, PSD; e(N) > 0;
e(sample) non-empty; N_clust >= 2; the CI column name carries the requested level();
and `datasignature` / `: sortedby` / variable count identical before and after the call.

---

## FINDINGS

### FZ-1  `contreatelasticitytype()` relabels a binary-mode ATT table as elasticities while changing nothing
`severity high / realistic true / confidence certain / surface contreatelasticitytype()`

With a **binary** treatment (no `contreat()`), `contreatelasticitytype(eyex|eydx|dyex)` is
accepted (rc 0), leaves every number bit-identical to the plain ATT, posts nothing to
`e(contreatelasticitytype)` — and rewrites the results header:

```
plain                      Main Results (Full Estimation Sample):
contreatelasticitytype(dyex)  ... : Marginal Treatment Effects Presented as Semi-elasticities (dyex)
contreatelasticitytype(eydx)  ... : Marginal Treatment Effects Presented as Semi-elasticities (eydx)
contreatelasticitytype(eyex)  ... : Marginal Treatment Effects Presented as Elasticities
```
att = 2.1109708, se = .00800688 in ALL FOUR cases, identical to the baseline run.

SCOPE, now measured: with `contreat(dose)` present the option DOES work -- the `main:ame`
row moves from -.0415657 (dydx / default) to -.03332325 (dyex), -.00079954 (eydx) and
-.00064271 (eyex), while the `main:att` row correctly stays at 1.5290969 throughout.
So the defect is exactly: inert in binary mode, but still rewriting the header.

Every other `contreat()` companion IS screened in that position (`lattice()`,
`contreatcontrols()`, `contreatwithin`, `makecfxplots` all error r(198) "requires
contreat()"), so the missing screen here is an inconsistency as well as a hazard: a
reader takes a level ATT away as an elasticity.
Repro: `repro/REPRO_fuzz_elasticity_label.do`

### FZ-2  every semi-elasticity form dies with a RAW Stata error when the prediction is not strictly positive
`severity medium / realistic true / confidence certain / surface semielasticity, contreatelasticitytype(eydx|eyex)`

On an outcome whose fitted values can reach or cross zero (my DGP: a mean-zero outcome),
`semielasticity` and `contreatelasticitytype(eydx|eyex)` fail on EVERY configuration tried
(13/13: i=unit and i=cohort, iid/robust/cluster, with and without `nofastcollapse`,
weighted and unweighted, with and without `lattice()`), with **no wooldid `Error:` line**:

* usual signature, r(459): `could not calculate numerical derivatives -- discontinuous
  region with missing values encountered`
* with `contreatcontrols(untreatedZero)`, r(3598): `variable __00001P not found` followed
  by a Mata traceback `st_data(): 3598 ... __marg_cr_eval() ... _deriv() ...` — a wooldid
  **tempvar name leaking into a raw error**
* under `permute()` the same root cause resurfaces as `Error: the observed comparison
  statistic for att could not be read from the rep-0 run`, which does not name the cause
* with `ccnoabsorb` added it becomes `inconsistent estimation sample levels 0 and 1 of
  factor __000004` — again a raw `margins` error carrying a wooldid tempvar

The DECISIVE control: rebuild the identical design with a strictly positive outcome
(min 49.2) and all seven variants return rc 0. So the feature works; what is missing is a
guard saying "a semi-elasticity needs a strictly positive prediction".
`contreatelasticitytype(dydx)` and `(dyex)` — the two forms that do not divide by the
prediction — never fail.
Repro: `repro/REPRO_fuzz_semielasticity.do`

### FZ-3  `level(99.9)` and `level(99.99)` are refused with a message that contradicts the input and wooldid's own range check
`severity medium / realistic true / confidence certain / surface level()`

`wooldid ..., level(99.9)` -> r(198) `level() can have at most two digits after the decimal
point`. So does `level(99.99)`. But `set level 99.9` and `set level 99.99` succeed in plain
Stata (c(level) comes back exactly 99.9 / 99.99), wooldid's own validator says "level() must
be between 10 and 99.99", and wooldid accepts `level(97.5)`, `level(90.25)` and `level(99.5)`.
The accepted set is exactly the binary-representable levels, so the level is being rendered
back to a string at full double precision (99.9 -> 99.900000000000006) before `set level`.
99.9% intervals are ordinary in multiplicity-adjusted reporting.
Repro: `repro/REPRO_fuzz_level999.do`

### FZ-4  `estimates replay` after wooldid fails r(198) and blames a command line the user never typed
`severity low / realistic true / confidence certain / surface replay`

`estimates store m1` then `estimates replay m1` prints
`Error: command line input should feature all of y i t ttre; too few variables specified.`
and exits r(198). Bare `wooldid` replay works; `estimates restore m1` + bare `wooldid`
works; `estimates table`, `estimates stats`, `esttab` and `coefplot` all work. The replay
gate is `if \`"\`0'"' == ""`, i.e. STRICTLY empty, so anything appended defeats it — a
trailing blank does too (`wooldid ` with a trailing space -> r(198), measured).
Repro: `repro/REPRO_fuzz_estreplay.do`

### FZ-5  `dots()` negative is refused under `bootstrap()` and silently accepted under `permute()`
`severity low / realistic false-ish / confidence certain / surface dots()`

`wooldid ..., cluster(cohort) permute(shuffle) permreps(25) dots(-5)` -> rc 0 (the engine's
`if \`dots' > 0` tests all fail and it behaves as `dots(0)`), while the same `dots(-5)` with
`bootstrap()` -> r(198) `dots() must be a nonnegative integer`. The code comments this
(ado:6265-6270) but nothing user-facing does. Measured in spec_R1 rows R043/R044.

### FZ-6  the requested event-study window moves the headline ATT by 5-11%
`severity medium / realistic true / confidence uncertain (mechanism certain; whether it is under-disclosed is the author's call) / surface esprelength()`

Same data, same estimand, same VCE; only the requested horizons differ:

```
canon: esprelength(1) espostlength(1)  ->  att = 1.165019  se = .113196
       esprelength(2) espostlength(3)  ->  att = 1.104531  se = .114076     (-5.2%)
unbal: esprelength(1) espostlength(1)  ->  att = 0.915300  se = .053022
       esprelength(2) espostlength(3)  ->  att = 1.029771  se = .086225     (+11.4%)
```
MECHANISM, measured (probe_eswin.log): it is the PRE window only, acting through the
reference-normalization pool. On ds_canon, `e(N)` stays 480 and `e(escomp_trtN)` stays 192
while `e(escomp_refN)` falls 120 -> 96 -> 72 -> 48 as `esprelength()` goes 0 -> 1 -> 2 -> 3,
and the att tracks it 1.0952323 -> 1.1650189 -> 1.1045313 -> 1.0732189. `espostlength(3)`
ALONE leaves the att bit-identical to the windowless run (1.0952323). This is `sthlp:1370`'s
rule (the pool is the rows with `relyr < -esprelength`) doing exactly what it says. That is a legitimate specification choice,
and `e(escomp_ref*)` exposes the composition — but the two options are named and described
as horizon *requests*, the `pretreatment` note at `sthlp:724` discusses a different nuance,
and nothing at the point of use tells a reader that asking for one more pre-horizon moves
the headline number by a tenth. Filed as an observation the author should rule on rather
than as an outright defect; my oracle arm confirms the numbers are *correct* for the sample
each window implies (the design truth recomputed on `e(sample)` matches to 1e-15 either way).

### FZ-7  (observation, not filed as a defect) reversing sort order moves analytic SEs in the 8th digit
`severity low / realistic true / confidence certain / surface reghdfe projection path`

In the metamorphic arm, `gsort -i -t` before an otherwise identical call leaves every point
estimate identical to 1e-16 but moves `se`, `t`, `p` and the interval by up to 2.6e-8
relative (P0076, P0052). `regtol(12)` and `regtol(13)` do not shrink it. Above the project's
1e-10 indifference threshold but far below anything substantive; recorded so that a future
bisect does not mistake it for a regression.


### FZ-8  `update` on a full command line runs the package installer and silently estimates nothing
`severity low / realistic true / confidence certain / surface update`

`wooldid y unit t ttre, att robust update` returns rc 0, prints
"checking wooldid consistency and verifying not already installed... all files already exist
and are up to date", and posts NOTHING to `e()` -- `e(cmd)` is empty afterwards. The
estimation request on the same line is discarded without a word. `update` reads like an
ordinary option (it sits in the same syntax line as `robust` and `att`), so a user who types
it expecting "re-estimate and update the stored results" gets a no-op plus a net operation.
Measured in spec_R2 row S028 (flagged by my `NO_ECMD` assertion).

---
## NEGATIVE RESULTS (probed hard, came back clean)

1. **The four estimand weightings are exact.** 160 zero-noise oracle runs: every `i`=unit run
   (113) reproduced all four closed-form design estimands to <= 2e-15, across random
   combinations of vce, controls, `fe()`, ES windows, `subgroups()`, `nofastcollapse`,
   `regtol()`, `level()`, `safety(off)`, weights, `if`, `permute()` and `bootstrap()`.
   Every deviation was an `i`=cohort run. A dedicated control arm on data with NO unit fixed
   effect AND no within-cohort effect heterogeneity (`ds_ex_hom`, 22 runs incl. `controls()`,
   `ccc_absorb()`, `coarsecohortcontrols()`, `fe()`, `timetrends`) came back exact to 2e-15,
   which confirms the coarse-`i` deviations are ordinary residual variance, not a defect.
2. **Every estimand is the weighted average of its own cells.** `probe_gtwt`, 84 configs:
   `sum(w_e * estimate)/sum(w_e)` over `e(gteffects)` reproduces the reported att/att_it/
   att_i/att_itime to 1e-15, every aggregation weight is nonnegative, and every estimand sits
   inside the convex hull of its post cells. `probe_hull`, 69 more configs: same, no exception.
   This is an EXACT check on noisy, realistic data - it needs no oracle DGP.
3. **Intervals really are at `level()`.** `probe_civ`, 94 configs / 1504 interval rows:
   `(ub-lb)/(2*se) == invttail(e(df_r), (1-level/100)/2)` to 1e-7 relative, at 80/90/95/99
   under iid / robust / cluster(region|unit|cohort) on six datasets. Zero mismatches; no row
   silently used a normal critical value.
4. **Permutation p-values sit exactly on their documented grids.** `probe_pgrid`, 17 configs:
   default two-sided add-one `2(m+1)/(n+1)`, `permrawp` `c/n`, `permabsp` `(m+1)/(n+1)` and
   `permside(right)` all integral on their own denominator; 0 off-grid out of 12 rows each.
5. **Engines replace inference, not estimates.** `probe_engine`, 61 comparisons: the observed
   point estimates under `permute(shuffle)` (+`permtstat`/`permanalyticse`/`permtci permabsp`/
   `permonsetcohorts`), `bootstrap(standard)` (+`boottstat`) and `bootstrap(bayes)` are
   identical to the analytic run's to <= 1.1e-16. No rc mismatches.
6. **Metamorphic invariance.** 125 pairs over 11 transforms (row shuffle, reverse sort, id
   relabel, time shift, y shift, y rescale, y negate, weight rescale, `if` vs `keep if`, a
   collinear control, row duplication): point estimates invariant to 1e-16, no rc differences,
   and the only movement is 1e-9..2.6e-8 in derived columns (FZ-7).
7. **Internal identities.** Adding estimands, adding `subgroups()`, adding diagnostic flags
   (`jointtests adversarial gteffects summarystats histogramcohorteffects`) and changing
   `level()` never move the other rows; `subgroups(sgall)` reproduces the main table exactly
   (12/12 rows, 4 datasets).
8. **The -1e100 placeholder behaves as documented.** `sthlp:1616-7` claims `test` drops the
   constraint and `lincom` propagates missing. Measured exactly that (`test` -> "Constraint 1
   dropped", `lincom`/`nlcom` -> se/t/p/CI missing, estimate still shown); a clean-coefficient
   control test still gives the right p.
9. **Section-7 item 11 resolves CLEAN on this build.** Under `fastnose`, wooldid DOES post
   `e(b_es_att)` (1x4) and `e(numpre_att)` (=2) and posts no `e(V)` (`e(properties)`="b"), so
   `sthlp:1325` listing `fastnose` under "the covariance is not returned but the coefficient
   vector still is" is accurate. `test`/`lincom` then refuse with r(321), which is correct.
10. **Guards fire.** 80 deliberately invalid option combinations: every one refused with an
   informative wooldid `Error:` line except the cases in FZ-1 and FZ-5. Across 300 random
   runs, only ONE failure out of ~83 was a raw Stata error with no wooldid message (FZ-2);
   the other raw-looking one (a 9-digit subgroup value overflowing Stata's 32-character name
   limit, `pretreatment_SG480000000_att_itime invalid name`) is preceded by an explicit
   wooldid Caution that predicts exactly that crash and names the remedy, so it is NOT filed.
11. **No state leakage.** `datasignature`, `: sortedby` and the variable count were identical
   before and after every one of the ~300 fuzz calls, successes and failures alike.
