# hv-contreat-sentinel — FINDINGS

**How to run any do-file here:** `cd` into this directory first, then
`"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do <file>.do`.
Each one builds its own panel; none touches `test/` or the project suite.

Verifier for recon **F5**: "contreat() can ship a 0 point estimate with a -1e100 variance,
and two exported copies disagree."

Build asserted in every probe: `which wooldid` ->
`/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`, `*! version 1.4 2026-08-19`.

---

## VERDICT ON F5: **TRUE, BUT DIFFERENT** — and the real damage is one step downstream

The mechanism reproduces exactly. The *harms* recon predicted mostly do not; a harm recon did
not predict does, and it is much worse than the one it did.

### What reproduces (p1_repro.do / p1_repro.log)
DGP: 100 units, 5 cohorts (4 treated onset 4/6/8/10, 1 never treated), T=1..12,
`dose = dose0_i * 1{t >= ttre_i}` — the canonical exposure-x-post continuous DiD.
`i = coh`, so every cohort-period cell holds 20 units with real dose variation.

In **5 of 5** configurations (untreatedZero/basic x iid/robust/cluster(unit)/cluster+fe(unit)):

* `e(b)` carries **exactly 0** on `pretreatment_ame{,_it,_i,_itime}` and `ES_ame{,_it}_pre4/3/2`;
* `e(V)` carries **-1.0000e+100** on those same 10 of 40 diagonal cells;
* `e(V_es_ame)` / `e(V_es_ame_it)` carry **exactly 0** on 3 of 6 diagonals — the same horizons;
* `mineig(e(V)) = -1e100`, `cholesky(makesymmetric(e(V)))` -> 820 missing cells.

**This is NOT the ruled AMESE R1 geometry.** R1 (cell containment, `cluster()` at or above `i`)
withholds `e(V_es_ame*)`. Here `e(V_es_ame)` **is** posted, `e(b_es_ame)` is posted, and the run
reproduces identically under plain iid, under `customvce(robust)` and under `cluster(unit)` —
all strictly finer than / orthogonal to the containment predicate. `amescrub = 0`.

### What recon predicted that is REFUTED (p4_consequence.do / .log)
* "the two exported copies disagree with each other" — they do, but **both are announced**. Two
  cautions fire on every one of the five configs:
  `"standard errors were not identified for 10 of 40 reported estimates..."` and
  `"e(V_es_ame) contains 3 of 6 horizons with zero or missing variance..."`.
* "zero-width esttab CIs / false-confidence lincom" from the `e(V_es_ame)` copy — **refuted.**
  `ereturn post e(b_es_ame) e(V_es_ame)` succeeds, and Stata then renders those rows as
  `(omitted)`, not as a zero-width interval: `lincom pre4` -> `est=0, se=0, lb=., ub=.`;
  `test pre4 pre3 pre2` -> all three constraints dropped, `chi2(0) = .`.
* "can a -1e100 reach test/lincom/esttab without a marker" — **no.**
  `esttab , se` prints `0.0000  (.)`; `esttab , ci` prints `0.0000  [.,.]`;
  `test ES_ame_pre4 ES_ame_pre3 ES_ame_pre2` -> `p=., df=0`, rc 0;
  `estat vce` runs; `estimates save` + `estimates use` preserves BOTH the 10 sentinels and
  `e(es_ndegen_ame) = 3`. The negative diagonal is itself the documented marker (`sthlp:1616-7`).
* **Realism gate (p2_realdose.do / .log).** With a dose that varies in *every* period (a genuine
  intensity index), **no sentinel fires in 10 of 10 configs** — `basic`/`byTime`/`byTreated`
  x iid/robust/cluster(unit), and `untreatedZero` too. So the trigger is precisely
  *pre-treatment dose identically zero*, which is exactly the premise
  `contreatcontrols(untreatedZero)` documents for itself ("the dose-zero-outside-treatment
  premise", sthlp:472). On every by-the-book `untreatedZero` run the ame pre-treatment
  event-study block is structurally degenerate.

### Which exported copy is "right"
Neither is wrong on its own terms, and the help documents each separately. `-1e100` in `e(V)`
means "this standard error is not identified" and renders missing on every Stata surface.
The exact `0` in `e(V_es_ame)` is `margins`' literal delta-method variance of a coefficient the
solver constrained to zero. They disagree because they answer different questions. Inside
Stata the disagreement is harmless (verified above). It becomes harmful the moment the matrix
is handed to an *external* consumer that reads a zero diagonal as information — which is
precisely the consumer wooldid's help routes users to.

---

## FINDING 1 (self-found) — HIGH · the advertised `honestdid` handoff returns a sensitivity analysis that is IDENTICAL AT EVERY M

Repro: `repro_ame_degenerate_preblock.do` (self-contained; builds its own panel).

On the same run, same data, same `mvec(0(0.5)1)`, `delta(rm)`, `parallel(0)`, under `cluster(unit)`:

```
--- honestdid AME  (b=e(b_es_ame), vcov=e(V_es_ame), numpre=e(numpre_ame)=3) ---
|    M    |   lb   |   ub   |
|       . |  0.447 |  1.684 | (Original)
|  0.0000 |  0.448 |  1.686 |
|  0.5000 |  0.448 |  1.686 |     <-- unchanged
|  1.0000 |  0.448 |  1.686 |     <-- unchanged
rc = 0, no error, no warning beyond honestdid's boilerplate M=0 note

--- honestdid ATT  (same run, att family) ---
|       . |  0.605 |  1.467 | (Original)
|  0.0000 |  0.604 |  1.467 |
|  0.5000 |  0.498 |  1.564 |     <-- widens, as it must
|  1.0000 |  0.348 |  1.723 |
```
Reproduced under **iid** (0.356/1.782 at every M) and **robust** as well. `delta(sd)` on the
ame family instead dies `r(3498)` ("HonestSensitivityResults(): function returned error").

**Why.** Delta^RM bounds the post-treatment parallel-trends violation by
`M * max|pre-period violation|`. wooldid hands honestdid three pre-treatment coefficients that
are exactly 0 with exactly-0 variance, so the feasible pre-period violation set is `{0}` and the
bound is `0` for every `M`. The run reports the ame effect as **completely robust to arbitrary
violations of parallel trends** when in fact it carries no pre-treatment information at all.

**What a user loses.** This is the failure mode Rambachan-Roth exists to prevent, delivered with
the opposite sign: maximally reassuring, silent, and rc 0. The help *does* state the precondition
("a covariance with no structurally degenerate rows (see `e(es_ndegen_type)`)") and
`e(es_ndegen_ame)=3` does survive `estimates save/use` — but nothing blocks, warns at the point
of use, or is visible in honestdid's own output. The att family on the identical call behaves
normally, so the output gives the user no contrast to notice.

## FINDING 2 (self-found) — MEDIUM · the `adversarial` ame p-value silently drops every pre-treatment horizon, and is the only pre-trend diagnostic on the run that does not disclose

Repro: `repro_ame_degenerate_preblock.do` (same file; also `p5b_adversarial.do`).

I reimplemented the minimum-distance line statistic independently from the help's own
description (GLS under a generalized inverse of `e(V_es_*)`; numerator df = rank(V) - 1) and
validated it: it reproduces wooldid's printed `att` p to 8 digits (0.01555611 vs printed 0.0156)
and its printed `ame` p exactly (0.06990167 vs printed 0.0699). Then, for each VCE:

| vce | W on the FULL 6x6 `e(V_es_ame)` | W on the POST-ONLY 3x3 submatrix | bitwise equal? | q |
|---|---|---|---|---|
| iid | `+1.5560e2ba7ef06X+002` | `+1.5560e2ba7ef06X+002` | **yes** | 2 (att: 5) |
| robust | `+1.02ca28d15d03cX-003` | `+1.02ca28d15d03cX-003` | **yes** | 2 (att: 5) |
| cluster(unit) | `+1.9987dcc78cab6X+003` | `+1.9987dcc78cab6X+003` | **yes** | 2 (att: 5) |

The 1e-12 eigenvalue truncation in the statistic removes exactly the three degenerate
pre-treatment directions, so the reported number is a test of *"the three post-treatment ame
coefficients lie on a line through the reference period"* — printed under the label
`"p-value for joint test of null that the ame coefficients and the reference period are on a
single line"`, in a section the option's own documentation calls **adversarial pre-trend
testing**.

**The disclosure gap is internal to the same run.** On this identical call:
* ame **joint pre-trend test** -> prints `.` (correctly withheld);
* **`pretrendpower()`** -> `"Caution: the ame pre-period variance-covariance block contains
  missing or non-positive variance entries; pretrendpower() diagnostics skipped for ame."`;
* ame **adversarial** -> prints `0.0699` with no caution, no printed `q`, and no `e()` record of
  the reduced rank. `e(adversarial_pval_ame)` is posted; `e(adversarial_line_ame)` is posted;
  nothing stores the numerator df or which horizons were used.

Two of three pre-trend diagnostics tell the user; the third silently degrades.
**The rule is stated; its firing is not.** The `adversarial` option's own entry
(`sthlp:750-754`) does say the numerator df equal "the rank of the event-study covariance
minus one", and `sthlp:389` names the 1e-12 truncation. What is missing is any signal that
it fired *on this run*: no printed q (the att line and the ame line look identical), no
`e()` record of the rank or of which horizons were used, and no caution — while its two
siblings caution on the same call. `e(adversarial_pval_ame)` and `e(adversarial_line_ame)`
are posted; neither carries the df.

## FINDING 3 (self-found) — LOW · no `e()` scalar records the `e(V)` sentinel count

`e(es_ndegen_ame)` records the `e(V_es_*)` zero-diagonal count and survives `estimates
save`/`use`. The partial-identification rescue's counterpart ("10 of 40 reported estimates")
exists only as a run-time display line (`ado:4079`). A replayed or restored estimate can
recover the fact only by scanning `e(V)`'s diagonal for negatives. Verified: after
`estimates use`, 10 sentinels are present and nothing names them.

---

### Realism, stated precisely (do not over-read)
The trigger is: `contreat()` + `esfixedbaseperiod` + at least one pre horizon + a dose that is
identically 0 in pre-treatment periods. That last condition is exactly the premise
`contreatcontrols(untreatedZero)` documents for itself — "Special option to formulate the
regression to work when a continuous treatment variable is always equal to 0 outside of
treated cohort-periods" (`sthlp:995-997`) — and the option's own entry says nothing about the
ame pre-treatment horizons becoming unidentified as a consequence. It also fires under
`contreatcontrols(basic)` on the same dose (P1 config B). What I did NOT show: that the
premise-satisfying case is the only trigger. The help's own `untreatedZero` example
(`sthlp:1551`) carries no ES options, and its honestdid example (`sthlp:1560-1565`) uses the
**att** family — so reaching Finding 1 requires the user to combine the two, which the Stored
Results entry invites by listing `ame` among the families `e(b_es_type)`/`e(V_es_type)` are
returned for.

## PROBED HARD, CAME BACK CLEAN
* `esttab` (se and ci), `ereturn display`, `lincom` (single, sum, and mixed with a
  well-identified row), `test` (pure-degenerate, mixed, and pure-att), `estat vce`,
  `estimates store/restore/save/use`: every one renders the sentinel as missing or drops the
  constraint, with rc 0 and no fabricated inference. The author's verification list holds.
* **Recon F5 lead 3 (`test` answers a different hypothesis) reproduced and is ruled.**
  `test ES_ame_pre4 ES_ame_post0` -> `p=0.0034, df=1` with `Constraint 1 dropped` printed —
  the one place a *scripted* `qui test` gets a finite p on a silently reduced df; and
  `test ES_ame_pre4 ES_ame_pre3 ES_ame_pre2` -> `p=., df=0`, rc 0. Both are the documented
  "test drops the constraint" behaviour (`sthlp:1616-1617`), so I do not file it. Finding 2
  is the same *mechanism* reaching a different consumer, where nothing is printed at all.
* `wooldid_vspl`'s splice: `e(wdid_vcvstate)=analytic-blocks`, `e(wdid_vblocks)=att att_it ame
  ame_it` — the ame block IS spliced and the sentinel diagonals ARE left verbatim, exactly as
  the code comment claims.
* Identical `e(V)` diagonals across `ES_ame_post0/1/2` (.1320351853) on iid data initially
  looked like the cell-containment "one number repeated" signature appearing with no
  `cluster()` at all — but it is legitimate here: `dose0_i` is time-invariant, so each cohort
  cell's within-cell dose variation is identical at every post period. With a time-varying dose
  (p2) the horizons' variances differ. Red herring, resolved.
