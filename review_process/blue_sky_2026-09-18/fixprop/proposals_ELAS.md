# Fix proposals — mechanism `ELAS` (elasticity / semi-elasticity denominators)

Findings covered: **FP-SEMIELAS-FE** (binary `semielasticity`) and **fp-contreat-A**
(`contreatelasticitytype(eydx|eyex)`).

**They are one defect with one fix site.** Both reach `margins` through the same prediction
plumbing and both fail for the same reason: on the default engine `predict, xb` after
`reghdfe` is `_cons + Xb`, and every absorbed effect (`i.i`, `i.t`, `fe()`, `ccc_absorb()`,
the absorbed dose ringer) is outside `e(b)` and so outside the prediction. `dydx`/`dyex`
**difference** the prediction, so the absorbed part cancels and they are exact — measured
exact. `eydx`/`eyex` **divide** by it, so they are not. One proposal (P1) fixes both; P2 is
the do-not-move-numbers fallback; P3 is a one-line companion worth doing under either ruling.

Everything below was measured on a patched copy of the pinned build. Workdir:
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-ELAS`
(patched build in `build_fix/`; probe do-files and logs `probe01`–`probe19`).

---

## The measurement that frames the ruling

Two exact-fit designs (`e(r2) == 1`, so `E[y|x] == y` and the truth is closed form) and one
noisy realistic panel. `customvce()`/`unconditionalse` write the fixed effects out as
explicit dummies, so their prediction *is* the fitted outcome; that arm is the reference.

| design / row | truth | default (reghdfe) | `customvce` | default error |
|---|---|---|---|---|
| binary exact-fit, `semielasticity` | 0.023688642186 | 0.019373047923 | 0.023688642186 | **−18.2 %** |
| binary exact-fit, `dydx` (control) | 1.453333333333 | 1.453333333333 | 1.453333333333 | 0 |
| contreat exact-fit, `eydx` | 0.019016464476 | 0.016409915568 | 0.019016464476 | **−13.7 %** |
| contreat exact-fit, `eyex` | 0.038068286220 | 0.032947192778 | 0.038068286220 | **−13.5 %** |
| contreat exact-fit, `dydx` / `dyex` (control) | 1.320000000000 / 2.658739710670 | identical | identical | 0 |
| noisy 40×6×10 panel, `semielasticity` | (no closed form) | 0.057985445692 | 0.076265354849 | **−24.0 %** |
| noisy panel, `eyex` | (no closed form) | 0.012265648125 | 0.017072811681 | **−28.2 %** |

`probe01_baseline.log`. The error is not a fixed rate: it is second-order in the absorbed
share of `y`'s level and ranges ~2 % to ~65 % across the designs this review has run.

**The reach is every ey\* row, not just `main_att`.** On the binary exact-fit design with
`esprelength(2) espostlength(3) subgroups(sg)` every one of the 8 non-zero rows differs from
its `customvce` twin (e.g. `ES_att_post0` 0.013378125946 vs 0.015477565352). On a noisy
`contreat + eyex + robust` run (`probe14_base.log` vs `probe14_fix.log`, 40 cohorts × 25
units × 15 periods) `main_ame` is 0.015327596504 against 0.018288888088 (**+19.3 %** after
the fix) and `ES_ame_pre3` **changes sign**, +0.000501404061 → −0.000195428748.

### A verification trap that cost this review a run — state it in the re-bless
`e(wdidmainresults)` row 2 is **`pretreatment_att`**, not the `ame` row. A timing probe read
`M[2,1]` on a `contreatelasticitytype(eyex)` run, saw it identical before and after, and
briefly concluded the contreat surface was untouched. It was reading a `dydx` row, which the
fix correctly leaves alone. Any re-bless comparison must key on **row names**
(`main_att`, `pretreatment_att`, `main_ame`, `pretreatment_ame`, `ES_att_*`, `ES_ame_*`),
never on row indices. `probe14_base.log` / `probe14_fix.log` dump them by name.

---

## P1 — PRIMARY: give `margins` the whole prediction on the default path

### The idea in one line
Reconstruct the sum of absorbed effects once, `d = y − xb − resid`, and point `margins`'
prediction at a wrapper that returns `xb + d`. That is the fitted outcome, so `eydx`/`eyex`
divide by the right thing and reproduce the explicit-dummy engine.

This is the arrangement the poisson path **already has**: `ado:1459` says *"ppmlhdfe requires
d for predict, which margins requires"* and `ado:1462/1468` pass `d` to `ppmlhdfe`, whose
`predict, mu/xbd` adds `e(d)` (verified in `ppmlhdfe_p.ado:66-70`). `reghdfe` never got the
twin. Note `ppmlhdfe_p`'s **`xb` does not add `d` either** — which is exactly why the
poisson-without-`poisexpresults` remap to `dydx` (BS-6) is safe, and why `poisexpresults`,
which uses the default `predict = mu`, is safe. Only the OLS path is exposed.

### The exact edits — six sites

**(1) `ado:1447-1449` → insert the gate immediately after the `vcespecifier` block.**

Now:
```stata
  local vcespecifier
  if "`robust'" == "robust" local vcespecifier vce(robust)
  if "`cluster'" != "" local vcespecifier cluster(`cluster')
```
Becomes (added lines only):
```stata
  * ---- ELAS-FIX (absorbed-effect denominator) ----------------------------------------
  * reghdfe's predict,xb is _cons + Xb: every absorbed effect (i.i, i.t, fe(),
  * ccc_absorb(), the absorbed dose ringer) is OUTSIDE e(b) and so outside the
  * prediction. dydx/dyex difference the prediction, so the absorbed part cancels and
  * they are exact. The ey* operators DIVIDE by it, so on the default path they report
  * the effect as a share of the partialled index rather than of the fitted outcome -
  * the quantity sthlp:243-245 and sthlp:905-907 define, and the quantity the regress
  * path already returns. Fix: ask reghdfe for the residual, reconstruct the absorbed
  * sum d = y - xb - resid once, and point margins' prediction at wooldid_pxbd (xb + d).
  * Gated on the operator that actually divides, so a dydx/dyex run pays nothing.
  local wdelas 0
  if "`poisson'" != "poisson" & "`unconditionalse'" == "" & "`customvce'" == "" {
    if "`semielasticity'" == "eydx" | "`contreatelasticitytype'" == "eydx" | "`contreatelasticitytype'" == "eyex" local wdelas 1
  }
  tempvar wdelasres wdelasd
  local wdelasopt
  if `wdelas' local wdelasopt residuals(`wdelasres')
```
The gate reads the **post-remap** operator locals, so poisson-without-`poisexpresults`
(already remapped to `dydx`/`dyex` at `ado:5034-5037`) never arms it, and the three
non-reghdfe engines are excluded explicitly.

**(2) `ado:1476` (the `reghdfe` call) — append one token.**

Now (tail of the line):
```stata
... `vcespecifier' tol(1e-`regtol')
```
Becomes:
```stata
... `vcespecifier' tol(1e-`regtol') `wdelasopt'
```
`residuals()` is empty unless the gate armed. `reghdfe 6.12.3` has no `d` option at all
(`reghdfe.ado:115` offers only `RESiduals(name)`/`RESiduals2`); a bare `d` copied from the
poisson line is refused r(198).

**(3) immediately after that `reghdfe` line, inside the same `else{` branch — build `d` and
register the wrapper.** Both `ereturn local`s go in **before** `qui est sto `wdfm''` at
`ado:1488`, so every later `est restore `wdfm'` (`ado:1797`, `3092`, `3132`, `3220`, `3432`,
`3478`) brings them back with the model.
```stata
      * ELAS-FIX: reconstruct the absorbed sum and register the prediction wrapper.
      * ORDER MATTERS - this predict must run BEFORE e(predict) is set, or it would
      * dispatch to the wrapper itself.
      if `wdelas' {
        qui predict double `wdelasd' if e(sample), d
        ereturn local wdid_dvar "`wdelasd'"
        ereturn local predict "wooldid_pxbd"
      }
```

**(4) NEW FILE `wooldid_pxbd.ado` (6 lines) — this cannot live inside `wooldid.ado`; see the
trap below.**
```stata
*! wooldid_pxbd : margins' prediction on the default reghdfe path -- xb + absorbed sum
program define wooldid_pxbd
	syntax newvarname [if] [in] [, *]
	tempvar wdpxb
	_predict double `wdpxb' `if' `in', xb
	gen double `varlist' = `wdpxb' + `e(wdid_dvar)' `if' `in'
end
```

**(5) `wooldid.pkg` — add one line beside `f wooldid.ado`:** `f wooldid_pxbd.ado`

**(6) a one-line preflight beside the gate**, so a stale install that copied only
`wooldid.ado` fails with a sentence instead of r(199) from inside `margins`:
```stata
  if `wdelas' {
    cap qui findfile wooldid_pxbd.ado
    if _rc != 0 {
      dis as error "Error: wooldid_pxbd.ado was not found on the adopath. semielasticity and contreatelasticitytype(eydx|eyex) need it on the default estimation path. Reinstall wooldid (the package is now two ado-files), or use unconditionalse or customvce()."
      error 601
    }
  }
```
`findfile` is one adopath walk, once per run, only when an ey\* operator was requested.

### THE TRAP (measured, and it cost two runs)
`wooldid_pxbd` **must be its own ado-file.** Defined as a subroutine inside `wooldid.ado` it
is callable from wooldid's own code — `probe09_dbg2.log`: `cap wooldid_pxbd` returns **rc
100** (`varlist required`), i.e. the command resolved and ran — but it is invisible to
`margins`' predict dispatcher, which fails `command wooldid_pxbd is unrecognized`, **rc 199**
(`probe05_dbg.log`). `prog list wooldid_pxbd` also cannot see it (rc 111). Ado subroutines are
scoped to the ado that loaded them and `margins`→`predict` is not in that scope. This is
exactly why `reghdfe` ships `reghdfe_p.ado` and `ppmlhdfe` ships `ppmlhdfe_p.ado` as separate
files.

### What it costs at runtime
**Nothing measurable, and nothing at all when the option is not requested.** The whole
arrangement is behind `if `wdelas''.

Wall clock, identical data and seeds, 40 cohorts × 15 periods, `esprelength(3)
espostlength(4)`, `robust` (`probe10_timing_base.log` / `probe10_timing_fix.log`):

| N | arm | pinned build | patched build | ratio |
|---|---|---|---|---|
| 15 000 | `att` plain (`dydx`, gate closed) | 7.673 s | 7.059 s | 0.92 |
| 15 000 | `att semielasticity` (gate armed) | 9.313 s | 9.047 s | **0.97** |
| 15 000 | `contreat + eyex` (gate armed) | 44.812 s | 40.325 s | **0.90** |
| 60 000 | `att` plain (`dydx`, gate closed) | 21.804 s | 22.789 s | 1.05 |
| 60 000 | `att semielasticity` (gate armed) | 30.450 s | 30.667 s | **1.01** |
| 60 000 | `contreat + eyex` (gate armed) | 140.380 s | 139.640 s | **0.99** |

Every ratio is inside the run-to-run noise of this machine (the closed-gate rows, which
execute byte-identical code, move 0.92 and 1.05). There is **no measurable cost**.

The reason it is free is the design choice that matters most here: **the wrapper keeps
`margins`' chain rule.** The prediction is still a function of the single linear index, so
`margins` differentiates it exactly as it does now. Measured in isolation on a 150-cell model
(`probe03_predprog.log`, timers 1/2/3):

| form | `margins, eydx(tf)` | correct? |
|---|---|---|
| today (`predict, xb`, no `d`) | **0.296 s** | no — 0.0580 vs 0.0763 |
| **P1 wrapper via `e(predict)`** | **0.302 s** | yes — 0.076265354849 |
| `expression(predict(xb) + d)` | **33.06 s** | yes, but forces `nochainrule` |

The `expression()` form is the obvious edit and it is **110× slower**; rejected below.

The only new per-run work is `residuals()` (one extra double column, 8 bytes × N — 0.5 MB at
60 k rows) and one `predict, d` (one `_predict` pass plus two `gen`s; 0.012 s at 2 400 rows,
`probe02_mech.log` timer 2). Both are O(N) once per run, against margins calls that are
already O(N × #cells) and already on the full data — the cell collapse is **already banned**
on exactly these paths (`ado:5086`, `ado:5121`), so nothing new is disabled.

### Does it move a reported number? **YES — this is the ruling.**
Every `eydx`/`eyex` point estimate **and standard error** on the default (reghdfe) engine
moves: `main_att`/`main_ame`, `pretreatment_*`, every ES horizon, every subgroup row and
subgroup contrast, `att_it`/`att_i`/`att_itime` twins, `e(histogramestimates)`, the `impvalues`
row, the `lattice()` rows and the cfx matrices. Measured magnitudes: **+22.3 %**
(0.019373 → 0.023689, binary exact-fit), **+15.9 %** (contreat `eydx`), **+15.5 %** (`eyex`),
**+31.5 %** (0.057985 → 0.076265, noisy realistic), **+39.2 %** (`eyex`, noisy), **+19.3 %**
(noisy `contreat + eyex` `main_ame`), and a **sign change** on `ES_ame_pre3`. SEs move with
them (noisy `eydx` se 0.001300962085 → 0.001745598763, +34 %).

**The post-fix SE is a *different estimand*, not a re-scaled one, and the author must rule on
that separately.** After P1 the default path reports se 0.001745598763 where `customvce`
reports 0.001822400665 — 4.2 % smaller, on a point estimate the two now agree on to 12
digits. That gap is **not** the reghdfe-vs-regress degrees-of-freedom difference, and my own
measurement rules that out: on the identical model the two engines' `dydx` standard errors
are **bit-identical**, 0.128352450346 on both (`probe16_fegap.log`, `L|dydx` rows), so no
common small-sample factor can open a 4.2 % gap on `eydx`. The cause is structural: the
wrapper hands `d` to `margins` as a **fixed constant**, so the delta-method Jacobian carries
no rows for the absorbed-effect coefficients, while the explicit-dummy path carries their
variance and their covariance with the cell coefficients. P1's interval is therefore a
*conditional-on-the-fixed-effects* interval. That is a defensible convention — it is what
treating the FEs as nuisance parameters means, and it is the same posture `reghdfe` takes
everywhere else — but it is a convention, it is 4 % narrower than the other engine's, and it
should be named in the help rather than discovered. I did not measure its coverage.

**Nothing else moves.** Verified bit-identical in `baseline.txt` vs `patched.txt` and
`probe14_base.log` vs `probe14_fix.log`: every `dydx` and `dyex` row on both engines, every
`customvce`/`unconditionalse` row, both poisson arms (`semielasticity` and
`poisexpresults semielasticity`, 0.242853462513 se 0.002254703421 on both builds), and — the
sharpest single datum — on a `contreat + eyex` run, where the wrapper **is** registered and
**is** used by every margins call, the whole att family comes back byte-identical *including
its standard errors*: `main_att` 6.101737830497 se 0.068978413268, `ES_att_post0`
4.553566435474 se 0.082701707575, and so on for all eleven att rows. The absorbed part
cancels in a contrast and the Jacobian is unchanged.

### The direction is *toward* the documented quantity, which is the argument for ruling it in
`sthlp:243-245` promises "the average percentage change ... between the two sets of predicted
outcomes"; `sthlp:905-907` defines `eyex` as "% change in y / % change in contreat". After the
fix the default-path answer **equals the closed-form truth on an exact-fit design** and
**equals the regress-path answer on realistic data**:

| | before | after | reference (`customvce`) |
|---|---|---|---|
| binary exact-fit `main_att` | 0.019373047923 | **0.023688642186** | 0.023688642186 (= truth, 12/12 digits) |
| all 8 ES + subgroup rows, same run | all differ | **all agree to 12/12 digits** | — |
| noisy `semielasticity` | 0.057985445692 | **0.076265354849** | 0.076265354849 (12/12 digits) |
| noisy `semielasticity`, `fe(i.unit)` | 0.057985445692 | **0.076268872945** | 0.076268872945 (12/12 digits) |
| noisy `semielasticity` `[aw=wt]` | 0.058257711419 | **0.076434677833** | 0.076434677833 (12/12 digits) |
| contreat exact-fit `eydx` | 0.016409915568 | 0.019016464560 | 0.019016464476 (rel 4.4e-9) |
| contreat exact-fit `eyex` | 0.032947192778 | 0.038068286387 | 0.038068286220 (rel 4.4e-9) |
| noisy `contreat + eyex + lattice(all 3)` `main_ame` | −0.193895614557 | −0.183863674753 | −0.183863674185 (rel 3.1e-9) |

**The ~4e-9 residual on the contreat arms, measured properly (this corrects an earlier
guess).** An earlier draft of this proposal said "the lever is `regtol()`, I have not swept
it." I have now swept it and **`regtol()` does not move it**: on the exact-fit contreat design
the relative error is 4.430e-09 at `regtol(8)`, `(10)`, `(12)` and `(14)` alike, while the
`customvce` twin sits at 6.5e-14 at every one (`probe14_base.log` / `probe14_fix.log`, part B).
Nor is it the `d` reconstruction: on a **matched** model `y − resid` reproduces the
explicit-dummy fitted value to **1.1e-15** (exact-fit design) and **6.4e-15** (noisy design),
tolerance-independent (`probe16_fegap.log`, `M|` rows). What is left is the only thing that
differs between the two engines on a `contreat()` run: wooldid puts the dose ringer in
`absorb()` on the reghdfe path (`contreatcontrolinserterABS`) and carries it as an explicit
regressor on the regress path (`contreatcontrolinserterFE`). Same span, different
parameterisation, so the two coefficient vectors agree only to solver precision and the
ratio inherits it. It shows up at ~3–4e-9 on every `contreat` arm and at exactly 0 (12/12
printed digits) on every binary arm. Above the project's 1e-10 threshold; on the record.

### `fe()` — a caveat this review raised and has now **retracted**
An intermediate battery showed the fixed default path at 0.076268872945 against `customvce`'s
0.076265354910 under `fe(unit)` — a 4.6e-5 gap that looked like reconstruction error. It is
not. With a **bare** `fe()` token the two engines fit **different models**, which
`sthlp:1051-1052` already documents ("a BARE token there is read as a CONTINUOUS variable"):
`e(r2)` is 0.9973888465 on the reghdfe side and 0.9970914908 on the regress side, so the
comparison was never valid. Written with the prefix the models match (`e(r2)` 0.9973888465
on both) and the fix is **exact**: `probe19_fix.log`, `fe(i.unit)`, robust 0.076268872945 vs
customvce 0.076268872945, 12/12 digits — against `probe19_base.log`, 0.057985445692 vs
0.076268872945. `fe()` absorbed effects flow into `d` correctly.

### Blast radius
- **Armed only when** the engine is reghdfe **and** an ey\* operator survived the BS-6 remap.
  A default run, a `dydx` run, a `dyex` run, every poisson run and every
  `customvce`/`unconditionalse` run take the identical code path they take today.
- **When armed, it is global to the fit** (it is `e(predict)`), so *every* margins call in the
  worker picks it up — the 40-odd sites at `ado:2030`–`3480`, `wooldid_hist` (`ado:1944`), the
  cfx block, `impvalues`' `over(cluster)` read, the `contrast(overjoint)` joint tests and
  `pwcompare`. That is a feature: it closes the surfaces the finding could only mark
  "INFERRED but not run", and no margins call can be missed.
- **Rows inside wooldid's sample but outside `e(sample)` do not poison anything** — the
  obvious objection to a `d` built `if e(sample)`, and it is measured clean
  (`probe14_*.log`, part D). (a) an outcome hole (`y` missing for one observation): fixed
  robust 0.076353014303 vs customvce 0.076353014303, 12/12 digits. (b) a never-treated unit
  observed in a single period, which reghdfe drops as a singleton so that `e(N)` is 2400 on
  the default path and 2401 on the regress path: fixed robust 0.076265354849 vs customvce
  0.076265354849, 12/12 digits. The dropped rows are never-treated and never enter a
  `targetflag==1` subpop, so the missing `d` there is never read.
- **No internal `predict` can be hijacked — verified, not assumed.** Setting `e(predict)`
  redirects *every* `predict` issued while the model is active, so the obvious hazard is a
  diagnostic elsewhere in the worker asking for `resid`/`xb`/`stdp` and silently receiving
  `xb + d`. `grep -n "predict" wooldid.ado` returns **no live `predict` command anywhere in
  the file**: every hit is either comment/display prose or margins' own `predict(xb)` /
  `expression(predict(xb))` option string at `ado:1846-1849`, and those are poisson-only,
  which the gate excludes. One latent hazard to note in the code comment: `wooldid_pxbd`
  declares `[, *]` and **ignores** whatever option it is handed, so it always returns
  `xb + d`. That is correct for every call that can reach it today, but a future edit that
  set `poismarginstext` off the poisson path would get `xb + d` for a `predict(xb)` request
  with no error. If the author wants belt and braces, the wrapper should `_assert` that no
  option other than `xb` was passed.
- **Not touched:** `e(V)` splicing, the AMESE R1 scrub, the collapse gate, the joint-test
  caution, the bootstrap/permute shells (the worker re-runs per replication, so the wrapper is
  rebuilt each time — `bootstrap(standard) reps(20)` verified, point 0.076265354849 matching
  the analytic run, se 0.001709855779 → 0.001907352859).
- **A behaviour change on pathological `y`, and it is a *convergence* not a *regression*.**
  Once the denominator is the fitted outcome rather than a grand-mean index, a `y` whose
  fitted values straddle zero fails instead of returning a number. Measured on a
  mean-centred outcome (`batt_base.txt`/`batt_fix.txt`, `Z|` rows): the default path goes
  from **r(3598)** to **r(459)**, which is exactly what `customvce(robust)` already returns
  on that data in *both* builds. So the fix does not introduce a failure — it replaces an
  obscure rc with the one the other engine already gives, and the case already failed. This
  is `sthlp:245`'s existing caveat ("tends to require ... a y-variable that is not sparsely
  distributed") finally biting where it should.

### Alternatives rejected
1. **`expression(predict(xb) + d)`** — correct, one token simpler, and **110× slower**
   (33.06 s vs 0.302 s on a 150-cell model) because `expression()` forces `nochainrule`, i.e.
   one interpreted pass over the sample per coefficient. It is also *less* accurate: its
   numerical Jacobian gives se 0.001745617064 against the wrapper's 0.001745598763. The ado's
   own comment at `ado:1822-1834` already prices this trade-off in the other direction for the
   poisson `dydx` sites. (I re-derived this route independently before finding the earlier
   measurement, and it is the trap a fresh reader will fall into: it is the obvious edit.)
2. **`margins ..., predict(xbd)` or `predict(d)`** — does not work and would fail silently in
   the direction of nonsense. `reghdfe_p.ado:50-58`: `xbd` is `y − resid` and `d` is
   `y − xb − resid`; neither responds to the counterfactual value of `targetflag`, so margins'
   contrast would be zero (`xbd`) or backwards (`d`).
3. **Refit with explicit dummies whenever an ey\* operator is asked for** — this is the
   att-rescue refit the author **removed on 2026-09-08**: a solve over hundreds of collinear
   dummy columns, "numerically far noisier than the partialled fit it replaced". It would cost
   the same and buy less, and it throws away the fast engine.
4. **Refuse ey\* on the reghdfe path and steer to `customvce()`** — breaks every existing
   `semielasticity` do-file, and the remedy changes the standard errors as well, so the user
   cannot get the documented point estimate on the default variance rule at all.
5. **Append a fixed-coefficient column to `e(b)`/`e(V)`** so `predict, xb` includes `d` —
   needs a full `ereturn post`, which would destroy the `e()` wooldid reads back (`e(df_r)`,
   `e(sample)`, `e(N_clust)`, `e(r2)`, reghdfe's DOF macros). Rejected as fragile, not wrong.
6. **Keep the wrapper inside `wooldid.ado`** — measured impossible (rc 199; see THE TRAP).

### The test that would prove it correct — no oracle needed
1. **Cross-engine identity (the strongest, and it needs no truth).** Same data, same
   specification, `robust` vs `customvce(robust)` and `cluster(v)` vs `customvce(cluster v)`:
   every `eydx`/`eyex` **point estimate** must agree to rtol 1e-10. The **standard errors
   must not** be asserted equal: P1's is conditional on the fixed effects and the
   explicit-dummy engine's is not (see above), so the right guard on the SE row is a
   direction-and-magnitude band (P1 narrower, measured 4.2 % here), not an identity. Today they do not agree at all.
   `probe01_baseline.do` / `probe04_patched.do` / `probe11_robust_*.do` are this test.
   **Write `fe()` with an explicit `i.` prefix in any such comparison**, or the engines are
   not fitting the same model (`sthlp:1051-1052`) and the test reports a false failure.
2. **Closed-form exact-fit truth.** Build `y` entirely inside the model's span (assert
   `e(r2) == 1` before trusting it) and compare against `mean[ln(yhat) − ln(yhat − δ_it)]` for
   the binary arm and `mean(β_it/y)`, `mean(β_it·x/y)` for the contreat arm.
   `probe01_baseline.do` and `gendata.do` build both.
3. **Metamorphic, no oracle at all:** scale the absorbed part of `y` by `s`. The finding's own
   graded probe (`fp-contreat/probe06_eydx_scale.do`) reads relative error 1.8e-9 at `s=0`
   (nothing to absorb) rising to 2.48e-2 at `s=3`; after P1 it must read ~1e-9 at **every** `s`.
4. **Regression guards that must not move:** `dydx`/`dyex` rows bit-identical *including SEs*;
   poisson and `poisexpresults` rows bit-identical; `customvce`/`unconditionalse` rows
   bit-identical. Verified (`baseline.txt` vs `patched.txt`, `probe14_*`).
5. **Surface battery — every default-path row must equal its `customvce` twin.** Run and
   measured here (`batt_base.txt` vs `batt_fix.txt`, `probe14/18/19`):

   | surface | before (default) | after (default) | `customvce` twin |
   |---|---|---|---|
   | `semielasticity` | 0.057985445692 | 0.076265354849 | 0.076265354849 ✓ |
   | `+ histogramcohorteffects` | 0.057985445692 | 0.076265354849 | 0.076265354849 ✓ |
   | `+ fastnose` | 0.057985445692 | 0.076265354849 | 0.076265354849 ✓ |
   | `+ esprelength(2) espostlength(3)` | 0.058889724980 | 0.078327559530 | 0.078327559530 ✓ |
   | `[aw=wt]` | 0.058257711419 | 0.076434677833 | 0.076434677833 ✓ |
   | `fe(i.unit)` | 0.057985445692 | 0.076268872945 | 0.076268872945 ✓ |
   | `+ impvalues cluster(unit)` | 0.056991718170 | 0.069523894293 | (refused on the regress path, r198 — no twin) |
   | `bootstrap(standard) reps(20)` | 0.057985445692 | 0.076265354849 | (= the analytic point) ✓ |
   | `contreat + eyex`, all 14 `e(b)` rows | all differ | all agree | to 3–7e-9 (see above) |
   | `contreat + eyex + lattice(all 3)` | −0.193895614557 | −0.183863674753 | −0.183863674185 (3.1e-9) |
   | mean-centred `y` (fitted values straddle 0) | r(3598) | r(459) | r(459) in both builds ✓ |

---

## P2 — FALLBACK if the author will not move the numbers: say what the number is

If P1 is declined, the defect must at least stop being silent. **This is strictly worse than
P1 and I would not choose it**, because it documents a number that is not the quantity the
same help file defines two paragraphs earlier — but it is cheap, moves nothing, and is better
than the present silence.

- **`sthlp:243-245`** (`semielasticity`), after "between the two sets of predicted outcomes",
  add: *"On the default estimation path the prediction `reghdfe` returns excludes the absorbed
  unit and time effects, so the percentage change is taken against the partialled index rather
  than against each observation's fitted outcome. `unconditionalse` and `customvce()`, which
  carry the fixed effects as explicit regressors, divide by the fitted outcome instead, and
  the two paths do not agree (measured 13–28 % on ordinary panels)."*
- **`sthlp:905-907`** (`contreatelasticitytype()`): the same sentence for `eydx`/`eyex`.
- **A runtime caution**, one `dis as text` gated exactly as P1's `wdelas`.
- Cost: zero. Moves no number. Blast radius: text only.
- Why I still prefer P1: the caution's own remedy (`customvce()`) changes the standard errors
  too, so the user cannot obtain the documented point estimate on the default variance rule at
  all. And the one message a default `semielasticity` run prints today is the collapse caution,
  whose text says "*the run is slower, and no answer changes*" — which a user will read as
  reassurance about precisely this.

---

## P3 — COMPANION (one string; do it under either ruling)

`ado:5140`, the `semielasticity` limb of the collapse caution, currently ends "...
`contreatelasticitytype(dydx)` and `(dyex)` are linear in the dose and stay eligible."

Under P1, append: *"Because they divide by the prediction, these operators also need the
whole fitted outcome: wooldid reconstructs the absorbed fixed effects and adds them back
before margins evaluates them (see semielasticity in the help)."*

Cost: zero (one string inside a block that only prints when the collapse is skipped). Moves
no number. This is the line that makes P1's number movement explicable to a user who has just
been told "no answer changes" and then sees a figure move 30 % between versions.

---

## What I did not test, stated plainly
- **`permute()` with `semielasticity`.** `probe18` tried `permute(standard)` and got r(198) —
  my option value was wrong, not a finding. `bootstrap(standard)` **was** run and passes.
  The worker re-runs per replication so the wrapper is rebuilt each time, but I did not
  demonstrate it on the permutation engine.
- **`makecfxplots` with an ey\* operator.** `probe18` returned r(198) from my own probe
  syntax. `lattice(all 3)` **was** run and passes, and both read the same
  `contreatelasticitytype` plumbing, so I expect it to be covered; I did not show it.
- **Very large panels.** Timing tops out at 60 000 rows / 40 cohorts. The added work is O(N)
  once per run and the margins work is unchanged in form, so the ratio should hold, but I have
  not measured at 10^6 rows.
- **Whether the ~3–7e-9 contreat residual can be removed** by also routing the dose ringer
  through `absorb()` on the regress path (or out of it on the reghdfe path). I diagnosed the
  parameterisation mismatch but did not test a change to it, and it is a separate edit.
- **The coverage of P1's conditional-on-the-FEs standard error.** I measured that it is
  4.2 % narrower than the explicit-dummy engine's and diagnosed why; I did not run a coverage
  simulation on either, and the project's own record already has the analytic SEs
  under-covering at `cluster() >= i`, so a simulation there would need the finer-cluster arm
  too.
- **Whether the author would rather define the semi-elasticity against the *untreated* fitted
  outcome** (`ln(yhat|1) − ln(yhat|0)`, what `margins, eydx` on a factor does and what both
  engines do) or against the observed `y`. P1 does not change that definition; it only makes
  the default engine compute the definition the package already uses.
