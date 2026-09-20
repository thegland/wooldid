# Fix proposals — mechanism `BAND` (finding SK-02)

Author: prop-BAND agent. Workdir `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-BAND`.
Everything below is measured on the pinned build `c1d8979` unless marked UNTESTED.

## Summary

| # | proposal | moves numbers | risk | status |
|---|---|---|---|---|
| **BAND-1** | Relative floor `se_r >= 1e-8 * se_obs_analytic` on the replication's own standard error, at the sites that share the screen — **four patched and run**: `boottp` (the printed table), the cell p, `bootpwt`, `bootband`; **one specified, untested**: `bootsgw`. Plus two cautions and two `e()` scalars | **YES, hugely**: printed `boottstat` CIs move ~12 orders of magnitude on affected runs | medium (needs a re-bless) | **implemented, run, before/after measured** |
| **BAND-1b** | The `e()` documentation and the two footer strings BAND-1 obliges | no | none | specified |
| **BAND-2** | *Preferred architecture:* one per-replication degeneracy verdict computed once, consumed by every studentized object — the only way to reach the joint / equality / adversarial statistics, which BAND-1 structurally cannot | YES, BAND-1's plus the joint p's | medium | kernel **sketched, not run**; the gap it closes is **measured** (joint p floored at 0.395 by 29 collapsed reps) |
| **BAND-3** | Help edit: `sthlp:528` promises the band "can never come out narrower than the interval it accompanies". It can, and does. Counterexamples measured | no | none | counterexamples measured |
| **BAND-4** | Three things NOT to fix, plus one correction to the finding's own record | no | none | argued |

**The headline is a re-scoping.** SK-02 was filed as a band-only, medium, "uncertain on
realism" defect. The measurement shows the same one-line screen is what prints a **95%
confidence interval of +/- 1e13 in the main results table**, on a fifteen-state panel
clustered at the state level, at `_rc == 0` with no caution. That is brief category 1.

---

## 0. What the measurement changed about the finding

SK-02 was filed as a **band-only, medium, "uncertain on realism"** defect: the
`simultaneousbands` critical value can reach 1e14. I re-ran its own arm A with
`bootsave()` and read the replication standard errors directly
(`p1_measure.do` -> `saveA.dta`, analysis `an1_A.txt`). Three things came back that the
original report could not see, and they change what needs fixing.

**(a) The degenerate class is a clean, separable population, not a tail.**
`log10( min_h se_rh / se_obs_h )` over the 293 listwise-complete replications of arm A:

| bucket | reps |
|---|---|
| [-15,-14) | 5 |
| [-14,-13) | 22 |
| [-13,-12) | 2 |
| **[-12, -3): EMPTY** | **0** |
| [-3,-2) | 16 |
| [-2,-1) | 41 |
| [-1,0) | 196 |
| [0,+1) | 11 |

Arm B (Gaussian, G=20) has the same shape: 3 reps at [-14,-12], then nothing until
[-2,-1). **A relative floor placed anywhere in 1e-12 .. 1e-3 selects exactly the same
replications.** That nine-order-wide empty gap is what makes a tolerance here principled
rather than arbitrary; it is the single most important fact in this document.

**(b) It is the denominator, and the mechanism is resample saturation.**
`max` over replications and horizons of `|b_rh - b_obs| / se_obs_h` = **16.58** — ordinary.
The blow-up replications have `se_r/se_obs ~ 1e-14` at *every* horizon at once. Crossing
the degeneracy against the resample multiplicities that `bootsave()` also writes
(`w_*` columns), by number of **distinct** clusters drawn:

| distinct clusters drawn | reps | of which collapsed |
|---|---|---|
| 6 | 1 | 1 |
| 7 | 10 | 3 |
| 8 | 44 | 6 |
| 9 | 86 | 12 |
| 10 | 94 | 6 |
| 11 | 41 | 1 |
| 12-13 | 23 | 0 |

`i` is `state` here, so every drawn cluster copy carries its own cohort x relative-period
dummy block. Few distinct clusters -> the model saturates -> residuals at machine zero ->
cluster-robust sandwich at ~1e-28 -> `se_r ~ 1e-14`, **positive**, so the shipped
`S[r,k] <= 0` screen passes it and `(b_r - b_obs)/se_r` runs to ~1e14.

**(c) THE FINDING IS MIS-SCOPED. This is not band-only, and the band is the least of it.**
Those 29 replications are ~10% of the distribution — four times the 2.5% tail — so they
reach the order statistics of **every printed pointwise percentile-t interval**. From the
arm-A run's own printed table (`p1_measure.log:175`), with **no caution and `_rc == 0`**:

```
             |   estimate  se_analy~c      t_pt         p     lb_95ci      ub_95ci
         att |  1.3761595   .48398569      2.84     0.332  -1.667e+13    2.129e+13
pretreat att |   .6075353   .18448361      3.29     0.495  -3.140e+12    5.216e+12
        pre1 |  .96054783   .24420566      3.93     0.481  -4.727e+12    1.389e+13
       post1 |  1.4192669   .46097059      3.08     0.312  -1.213e+13    3.367e+13
```

Every reported interval is +/- 1e12-1e13 and every p-value is inflated to the middle of
the range at |t| ~ 3-4. This is the reviewer's fourth verdict confirmed, on realistic
input, at severity **critical**, not medium. (One correction to that verdict: the p is
**not** 0 — F033's add-one floor holds — it is *too large*, 0.31-0.50, which is the
worse failure mode because it silently fails to reject.)

**The fix site the brief points at (`ado:16830`, the band) is one of five reachable instances of the
same one-line screen** — `wooldid_boottp` (`ado:14379`), the cell p (`15067`),
`wooldid_bootpwt` (`15493`), `wooldid_bootband` (`16830`) and `wooldid_bootsgw` (`15380+`).
I patched and ran the first four; the fifth is specified and untested (see edit (5)). A
sixth site, `wooldid_bootjwr` (`15238`), carries the same failure but cannot be reached this
way at all — that is what BAND-2 is for. The band is simply where it is loud enough to
notice.

**Is this pathological data?** No. Arm A is skewed (log-normal errors, sigma 1.3), but arm
B is plain **Gaussian at G=20** and still collapses 3 of 299 replications — the mechanism
fires on well-behaved data, it just sits below the 2.5% tail there and so does not reach
the printed order statistics. What arm A adds is not pathology but *rate*: at G=15 the
collapse rate crosses 2.5% and the printed interval goes with it. Fifteen states clustered
at the state level with `i` = state is an ordinary applied setup, and `sthlp:403` calls
clustering on a variable that nests the cohort variable — "`cluster(state)` where the
cohorts…" — "the common route"; reviewer verdict 2 goes further and calls it the help's
recommended setup, which I have not verified and do not need. (Both arms
fire the "reference period forced" caution; that is the panel being nine periods long, not
a defect trigger.)

---

## Proposal BAND-1 — a relative floor on the replication's own standard error

**Status: implemented and run.** Patch `probes/prop-BAND/patch.diff` (65 added lines, 10
touched), applied to `probes/prop-BAND/patched/wooldid.ado`; before/after runs
`run_pinned.do` / `run_patched.do` in the same directory.

### The rule

> A replication is admitted to a studentized object at column *k* only if its own standard
> error there is at least **1e-8 times the observed analytic standard error of the same
> column**. The observed analytic standard error is already required positive and finite at
> every one of these sites — it is the denominator of `t_obs` and the endpoint scale — so
> the floor costs no new quantity and is scale-free by construction.

Why 1e-8 and not something else: it is the package's existing relative tolerance
(`wooldid_vspl`'s diagonal agreement, the `1e-8*max(ev)` eigenvalue rank guards), and it
sits in the middle of the nine-order empty gap measured in §0(a). 1e-12 and 1e-3 give
**bit-identical** results on both arms. This is the opposite of a tuned constant.

### The edits

All five are the same screen. Line numbers are on the pinned build.

**(1) `ado:14379`, `wooldid_boottp` — the PRINTED pointwise percentile-t p and interval.**

```mata
    to = ob[1, k] / so
    tr = J(0, 1, .)
    for (r = 1; r <= R; r++) {
      if (A[r, k] >= . | S[r, k] >= . | S[r, k] <= 0) continue
      tr = tr \ ((A[r, k] - ob[1, k]) / S[r, k])
    }
```
becomes
```mata
    to = ob[1, k] / so
    tol = 1e-8 * so
    nde = 0
    tr = J(0, 1, .)
    for (r = 1; r <= R; r++) {
      if (A[r, k] >= . | S[r, k] >= . | S[r, k] <= 0) continue
      if (S[r, k] < tol) {
        nde++
        continue
      }
      tr = tr \ ((A[r, k] - ob[1, k]) / S[r, k])
    }
    if (nde > ndeg) ndeg = nde
```
with `tol, nde, ndeg` added to the `real scalar` declaration, `ndeg = 0` beside
`al = (100 - lvl)/200`, the signature gaining `| string scalar degl`, and
`if (args() >= 9) st_local(degl, strofreal(ndeg))` after `st_matrix(outname, out)`.

**(2) `ado:15067`, the studentized CELL p (`gteffects` / `e(cfx*)` cells).**
```mata
      use = selectindex((A[., k] :< .) :& (S[., k] :< .) :& (S[., k] :> 0))
```
becomes
```mata
      use = selectindex((A[., k] :< .) :& (S[., k] :< .) :& (S[., k] :>= 1e-8*so))
```
(`so` is the observed analytic cell standard error, required positive eight lines above.)

**(3) `ado:15493`, `wooldid_bootpwt` — the pairwise subgroup contrasts.**
```mata
  ok = (rowmissing(A) :== 0) :& (P :< .) :& (P :> 0)
```
becomes
```mata
  ok = (rowmissing(A) :== 0) :& (P :< .) :& (P :>= 1e-8*seobs)
```
(`seobs` is required positive in the function's first guard.)

**(4) `ado:16830`, `wooldid_bootband` — SK-02 as filed.** This one has an **ordering
wrinkle that will bite anyone who writes the obvious edit**: `AO`'s own positivity is not
established until the `sd[j] = AO[1, sel[j]]` block *below* the completeness loop. Write
`S[r,sel[j]] < 1e-8*AO[1,sel[j]]` unguarded and, on a family whose observed analytic
standard error is missing at some horizon, `S < 1e-8*.` is true for every replication,
every row drops, and the family is refused with the **wrong** reason ("fewer than two
replications are complete over every horizon") instead of the right one ("the observed
analytic standard error at horizon H is missing or non-positive, so it cannot
studentize"). So:

```mata
  idx = J(0, 1, .)
  ndeg = 0
  for (r = 1; r <= R; r++) {
    mx = 1
    dg = 0
    for (j = 1; j <= m; j++) {
      if (A[r, sel[j]] >= .) mx = 0
      if (stud) {
        if (S[r, sel[j]] >= . | S[r, sel[j]] <= 0) mx = 0
        else if (AO[1, sel[j]] > 0 & AO[1, sel[j]] < . &
                 S[r, sel[j]] < 1e-8*AO[1, sel[j]]) dg = 1
      }
    }
    if (mx & dg) ndeg++
    if (mx & !dg) idx = idx \ r
  }
  if (args() >= 12) st_local(degl, strofreal(ndeg))
  n = rows(idx)
```
`dg` is kept separate from `mx` so the count is "replications that were otherwise complete
and were dropped **only** by the floor" — the number the caution needs. The signature
gains `| string scalar degl`.

**(5) `ado:15380+`, `wooldid_bootsgw` — the across-subgroup equality W. UNTESTED; include
it or say why not.** This is the fifth instance of the screen and, unlike `wooldid_bootjwr`
(BAND-2's case), it **is** reachable by a per-column floor: its call at `ado:7129` passes it
both `aseobs` and `pwobs0`, so the observed scales are in hand as `ase` and `pwo`. Its own
`min(ev) < 1e-8*max(ev)` eigenvalue guard does **not** shield it — that guard runs on
`R :/ (sc' * sc)`, a correlation matrix, so it is scale-invariant and a uniformly collapsed
`Sr` sails through exactly as in `bootjwr`. The edits are
`if (min(sr) <= 0) continue` -> also `continue` when `sr[a] < 1e-8*ase[idx[a]]` for any `a`,
and `if (P[r, j] >= . | P[r, j] <= 0)` -> `| P[r, j] < 1e-8*pwo[j]`.
**I did not apply or test these**: no arm in this probe used `subgroups()`, so I would be
shipping an edit I never ran. Flagged rather than hidden — if BAND-2 is taken instead, its
row mask covers this site for free.

**(6) The two cautions and the two `e()` scalars** (full text in `patch.diff`):
- in `wooldid_bootwire`, beside the existing `wdbtmiss`/`wdbtlow` pair at `ado:10216-17`:
  a caution naming `up to N of the DONE completed replications` dropped for a collapsed
  standard error, what that means (a saturated resample), and that the design is one the
  pairs cluster bootstrap is not reliable on. Plus `e(boot_ndegse)`.
- in `wooldid_esband`, ahead of the existing D-7 "banded some families and not others"
  caution: the per-family list of discarded counts, plus `e(esband_ndeg_<fam>)`
  (length-checked at 32 characters like its `esband_crit_`/`esband_n_` siblings).

There is deliberately **no new refusal**. The object is computable on the surviving
replications and the existing `wdbtlow` caution already covers the case where the loss is
catastrophic (>50% of a column).

### Runtime cost — this is a free guard

Nothing is added that is not already being done.

- Edits (2) and (3) are colon comparisons against a scalar replacing colon comparisons
  against a scalar. `1e-8*so` / `1e-8*seobs` is one multiply computed **once**, outside
  the vector operation. **Element count unchanged; cost identical.**
- Edit (1) hoists `tol = 1e-8*so` out of the replication loop and adds **one double
  comparison per (replication, column) cell** to a loop that already visits that cell and
  already does a subtraction and a division there.
- Edit (4) adds one comparison and two array reads per (replication, horizon), with
  horizons in the tens.
- No new pass over the data, no new matrix allocation, no new `e()` matrix, no sort.

Worst realistic accumulator, `reps(2000)` x 500 coefficient columns: one million extra
double comparisons, under a millisecond, against 2000 full model refits with variance
estimation. The added work is **not measurable**. This is the shape of guard the author
has said he will take: it rides inside an existing loop rather than adding a pass.

### Does it move a reported number? YES — and this is the part that needs a ruling

It moves, on any `boottstat` run where at least one replication's own standard error
collapses: the printed `p`, `lb_#ci`, `ub_#ci` on every coefficient row; `e(boot_tstud)`;
the studentized cell p / interval in `e(gteffects)` and `e(cfx*)`; the pairwise subgroup
contrast rows; `e(esband_*)` and `e(esband_crit_*)`; and the dashed band on the plots.

**Measured, arm A (G=15 states, `i`=state, `cluster(state)`, log-normal errors, reps 299),
live before/after — `p1_measure.log` vs `run_patched.log`:**

| row | p before | p after | 95% CI before | 95% CI after |
|---|---|---|---|---|
| `att` | 0.332 | **0.263** | [-1.667e+13, 2.129e+13] | **[-4.153, 34.588]** |
| `pretreatment att` | 0.495 | 0.474 | [-3.140e+12, 5.216e+12] | [-2.090, 6.260] |
| `ES_att pre3` | 0.932 | 0.913 | [-8.710e+12, 6.946e+12] | [-127.44, 2.653] |
| `ES_att pre2` | 0.617 | 0.609 | [-2.813e+12, 1.738e+13] | [-5.353, 8.176] |
| `ES_att pre1` | 0.481 | 0.466 | [-4.727e+12, 1.389e+13] | [-4.976, 16.587] |
| `ES_att post0` | 0.427 | 0.368 | [-1.989e+13, 1.212e+13] | [-4.121, 17.573] |
| `ES_att post1` | 0.312 | 0.248 | [-1.213e+13, 3.367e+13] | [-2.097, 9.178] |

Band: `e(esband_crit_att)` **1.17727e+14 -> 369.99906**, `e(esband_n_att)` 293 -> 264,
new `e(esband_ndeg_att)` = 29, new `e(boot_ndegse)` = 29. Both new cautions fire.

The intervals move by **twelve orders of magnitude**; the p-values by 0.02-0.07. A
re-bless of any `boottstat` baseline whose data can produce a saturated resample is
required. **Note that 370 is still a vacuous critical value.** The floor makes the object
finite, reproducible and honestly reported; it does **not** make a 15-cluster percentile-t
well behaved, and nothing can — see §"What this does NOT fix".

### Blast radius

**Touched:** everything studentized under `boottstat` — the printed coefficient table's
`p`/`lb`/`ub` and `e(boot_tstud)`; the studentized cell p and interval in `e(gteffects)`
and `e(cfx*)`; the pairwise subgroup contrasts; `e(esband_*)` / `e(esband_crit_*)` /
`e(esband_n_*)` and the dashed plot layer.

**Untouched, by construction:** every point estimate; the `se_analytic` column (it is the
observed analytic standard error, the floor's own reference); `e(V)` and `e(df_r)` (the
standard-deviation family, so `test`/`lincom`/`esttab`/`coefplot`/`honestdid` exports are
bit-identical); the percentile / bc / bca families (they never touch `se_r`); every run
without `boottstat`; every `permute()` run; the analytic band engine `wooldid_esbanda`.

**Not reached, and this is the gap:** the **joint** pre-trend and **adversarial** line test
statistics under `boottstat` — `wooldid_bootjwr` runs inside a replication and cannot see
the observed analytic scale, so no per-column floor can reach it. See BAND-2. The
**across-subgroup equality** W (`wooldid_bootsgw`) *is* reachable — it receives `aseobs` —
but edit (5) above is untested, so as shipped-and-run BAND-1 does not cover it either.

### Alternatives rejected

1. **Revert the band's denominator to the coefficient-draw SD.** Rejected twice over:
   PS-45 (author ruling 2026-09-10) deliberately moved *away* from that because the band
   could then come out narrower than the pointwise interval; and it leaves the pointwise
   percentile-t — the bigger half of the defect — completely untouched.
2. **Winsorize or cap `t_r`.** Rejected: no principled cap, and it moves every number on
   every run. Dropping is right *because the replication's variance estimate is invalid*,
   not because the statistic is large.
3. **An absolute floor (`se_r > 1e-12`).** Rejected: not scale-free. A model whose outcome
   is measured in units where legitimate standard errors are ~1e-9 would lose every
   replication. The package's own IM guard (CORE-06) is likewise relative, deliberately.
4. **Withhold the object, or refuse the run.** Rejected: the object is computable on the
   surviving replications, withholding hides the diagnosis, and the package already has a
   graded response for this shape (`wdbtlow` at 50% loss). A caution plus a count is the
   proportionate answer and matches `e(boot_failed)`'s existing convention.
5. **Screen upstream, at harvest in `wooldid_bootgrab`, on the replication's own `e(r2)`
   or residual df.** The most tempting alternative, because "this replication fit exactly"
   is a property of the replication rather than of a column. Rejected: `bootgrab` cannot
   see the observed analytic standard error, so it has no scale to compare anything to; it
   would have to read engine-specific `e()` scalars whose names and meanings differ across
   `reghdfe` / `ppmlhdfe` / `regress` / `poisson`; and it would miss the near-saturated
   replication that is not exactly saturated. The floor catches the symptom directly and is
   engine-agnostic.
6. **A robust reference scale (median over replications of `se_rh`) instead of the observed
   analytic one.** Rejected: needs a second pass and a sort, and it would put the floor on
   a different scale from the endpoints, which use the observed analytic SE. With a
   nine-order gap a robust scale buys nothing.
7. **Nothing; document it.** Defensible for the band alone (opt-in, not printed in the
   table). **Not** defensible for the pointwise path: a printed 95% interval of +/- 1e13 at
   `_rc == 0` with no caution is brief category 1, a silently wrong number on realistic
   input.

### Tests that would prove it correct

**T1 — tolerance invariance (metamorphic, no oracle). ALREADY RUN.**
The floor's constant must not matter. Recompute every studentized object at
tolerance 1e-12, 1e-10, 1e-8, 1e-6, 1e-4 and 1e-3 on the same replication draws. All six
give **identical** retained sets and identical critical values on both arms
(`an1_A.txt`, `an1_B.txt`: arm A n=264, c=369.999 at all six; arm B n=292, c=12.3031 at all
six). This is the test that shows 1e-8 is not a tuned constant, and it is the one I would
put in the suite, because it is cheap and it fails loudly if anyone ever "improves" the
tolerance into the populated region.

**T2 — the no-op invariant.** On a design with enough clusters that no replication
saturates, every printed number must be **bit-identical** to the pre-patch build in `%21x`.
(`run_pinnedC.do` vs the arm C block of `run_patched.do`; result in §Results below.)

**T3 — the sup-band dominance invariant (algebraic, no oracle, and it is what forces the
five edits to ship together).** For a studentized sup band,
`c = q_lev(max_h |t_rh|) >= q_lev(|t_rh|)` for **every** horizon h, computed on the
**same** retained replications — the quantile of a maximum dominates the quantile of each
component, exactly. Ship the floor in `wooldid_bootband` but not in `wooldid_boottp` (or
vice versa) and the two objects are built on different retained sets and this identity
**fails**. It is a one-line assertion over `e(bootdist)` / `e(bootdist_se)` and it needs no
truth. **RUN, offline, on the arm A and arm B draws:**

| configuration | arm A | arm B |
|---|---|---|
| shipped build (no floor either side) | holds | holds |
| BAND-1, floor on both sides | holds | holds |
| band floored, pointwise **not** (the mis-ship) | **FAILS at all five horizons**, `c` = 370 against `q95(|t_h|)` = 3.5e13-5.8e13, a factor of 1e11 | holds |

So T3 is not a formality: it is the test that catches the exact partial-application mistake
a reviewer reading only the SK-02 brief would make, and it catches it by eleven orders of
magnitude. It costs one quantile per horizon over matrices already in `e()`.

**T4 — coverage, and the reason the current build cannot pass it.** 500 draws at G=15 with
a known `att`; measure the realized coverage of the printed `boottstat` interval. The
pre-patch answer **would be** 100.0% at nominal 95% — not because the interval is
conservative but because it is +/- 1e13 and covers every conceivable truth, so that
prediction is arithmetic, not a measurement. A 100% coverage rate is not
a pass, it is the signature of a degenerate interval. NOT RUN (500 x 299 refits is far
outside this budget); flagged as the acceptance test for the ruling.

---

## Proposal BAND-2 — one replication verdict instead of per-column screens (PREFERRED)

**Status: UNTESTED — code sketched, not run. It is the design I would ship if the author
will take the larger diff.**

BAND-1's weakness is that it screens **per column**, so it can only reach objects whose
kernel can see the observed analytic standard error for that column. There is one family
it therefore **cannot** reach: the joint pre-trend, across-subgroup equality and
adversarial statistics under `boottstat`. `wooldid_bootjwr` (`ado:15238`) runs *inside* a
replication:

```mata
    R = (R + R')/2
    sc = sqrt(diagonal(R))'
    if (min(sc) <= 0 | missing(sc)) continue
    R = R :/ (sc' * sc)                                      // -> a CORRELATION matrix
    ev = symeigenvalues(R)
    if (min(ev) < -1e-8*max(abs(ev)) | max(ev) <= 0 | min(ev) < 1e-8*max(ev)) continue
    Si = invsym(R) :/ (sc' * sc)                             // scale restored
    W = d' * Si * d
```

Every guard there is **scale-invariant**: a uniformly collapsed `V_r` (all of `sc` at
~1e-14, which is exactly what §0(b) measured) has a perfectly well-conditioned correlation
matrix and sails through, producing `W_r ~ 1e28`. The joint p is
`(1 + #{W_r >= W_obs}) / (n + 1)`, so **every collapsed replication pushes the p up**: with
29 of 294 collapsed, the joint pre-trend p acquires an effective floor near 0.10 and the
test cannot reject at any conventional level. That is the *conservative* direction, which
is why nothing has ever noticed it — but it means a real pre-trend can be masked. And
`bootjwr` cannot see the observed analytic scale, so no per-column floor can fix it.

### The design

Compute the degeneracy verdict **once**, per replication, and let every studentized object
exclude the same rows.

```mata
// SK-02: one pass over the accumulated replication standard errors against the observed
// analytic ones. Returns the number of degenerate replications and writes an R x 1 0/1
// column. Called once in wooldid_bootwire, right where `wdbts' is established.
real scalar wooldid_bootdeg(string scalar seaccname, string scalar aseobsname,
                            string scalar exclname, real scalar tol, string scalar outn)
{
  real matrix S, ase, ex, D
  real scalar R, K, r, k, nd, nk, ndeg
  S = st_matrix(seaccname); ase = st_matrix(aseobsname); ex = st_matrix(exclname)
  R = rows(S); K = cols(S)
  D = J(R, 1, 0); ndeg = 0
  for (r = 1; r <= R; r++) {
    nd = 0; nk = 0
    for (k = 1; k <= K; k++) {
      if (ex[1, k] != 0) continue
      if (ase[1, k] >= . | ase[1, k] <= 0) continue
      if (S[r, k] >= . | S[r, k] <= 0) continue
      nk++
      if (S[r, k] < tol*ase[1, k]) nd++
    }
    // PERVASIVE, not any-one-column. The measured collapse takes a replication's WHOLE
    // covariance at once (se_r/se_obs ~ 1e-14 at all five horizons of arm A
    // simultaneously), while a single column can lose its standard error for ordinary
    // reasons - a resample that lost dose variation in an ame column, a horizon this
    // draw cannot identify. So the verdict is a MAJORITY of the columns that had
    // anything to say, which is what separates "this replication's variance estimate
    // collapsed" from "this column is awkward in this draw".
    if (nk > 0 & nd*2 >= nk) { D[r] = 1; ndeg++ }
  }
  st_matrix(outn, D)
  return(ndeg)
}
```

It needs its own `cap mata: mata drop wooldid_bootdeg()` line in the drop list beside the
other kernels at `ado:14230-14270`; a new Mata function without one behaves differently
under `cleanmatrices` (SI-17/SI-18).

Consumers, each gaining one optional `degn` argument and one `& D[r] == 0` (or
`:& (D :== 0)`) in its existing screen: `wooldid_boottp`, the cell-p kernel,
`wooldid_bootpwt`, `wooldid_bootband`, and — the one BAND-1 cannot reach — the joint /
equality / adversarial add-one shares in `wooldid_bootjnt` (`ado:6905`, `6931`), where
`select(st_matrix("`jwacc'")[., `jk'], ...)` gains the same row mask.

### Why it is better than BAND-1

- It reaches the joint / equality / adversarial family, which BAND-1 structurally cannot.
- **Every studentized object is then built on one retained set.** That makes the sup-band
  dominance invariant (T3) hold by construction instead of by coincidence, and it removes
  a whole class of "these two objects used different masks" defects rather than one.
- One rule, one count, one caution, one `e()` scalar, one place to change the tolerance.
- It is not more expensive: one `R x K` pass replaces the per-column re-derivation, and the
  accumulator is already read whole by every kernel.

### Cost, number-moving, blast radius

Cost: one `R x K` double comparison pass, done once. At `reps(2000)` x 500 columns that is
one million comparisons against 2000 refits — not measurable. Number-moving: the same rows
as BAND-1 on the coefficient side (so the same table as above), **plus** the joint,
equality and adversarial p-values, which BAND-1 leaves alone. Blast radius: BAND-1's, plus
`e(wdidjointtests)`, `e(boot_jointW)`-derived p's and the adversarial line tests.

### The evidence — measured, not inferred

`run_joint.do` (arm A + `jointtests` + `bootsave`, pinned build) confirms it exactly:
**29 of the 75 usable joint-test replications carry `W_r > 1e10` (max 1.29e30), and those
29 are precisely the 29 replications whose standard error collapsed** — a 1:1 overlap.
The joint pre-trend p on that run therefore has a hard floor of `(1 + 29)/(75 + 1)` =
**0.395**. Full numbers in §Results.

**What I still did not do:** I did not implement or run BAND-2 itself. The kernel above is
a sketch; the fix-site analysis and the evidence are solid, the code is not tested.

---

## Proposal BAND-3 — a help edit: `sthlp:528` promises a guarantee the construction does not give

**Status: measured counterexamples in hand. Help-only. Moves no number. Risk: none.**

The `simultaneousbands` entry, `sthlp:528` (the sentence occurs exactly once in the file,
so the edit anchors cleanly), says of the `boottstat` band:

> "...a genuine percentile-{it:t}, on the same scale as the pointwise intervals printed
> beside it, **so the band can never come out narrower than the interval it accompanies**"

That clause is false, and it contradicts a sentence eight lines above it in the same
paragraph, which gets it right for the other display families:

> "So where the printed pointwise interval is a rank interval and not symmetric — under
> {cmd:bootstrap()} with the default {cmd:bootdisplay(percentile)}, or with
> {cmd:bootdisplay(bc)} or {cmd:bootdisplay(bca)} — the band does not nest it endpoint for
> endpoint: one of the band's two ends can sit inside the printed interval, even where the
> band is the wider of the two."

The `boottstat` pointwise interval is **also** a rank interval and **also** not symmetric —
the run's own caution says so: *"the interval is equal-tailed from the signed t-quantiles
and is not symmetric about the estimate"*. Excluding `boottstat` from that sentence and
then promising the opposite for it is the error.

**Why it is false.** The pointwise interval is equal-tailed from *signed* quantiles,
`[b - q97.5*se, b - q2.5*se]`; the band is symmetric, `b +/- q_lev(max_h |t_h|)*se`. On a
skewed studentized distribution `q95(|t_h|)` can be smaller than `|q2.5(t_h)|`: if most of
the exceedance mass sits in one tail, `P(|t| > |q2.5|)` can be below 0.05. The sup over
horizons pushes `c` up but is not guaranteed to close the gap.

**Measured counterexamples** (both from `an1_A.txt` / `an1_B.txt`, quantities in `t` units
so they are directly comparable):

| arm | horizon | pointwise half-widths (-qlo, qhi) | band `c` | |
|---|---|---|---|---|
| B (Gaussian, G=20), post-patch | `post1` | (12.995, 2.249) | **12.303** | band narrower |
| A (log-normal, G=15), post-patch | `pre3` | (13.687, **651.63**) | **370.0** | band narrower |
| B, **shipped build** | `post1` | (13.705, 2.318) | 13.799 | holds by **0.7%** |

So on the shipped build the claim survives on my arms only by a 0.7% margin, and it fails
outright once the collapsed replications are removed. It is not a property of the
construction; it is luck.

**The edit.** In `wooldid.sthlp`, `simultaneousbands` entry, replace

```
a genuine percentile-{it:t}, on the same scale as the pointwise intervals printed beside
it, so the band can never come out narrower than the interval it accompanies;
{it:e(esband_engine)} reads {it:bootstrap-studentized} there.
```
with
```
a genuine percentile-{it:t}, on the same scale as the pointwise intervals printed beside
it. The scale is shared, but nesting is not: that pointwise interval is equal-tailed from
the signed {it:t}-quantiles while the band is symmetric about the estimate, so on a skewed
studentized distribution one of the band's two ends can sit inside the printed interval -
exactly as it can under {cmd:bootdisplay(percentile|bc|bca)}, and for the same reason.
{it:e(esband_engine)} reads {it:bootstrap-studentized} there.
```

This also settles finding **BT-01** ("band narrower than the pointwise percentile-t"),
which the reviewer graded PARTIAL/medium on the grounds that the narrowness it saw was
SK-02's inflation. It is not: BT-01's mechanism is independent and survives the fix.

---

## BAND-4 — three things that should NOT be fixed, and one correction to the record

### (a) The 13.8 critical value at G=20 with Gaussian errors is not a defect. NO CHANGE.

The sweep cell that made SK-02 look like it fires on ordinary data (`r08_bandsweep`
BCELL3, `e(esband_crit_att) = 13.80` at G=20 with plain Gaussian errors) is **not** this
mechanism. Measured: arm B has only 3 collapsed replications out of 299 — 1%, well inside
the 2.5% tail, so it touches neither the pointwise intervals (`ES_att_pre3` prints
`[-2.31, 1.45]`, entirely sane) nor materially the band. The replication that actually
*sets* `c = 13.799` has standard-error ratios 0.19-0.57 and `|db|/se_obs` up to 3.1 — an
ordinary, legitimate few-cluster draw. Under the floor `c` goes 13.799 -> 12.303, an 11%
move and the same order of magnitude.

So reviewer verdict 2 ("a band 4-5 times wider than the one wooldid's other engine would
draw on the same data") is measuring the **engine gap** — a percentile-t sup over five
horizons at G=20 against an analytic-simulated band from a normal-theory covariance — and
that gap is the PS-45 ruling working as designed. The percentile-t is wider because the
denominator really is unstable at 20 clusters; that is the information, not the bug. If
anything is owed here it is a sentence in the help (BAND-3 already adds one), not code.

### (b) "A run that produced no coefficients still exits rc 0" is a different mechanism.

Two of the four reviewer verdicts filed under SK-02 are about a run that estimates nothing
exiting `_rc == 0` with `e(cmd) = wooldid`, so a `capture` guard passes and the next
`esttab` / `_b[main_att]` dies on "matrix e(b) not found". That is real, but it is not the
BAND mechanism and shares no fix site with it. Out of scope here; it belongs with whichever
mechanism owns the estimability-failure exit path.

### (c) Correction to the record: the p-value is not 0, it is too LARGE.

Reviewer verdict 4 says a user "gets a printed 95% confidence interval of plus or minus
1e13 **and a p-value of 0**". The interval half is exactly right and I reproduce it. The
p-value half does not reproduce and cannot: F033's add-one convention
`2*(1 + min(m0,m1))/(n+1)` floors the studentized coefficient p at `2/(n+1)`, and the run's
own footer says so. What actually happens is worse to read and easier to miss: the 29
garbage replications land on **both** sides of `t_obs`, so they inflate `min(m0, m1)` and
the p-values come out **too large** — `att` p = 0.332 at `t_pt` = 2.84, `ES_att pre1`
p = 0.481 at `t_pt` = 3.93. Inverting the add-one formula on the `att` row:
`min(m0, m1)` = 48 before and 34 after, so 14 of the 29 removed replications sat on the
minority side of `t_obs` and 15 on the majority side. That split is the expected one:
`t_r = (b_r - b_obs)/se_r` takes its **sign** from the coefficient, and the coefficients in
these replications are perfectly sensible (`|db|/se_obs` of 2-11); only the **magnitude** is
roundoff. So the collapsed replications scatter symmetrically about `b_obs` and land on
both sides — which is precisely why they inflate `min(m0, m1)` rather than shifting the p
in one direction. A p of 0 would at least look wrong. A p of 0.33 beside a t of
2.84 looks like an honest non-rejection, and it is the failure mode a user will act on.

Worth recording because it is the only part of the finding that the measurement refutes.

---

## What this does NOT fix, and should be said out loud

After BAND-1, arm A's band critical value is **370** and `ES_att pre3`'s printed interval is
`[-127.4, 2.65]`. Those are finite, reproducible and honestly cautioned — and still
useless. The residual heavy tail is the 16 replications at `se_r/se_obs` in 1e-3..1e-2:
not numerically dead, just badly conditioned. That is the pairs cluster bootstrap being
badly sized at 15 clusters, which `sthlp` already states at length ("**Neither mode is a
remedy for few treated clusters** ... the few-cluster literature finds the pairs cluster
bootstrap badly sized there").

So the honest claim for BAND-1 is: **it converts a number that was determined by roundoff
into a number that is determined by the data.** It does not make a 15-cluster percentile-t
trustworthy, and no guard can. The new caution says this in as many words, and deliberately
points at `permute()` and at more clusters rather than at another option.

---

## BAND-1b — the documentation and consistency work BAND-1 obliges (ships with it)

BAND-1 adds two `e()` scalars and changes which replications the printed numbers come
from. Three consequences, each of which is a doc-vs-code divergence by this project's own
standard if left undone:

1. **The `e()` list needs `e(boot_ndegse)` and `e(esband_ndeg_`*`)`**, and the `boottstat`
   and `simultaneousbands` entries each need one sentence stating the rule: *a replication
   whose own standard error at a column is below 1e-8 of the observed analytic one is
   excluded from that column's studentized distribution, because its variance estimate has
   collapsed rather than merely being small; the count is in `e(boot_ndegse)` and, per
   band family, in `e(esband_ndeg_`*`)`.*

2. **`e(bootdist)`, `e(bootdist_se)` and `bootsave()` still hold the excluded rows.** The
   accumulator is deliberately untouched — that is right, because the raw draws are the
   evidence and a user must be able to see the collapse for themselves (it is how I
   diagnosed this). But the footer currently invites reproduction from those matrices, and
   after BAND-1 a by-hand reproduction that does not apply the floor gets the **old**
   +/-1e13. One sentence is needed: *the replication draws in `e(bootdist_se)` are the raw
   harvest and include the excluded rows; reproducing the printed interval by hand requires
   applying the same floor.* The T3 assertion must apply the floor for the same reason.

3. **The existing `wdbtlow` caution's reason clause becomes partly false.** It says
   replications were dropped "because the remaining replications produced no usable
   standard error for that coefficient" — after BAND-1 some are dropped because the
   standard error was usable-looking and collapsed. Amend to "...produced no usable
   standard error for that coefficient, or produced one that had collapsed (see the
   caution above)". One string.

4. **A posting convention to state or change.** As patched, `e(boot_ndegse)` and
   `e(esband_ndeg_`*`)` are posted **only when the count is non-zero** — so a script that
   tests `e(boot_ndegse) == 0` gets missing on a clean run, not 0. That follows
   `e(esband_ndeg_`*`)`'s band-family siblings (`esband_crit_`*`) which are posted only for
   families that banded, but it is the **opposite** of `e(boot_failed)`, which always posts.
   My preference is to switch `e(boot_ndegse)` to always-post (it is a run-level scalar like
   `boot_failed`, and always-post is what a script wants) and leave the per-family
   `esband_ndeg_`*` conditional. Either way it must be documented, because a silent absence
   is exactly the shape of defect this section exists to prevent.

Cost: zero runtime. Moves no number. Risk: none. But it is not optional — undocumented
`e()` and a footer that no longer describes the arithmetic are exactly the class the
2026-09-10/11 documentation review was run to remove.

---

## Results — what was actually run

All files in `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-BAND/`.
The patched build is `patched/wooldid.ado` (`patch.diff` against the pinned build); every
run asserts `which wooldid` and the logs show the patched path.

| run | do-file | log | what it shows |
|---|---|---|---|
| measurement | `p1_measure.do` | `p1_measure.log`, `saveA.dta`, `saveB.dta` | pinned build, arms A and B with `bootsave()`; the printed +/-1e13 table; the raw replication standard errors |
| offline analysis | `an1.py` | `an1_A.txt`, `an1_B.txt` | the bimodality, the resample-saturation cross-tab, the tolerance-invariance sweep, the T3 dominance check |
| patched | `run_patched.do` | `run_patched.log` | arms A, B, C on the patched build |
| pinned control | `run_pinnedC.do` | `run_pinnedC.log` | arm C on the pinned build, for the no-op invariant |
| joint probe | `run_joint.do` | `run_joint.log`, `saveJ.dta` | arm A + `jointtests`, pinned, for BAND-2's evidence |

**Arm A (G=15, log-normal, `cluster(state)`, `i`=state, reps 299)** — live before/after in
the number-moving table above. `e(esband_crit_att)` **1.17727e+14 -> 369.99906**
(`%21x` `+1.71ffc2445b50bX+008`), `e(esband_n_att)` 293 -> 264, `e(esband_ndeg_att)` = 29,
`e(boot_ndegse)` = 29, both new cautions printed. The live numbers reproduce the offline
recomputation from `saveA.dta` to every digit printed, which is the cross-check that the
patch implements the rule I analysed.

**Arm B (G=20, Gaussian, reps 299)** — `e(esband_crit_att)` **13.799023 -> 12.303097**
(`%21x` `+1.89b2f8495cdc0X+003`); 3 replications dropped. Again matches the offline
recomputation exactly. This arm is the evidence for BAND-4(a): the mechanism fires here
too, but at 1% it is below the 2.5% tail, so the printed pointwise intervals were already
sane and the band's 13.8 was a legitimate few-cluster tail, not this defect.

**Arm C (G=30, Gaussian, reps 299)** — the no-op invariant, **PASSES**. All 299
replications completed and all 299 were retained by the band (`e(esband_n_att)` = 299);
**zero** were dropped by the floor, so `e(boot_ndegse)` and `e(esband_ndeg_att)` are not
posted and neither caution fires. Pinned vs patched:

```
pinned  : C_CRIT=+1.423c63fea5eb7X+002  dec=5.0349359   C_NBAND=299 NDEG=. NDEGSE=.
patched : C_CRIT=+1.423c63fea5eb7X+002  dec=5.0349359   C_NBAND=299 NDEG=. NDEGSE=.
e(wdidmainresults), all 35 %21x lines: BIT-IDENTICAL
```

Not "agrees to rtol 1e-10" — **bit-identical in `%21x`**, which is what a floor that fires
on nothing must produce. It also incidentally shows what a healthy `boottstat` sup band
looks like: `c` = 5.03 over five horizons at 30 clusters, which is the number arm A should
have been in the neighbourhood of and was not.
(`cmpC.sh run_pinnedC.log run_patched.log` reproduces the comparison.)

**Joint probe (arm A + `jointtests`, pinned build)** — **BAND-2's mechanism is MEASURED,
not inferred.** `run_joint.log`, `saveJ.dta`:

```
jointW_att : n = 75 usable replications,  max W_r = 1.28979e+30,  #W_r > 1e10 = 29
reps whose pre-horizon standard error collapsed (se_r < 1e-8*se_obs) = 29
overlap between the two sets = 29 of 29        <- EXACTLY the same replications
printed: "p-value for joint test that all pre-treatment att horizon estimates are 0: 0.9474"
```

The 29 collapsed replications are **precisely** the 29 whose joint-test statistic exceeds
1e10, and they carry `W_r` up to 1.3e30. Since the joint p is the add-one share
`(1 + #{W_r >= W_obs}) / (n + 1)` and every one of those 29 exceeds any finite `W_obs`,
this run's joint pre-trend p has a **hard floor of 30/76 = 0.395**: no pre-trend, of any
magnitude, could have produced a joint p below 0.395 on this data. Removing them drops the
usable count to 46 and the floor to 1/47 = 0.021.

That is the conservative direction, which is why it has never shown up — an over-large
joint p reads as "no pre-trend", and a masked pre-trend leaves no trace. BAND-1 does
**not** fix it (`wooldid_bootjwr` cannot see the observed analytic scale); BAND-2 does.

Incidentally: only **75 of the 295 completed replications** contributed a joint statistic
at all, the rest being rejected by `bootjwr`'s own guards. I did not chase that, but a 75%
loss on the joint family is worth someone's attention independently.

### Honest list of what I did NOT do

- **No coverage simulation** (T4). It is the acceptance test and it is outside this budget.
- **BAND-2 is not implemented or run.** Its *mechanism* is now measured (the joint probe
  above), but the proposed kernel is a sketch that has never been compiled.
- **Edit (2), the studentized cell p, is patched but never exercised.** No arm requested
  `gteffects` or `makecfxplots`, so the Mata compiled but the branch never ran.
- **`subgroups()` was never in the picture.** The pairwise-contrast edit (3) and the
  subgroup band families are patched but untested — no arm had `subgroups()` — and edit
  (5), `wooldid_bootsgw`, is specified but deliberately NOT applied for the same reason.
- **`poisson` / `ppmlhdfe` untested.** `boottstat` is refused with `poisson` at parse, so
  the studentized path should be unreachable there, but I did not confirm it.
- **One more site with the same screen shape, reasoned but not probed:** `ado:230`
  (`permtci`'s quantile, `(sacc :> 0)`) and `wooldid_permstat`
  (`if (A[r,k] < . & S[r,k] < . & S[r,k] > 0)`). `permute()` permutes onsets rather than
  resampling clusters with replacement, so no cluster is duplicated and the
  resample-saturation mechanism should not arise the same way. I did not test it, and it
  is the obvious next probe if the author takes BAND-1.

---

## Open questions for the coordinator / author

1. **BAND-1 or BAND-2?** BAND-1 is tested and ships today; BAND-2 is the right shape and
   is the only one that reaches the joint / equality / adversarial family, where the defect
   is *measured* (a joint pre-trend p floored at 0.395 by 29 collapsed replications). My
   recommendation is BAND-2, with BAND-1's patch as the fallback if the larger diff is not
   wanted this cycle — the two are compatible, since BAND-2 replaces BAND-1's per-column
   test with a row mask at the same sites, and reaches the two BAND-1 cannot.
2. **The re-bless.** Every `boottstat` baseline whose data can produce a saturated resample
   moves, and it moves by twelve orders of magnitude. The author has to rule. My view: the
   current numbers are not numbers, so this is a correction, not a regression — but it is
   his call and it is the one thing in this document that cannot be decided by measurement.
3. **The acceptance test I could not run:** 500-draw coverage of the printed `boottstat`
   interval at G=15, before and after. Worth the compute before ruling, because it is the
   only thing that turns "the interval is finite now" into "the interval is an interval".
4. **The one site I reasoned about and did not probe:** `ado:230` (`permtci`'s quantile) and
   `wooldid_permstat` carry the same `S > 0` screen. `permute()` permutes onsets instead of
   resampling clusters with replacement, so no cluster is duplicated and the
   resample-saturation mechanism should not arise — but that is an argument, not a
   measurement, and it is the obvious next probe.
5. **Out of scope here:** two of SK-02's four reviewer verdicts are about a run that
   estimates nothing exiting `_rc == 0` with `e(cmd) = wooldid`. Real, but a different
   mechanism with no shared fix site; it needs an owner elsewhere.
6. **Independently worth a look:** only 75 of 295 completed replications contributed a
   joint-test statistic at all on arm A, the rest rejected by `wooldid_bootjwr`'s own
   guards. A 75% loss on the joint family is not the BAND mechanism and I did not chase it.
