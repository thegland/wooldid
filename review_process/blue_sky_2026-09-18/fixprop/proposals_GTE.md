# Fix proposals — mechanism `GTE` (gteffects cell variances and weights)

Pinned build `.../bluesky/build/wooldid.ado` (sha c1d8979, md5 `444892481be8d50ec962a102d6abef51`).
Every `ado:` / `sthlp:` line number below is that build's.

**Everything in this document was patched into a private copy of the ado and RUN.**

| build | contents | logs |
|---|---|---|
| `base` | pinned, unmodified | `probes/prop-GTE/run_base.log`, `run2_base.log` |
| `p1` | GTE-1 + GTE-2 only | `run_p1.log` |
| `p2` | p1 + GTE-4 (first form) | `run_p2.log`, `run2_p2.log` |
| `p3` | p2 + GTE-3 (first form) | `run_p3.log` |
| **`pf`** | **GTE-1 + GTE-2 + GTE-3 + GTE-4, the shipping candidate** | **`run_pf.log`, `run2_pf.log`** |
| — | cost of the added passes, 2,000,000 rows | `cost.log` |
| — | 80-draw Monte Carlo that fixes the predicate | `a1_vcearms.log` |

`which wooldid` is printed at the head of every log and resolves to the intended directory in
each. Patchers (line-anchored; they assert the original text before replacing it):
`probes/prop-GTE/patch.py` (p1/p2), `patchf.py` (pf). Probe bodies `body.do`, `body2.do`;
drivers `run_{base,p1,p2,p3,pf}.do`, `run2_{base,p2,pf}.do`.

**Headline no-change control.** Every `e(gteffects)` cell of every arm that does not contain the
defect is **bit-identical** between `base` and `pf` — estimate, se, lb and ub each dumped in
`%21x` per row and compared:

```
P1B cluster(state) below i      39 rows   est 0  se 0  lb 0  ub 0
P1C robust, cohort-level i      39 rows   est 0  se 0  lb 0  ub 0
P2A cluster(st), unit i        109 rows   est 0  se 0  lb 0  ub 0
P2B robust,      unit i        109 rows   est 0  se 0  lb 0  ub 0
P2C iid,         unit i        109 rows   est 0  se 0  lb 0  ub 0
C1  att att_i att_it att_itime  57 rows   est 0  se 0  lb 0  ub 0
C2  poisson                     45 rows   est 0  se 0  lb 0  ub 0
C4  bootstrap(standard), cluster below i  30 rows  est 0  se 0  lb 0  ub 0
```
**The `estimate` column never moves anywhere, on any arm, including the defect arms.** Only
`se`, `p`, `lb` and `ub` move, and only on rows the floor identifies.

---

# NEW FINDING, not in either brief entry — read this first

`run2_base.log` arm **C5**: the SKEP-NEW-2 data with
`wooldid y cohid year ttre, att gteffects cluster(cohid) bootstrap(standard) reps(40) boottstat`.
On the **pinned build**, 24 of the 30 `e(gteffects)` cells are perfectly healthy — and every one
of them is returned with a studentized 95% interval endpoint of order 1e13:

```
   i relyr        estimate        se        lb(base)      ub(base)  |   lb(pf)    ub(pf)
   2     0        1.868828    0.14717      -5.367e+12       2.046   |   0.8159     2.046
   2     1        1.191397    0.18447           1.115    8.863e+13  |    1.115     1.511
   2     2        1.921828    0.02045       -9.500e+10       1.941  |    1.627     1.941
   3     0        1.758141    0.06394           1.725    1.643e+12  |    1.725     1.791
   4     0        1.629200   3.369e-15           1.629       1.629  |  (missing)
```
`grep -c Caution` over the base C5 region returns four cautions, none of which is about this
(few clusters, failed replications, and the standard boottstat provenance note).

**Mechanism, confirmed.** The replications pass `gteffects` through `repopts` (`ado:342`), so each
replication produces its own `e(gteffects)` and `wooldid_bootcellkeys` harvests the **analytic**
cell standard error from it (`ado:9155`) to studentize by. In the replications whose covariance is
rank-deficient, that analytic standard error is residue (~1e-15), the replication's t-ratio is of
order 1e15, and that value sets the studentized quantile for **every** cell. The confirmation is
that under `pf` the package's own existing caution at `ado:6705` fires and names it —
"Caution: 6 gteffects/cfx cell(s) have no analytic standard error to studentize by ... left blank"
(`run2_pf.log:1223`) — and it does **not** appear anywhere in `run2_base.log`.

So GTE-1 is not only a cosmetic repair of a table nobody may read. It removes a silently wrong
number from a documented inference path (`sthlp` gteffects block: "the standard error is ... its
analytic standard error under {cmd:boottstat}"). I would rate it **critical** and flag it for the
coordinator separately from SKEP-NEW-2.

---

# GTE-1 + GTE-2 — put a relative floor on the gteffects cell variance, count it, say it
*covers:* **SKEP-NEW-2**, and the new boottstat finding above. One fix site.

## The defect on the pinned build

`run_base.log` arm P1A — 40 states, 20 years, three adoption waves plus 10 never-treated,
`i` = the onset cohort, `cluster(cohid)`:

```
P1A|rows=39|semiss=0|sezero=0|setiny=9|minse= 1.035e-15|maxse= 3.369e-01
```
Nine of 39 cells carry `se` of order 1e-15, a `p` down to 1.8e-45 and an interval whose endpoints
equal the point estimate. The `sthlp:1148` sentence says those cells carry a **missing** `se`.
The guard at `ado:13497` is `se = (v < . & v > 0) ? sqrt(v) : .` — an **exact**-zero test, and the
rank-deficient sandwich leaves residue of order 1e-30, not zero.

## The edit

### 1a. `wooldid.ado:13497`, inside `mata: wooldid_gtx()`

**Now:**
```mata
        se = (v < . & v > 0) ? sqrt(v) : .
```

**Becomes:**
```mata
        se = (v < . & v > 0) ? sqrt(v) : .
        // GTE-1: the help's "degenerate zero variances are always reported as missing,
        // never zero" needs a FLOOR, not an exact-zero test. A rank-deficient robust or
        // cluster-robust stage-1 sandwich leaves cancellation RESIDUE (~1e-30) on the
        // affected diagonal; it sqrt()s to ~1e-15 and prints as an ordinary row with a
        // p-value of order 1e-45 and a zero-width interval. Under boottstat the same
        // residue is what each replication studentizes by, which sends the studentized
        // quantile - and so the interval of every OTHER cell - to order 1e13.
        // Two SCALE-FREE floors (an absolute floor is not equivariant in the units of y):
        //   (a) se < 1e-10*|est| - the project's own IM degeneracy rule, ado:2005;
        //   (b) v <= 1e-12*max(cell variances) - the rank tolerance the eigenvalue
        //       screens use (ado:14136/14150/14189), which also catches a residue cell
        //       whose point estimate is itself ~0.
        // (a) is load-bearing on its own: when EVERY cell is residue (a two-cluster
        // design) the anchor in (b) is residue too and (b) sees nothing.
        if (se < . & (se < 1e-10*abs(est) | v <= 1e-12*gtvmx)) {
            se   = .
            gtnd = gtnd + 1
        }
```

### 1b. `wooldid.ado:13490`, the anchor for floor (b), immediately before the cell loop

**Now:**
```mata
    out = J(n, 15, .)
```

**Becomes:**
```mata
    out = J(n, 15, .)
    // GTE-1: scale anchor for the degeneracy floor - the largest CELL variance on the
    // stage-1 diagonal. One pass over the n cell columns (n is the table's own row
    // count), no pass over the data.
    gtvmx = 0
    gtnd  = 0
    for (r = 1; r <= n; r++) {
        v = sed[1, cells[r, 3]]
        if (v < . & v > gtvmx) gtvmx = v
    }
```

### 1c. `wooldid.ado:13473`, the declaration

**Now:** `    real scalar j, r, q, df, g, tt, est, v, se, z, p, lb, ub, crit, n`
**Becomes:** `    real scalar j, r, q, df, g, tt, est, v, se, z, p, lb, ub, crit, n, gtvmx, gtnd`

### 1d. two `st_local` lines so the caller can see the count

After `ado:13484` (`st_local("gtn", strofreal(n))`) add `st_local("gtndegen", "0")` — it must be
set before the `n == 0` early return. Immediately before `out = sort(out, (1, 4))` add
`st_local("gtndegen", strofreal(gtnd))`. And at `ado:1767` (`if _rc == 0 {` in the gteffects
block) add `local gtndegen 0` as the first line inside the brace, so a failed `mkmat` leaves the
local defined.

### 2. `wooldid.ado:4169-4173`, the posting site — the count, the caution, `e(gt_ndegen)`

**Now:**
```stata
  if "`gteffects'" != "" {
    if `gtn' > 0 {
      cap ereturn matrix gteffects = `gtm'
      if _rc != 0 dis as text "Caution: e(gteffects) could not be posted (rc `=_rc')."
    }
```

**Becomes** (the GTE-4 block below inserts into the same brace; this is the GTE-1/2 half):
```stata
  if "`gteffects'" != "" {
    if `gtn' > 0 {
      cap ereturn matrix gteffects = `gtm'
      local gtprc = _rc
      if `gtprc' != 0 dis as text "Caution: e(gteffects) could not be posted (rc `gtprc')."
      * GTE-1/2: report the degeneracy floor, in the shape NFES-02 uses just above for
      * e(V_es_*). Costs one comparison on a local already in hand. Under a plain
      * bootstrap()/permute() the cell inference columns are REPLACED by replication
      * quantities (:6609, :8121), so the scrubbed analytic value is not what the user
      * reads and the caution would be false; e(gt_ndegen) still records it. Under
      * boottstat the scrubbed value IS the reported one and is studentized by, so the
      * caution says so in those terms.
      if `gtprc' == 0 & `gtndegen' > 0 {
        cap ereturn scalar gt_ndegen = `gtndegen'
        if "`bootstrap'`permute'" == "" dis as text "Caution: `gtndegen' of `gtn' e(gteffects) cell(s) have a stage-1 variance that is cancellation residue rather than an estimated variance (the robust covariance is rank-deficient in that cell's direction). Their se, p and interval columns are reported as missing rather than as the standard error of order 1e-15 the residue implies; the point estimates are unaffected. bootstrap() or permute() rebuild the cell inference from replications, and a cluster() strictly finer than the cohort-by-period cells avoids the degeneracy."
        else if "`boottstat'" != "" dis as text "Caution: `gtndegen' of `gtn' e(gteffects) cell(s) have a stage-1 variance that is cancellation residue rather than an estimated variance. Under boottstat that analytic standard error is the reported one and is what every replication studentizes by, so those cells report a missing se, p and interval rather than a standard error of order 1e-15; the point estimates are unaffected."
      }
    }
```
Note `local gtprc = _rc` replaces the inline `` `=_rc' `` so the second test cannot be disturbed;
the rendered text of the existing message is unchanged.

### 3. `wooldid.ado:11628`, the note printed above the table in `wooldid_display`

Append one line after the existing analytic note, gated so it never contradicts the bootstrap or
permute rewiring:
```stata
    if "`e(gt_ndegen)'" != "" & "`e(gt_ndegen)'" != "0" & "`e(boot_method)'`e(perm_method)'" == "" dis as text "Note: `=e(gt_ndegen)' cell(s) below carry a degenerate stage-1 variance (cancellation residue, not an estimated variance); their se, p and interval columns are missing by design."
```
It is existence-gated on `e()`, so replay reproduces it and older stores print nothing.

## Measured before/after

| arm | base | pf |
|---|---|---|
| P1A `i`=cohort, `cluster(cohid)` | 39 rows, **9 with se 1.0e-15 … min p 1.8e-45** | 39 rows, **9 missing**, min real se 2.188e-03, `e(gt_ndegen)=9` |
| P5 two-cluster design | 17 rows, **all 17 se in [9.5e-16, 2.4e-14]** | **all 17 missing**, `e(gt_ndegen)=17` |
| C5 `boottstat` | 6 residue cells **and 24 healthy cells with ±1e13 interval endpoints** | 6 missing; healthy intervals sane; `ado:6705` caution now names the blanking |
| P1B, P1C, P2A/B/C, C1, C2, C4 | — | **bit-identical** |

## Cost

**Zero measurable.** Both floors are Mata arithmetic on objects already in registers. Floor (b)'s
anchor is one extra loop over `n` — the *row count of the table*, tens to a few thousand — with no
pass over the data. On the 2,000,000-row panel the gteffects block's six existing `gegen` passes
cost 0.141 s **each**; this adds nothing to any of them.

## Does it move a reported number?

**Yes, and only on rows the floor identifies.** `e(gteffects)` columns `se`, `p`, `lb_*ci`,
`ub_*ci` go from residue to missing on degenerate cells. Under `boottstat` it additionally moves
the `lb`/`ub` of *non-degenerate* cells, because dropping replications whose studentizing standard
error is residue changes the studentized quantile — measured 25 `lb` and 11 `ub` of 30 rows on C5.
That movement is a repair (from ±1e13 to sane), but it is a movement and needs a re-bless.
The `estimate` column, `e(b)`, `e(V)` and every aggregate row are untouched on every arm tested.

## Blast radius

`e(gteffects)` and everything downstream of it: the printed cohort-by-relative-period table
(`ado:11616`), the bootstrap cell harvest (`ado:9152-9155`, `wooldid_bootcellsx` at `ado:6609`),
the permute blanking (`ado:8121`, which blanks these columns anyway). Nothing reads
`e(gteffects)` back into an aggregate. `makecfxplots` / `e(cfx*)` are continuous-treatment
surfaces and `gteffects` is refused with `contreat()` (`ado:4695`), so they cannot interact.

## Alternatives rejected

- **A single floor relative to the largest cell variance** (`v <= 1e-18*max(cell v)`) — this is
  what an earlier agent in this directory built (`probes/prop-GTE/patched/`, `NOTES.md`). It is
  **measurably insufficient**: on the two-cluster design (`probe2.log`, K=1) *every* cell variance
  is residue (largest 6.0e-28), so the anchor is residue too and 17 of 17 degenerate rows escaped
  — `K1_ndegen=` came back empty while all 17 cells printed se ~1e-15. Clause (a), relative to the
  cell's own estimate, is the one that catches that case. I keep (b) anyway for the cell whose
  point estimate is itself ~0, where (a) cannot fire.
- **An absolute floor** (`v <= 1e-20`) — rejected: not equivariant in the units of `y`. Rescale
  `y` from dollars to millions and a legitimate variance falls under any absolute floor.
- **Dropping the degenerate rows from the table** — rejected: the help promises the row with a
  missing `se`, the point estimate is still correct and still carries the weight columns the
  documented weighted-mean identity needs, and dropping rows would break that identity.
- **Rescaling the residue into a usable standard error** — rejected, and `sthlp:255` already rules
  it out in general: "at that level the analytic covariance is the wrong object rather than a
  mis-scaled one".

## The test that proves it

Two invariants, no oracle:

1. **Metamorphic, scale.** Run any design twice with `y` and with `1000*y`. Every cell's
   `se/|estimate|` ratio is invariant, so the set of cells the floor scrubs must be **identical**
   across the two runs. An absolute floor fails this; both of the proposed floors pass it by
   construction. (Not run — proposed as the regression test.)
2. **Cross-engine identity, already run.** On a design where `cluster()` is strictly finer than
   `i`, `e(gteffects)` must be bit-identical between the pinned build and the patched one. Done:
   eight arms, 537 rows, `est/se/lb/ub` all zero differences (table at the top).
3. **The closed-form degenerate design, already run.** A cohort whose every observation lies in
   one cluster has, by construction, no within-cluster variation to identify its cell coefficients:
   every one of its cells must return a missing `se`. `run_pf.log` P1A and P5.

---

# GTE-3 — the all-zero-VCV screen tests for exact zero too
*covers:* the escape case that GTE-1 clause (a) catches at the cell level but nothing catches at
the run level. Touches a **ruled** message (AMESE R4, shipped 2026-08-24; the refit removed
2026-09-08), so it needs the author.

## The defect

`run_base.log` arm P5: 40 states, one treated wave, `i` = cohort, `cluster(cohid)` — **two
clusters**. The whole stage-1 covariance is cancellation residue: `probe2.log` reports the largest
cell variance at 6.0e-28 and the largest non-cell variance at 2.8e-30. The run prints

```
P5_ATT b= 2.31227e+00  se= 1.60299e-15  p= 4.41339e-16
```

`wooldid_amescreen, mode(vzero)` (`ado:12978`) exists precisely to catch this and reads
`max(abs(e(V))) == 0`. Residue is not zero, so it stays silent and neither of its two cautions
prints. The run is not *wholly* silent — the `sescan` belt-and-braces guard does print
"8 reported standard errors appear to be identical under a cluster-type VCE" — but nothing says
that every standard error in the run is cancellation residue.

## The edit — `wooldid.ado:12978`

**Now:**
```stata
      mata: st_local("vzero", strofreal(max(abs(st_matrix(st_local("vchk")))) == 0))
```

**Becomes:**
```stata
      * GTE-3: a FLOOR, not an exact-zero test. A cluster-robust sandwich that is
      * rank-deficient in EVERY direction returns cancellation residue, not zeros -
      * measured max|e(V)| = 6.0e-28 on a two-cluster binary design in which every
      * reported standard error was ~1e-15 and this screen stayed silent. The floor is
      * relative to the fit's own largest coefficient, so it is equivariant in the units
      * of y (an absolute floor is not); 1e-12 on the standard-error scale leaves about
      * twelve orders between it and the smallest standard error a real fit produces.
      mata: __wdvzb = abs(st_matrix("e(b)")); __wdvzb = select(__wdvzb, __wdvzb :< .); st_local("vzero", strofreal(max(abs(st_matrix(st_local("vchk")))) <= (1e-12 * (cols(__wdvzb) > 0 ? max(__wdvzb) : 0))^2))
      cap mata: mata drop __wdvzb
```
The `select` is not decoration: `max()` over a vector containing a missing returns missing in
Mata, and `x <= .` is true for every finite `x`, so an unfiltered `e(b)` would make this screen
fire on **every** run with an omitted coefficient.

Because the screen now fires where the covariance is residue rather than zeros, the two caution
texts it controls become false as written ("every standard error … is missing" — they are 1e-15,
not missing). Both are amended in the same patch, substring-exact:

- `ado:12984`: `Caution: estimation returned an all-zero variance-covariance matrix. This usually reflects a degenerate clustering` → `Caution: estimation returned a variance-covariance matrix that is all zeros or cancellation residue (no standard error in it exceeds 1e-12 of the largest coefficient), so none of the standard errors, test statistics, p-values or intervals below carries information - a residue covariance prints standard errors of order 1e-15 and p-values of order 1e-45 rather than missing values. This usually reflects a degenerate clustering`
- `ado:1567`: `… is missing.` → `… is missing or, under residue, is a value of order 1e-15 with a p-value of order 1e-45 - read none of them.`

## Measured before/after

`run_base.log` / `run_p1.log` / `run_p2.log`: neither vzero caution ever prints on P5.
`run_p3.log` and `run_pf.log`: both print (`run_pf.log:2704` and `:2712`). No other arm in any
body changes: the screen fires on P5 and nowhere else across the eleven arms.

## Cost

Two Mata expressions on `e(b)` and `e(V)`, both already in memory, once per run. Unmeasurable.
`e(b)` for a large design is a few thousand elements.

## Does it move a reported number?

**No.** `r(vzero)` controls only which cautions print. The residue standard errors are still
reported exactly as before. That is the honest limitation: GTE-3 makes the run *say* its
covariance is worthless; it does not scrub it.

## Alternatives rejected

- **Scrub every standard error in the run when the screen fires** — the R1 machinery to do it
  exists. Rejected *for this proposal*, not on the merits: it turns a caution-only patch into one
  that blanks every inference column of a whole run, which is a different order of ruling. Offered
  as variant B if the author wants it.
- **Anchor the floor on `e(rmse)`** (the physically right scale) — rejected because `e(rmse)` does
  not exist on the `ppmlhdfe` path, and the fallback would make the screen behave differently
  under `poisson`.
- **An absolute floor** — same objection as GTE-1: not equivariant in the units of `y`.

## The honest gap

At K=1 the anchor of GTE-1 floor (b) is itself residue, so a cell whose point estimate is smaller
than about 1e-5 would slip floor (a) as well and keep its residue standard error. I could not
construct that on realistic data (every cell in the measured K=1 design has an estimate of order 1),
but the gap is real and GTE-3 only cautions about it, so say so rather than claim coverage.

---

# GTE-4 — say that a contained cell's standard error does not cover
*covers:* **bbr-N1-gteffects-cell-ci-covers-45pct**. Caution and `e()` marker only; no number moves.

## The measurement that changes the fix

`a1_vcearms.log`, 80 draws, 30 states x 10 years, onsets 2004/2006/2008, 9 never-treated,
`y = a_i + 0.03*(year-2000) + 1.0*post + 0.7*N(0,1)`. For the cell `(i = 1, relyr = 0)`, whose true
effect is exactly 1:

| arm | `i` | vce | cell se/sd | cell 95% coverage | att se/sd |
|---|---|---|---|---|---|
| 1 | state | `cluster(st)` | 0.291 | **0.438** | 0.865 |
| 2 | state | **`robust`** | **0.636** | **0.700** | 1.071 |
| 3 | state | iid (default) | 1.095 | 0.963 | 0.955 |
| 4 | cohort | `cluster(st)` | 0.919 | 0.938 | 0.944 |

Arm 1 reproduces bbr-N1 (it measured 0.287 and 0.450 on different seeds). **Arm 2 is the point.**
It has no `cluster()` at all, so the finding's own framing ("confined to `cluster()` at-or-above
`i`") and the `known_bydesign` §1 entry ("Analytic SEs are understated when `cluster()` is AT or
ABOVE the level of `i`") do not reach it — and its cell intervals still cover 70% at nominal 95%,
while its `att` row is fine (1.071). **I am filing that as a contradiction of a ruled entry**, per
the coordinator note: the cell-level failure is not a clustering-level phenomenon. It is a
saturation phenomenon. A treated cell whose observations are exhausted by the cell's own dummy has
residuals identically zero there, so it contributes **nothing** to the sandwich meat, whatever the
sandwich is clustered on. Arm 4, where the cell holds eight independent states, is fine.

That is why a nesting-based screen — the `jtcoarse` route at `ado:1414`, which is what an earlier
agent in this directory used — is the wrong instrument: it is silent on arm 2.

## The right predicate already exists in the file, and cannot currently reach these designs

`wooldid_amescreen, mode(detect)` (`ado:12864-12878`) computes exactly this:
```stata
qui gegen double `amecl' = nunique(`cv') if `esamp' & `tflag' == 1, by(`icell' `tcell')
qui gegen byte `ametag' = tag(`icell' `tcell') if `esamp' & `tflag' == 1 & `amecl' == 1
qui count if `ametag' == 1
```
It is gated at `ado:1534` to `"`contreat'" != ""`, and `gteffects` is **refused** with `contreat()`
(`ado:4695`). The two are mutually exclusive by construction: **no `gteffects` run has ever met the
package's own cell-containment screen.** That is the coherence gap in one sentence.

## The edit — `wooldid.ado:1764`, immediately before the `gttag` pass in the gteffects block

**Now:**
```stata
    qui gegen byte `gttag' = tag(`_i' `_t') if `gtin'
```

**Becomes:**
```stata
    * GTE-4: cell-containment disclosure for the BINARY designs. The AMESE R1 screen
    * (:1533) computes exactly this predicate but is gated to contreat(), and gteffects is
    * refused with contreat() (:4695), so no gteffects run has ever met it. A treated cell
    * that contributes no independent residual unit to the sandwich meat gets a cell
    * standard error far below the estimate's true sampling sd, and a cell interval that
    * does not cover at its nominal level (the measured factors are in the help).
    * SKIPPED INSIDE A REPLICATION (`wdrepmode'): the replications pass `gteffects' through
    * repopts, so an ungated screen would cost these passes reps+1 times.
    * Cost when it runs: one gegen nunique + one gegen tag + count per cluster variable,
    * every one restricted to `gtin' (TREATED rows only) - measured 0.238s on a 2,000,000
    * row / 805,000 treated panel, against 0.141s for each of the SIX by-(i,t) gegen passes
    * this block already runs. Under a robust VCE it is one tag + count on a variable
    * already built (0.128s). Under an iid VCE there is no pass at all.
    local gtcontain 0
    local gtcmode
    local gtcvce
    if "`fastnose'" == "" & "`wdrepmode'" == "" {
      local gtcclus
      if "`cluster'" != "" local gtcclus `cluster'
      local gtccvw = word(`"`customvce'"', 1)
      if "`gtcclus'" == "" & strlen("`gtccvw'") >= 2 & strpos("cluster", "`gtccvw'") == 1 {
        local gtcclus : subinstr local customvce "`gtccvw'" "", word
      }
      if "`gtcclus'" != "" {
        local gtcmode cluster
        local gtcvce cluster-robust
        tempvar gtccl gtcct
        foreach gtcv of local gtcclus {
          cap drop `gtccl'
          cap drop `gtcct'
          qui gegen double `gtccl' = nunique(`gtcv') if `gtin', by(`_i' `_t')
          qui gegen byte `gtcct' = tag(`_i' `_t') if `gtin' & `gtccl' == 1
          qui count if `gtcct' == 1
          if r(N) > `gtcontain' local gtcontain = r(N)
        }
        cap drop `gtccl'
        cap drop `gtcct'
      }
      else if "`robust'" != "" | "`weight'" == "pweight" {
        local gtcmode observation
        local gtcvce heteroskedasticity-robust
        tempvar gtcct2
        qui gegen byte `gtcct2' = tag(`_i' `_t') if `gtin' & `gtnob' <= 1
        qui count if `gtcct2' == 1
        local gtcontain = r(N)
        cap drop `gtcct2'
      }
    }
    qui gegen byte `gttag' = tag(`_i' `_t') if `gtin'
```
The robust branch costs nothing new: `gtnob` is the per-cell observation count the block already
builds at `ado:1756` and returns as the `nobs` column. `weight` is a live local in
`wooldid_worker` (44 references), so the pweight arm of the test is real, and `gtcmode` carries the
wording so a pweight run is not told "cluster".

And in the posting brace at `ado:4169` (same brace as GTE-2, printed first):
```stata
      * GTE-4: the containment disclosure. `gtcontain' is 0 unless the screen above ran.
      if `gtprc' == 0 & `gtcontain' > 0 {
        dis as text "Caution: `gtcontain' of `gtn' e(gteffects) cell(s) lie wholly inside a single `gtcmode', so they contribute no independent residual unit to the `gtcvce' variance: their se, p and interval columns are understated by a design-dependent factor and the cell interval does not cover at its nominal level. The aggregate rows are a separate question - see the Standard Errors section of the help. The remedies are a cluster() strictly finer than the cohort-by-period cells, or bootstrap() or permute(), which rebuild the cell inference from replications."
        cap ereturn scalar gt_ncontained = `gtcontain'
      }
```
The measured factors deliberately stay out of the runtime text and go in the help (GTE-5b): they
come from one 80-draw Monte Carlo on one design, and the package's register puts measured ranges
in the help and mechanism-plus-remedy in the caution.

## Measured before/after

`run_pf.log`, which arm fires and with which wording:
```
P1A  i=cohort cluster(cohid)  -> "39 of 39 ... lie wholly inside a single cluster"      NCONT=39
P2A  i=state  cluster(st)     -> "109 of 109 ... lie wholly inside a single cluster"    NCONT=109
P2B  i=state  robust          -> "109 of 109 ... lie wholly inside a single observation" NCONT=109
P1B  i=cohort cluster(state)  -> silent                                                  NCONT=.
P1C  i=cohort robust          -> silent                                                  NCONT=.
P2C  i=state  iid             -> silent                                                  NCONT=.
C1   cluster below i          -> silent        C8  i=cohort robust  -> silent
```
Exactly the arms A1 shows are broken, and no others. **No `e(gteffects)` number moves** —
`run_base.log` vs `run_pf.log` on P2A/P2B is zero differences in all four columns.

**The `wdrepmode` gate is verified, not assumed.** `grep -c "lie wholly inside a single"
run2_pf.log` returns **3** — one per top-level call (C5, C6, C7) across two 40-replication
bootstrap runs. Ungated it would have been 3 + 80.

## Cost

`cost.log`, 2,000,000 rows of which 805,000 treated, `set processors 2`:

| pass | s |
|---|---|
| the block's existing `gegen total() if gtin, by(i t)` (it runs **six** of these) | 0.141 |
| NEW `gegen nunique(clus) if gtin, by(i t)` | 0.110 |
| NEW `gegen tag(i t) if gtin & ncl==1` + `count` | 0.128 |
| NEW robust branch: `tag(i t) if gtin & nobs<=1` + `count` | 0.128 |
| the REJECTED `jtcoarse` route (two whole-sample tags + counts) | 0.185 |

**0.238 s per 2,000,000 rows** for the cluster branch, 0.128 s for the robust branch, **nothing at
all under an iid VCE, under `fastnose`, or inside a replication.** Against ~0.85 s the gteffects
block already spends on its own six `gegen` passes, and against a `reghdfe` fit on the same data.
It is also cheaper than the nesting route it replaces, and it is only paid when `gteffects` is
asked for.

## Does it move a reported number?

**No.** One caution line, one new `e(gt_ncontained)` scalar. Verified bit-identical on every arm.

## Blast radius

The estimation-time output of any `gteffects` run whose cells are contained — which includes the
package's central use case, `i` = state with `cluster(state)`, where it will fire on **every** cell.
That is the point, and it is also the reason it may need a ruling: PS-64 is the precedent (the
author chose "compute it with a named caution" over refusal on exactly this geometry, for
`jointtests`), but a caution that fires on the commonest specification is a UX decision.
If the author prefers, the same predicate supports two weaker variants at identical cost:
post `e(gt_ncontained)` silently and put the warning in the help only; or print only when the
contained cells are a strict subset of the table (i.e. suppress the all-cells case, where the user
has effectively chosen a cohort-level design).

## Alternatives rejected

- **The `jtcoarse` nesting screen at `ado:1414`, widened to `gteffects`.** Rejected on measurement:
  A1 arm 2 (`robust`, no `cluster()`) undercovers at 70% and the nesting test cannot see it. It is
  also more expensive (0.185 s of whole-sample passes against 0.128 s of treated-rows-only passes
  for the case it does cover) and it answers the wrong question — nesting of `i` in the cluster
  variable, rather than containment of a cell in one independent unit.
- **Scrubbing the contained cells' `se`/`p`/interval**, the way AMESE R1 scrubs the ame columns on
  the identical predicate. Consistent, and I nearly proposed it — rejected because on `i` = state
  with `cluster(state)` it blanks *every* cell of the table, the estimate is only understated (by
  1/0.29 and 1/0.64), not meaningless, and the author has ruled twice against removing numbers
  where a caution will do (PS-64; the T4 refusal of 2026-09-15). It stays the right variant B.
- **Correcting the standard error.** `sthlp:255` already rules it out and the arithmetic agrees:
  the missing variance component is the treated cell's own idiosyncratic variance, which the
  saturated dummy sets to zero residual and nothing in the sandwich can recover.
- **Widening the predicate to `ncl <= 2` or `<= 3`.** Not measured, so not proposed. A1 gives
  `ncl = 1` catastrophic (0.29, 0.64) and `ncl = 8` fine (0.92); everything between is unmeasured.

## The test that proves it

An invariant with no oracle: **the containment count must equal the number of cells whose
distinct-cluster count is 1**, and that is computable independently of wooldid from the raw data —
`gegen nunique(clus), by(i t)` on treated rows, compared with `e(gt_ncontained)`. Run on P1A/P2A/P2B
implicitly (39/109/109 all-cell cases); the discriminating case is a design where only *some* cells
are contained, which I did not build. **Untested on a partial-containment design.**

---

# GTE-5 — the help edits
*covers:* the documented-vs-actual halves of both findings. No code, no cost, no number moves.

### 5a. `sthlp:1148` — the sentence SKEP-NEW-2 quotes

**Now:**
```
{it:se} column is missing (degenerate zero variances are always reported as missing, never zero);
```
**Becomes:**
```
{it:se} column is missing. A degenerate cell variance is always reported as missing, never as a standard error, and the test is a relative floor rather than an exact zero: a robust or cluster-robust stage-1 covariance that is rank-deficient in a cell's direction leaves cancellation residue of order 1e-30 rather than a zero, so a cell whose standard error falls below 1e-10 of its own estimate, or whose variance falls below 1e-12 of the largest cell variance in the fit, is reported with a missing {it:se}, {it:p} and interval. The count is in {it:e(gt_ndegen)} and a caution names it.
```

### 5b. `sthlp:255` — "two things follow" is three, and the third is not about clustering

**Now (substring):**
```
and two things follow that the table does not announce: the pointwise standard errors on the {it:att} rows are understated (measured at roughly 0.66 to 0.75 of the reference on the diagnostic fixtures), and the joint pre-trend tests are oversized, which is why those carry a caution of their own (see {cmd:jointtests}).
```
**Becomes:**
```
and three things follow: the pointwise standard errors on the {it:att} rows are understated (measured at roughly 0.66 to 0.75 of the reference on the diagnostic fixtures) and the table does not announce it; the joint pre-trend tests are oversized, which is why those carry a caution of their own (see {cmd:jointtests}); and the {cmd:gteffects} cells fail far harder than the {it:att} rows do, because a cell has no aggregation to average the error over. Measured on a 30-state panel with {it:i} = the unit and {cmd:cluster()} on it, a cell standard error was 0.29 of the cell estimate's true sampling standard deviation and the cell's own 95% interval covered its true effect 44% of the time, against 0.92 and 94% for the same cell with {it:i} at the cohort and {cmd:cluster()} finer than it. That third one is not confined to clustering: {cmd:robust} with a unit-level {it:i} leaves each treated cell holding a single observation whose residual its own dummy drives to zero, and the same panel then gives 0.64 and 70%. A run that asks for {cmd:gteffects} counts the affected cells, prints a caution and returns the count in {it:e(gt_ncontained)}.
```

### 5c. `sthlp:1146` — the gteffects "Notes:" list, cross-reference

**Now (substring):**
```
standard errors are the regression's coefficient-level (conditional) VCE — {it:unconditionalse} and {it:impvalues} do not alter them;
```
**Becomes:**
```
standard errors are the regression's coefficient-level (conditional) VCE — {it:unconditionalse} and {it:impvalues} do not alter them, and where a treated cell lies wholly inside one cluster (or, under {cmd:robust}, inside one observation) they are understated and the cell interval does not cover at its nominal level, which {help wooldid##standarderrors:What level to cluster at} measures and {it:e(gt_ncontained)} counts;
```

Register note: the surrounding prose runs long sentences with em-dash asides and parenthetical
measurements, and `sthlp:255` already carries a measured range in exactly this shape, so 5b is
written to that register rather than to a shorter one. It should still be measured against
`2a2da17` on caps, " - " per 1000 words and asides before it ships.

---

# What should NOT be fixed

- **The `att` row's residue standard error at K=1.** GTE-3 makes the run say the covariance is
  worthless; it does not scrub the number. Scrubbing every inference column of a whole run is a
  bigger ruling than this evidence supports, and the two-cluster design that produces it is thin
  enough that a loud caution is the proportionate answer.
- **The plain-`bootstrap()` cell standard errors on the same design.** `run2_pf.log` arm C6:
  `bootstrap(standard) reps(40) cluster(cohid)` without `boottstat` still returns 6 cells with
  `se = 3.9e-15` — but that number is the **bootstrap standard deviation** of the cell across
  replications, not the analytic residue, and a cell that genuinely does not move across resamples
  is a real (degenerate) bootstrap result, not a floating-point artefact. Flooring it would need
  its own measurement of what a legitimate near-zero bootstrap sd looks like. **Reported, not
  proposed.**
- **`ncl <= 2` or `<= 3` in GTE-4.** Unmeasured; do not widen a predicate on intuition.

---

# Re-bless scope

The author can find every battery configuration this touches with one predicate:

> a config specifying `gteffects` **and** (`cluster()` at or above `i`, **or** `robust` with a
> unit-level `i`, **or** fewer than about three clusters, **or** `boottstat`).

Configs outside it are bit-identical — eight arms and 537 `e(gteffects)` rows verified, plus
`e(wdidmainresults)` on C1. I did **not** run the project battery (the brief puts `test/suite/**`
off limits) and did not count how many configs match; the predicate is offered so the author can.

New `e()` members: `e(gt_ndegen)` and `e(gt_ncontained)`, both `cap ereturn scalar`, both absent
when zero, both existence-gated at every read.

---

# Untested / open

- **Multi-way `cluster()`**: GTE-4 loops the cluster variables and keeps the largest count. Not
  measured, and the right aggregation for a CGM multiway meat is arguable — one degenerate
  dimension does not necessarily kill the sum.
- **`customvce(cluster ...)`**: parsed by the same three lines R1 and the `jtcoarse` screen use, so
  it should behave identically. Not measured.
- **A partial-containment design** (some cells contained, some not) — the discriminating test for
  GTE-4's count. Not built.
- **`permute()`**: `e(gteffects)` inference columns are blanked outright there (`ado:8121`), so
  GTE-1 and GTE-2 are inert. Not measured.
- **`fastnose`**: the help says the regression's own VCE is intact under it, so GTE-1's floor
  applies and GTE-4 is skipped by its `fastnose` gate. My control arm for it (`body.do` P4) errored
  r(459) on the panel I built, so this is asserted from the code, not measured.
- **The scale-metamorphic regression test for GTE-1** (`y` vs `1000*y` must scrub the same rows) —
  proposed, not run.
