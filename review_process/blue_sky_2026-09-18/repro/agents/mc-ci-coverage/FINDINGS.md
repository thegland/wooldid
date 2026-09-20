# FINDINGS — mc-ci-coverage
## target: coverage of every confidence-interval variant wooldid can report
## method: Monte Carlo (interval coverage) + interval-shape sanity

Build: the pinned `.../bluesky/build`, loaded with `adopath ++` in every do-file (`which wooldid`
is echoed in every log). No file under `test/` was read or run; every dataset is generated inside
the probe. Pre-registered decision rule, fixed before any coverage number was looked at:
**flag when |coverage - nominal| > 3*sqrt(p(1-p)/R)**.

**Harness power check (`mc_cheap.do` arm A, R=250).** Analytic interval with `cluster()` AT `i`
(i=cg, cluster(cg), G=4): coverage **0.544** at nominal 0.95, reported se / true sd = **0.280**.
So the harness detects undercoverage when it is there; a clean number below is a real negative.
(This is the already-ruled understatement of known_bydesign s1, measured at G=4. Control only.)

---

# MC1 — CRITICAL — the `impvalues` (Ibragimov-Muller) interval covers 5% at nominal 95%

**surface:** `impvalues` — every reported `se`, `t`, `p`, `lb_*ci`, `ub_*ci`
**severity:** critical **confidence:** certain **realistic:** yes
**seeded:** partly — the brief listed `impvalues` among the interval machines to measure at
`i = cohort`. The mechanism diagnosis, the geometry rule and the permute/bootstrap-vs-impvalues
guard contrast below are mine.

## Claim
Under `impvalues` the G "per-cluster estimates" whose dispersion becomes the standard error are
not G separate fits. They are `margins ... over(cluster)` evaluated on the ONE pooled model
(`wooldid.ado:2042` the margins call, `2047-2048` the column grab, `2054` the `ttest` whose
`r(se)` becomes the reported standard error): deterministic re-weightings
`theta_g = sum_c w_{g,c} * beta_c` of a single estimated coefficient vector. The identical
construction supplies every `impvalues` row — the pretreatment placebo (`2100-2110`), each
event-study horizon (`2219-2229`) and the continuous-treatment ame rows (`2320`, `2375`,
`2715`) — so the whole `impvalues` inference surface is affected, not just `att`. Wherever **several clusters share the same levels of `i`** the G reads take only as many distinct
values as the clusters have distinct weight profiles over the `i`-by-relative-period cells —
three, one per treated cohort, when `cluster(unit)` sits inside `i = cg`; exactly one when every
cluster has the same cohort composition. Their spread is then set by **cross-cohort effect
heterogeneity**, not by sampling error, and the reported interval is far too narrow. IM is
allowed to pick up true cluster heterogeneity — what no version of IM permits is the sampling
component being **absent entirely**, which is what a deterministic re-weighting of one fit gives.

## Repro
`repro_MC1_impvalues_ci.do` (self-contained, ~2 min, two blocks).

## Observed

**(a) The mechanism, exactly.** `probeD.do` / `probeD.log`: 100 units, `i = cg` (3 treated
cohorts + never-treated), `cluster(reg25)` = 25 regions that cut across the cohorts, each holding
exactly one unit of each cohort. Idiosyncratic sd 0.60 plus a unit shock of sd 0.50 — `sd(y)` is
about 1.4. wooldid blanks **every** se/t/p/lb/ub for all 6 reported rows and prints

> Caution: for one or more Ibragimov-Muller estimates the per-cluster estimates were
> **numerically identical**, so the cross-cluster (IM) standard error is not identified ...
> **A collapsed IM standard error usually reflects deterministic or noiseless data - please
> verify your input if this was not expected.**

Twenty-five independent fits on noisy data cannot agree to the last bit. Twenty-five
re-weightings of ONE coefficient vector by IDENTICAL weights must — which is exactly what these
clusters have, since each holds one unit of each cohort. The analytic clustered run on
the same data returns an ordinary se of 0.0497. The caution sends the user to look for a data
problem that does not exist.

**(b) The consequence, at the pairing the help RECOMMENDS.** `mc_imp2.do` / `mc_imp2.txt`:
`i = cg`, `cluster(unit)` (the cluster is FINER than `i` — "an `i` set at the treatment-onset
level can therefore sit above the randomization unit, which is the pairing ... recommends"),
40 units, effects 1.00/1.05/0.95, idio sd 0.60 + unit sd 0.40:

R = 200:

| arm | mean reported se | true sd of the point estimate | 95% coverage |
|---|---|---|---|
| `impvalues` | 0.0211 | 0.1239 | **0.290** (3 SE = 0.046 -> FLAGGED) |
| analytic `cluster(unit)` on the same draws | 0.1377 | 0.1252 | 0.955 |

The reported IM standard error is **5.9x too small**.

**(c) The consequence, at a cross-cutting cluster.** `mc_imp.do` / `mc_imp.txt`, R=250,
100 units, 25 regions cutting across 4 cohorts, effects 1.00/1.02/0.98:

| arm | mean reported se | true sd | mean(b) | target | 95% coverage |
|---|---|---|---|---|---|
| `impvalues` | 0.00257 | 0.0857 | 0.9983 | 1.00142 | **0.048** |
| analytic `cluster(reg)` | 0.0889 | 0.0859 | 0.9984 | 1.00154 | 0.968 |

The point estimate is unbiased for its own target (0.9983 against 1.00142, 0.36 MC standard
errors) and has sampling sd 0.0857, so ANY interval claiming 95% coverage of ANY fixed target
must be at least about 0.34 wide. The reported interval is 0.011 wide.

**(d) At realistic noise the IM se is far below the true sd, and it is not a noise quantity.**
`probeF.do` / `probeF.log` sweeps the idiosyncratic sd 0.60 -> 0.06 -> 0.006. The analytic
clustered se tracks it exactly (0.0497 -> 0.00497 -> 0.000497); the IM se barely moves
(0.1205 -> 0.1028 -> 0.1011, and 0.0398 -> 0.0432 -> 0.0436 in the unbalanced arm). Textbook IM
with heterogeneous cluster parameters also has a floor, so this is corroboration rather than
proof — the proof is (a). What it does establish is that at realistic noise the IM se is set by
something other than the sampling error of the estimate.

## Expected
An interval reported at `level(95)` covers its target about 95% of the time, and its standard
error is within sampling error of the estimator's sd. The help sells `impvalues` as the valid
alternative where the analytic path is not ("Valid ame inference under such designs: ...
`impvalues` - cross-cluster dispersion inference, which remains valid").

## Scope (measured, `mc_cheap.do` arm C, R=250)
`impvalues` is approximately right when the cluster **nests** `i`: `i` = unit, `cluster(reg)`
with each region holding 3 whole units gives coverage **0.916** and se/sd **0.783** — not flagged
by the 3-SE rule, though still on the narrow side. The failure is confined to a cluster that is
**finer than `i` or crossed with it**.

## Why this reads as an oversight rather than a ruled trade-off
`permute()` REFUSES a cluster that does not nest `i`, with a named r(459)
("permute() requires the cluster() variable to nest i (every i wholly inside one cluster)"),
and `bootstrap()` carries an explicit two-regime rule for the same geometry. `impvalues` has
**no guard at all** — no error, no caution, no help sentence — and silently reports a standard
error an order of magnitude too small.

## Not the already-ruled T4 item
known_bydesign s1 rules that `impvalues` reports the cluster-equally-weighted mean and that
**non-`att`** estimands print under **cross-cutting** clusters. Every run above uses `att` alone,
the headline geometry (b) is a cluster **nested inside** `i` rather than crossed with it, and the
defect is in the standard error, not in which estimand is reported.

---

# MC2 — MEDIUM — the analytic sup-t simultaneous band uses a NORMAL critical value while the table beside it uses t(G-1), so the band under-covers, badly at few clusters

**surface:** `simultaneousbands` (analytic engine, `e(esband_engine) == analytic-simulated`)
**severity:** medium **realistic:** yes
**confidence — split deliberately:**
* that the band under-covers — **certain**. 0.752 at nominal 0.90 and 0.843 at nominal 0.95,
  R=400, about 9 binomial SE in both cases, with the pointwise intervals on the same runs
  correctly calibrated.
* that the cause is the z-versus-t reference — **probable**. The z->t rescale buys about 12%
  extra width; whether that closes the whole 10.7-point deficit, or whether the other half comes
  from taking the sup over an ESTIMATED correlation matrix as if it were known, is not
  established here.
**seeded:** no — the seed asked whether a band can be NARROWER than the pointwise interval (it
cannot; see the clean result below). Band COVERAGE is mine.

## Claim
`wooldid_esbanda` (`wooldid.ado:16695-16746`) forms the simultaneous critical value as the
`level()` quantile of `max_j |Z_j|` over 10,000 draws of `Z = L * rnormal(m, nd, 0, 1)` — a
**standard-normal** sup quantile taken over the ESTIMATED correlation matrix treated as known. The pointwise intervals printed on the very same table use
`invttail(G-1)` (verified exactly in `probeA.log`: implied critical value 1.980099876 at
level(95), G=40, = `invttail(39,.025)` to 9 digits). The band is therefore too narrow by roughly
`t(G-1)/z`, and its measured simultaneous coverage falls below `level()` — by 2 points at G=40
and by 8-10 points at G=12.

## Repro
`repro_MC2_band_coverage.do` (self-contained, ~4 min). Its own run at level(90), G=12, R=200:

```
  SIMULTANEOUS band covers all 6 horizons :  0.780   (nominal 0.90)
  mean band critical value                :   2.20272
  invttail(11,(100-90)/200)               :   1.79588
  invnormal(1-(100-90)/200)               :   1.64485
```
The band's sup multiplier over 6 correlated horizons is 2.20272/1.64485 = 1.3392 relative to a
z reference. Applied to the t(11) reference the same multiplier gives 2.4047 — a 9% wider band.
Also `mc_band10.do` (G=12, R=316) and `mc_band.do` (G=40, R=400).

## Observed
`mc_band10.txt`, i=cg, `cluster(unit)` G=12 units, 10 periods, 6 horizons, R=400 (final):

| nominal | band simultaneous coverage | 3 binomial SE | flag | mean critical value |
|---|---|---|---|---|
| 0.90 | **0.752** | 0.045 | **FLAG** (14 points of deficit, 9 SE) | 2.2031 |
| 0.95 | **0.843** | 0.033 | **FLAG** (10 points of deficit, 9 SE) | 2.4886 |

On the very same runs the analytic POINTWISE intervals are correctly calibrated per horizon
(0.877-0.951 at nominal 0.90; 0.944-0.975 at nominal 0.95), so the covariance and the pointwise
machinery are fine — it is the band's critical value.

At level(95), G=12: `invttail(11,.025) = 2.2010`, so the band's 2.4884 is only 1.13x the
pointwise critical value. A sup-t critical value built on the same t(11) reference would be about
`2.4884 * 2.2010/1.9600 = 2.794`, i.e. **12% wider**.

`mc_band.txt`, the same measurement at G=40, R=400: 0.875 at nominal 0.90 and 0.930 at nominal
0.95 — 1.7 and 1.8 binomial SE low, below nominal at both levels in the same direction, and
`t(39)/z = 1.010` predicts exactly a deficit that small. The size of the miss tracking `t(G-1)/z`
across two cluster counts is what makes the mechanism more than a guess.

## Expected
The help: "a simultaneous band covers **every** horizon of the family at once with that
probability" at `level()`. It does not.

## Realistic
Yes, and the failure is worst exactly where simultaneous bands are most wanted — a
state-policy panel with a dozen clusters.


---

# MC3 — LOW — requesting event-study horizons silently moves the headline `att` and inflates its standard error by up to 60%

**surface:** `esprelength()` / `espostlength()` interaction with the reported `att` interval
**severity:** low **confidence:** certain **realistic:** yes **seeded:** no

**Observed** (`probeI.do` / `probeI.log`, ONE dataset, 40 units, i=cg, `cluster(unit)`).
The load-bearing observation a skeptic should check first: **`e(N)=320` and `e(df_r)=39` are
identical in all four runs**, so no row and no degree of freedom changed — only the reference pool
did.

| options | `att` | `se` |
|---|---|---|
| (none) | 1.093400346 | 0.132401198 |
| `esprelength(1) espostlength(1)` | 1.117099722 | 0.158879296 |
| `esprelength(2) espostlength(2)` | 1.132362929 | 0.173433656 |
| `esprelength(3) espostlength(3)` | 1.065211265 | 0.211200984 |

**Mechanism, and it IS documented — but in the wrong place.** The `e(escomp_refN)` stored-result
entry says the reference-normalization pool is "the union of treated rows earlier than the
requested pre-window (`relyr < -esprelength`) and each treated unit's forced first in-sample
observation". So asking for pre-horizons removes those periods from the reference pool. The
`esprelength()` option entry itself says nothing about it, and nothing is printed at run time.

**Why it still matters.** A user who adds `esprelength(3)` only to draw a plot gets a different
ATT and a 60%-wider confidence interval, with no message. All four estimates remain consistent
for the same truth (1.0 in this DGP), so this is efficiency, not bias — hence low severity — but
it is exactly the kind of thing that makes two runs of "the same" specification disagree.

---

# PROBED HARD AND CLEAN (negative results)

Every one of these would have shown up if broken; the harness's power is established by the
`cluster(cg)` control (coverage 0.544 at nominal 0.95).

1. **Analytic intervals honour `level()` exactly.** `probeA.log`: at 90 / 95 / 99 / 87.5, on all
   20 rows (att family, pretreatment, ES horizons), the implied critical value is
   `invttail(G-1,(100-L)/200)` to nine digits (1.657759285 / 1.980099876 / 2.617775976 /
   1.545006301); `e(level)` is right; columns are named `lb_87p5ci`. Zero lb>=ub, zero
   b-outside-CI, zero half intervals.
2. **Analytic interval coverage at a sound geometry.** `mc_analytic.txt`, R=500, i=cg,
   `cluster(unit)` G=40: coverage 0.946-0.956 at nominal 95 and 0.880-0.904 at nominal 90 on all
   16 parameters; mean(se) = sd(b) to three digits. Nothing exceeds 3 binomial SE.
3. **All four bootstrap interval families reconstruct exactly.** `mcboot_*.txt`: an independent
   re-implementation of normal / percentile / bias-corrected / BCa reproduced
   `e(boot_ci_normal|percentile|bc|bca)` with a maximum absolute difference of **0.0** on every
   draw (the `NATIVEDIFF` column), using `invt(G-1)` for the normal family. So `level()` reaches
   all four bootstrap families correctly.
4. **Bootstrap interval nesting.** 0 violations of `CI(80) subset CI(90) subset CI(95) subset
   CI(99)` in 2582 ordered pairs, BCa included (`nest_boot.py`).
5. **Point estimate inside its own interval.** 0 violations in 3448 (family x level x draw)
   checks across normal / percentile / bc / bca. The seeded worry that bc/bca can exclude the
   estimate did not fire on realistic data.
6. **A simultaneous band is never NARROWER than the pointwise interval beside it.** The seeded
   suspicion at few clusters: `probeE2.log`, 12 clusters with one asymmetric contaminated unit,
   `boottstat` + `simultaneousbands`, levels 99 and 95, 30 horizon comparisons — the band was
   wider in every one (ratio 1.78 to 8.02). Under the analytic engine (`probeB.log`) the band's
   critical value exceeded the pointwise one at 90/95/99 as well. The help's "under `boottstat`
   the band can never come out narrower" survives. (Its coverage does not — see MC2.)
7. **`permute()` interval shape.** `probeH.log`: at `cluster(cg)` (G=4) and at
   `cluster(unit) + permonsetcohorts` (G=40, `e(perm_irebuilt)=1`, the 2026-09-16 feature), 12
   rows each at levels 90/95/99/87.5 — lb<ub everywhere, estimate inside the interval everywhere,
   no half intervals, CI(90) strictly inside CI(99). `probeI.log` maps the objects:
   `e(perm_lb_*)`/`e(perm_ub_*)` are the permutation-distribution QUANTILES, the printed interval
   is `e(perm_ci)` = `[b - q_hi*se, b - q_lo*se]`, and `e(perm_tlb/tub_*)` is `permtci`
   (symmetric about `b`).
8. **BCa half-interval guard.** The `1 - a*w <= 0` guard in `wooldid_bootbca` can in principle set
   `ci[1,k]` and then `continue` before `ci[2,k]`, leaving an interval with a lower bound and no
   upper bound. It never fired: 0 half intervals in ~3500 BCa columns, including 12-cluster runs
   with a contaminated unit. Not reachable on realistic data.
9. **`unconditionalse`.** `mc_cheap.txt` arm B, R=250: coverage 0.980 at nominal 0.95, se/sd
   1.230 — NOT flagged, and over-covering. `vce(unconditional)` prices in sampling of the
   covariate distribution that my fixed-design DGP does not have, so this is the expected
   direction and not a defect.
10. **`impvalues` when the cluster NESTS `i`.** `mc_cheap.txt` arm C, R=250: coverage 0.916,
    se/sd 0.783 — not flagged. Reported as the scope boundary of MC1.
11. **The two percentile-t endpoint constructions are right, by code read.**
    `wooldid_bootcellwire`'s studentized branch sets `M[i,lbcol] = b - hi*so` and
    `M[i,ubcol] = b - lo*so` (`ado:15092-15093`) and `wooldid_boottp` takes equal-tailed SIGNED
    quantiles of `t_r = (b_r - b_obs)/se_r` and the caller places `[b - q_hi*se, b - q_lo*se]`
    (`ado:14332-14405`). Both have the reversal that a percentile-t needs; neither is the
    `b +/- q*se` mistake. `wooldid_bootbca`'s acceleration is Efron (1987)
    `a = sum(d^3)/(6*[sum(d^2)]^1.5)` over a delete-one-CLUSTER jackknife, as documented.
12. **Shape scans on every other interval surface.** `e(wdidsgresults_0/_1)`, `e(pwsgs)`,
    `e(gteffects)`, `poisson`, `poisson + poisexpresults`, `semielasticity`, at level 95 and
    across 90/99 nesting — zero lb>=ub, zero b-outside-CI, zero half intervals
    (`probeB.log`, `probeD.log`). Repeated UNDER `bootstrap(standard)` with `subgroups(sg)` and
    `gteffects` (`probeJ.log`, levels 90 and 95): `e(wdidmainresults)`, `e(wdidsgresults_0/_1)`,
    `e(pwsgs)`, `e(gteffects)` (18 cells) and all 36 columns of
    `e(boot_ci_normal|percentile|bc)` — zero lb>=ub, zero half intervals, zero
    point-estimate-outside-interval, and the level label follows into `e(gteffects)`
    (`lb_90ci` / `lb_95ci`).

13. **`impvalues` at a cross-cutting cluster is not silently wrong about its LEVEL**: columns are
    named `lb_87p5ci` at `level(87.5)` and `e(level)` is right, `e(df_r)` is the IM `G-1` = 24
    at G=25 (`probeD.log`). The level plumbing is fine; only the standard error is not.

---

# UNDERPOWERED / STILL RUNNING AT WRITE-UP (honest status)

* **Bootstrap coverage** (`mcboot_1.txt`, `mcboot_2.txt`): 36 Monte Carlo draws at
  reps(199), i=cg, `cluster(unit)` G=40, `att` + `att_i` — this is the §7 "finer-than-`i`
  bootstrap regime" arm (units kept whole, drawn copies arrive as duplicate rows), whose
  assessment suite "has not run". 3 binomial SE at R=36 is 0.11, so coverage there can only
  exclude gross failure. What it does show: coverage at nominal 95 is at or above nominal for
  every family and every parameter (no family ever under-covers), and the bootstrap draw
  standard deviation runs about 9% above the analytic standard error that `mc_analytic.txt`
  showed to be correctly calibrated at R=500 (mean bootstrap sd 0.1997 for `main_att` against
  0.1830). Direction: mildly conservative, consistent with the §7 note that
  `bootstrap(standard)` over-covers. **Not a finding, and not a clean bill of health either.**
* **`boottstat` percentile-t coverage** (`mctstat_1.txt`) and **`permute()` interval coverage**
  (`mcperm_1.txt`, i=cg + `cluster(unit)` + `permonsetcohorts`) reached only single-digit R
  under machine contention. Both are shape-clean so far (0 lb>=ub, 0 point-estimate-outside-CI)
  but carry no coverage verdict. Their do-files are in the workdir and simply need more draws.
