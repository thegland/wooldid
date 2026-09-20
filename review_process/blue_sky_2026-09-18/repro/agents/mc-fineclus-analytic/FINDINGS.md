# mc-fineclus-analytic — FINDINGS

Agent method: MONTE CARLO. Target: does the documented remedy for the analytic-SE
undercoverage ("choose `i` COARSE, cluster FINER"; `build/wooldid.sthlp:251`) actually work?

## FLAG RULE, fixed before any result was inspected (see NOTES.md)
An arm FAILS if, in draws where BOTH external controls (csdid, did_imputation) sit inside
the bands: median over horizons of mean(reported SE)/empirical SD < 0.90, OR 95% coverage
< 0.93. Same two thresholds applied separately to the headline `att` row.
At ~290 draws a coverage estimate carries a +/- 1.3pp standard error.

---

## HEADLINE (this is the answer to the assigned question)

**THE REMEDY WORKS on the analytic path — at 10 clusters per treated cohort-by-period cell.**
Below that it degrades, and the degradation is measured and mechanised in FINDING MC-FC-1.
DGP-S, 283 draws, 40 states x 5 units x 12
periods, onset at the state level, dependence at the state level (state RE + state-by-time
shock), cluster(state) the correct level. Every arm on the SAME draws, and `b` is EXACTLY
identical across every arm (max mreldif over 283 draws = 0).

| arm | geometry | med SE/SD | med cov95 | att SE/SD | att cov95 | verdict |
|---|---|---|---|---|---|---|
| S2  | i=state, cluster(state) — **AT i** | 0.559 | 0.742 | 0.610 | 0.763 | **FAILS** |
| S2p | same, pooled `esprelength(0)`      | 0.578 | 0.739 | 0.707 | 0.841 | **FAILS** |
| S3  | i=coh, cluster(state) — **FINER**  | 0.933 | 0.943 | 0.996 | 0.947 | **PASSES** |
| S3p | same, pooled `esprelength(0)`      | 1.017 | 0.958 | 1.067 | 0.968 | **PASSES** |
| S4  | i=coh, cluster(unit) — below the DGP's dependence | 0.532 | 0.707 | 0.574 | 0.721 | expected (see below) |
| CS  | csdid long2 (control)              | 0.982 | 0.936 | — | — | in band |
| DI  | did_imputation (control)           | 0.987 | 0.943 | — | — | in band |

Both controls are inside the bands, so the gate APPLIES here (it does not abstain as it
did in the project's own 2026-09-14 run).

S4 was logged **in advance** as expected-to-undercover: cluster(unit) sits BELOW the level
at which this DGP's dependence lives, so its failure is the DGP's, not wooldid's. It is
reported because it is a useful positive result in its own right — wooldid tracks whatever
cluster level it is handed rather than silently substituting one.

Repro: `mc_fineclus.do` (driver), `analyze.py` (analysis), full table `RESULTS_DGP_S.txt`.
Self-contained single-draw demonstration: `repro_remedy_works.do`.

## HEAD-TO-HEAD AT THE GEOMETRY WHERE THE DEFECT WAS DOCUMENTED (DGP-P, 142 draws)

40 INDEPENDENT units, T = 12, units 1-10 adopt at t=5, 11-20 at t=8, 21-40 never; unit RE
N(0,1) plus AR(1) rho .5 idiosyncratic errors. Units are independent, so cluster(unit) is
exactly the right level and any undercoverage is the estimator's. All arms have 40 clusters
and identical points (max |b(i=unit) - b(i=coh)| = 4.3e-14 over 141 draws).

| arm | geometry | med SE/SD | med cov95 | att SE/SD | att cov95 |
|---|---|---|---|---|---|
| P1  | i=unit, cluster(unit) — **cluster == i** | 0.739 | 0.852 | 0.816 | 0.866 |
| P1p | same, pooled `esprelength(0)`            | 0.644 | 0.796 | 0.734 | 0.859 |
| P2  | i=coh,  cluster(unit) — **FINER**        | 1.002 | 0.951 | 1.024 | 0.965 |
| P2p | same, pooled `esprelength(0)`            | 0.971 | 0.957 | 0.971 | 0.965 |
| P3  | i=icoh, cluster(icoh) — **cluster == i** | 0.556 | 0.738 | 0.618 | 0.738 |
| DI  | did_imputation, cluster(unit) (control)  | 0.962 | 0.929 | — | — |

**The positive control reproduces the seeded result.** P1 covers at 0.852 and P1p at 0.796
against the project's own 500-draw figures of 0.838 (esfixedbaseperiod, i=unit) and 0.794
(pooled, i=unit). P1p's SE/SD of 0.644 sits on top of the project's 0.661. The harness is
measuring the same thing the project measured.

**On the same draws, the remedy restores coverage to nominal.** P2 0.951 / P2p 0.957.

**Honest caveat on the gate in this wave.** csdid returns r(198) on all 141 draws when handed
`cluster(unit)` with `ivar(unit)` — the message was captured, so the reason is inferred
rather than read, but it is presumably that cluster() coincides with ivar(), so only did_imputation survives as a control, at a median
covT of 0.929 — 0.1pp under the 0.93 gate, i.e. about 1.1 standard errors of MC noise from
0.95. On a strict reading of "both controls inside band" the gate abstains in this wave.
The magnitudes are not close (0.852 / 0.796 against 0.951 / 0.957), and in the DGP-S wave
both controls ARE inside the band and the gate fires cleanly. Abstention is not exoneration.

### A sharpening of the project's own 2026-09-14 record (context, not a package defect)
`review_process/se_coverage_2026-09-14/se_cov_mc.do:56` declares the unit id as the variable
named `i`, and lines 111 and 145 then run `wooldid y icoh t ttre, att cluster(i)`. So the arm
that README calls "i = the cohort id" was ALREADY clustering at the UNIT, i.e. STRICTLY FINER
THAN `i` — it was the remedy arm all along, and its 0.934-0.953 coverage is the remedy
working. My P2 reproduces it at 0.951 on independent code. The README's conclusion, "The
defect is specific to i = the unit id", therefore reads as a claim about the CHOICE OF i when
the operative variable is the RELATION between cluster() and i. My P3 separates them: i is
coarse there (the onset cohort for treated units) but cluster() sits AT i, and coverage is
0.738. Choosing a coarse i does not fix anything by itself; clustering finer than i does.

## DETERMINISTIC CHECKS AROUND THE MACHINERY — ALL CLEAN (`probe_determ.do` / `.log`)
On one DGP-S draw, with `i` = the onset cohort and `cluster()` strictly finer:
* `nofastcollapse` moves nothing: mreldif(b) = 2.6e-16, mreldif(V) = 6.9e-18. The default
  cell collapse is exact under a finer-than-`i` cluster (the collapse is a margins-side
  device and the clustered VCE is formed on the full estimation sample).
* The point estimate is **bit-identical** across `robust`, `cluster(state)` and
  `cluster(unit)`: att = +1.9416320a9e57aX-001 in all three. cluster() does not move b.
* The `jointtests` degeneracy caution fires in **exactly** the nested geometry
  (i=state, cluster(state)) and does NOT fire for cluster-finer-than-i or for a
  cross-cutting cluster. The predicate at `ado:1421-1432` (count distinct i vs distinct
  (i,cluster) pairs) is the right test and is correctly applied.
* `e(N_clust)`/`e(df_r)` are right (200/199 for cluster(unit)); two-way
  `cluster(state region)` runs and reports the min-G convention (8 clusters, df 7).
* A cluster variable with missing values drops the affected rows and says so in the
  header ("N = 2220 (180 obs dropped from initial sample)"), N_clust 37 of 40.

## CROSS-CUTTING cluster() — CLEAN (DGP-X, 200 draws, `RESULTS_DGP_X.txt`)
Same 40 states x 5 units x 12 periods and the same state-level dependence, but onset is
assigned INSIDE the state (unit 1 adopts at t=5, unit 2 at t=8, units 3-5 never), so every
state holds all three cohorts and `state` neither nests nor is nested in `coh`.

| arm | geometry | med SE/SD | med cov95 | att SE/SD | att cov95 |
|---|---|---|---|---|---|
| X3  | i=coh, cluster(state) — **CROSS-CUTTING** | 0.999 | 0.950 | 0.973 | 0.955 |
| X3p | same, pooled                              | 0.959 | 0.945 | 0.959 | 0.945 |
| X4  | i=coh, cluster(unit)                      | 1.465 | 1.000 | 1.365 | 0.995 |
| CS  | csdid (control)                           | 0.978 | 0.933 | — | — |
| DI  | did_imputation (control)                  | 0.966 | 0.940 | — | — |

The cross-cutting case is **not** silently wrong: X3 sits on nominal, on the same draws where
both controls do. `att_i` and `att_it` behave the same (0.99-1.01, coverage 0.95-0.965), so
the cell-weighted estimands are not specially broken by a cross-cutting cluster on the
analytic path either. X4's over-coverage is the DGP's (with onset varying inside the state,
the state shock cancels across treated and control units in the same state; cluster(unit)
cannot see that cancellation and adds variance the estimator does not have) — another
demonstration that wooldid reports the level it was handed rather than substituting one.

## RAGGED, UNEQUAL-SIZED ADMINISTRATIVE PANEL — CLEAN (DGP-U, 92 draws, `RESULTS_DGP_U.txt`)
State sizes unequal and fixed (M_g = 2 + mod(3g,9), so 2 to 10 units), panel unbalanced on a
fixed deterministic pattern (~6% of unit-periods dropped), onset at the state level.

| arm | geometry | med SE/SD | med cov95 | att SE/SD | att cov95 |
|---|---|---|---|---|---|
| U2  | i=state, cluster(state) — AT i            | 0.589 | 0.739 | 0.589 | 0.772 |
| U3  | i=coh,   cluster(state) — **FINER**       | 0.970 | 0.935 | 0.957 | 0.946 |
| U3p | same, pooled                              | 0.970 | 0.935 | 0.973 | 0.957 |
| U4  | i=coh,   cluster(unit)                    | 0.513 | 0.685 | 0.513 | 0.685 |
| CS  | csdid (control)                           | 0.975 | 0.934 | — | — |
| DI  | did_imputation (control)                  | 0.949 | 0.923 | — | — |

The remedy survives unequal cluster sizes and an unbalanced panel. `att_i` (0.982 / 0.935)
and `att_it` (0.954 / 0.946) hold up too. The ragged design moves the estimand itself (the
true att is 0.9060 rather than 0.9154), and the zero-noise truth tracks that correctly:
mean(b) - truth is +0.005, 0.04 of a sampling SD, and csdid and did_imputation carry the
same small bias, so it is the estimand's composition, not a defect.

---

# FINDING MC-FC-1 — the understatement is a GRADIENT in clusters-per-treated-cell, not a binary in the nesting relation; the finer-than-`i` remedy is incomplete when a treated cohort holds few clusters

**surface**: `cluster()` finer than `i` on the default analytic path (`sthlp:251`).
**severity**: medium-high  **confidence**: probable  **realistic**: yes  **seeded**: no

**Claim.** `sthlp:251` presents the analytic-SE problem as a property of the NESTING
RELATION — "At or above the level of `i` ... the pointwise standard errors on the att rows
are understated ... the remedies are a finer `cluster()`". Measured, the shortfall is
CONTINUOUS in K, the number of clusters inside a treated cohort-by-period cell, and a finer
`cluster()` only moves K off 1. With `cluster()` strictly finer than `i` throughout — the
documented remedy, unchanged — wooldid's mean reported SE over the empirical sampling SD is
0.912 at K = 10, 0.884 at K = 5 and 0.845 at K = 3, against a nominal 1.0. Coverage of the
95% interval falls with it: ~0.93 at K = 10, ~0.89 at K = 5, ~0.87 at K = 3.

**The mechanism, and why this is a documentation gap rather than a coding defect.** The
measured ratios track `sqrt((K-1)/K)` — the residual downward bias of a CR1 cluster-robust
variance for an effect identified from K clusters — across seven independently generated
designs:

| design | K | sqrt((K-1)/K) | measured median SE/SD |
|---|---|---|---|
| ladder K=3, 150 draws        |  3 | 0.816 | 0.845 |
| DGP-F (5 waves x 3), 200 dr  |  3 | 0.816 | 0.865 |
| ladder K=5, 150 draws        |  5 | 0.894 | 0.884 |
| ladder K=10, 150 draws       | 10 | 0.949 | 0.912 |
| DGP-S (2 waves x 10), 283 dr | 10 | 0.949 | 0.933 |
| DGP-U ragged, 92 draws       | 10 | 0.949 | 0.970 |
| DGP-X cross-cutting, 200 dr  | 40 | 0.987 | 0.999 |

reghdfe applies a global G/(G-1) correction with G = 40 clusters; nothing corrects the K
clusters that actually identify a given treated cohort-by-period cell. `cluster() == i` is
the K = 1 limit of the same formula, where the prediction is 0 and the measurement is the
0.18-0.71 collapse the help already documents. So this is an INHERITED property of CR1
cluster-robust inference, not a wooldid coding error — but the help presents the remedy as
a fix for a binary condition when what it actually does is move K from 1 to whatever the
design supplies, and it gives the user no way to see what K they have.

**The ladder.** 40 states x 4 units x 12 periods, two adoption waves at t=5 and t=8 with K
states each and 40-2K never-adopting states; state RE + state-by-year shock + unit RE +
AR(1) errors; `cluster(state)` (40 clusters) in every rung; `i` = the onset cohort in every
rung. 150 draws per rung. Only K changes.

| K | hor | SD | SE csdid | SE wooldid | wd/cs | SE/SD cs | SE/SD wd | cov cs | cov wd |
|---|---|---|---|---|---|---|---|---|---|
| 3 | post0 | .2306 | .2573 | .2160 | **.839** | 1.116 | .937 | .947 | .893 |
| 3 | post1 | .2784 | .2623 | .2184 | **.833** | .942 | .785 | .900 | .853 |
| 3 | post2 | .2596 | .2465 | .2011 | **.816** | .949 | .775 | .920 | .847 |
| 3 | post3 | .2461 | .2630 | .2228 | **.847** | 1.068 | .905 | .913 | .900 |
| 5 | post0 | .1946 | .2151 | .1826 | .849 | 1.105 | .938 | .940 | .913 |
| 5 | post1 | .2237 | .2252 | .1956 | .869 | 1.007 | .875 | .947 | .847 |
| 5 | post2 | .2069 | .2159 | .1848 | .856 | 1.044 | .893 | .940 | .900 |
| 5 | post3 | .2372 | .2224 | .1897 | .853 | .938 | .800 | .927 | .900 |
| 10 | post0 | .1677 | .1738 | .1601 | .921 | 1.037 | .955 | .953 | .947 |
| 10 | post1 | .1967 | .1771 | .1637 | .924 | .900 | .832 | .907 | .887 |
| 10 | post2 | .1704 | .1653 | .1490 | .901 | .970 | .875 | .920 | .900 |
| 10 | post3 | .1774 | .1820 | .1683 | .925 | 1.026 | .948 | .927 | .927 |

Corroboration from the oracle: csdid's SE/SD stays at 0.94-1.12 across the whole ladder
while wooldid's falls to 0.775, and the paired SE ratio wd/cs moves monotonically
0.834 -> 0.857 -> 0.918 as K goes 3 -> 5 -> 10. (Read the paired ratio as corroboration
only: csdid runs mildly conservative at small K, 1.05-1.12, which inflates the gap.)
The SD column is shared because the point estimates are the same number in both packages
(DGP-S confirms this to four decimals horizon by horizon over 283 draws), so this is a
PAIRED comparison of reported uncertainty: the only quantity with MC noise is each mean SE,
whose standard error is ~1-2% at 150 draws. The wd/cs column moves monotonically
0.834 -> 0.857 -> 0.918 as K goes 3 -> 5 -> 10, and the 283-draw DGP-S run at 10 clusters
per cohort independently gives 0.912 (.1578/.1727, .1622/.1761, .1484/.1652, .1643/.1792).
A fifth point: DGP-F (five treated waves of 3 states each, 200 draws) gives 0.798-0.830.

**Where it crosses the pre-registered gate.** DGP-F, 200 draws, i=coh + cluster(state),
5 adoption waves of 3 states each — a wholly ordinary staggered-adoption design:
wooldid SE/SD 0.793-0.870 and coverage 0.860-0.915 (median 0.865 / 0.905), while csdid on
identical points sits at 0.979-1.062 and 0.915-0.950. wooldid fails both halves of the gate;
csdid passes the SE/SD half. This is the "finer-than-`i` arm that FAILS coverage" that the
review's own open-questions list names as genuinely new.

**Honest scope.** The gate's "both controls inside band" clause does not strictly fire:
did_imputation also undercovers on these designs (SE/SD 0.76-0.85), and in the two-wave K=3
rung csdid's own coverage drops to 0.90-0.947. So at three clusters per treated cell the
design is hard for everyone and wooldid is not an outlier in coverage. What is NOT shared
is the SE level: csdid's SE/SD stays at 0.94-1.12 across the whole ladder while wooldid's
falls to 0.775.

**Consequence for users.** The help's remedy sentence reads as sufficient. On a design with
three adopting states per wave it leaves a 17% SE shortfall and ~0.87 coverage, with no
runtime signal. A quantitative form of the help's own qualitative clause ("so that a treated
cohort-by-period cell contains several clusters") would be the fix — e.g. a caution when
the smallest treated cohort-by-period cell holds fewer than some number of clusters.

**Repro**: `repro_fewclus_gap.do` (self-contained, ~8 min; prints the table above for
K = 3 and K = 10). Full MC output `RESULTS_LADDER_K{3,5,10}.txt`, `RESULTS_DGP_F.txt`,
`RESULTS_DGP_S.txt`; drivers `mc_wave3.do`, `mc_wave2.do`, `mc_fineclus.do`.

---

# FINDING MC-FC-2 — `sthlp:251`'s "roughly 0.66 to 0.75" is not the size of the at-or-above-`i` understatement on realistic data; it runs to 0.18

**surface**: `cluster()` at or above `i` on the default analytic path (`sthlp:251`).
**severity**: medium  **confidence**: certain (on the measurement)  **realistic**: yes  **seeded**: partly (the effect is known; the magnitude is not)

**Claim.** The help quantifies the at-or-above-`i` understatement as "measured at roughly
0.66 to 0.75 of the reference on the diagnostic fixtures". On two ordinary state-policy
DGPs the measured ratio of mean reported SE to empirical sampling SD is far below that band,
and the resulting 95% coverage falls to 0.29.

**Measured** (all with cluster() AT the level of `i`, both controls covering normally in the
same draws in DGP-S):

| DGP | arm | att SE/SD | att cov | worst ES row SE/SD | worst ES row cov |
|---|---|---|---|---|---|
| DGP-S, 283 draws | i=state, cluster(state) | 0.610 | 0.763 | 0.430 (post2) | 0.580 |
| DGP-S, 283 draws | same, pooled `esprelength(0)` | 0.707 | 0.841 | 0.529 (post2) | 0.700 |
| DGP-F, 200 draws | i=state, cluster(state) | 0.384 | 0.595 | **0.175** (post0) | **0.294** |
| DGP-U, 92 draws  | i=state, cluster(state) | 0.589 | 0.772 | 0.452 (post2) | 0.663 |
| DGP-P, 142 draws | i=icoh,  cluster(icoh)  | 0.618 | 0.738 | 0.470 (pre4)  | 0.667 |

At 0.175 the reported standard error is 5.7 times too small and the nominal 95% interval
covers 29% of the time. The help's band is scoped to "the diagnostic fixtures", so nothing
in it is false; but a user reads 0.66-0.75 as the size of the problem and calibrates to
"my t-statistics are ~1.4x too large", when on an ordinary staggered state-policy panel the
factor is 2 to 6. The band is not a floor and the help does not say so.

Why my DGPs are harsher than the fixtures: they contain a **state-by-year shock**, the single
most standard component in applied DiD variance (Bertrand-Duflo-Mullainathan), and that is
precisely the component a saturated `i`-by-relative-period dummy absorbs entirely. Remove it
and the understatement shrinks — DGP-P, which has no cluster-level shock at all, gives
0.739 / 0.852 on the i=unit arm, squarely inside the documented band and matching the
project's own 500-draw numbers (0.838 / 0.794). So the documented band describes a DGP
without cluster-level shocks; adding the most common one doubles the damage.

**Repro**: `repro_remedy_works.do` (self-contained, ~6 min) prints the at-i arm and the
remedy arm side by side on DGP-S; `RESULTS_DGP_S.txt`, `RESULTS_DGP_F.txt` for the rest.

---

## REPRO VERIFICATION
`repro_remedy_works.do` ran to completion (`repro_remedy_works.log`, sentinel
`REPRO_RRW_DONE`), 150 draws, and printed:

    Largest |b(i=state) - b(i=coh)| over 150 draws:  6.37268e-14
    i = state, cluster(state)   [cluster AT i, the documented-bad regime]
       n = 150   empirical SD = 0.1343   mean reported SE = 0.0759
       mean(SE)/SD =  0.565     95% coverage of the true att =  0.693
    i = coh,   cluster(state)   [cluster FINER than i, THE REMEDY]
       n = 150   empirical SD = 0.1343   mean reported SE = 0.1266
       mean(SE)/SD =  0.942     95% coverage of the true att =  0.940
    csdid Tp0 (external control, cluster(state)):  mean(SE)/SD =  1.082

which reproduces the 283-draw DGP-S numbers (at-i 0.610/0.763, remedy 0.996/0.947) on an
independent seed stream. Note the shared empirical SD of 0.1343 in both arms: the two arms
are the same point estimate, so the only thing that differs is the reported uncertainty.

`repro_fewclus_gap.do` ran to completion (`repro_fewclus_gap.log`, sentinel
`REPRO_RFG_DONE`), 150 draws per rung on an independent seed stream, and printed:

      K  hor       SD   SE_csdid  SE_wooldid   wd/cs   SE/SD_cs  SE/SD_wd  cov_cs  cov_wd
      3  post0   0.2285    0.2543     0.2086   0.820      1.113    0.913   0.947   0.900
      3  post1   0.2821    0.2608     0.2122   0.814      0.924    0.752   0.880   0.847
      3  post2   0.2542    0.2636     0.2110   0.800      1.037    0.830   0.927   0.873
      3  post3   0.2656    0.2744     0.2251   0.820      1.033    0.848   0.920   0.847
     10  post0   0.1602    0.1722     0.1590   0.923      1.075    0.992   0.967   0.967
     10  post1   0.1609    0.1717     0.1587   0.924      1.067    0.986   0.947   0.927
     10  post2   0.1627    0.1697     0.1521   0.896      1.043    0.935   0.933   0.933
     10  post3   0.1759    0.1775     0.1639   0.923      1.009    0.931   0.953   0.947

Median wooldid SE/SD 0.839 at K=3 against sqrt(2/3)=0.816, and 0.961 at K=10 against
sqrt(9/10)=0.949. Median wooldid coverage 0.860 at K=3 (csdid 0.923) and 0.940 at K=10
(csdid 0.950). Same monotone pattern as the 150-draw ladder on a different seed stream.

---

## PROBED HARD AND CLEAN (negative results — what would have shown up)
1. **Cross-cutting `cluster()`** (DGP-X, 200 draws): a `cluster()` that neither nests nor is
   nested in `i` covers at 0.950 with SE/SD 0.999. Had the two-regime machinery mishandled
   the non-nesting case on the analytic path, SE/SD would have moved off 1 the way the at-i
   arms do. It does not.
2. **Ragged, unequal-sized panels** (DGP-U, 92 draws): state sizes 2-10, ~6% of unit-periods
   missing on a fixed pattern. Remedy arm 0.970 / 0.935; controls 0.975 / 0.934 and
   0.949 / 0.923. The estimand itself moves (true att 0.9060, not 0.9154) and the zero-noise
   truth tracks it; bias +0.005 = 0.04 SD, and csdid/did_imputation carry the same.
3. **Point-estimate invariance to `cluster()`**: max `mreldif(b)` over 283 (DGP-S), 200
   (DGP-X), 200 (DGP-F) and 91 (DGP-U) draws is **exactly 0** in every wave; across `i` as
   well (i=unit vs i=coh) it is 4.3e-14 over 141 draws and 6.4e-14 in the repro. Single-draw
   hex: att = +1.9416320a9e57aX-001 under `robust`, `cluster(state)` and `cluster(unit)`
   alike. `cluster()` does not move `b`.
4. **`nofastcollapse` under a finer-than-`i` cluster**: mreldif(b) 2.6e-16, mreldif(V)
   6.9e-18. The default cell collapse is exact here. Had the collapse dropped the cluster
   identity, V would have moved by percent, not by 1e-18.
5. **The `jointtests` degeneracy caution fires exactly where it should**: present for
   i=state + cluster(state), absent for cluster-finer-than-`i` and absent for a cross-cutting
   cluster. The predicate (`ado:1421-1432`, distinct `i` vs distinct `(i,cluster)` pairs) is
   the correct nesting test and is correctly applied to all three geometries.
6. **`e(N_clust)` / `e(df_r)`**: 200/199 for `cluster(unit)` on 200 units, 40/39 for
   `cluster(state)`; two-way `cluster(state region)` runs and reports the min-G convention
   (8 clusters, df 7). A cluster variable with missing values drops those rows and discloses
   the count in the header (N = 2220, 180 dropped, N_clust 37 of 40).
7. **The cell-weighted estimands `att_i` and `att_it` are not specially broken under a
   finer-than-`i` or cross-cutting cluster** on the analytic path: DGP-X 1.012/0.965 and
   0.973/0.955; DGP-U 0.982/0.935 and 0.954/0.946. They degrade in lockstep with `att` on
   the at-i arms and on the small-K designs, never differently. (`att_itime` untested.)
8. **The seeded positive control reproduces.** i=unit + cluster(i) covers at 0.852
   (`esfixedbaseperiod`) and 0.796 (pooled) against the project's own 500-draw 0.838 and
   0.794, with SE/SD 0.644 against their 0.661. If my harness were measuring something else,
   this is where it would have shown.
