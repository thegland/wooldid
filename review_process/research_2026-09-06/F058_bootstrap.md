# F058 — extended calibration simulations: the BOOTSTRAP surfaces

Author ruling 2026-09-06 ("run them"). Package tree: branch `bugfixes`, HEAD `01a7f7a`,
read-only. New files, all under `test/suite/mc/`:

| file | what it is |
|---|---|
| `mc_bootse2_dgp.doh` | the two designs, shared by the driver and the verification script |
| `mc_bootse2.do` | the driver: one CSV row per arm run, appended as it completes |
| `run_bootse2.sh` | shard launcher with the instance-cap wait |
| `sched_bootse2.sh` | wall-clock time-slicer (see "What the budget bought") |
| `analyze_bootse2.py` | the tables below |
| `smoke2/smoke2.do`, `smoke2/verify2.do` | structure/timing probe and the design checks |
| `out2/bootstrap/mcb2_<cfg>.csv` | the raw rows |
| `out2/bootstrap/tables.md` | the generated tables (reproduced below) |

Re-run or top up with, e.g.,
`./run_bootse2.sh s3 AN:61:200` — the driver appends and the analyzer globs, so more
rows simply sharpen the same cells.

## The two designs

**Design A** is `mc_bootse.do`'s design A unchanged: 24 clusters of 2 units, T = 6,
clusters strictly above the unit, cohorts `ttre` in {3, 5, never} by `mod(cl,3)`, a
cluster-by-period shock so a clustered variance is actually needed, within-unit AR(1),
`cluster(cl)`. The treatment effect is a constant on every treated post cell, so the
true `att` is that constant **under any event-study window** — verified noiselessly at
2.0000000000000 with `esprelength(2) espostlength(1)`.

**Design B is new and is the unequal-cluster design the ruling asked for**: 12 states,
state *s* holding *s* units (78 units), T = 8, odd states treated at t = 5, even states
never treated, τ_s = 0.1·n_s, `cluster(state)`. Six treated clusters of sizes 1, 3, 5,
7, 9, 11 — the cluster sizes span 12× and the treatment effect is perfectly correlated
with cluster size, which is the case a cluster bootstrap is most exposed on. Because
every treated unit carries the same four post periods, the unit-weighted mean of τ, the
observation-weighted `att` and the cell-weighted `att_it` all coincide exactly:

    truth = 0.1 · Σ_{s odd} s² / Σ_{s odd} s = 0.1 · 286/36 = 0.7944444444…

confirmed to 16 digits by a noiseless run reporting `main_att` = `main_att_it` =
0.79444444.

## Arms

`reps(199)` everywhere, `rseed()` set per dataset per arm, `dots(0)`, `safety(off)`,
level left at the **default 95**.

| cfg | design | DGP | arms run per dataset |
|---|---|---|---|
| `AN` | A | τ = 0 | analytic; `bootstrap(standard)`; `bootstrap(bayes)` — all three with `esprelength(2) espostlength(1) jointtests` |
| `AT` | A | τ = 0 | analytic; `bootstrap(standard) boottstat` with the same window and `jointtests` |
| `BA` | B | τ_s = 0.1·n_s | analytic; `bootstrap(standard) bootdisplay(bca)`; `bootstrap(bayes)` |
| `BN` | B | τ = 0 | the same three |
| `BI` | B | τ_s = 0.1·n_s | analytic `att_it`; `att_it bootstrap(standard)` — the physical-resample re-identification leg |

**Design A is run at τ = 0, and that is not a substitution — it is the same run.**
The ruling asked for design A with a true `att` of 2. Because the treatment indicator
lies in the column space of the saturated cohort×period design, shifting `y` by τ·D
shifts the estimate, every replication and therefore every interval endpoint by exactly
τ. That was checked, not assumed (see "Design checks", item 2: max deviation 3.1e-15).
So the τ = 0 rows report coverage of the true `att` *and*, endpoint for endpoint, what
the τ = 2 design would report — and they additionally deliver the size arms and, through
the level-independent dual p-values, coverage at 90 as well as 95. Design B's τ_s varies
by state, so no such identity holds there and its alternative and null are separate
configs.

Two economies are worth recording because they cost nothing and halve the bill:

* **`bootdisplay(bca)` is a superset run.** Every family is computed whichever is
  displayed, so one `bootdisplay(bca)` fit yields the percentile, bias-corrected, Wald
  *and* BCa intervals and all four p-values — the BCa jackknife (12 extra fits on design
  B) is the only marginal cost. Arms (2) and (4) of the ruling are therefore one run.
* **`bootstrap(standard)` with the event-study window is also a superset run**: the same
  199 replications deliver the main-`att` interval families *and* the fixed-metric
  joint pre-trend test. Arms (2) and (6) are one run.

One thing the event-study window is **not** free of: it changes the estimator. On the
same design-A dataset the plain fit gave `main_att` = 1.7817 (se .218) and the
`esprelength(2) espostlength(1)` fit gave 1.4835 (se .349) — the windowed run estimates
the att on the reported window only. The estimand (and so the truth) is unchanged, but
every comparison below uses each arm's **own** estimates as the sampling-SD denominator
rather than a shared one, and the analytic comparator in each config carries the same
window as the bootstrap arms it sits beside.

## Levels: why 95 is measured and 90 is derived

`wooldid_bootci` computes the percentile p-value as the exact inversion of the
percentile interval, and the bias-corrected p as the level at which a BC endpoint lands
on the null; both are **level-independent** and are documented as duals of their
intervals. So under a τ = 0 DGP the rejection rate at α *is* the non-coverage rate of the
truth at level 100(1−α), and one set of replications delivers coverage at 95 (from the
intervals) and at 90 (from the p-values) at once. The Wald family is recomputed
directly as b ± t_{.95}(df)·se, with df recovered from the reported endpoints.

The one cell no such trick reaches is percentile/BC/BCa at **90** under the
*heterogeneous alternative* (config `BA`), where the p tests 0 and the truth is 0.7944.
It is left blank rather than bought with a second set of replications.

## What the budget bought — read this before reading the tables

This is the part of the result that matters most for what to do next, so it goes first.

A single `wooldid` refit on these designs costs about **0.3 s** (design A, 288
observations) to **0.75 s** (design B with the event-study window, 624 observations) on
this machine while six Stata instances share its 12 cores. A `reps(199)` arm is
therefore ~60 s, a `boottstat` arm ~205 s, and the ruling's own sizing — S = 200 for the
cheap arms and 100 for `bca`/`boottstat`/re-id, across both designs, alternative and
null — comes to roughly **190,000 core-seconds ≈ 53 core-hours ≈ 27 hours of wall clock
at the two-shard cap**. Against a ~90-minute ceiling that is a factor of about 18.

Two ways to spend 90 minutes were available: keep `reps(199)` and cut S, or cut reps and
buy S. **`reps(199)` was kept**, on the grounds that the 199-replication percentile tail
is the object being calibrated (a 99-rep percentile interval is a different estimator,
not a cheaper measurement of the same one) and that a low-S run at 199 reps can be
topped up into the full result later by raising `LAST`, whereas a 99-rep run cannot.
The cost of that choice is precision: **the n per cell below is in the low tens, so a
Monte-Carlo standard error of 3–8 percentage points is normal and only gross
miscalibration is detectable.** Cells are flagged at |observed − nominal| > 2 MCSE
computed *under correct calibration*, √(p₀(1−p₀)/n), and no cell with n < 15 is flagged
at all.

Two of the ruling's surfaces got no slice and are **not measured**: design B's
fixed-metric joint pre-trend test (`BJ`) and design B's `boottstat` (`BT`). Both drivers
exist and are one command away.

*Aside for the orchestrator:* the prescribed instance-cap probe,
`pgrep -fl stata-mp | grep -c " -e do "`, **over-counts badly** — sibling agents' own
polling shells contain that pattern in their command lines and match it. It reported 9–12
instances while only 3–4 Stata processes were actually running, and it would have
deadlocked this part indefinitely. The scripts here count
`^[0-9]* /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do ` instead.

## Design checks run before the simulation (all passed)

1. **Windowed truth.** Noiseless fits with `esprelength(2) espostlength(1)`: design A's
   windowed `att` = 2.0000000000000 (16 digits), design B's windowed `att` and `att_it`
   both 0.7944444444444447 against the exact 0.7944444444444444. The event-study window
   changes the estimator but not the estimand.
2. **Shift equivariance on design A — this is what lets the τ = 0 run stand for the
   τ = 2 run the ruling asked for.** Same dataset, same `rseed`, `y` versus `y − 2·D`:
   `b`, `se`, and the percentile, bias-corrected and Wald endpoints all moved by exactly
   2 (max deviation 3.1e-15), and the joint pre-trend p-value was byte-identical
   (0.47 both ways). D lies in the column space of the saturated cohort×period design,
   so every replication of the resample shifts by exactly τ; the whole bootstrap
   distribution translates. **Consequence:** design A's `AN` rows report coverage of the
   true `att` at τ = 0, and by this identity they are *the same numbers* the ruling's
   τ = 2 design A would produce. Nothing is inferred loosely here — the endpoints were
   checked, not argued.
3. `att att_it` together with an event-study window returns rc = 0 (an earlier r(198)
   in the scratch script was a degenerate y ≡ 0 dataset, not this option combination —
   worth one line only because "invalid syntax" is an odd message for a constant
   outcome, and it fires only on a hopeless dataset).

## What was actually run

Compute ran 16:17–17:47 (**90 minutes exactly**, two shards, `reps(199)` throughout,
level at the default 95). A wall-clock time-slicer (`sched_bootse2.sh`) retired each
config at a set moment so the first one did not eat the budget; because the driver
appends a row per completed arm, a kill costs at most the arm in flight.

| cfg | datasets | why it stopped there |
|---|---|---|
| `AN` (design A, std + bayes + analytic, ES window, jointtests) | **25** (24 with a completed `std` arm, 23 with `bayes`) | time-sliced at 17:12, topped up 17:35–17:47 |
| `BA` (design B, std+BCa + bayes + analytic) | **15** (14 with `bayes`) | time-sliced at 17:10 |
| `BI` (design B, `att_it` re-identification leg + analytic) | **29** (28 with the bootstrap arm) | ran to the deadline |
| `AT` (design A, `boottstat` + analytic) | **10** | hit its own `LAST=10` at ~17:35 |
| `BN` (design B, τ = 0 size arms) | **0** | dropped: at ~200 s a dataset the remaining slice bought ~6, against ~28 for `BI` at ~63 s |
| `BJ`, `BT` (design B joint test; design B `boottstat`) | **0** | no slice |

Zero `rc != 0` in all **193** arm runs. All bootstrap arms completed 199 of 199 replications
except the `att_it` leg, where 3 of 28 runs lost 1–2 replications (mean
`e(boot_N_V)` 198.9) — replications that drew no treated cluster, which is the
documented behaviour. `e(V_analytic)` reproduced the separate analytic run's standard
error to **0.0 absolute difference on all 24 design-A datasets**.

## Run inventory

| cfg | arm | n rows | rc!=0 | mean secs | reps completed | reps failed | boot_N_V |
|---|---|---|---|---|---|---|---|
| AN | bayes | 23 | 0 | 114.8 | 199.0 | 0.0 | 199.0 |
| AN | plain | 25 | 0 | 1.1 | - | - | - |
| AN | std | 24 | 0 | 63.8 | 199.0 | 0.0 | 199.0 |
| AT | bts | 10 | 0 | 111.3 | 199.0 | 0.0 | 199.0 |
| AT | plain | 10 | 0 | 0.8 | - | - | - |
| BA | bayes | 14 | 0 | 136.2 | 199.0 | 0.0 | 199.0 |
| BA | plain | 15 | 0 | 0.7 | - | - | - |
| BA | stdbca | 15 | 0 | 83.2 | 199.0 | 0.0 | 199.0 |
| BI | attit | 28 | 0 | 76.3 | 198.9 | 0.1 | 198.9 |
| BI | plainit | 29 | 0 | 0.4 | - | - | - |

## Point estimate and the standard-error ratio

(ratio = mean reported SE / SD of the estimate across datasets; 1.00 = calibrated. The SD denominator is that arm's OWN estimates, because the event-study window changes the estimator.)

| cfg | arm | n | truth | mean b | bias | bias t | SD(b) | mean SE | ratio (SE) |
|---|---|---|---|---|---|---|---|---|---|
| AN | bayes | 23 | 0.00000 | 0.1193 | +0.1193 | +1.76 | 0.3244 | 0.3081 | 0.950 (0.148) |
| AN | plain | 25 | 0.00000 | 0.1189 | +0.1189 | +1.76 | 0.3370 | 0.3450 | 1.024 (0.157) |
| AN | std | 24 | 0.00000 | 0.1002 | +0.1002 | +1.48 | 0.3308 | 0.3453 | 1.044 (0.160) |
| AT | bts | 10 | 0.00000 | 0.2118 | +0.2118 | +1.66 | 0.4035 | 0.3324 *coefficient-draw SD | 0.824 (0.202) |
| AT | plain | 10 | 0.00000 | 0.2118 | +0.2118 | +1.66 | 0.4035 | 0.3368 | 0.835 (0.211) |
| BA | bayes | 14 | 0.79444 | 0.8066 | +0.0121 | +0.14 | 0.3248 | 0.3030 | 0.933 (0.199) |
| BA | plain | 15 | 0.79444 | 0.7821 | -0.0123 | -0.15 | 0.3270 | 0.2400 | 0.734 (0.146) |
| BA | stdbca | 15 | 0.79444 | 0.7821 | -0.0123 | -0.15 | 0.3270 | 0.3612 | 1.105 (0.228) |
| BI | attit | 28 | 0.79444 | 0.7695 | -0.0250 | -0.45 | 0.2947 | 0.3417 | 1.159 (0.172) |
| BI | plainit | 29 | 0.79444 | 0.7636 | -0.0308 | -0.57 | 0.2911 | 0.2045 | 0.702 (0.104) |

## Interval coverage of the true att, nominal 95 (percent, MCSE in brackets)

| cfg | arm | percentile | bias-corrected | BCa | Wald | displayed |
|---|---|---|---|---|---|---|
| AN | bayes | 87.0 (7.0) n=23 | 87.0 (7.0) n=23 | - | 87.0 (7.0) n=23 | 87.0 (7.0) n=23 |
| AN | plain | - | - | - | - | 96.0 (3.9) n=25 |
| AN | std | 87.5 (6.8) n=24 | 87.5 (6.8) n=24 | - | 91.7 (5.6) n=24 | 87.5 (6.8) n=24 |
| AT | bts | 70.0 (14.5) n=10 | 70.0 (14.5) n=10 | - | 70.0 (14.5) n=10 | 100.0 (0.0) n=10 |
| AT | plain | - | - | - | - | 90.0 (9.5) n=10 |
| BA | bayes | 92.9 (6.9) n=14 | 85.7 (9.4) n=14 | - | 92.9 (6.9) n=14 | 92.9 (6.9) n=14 |
| BA | plain | - | - | - | - | 86.7 (8.8) n=15 |
| BA | stdbca | 93.3 (6.4) n=15 | 93.3 (6.4) n=15 | 86.7 (8.8) n=15 | 100.0 (0.0) n=15 | 86.7 (8.8) n=15 |
| BI | attit | 96.4 (3.5) n=28 | 92.9 (4.9) n=28 | - | 100.0 (0.0) n=28 | 96.4 (3.5) n=28 |
| BI | plainit | - | - | - | - | **82.8 (7.0) n=29** |

## Interval coverage of the true att, nominal 90 (percent, MCSE in brackets)

(percentile / BC / BCa taken from the level-independent dual p-values, so this row exists only where the DGP is the null the p tests; Wald recomputed as b +/- z*se.)

| cfg | arm | percentile | bias-corrected | BCa | Wald |
|---|---|---|---|---|---|
| AN | bayes | 87.0 (7.0) n=23 | 87.0 (7.0) n=23 | - | 87.0 (7.0) n=23 |
| AN | plain | - | - | - | 88.0 (6.5) n=25 |
| AN | std | 87.5 (6.8) n=24 | 87.5 (6.8) n=24 | - | 87.5 (6.8) n=24 |
| AT | bts | 70.0 (14.5) n=10 | 70.0 (14.5) n=10 | - | - |
| AT | plain | - | - | - | 70.0 (14.5) n=10 |
| BA | bayes | n/a | n/a | n/a | 85.7 (9.4) n=14 |
| BA | plain | n/a | n/a | n/a | 86.7 (8.8) n=15 |
| BA | stdbca | n/a | n/a | n/a | 93.3 (6.4) n=15 |
| BI | attit | n/a | n/a | n/a | 92.9 (4.9) n=28 |
| BI | plainit | n/a | n/a | n/a | 82.8 (7.0) n=29 |

## Mean interval width by family (same datasets, so widths are directly comparable)

| cfg | arm | n | percentile | bias-corrected | BCa | Wald | displayed |
|---|---|---|---|---|---|---|---|
| AN | bayes | 23 | 1.2410 | 1.2056 | - | 1.2746 | 1.2410 |
| AN | plain | 25 | - | - | - | - | 1.4275 |
| AN | std | 24 | 1.3589 | 1.3355 | - | 1.4285 | 1.3589 |
| AT | bts | 10 | 1.3170 | 1.2846 | - | 1.3752 | 1.7543 |
| AT | plain | 10 | - | - | - | - | 1.3933 |
| BA | bayes | 14 | 1.1849 | 1.1688 | - | 1.3340 | 1.1849 |
| BA | plain | 15 | - | - | - | - | 1.0566 |
| BA | stdbca | 15 | 1.4299 | 1.3506 | 1.3440 | 1.5899 | 1.3440 |
| BI | attit | 28 | 1.3532 | 1.3176 | - | 1.5040 | 1.3532 |
| BI | plainit | 29 | - | - | - | - | 0.9001 |

## Paired standard-error ratios (per dataset, so free of sampling noise in b)

(each bootstrap arm's reported SE divided by the SAME dataset's analytic clustered SE, from the analytic arm run in the same config with the same event-study window.)

| cfg | arm | n | mean se_boot/se_analytic (SE) | min | max |
|---|---|---|---|---|---|
| AN | bayes | 23 | 0.942 (0.042) | 0.668 | 1.531 |
| AN | std | 24 | 1.039 (0.039) | 0.759 | 1.632 |
| AT | bts | 10 | 1.043 (0.091) *coefficient-draw SD from the hidden Wald table, not the displayed studentized se | 0.699 | 1.514 |
| BA | bayes | 14 | 1.276 (0.065) | 0.908 | 1.608 |
| BA | stdbca | 15 | 1.526 (0.102) | 0.932 | 2.274 |
| BI | attit | 28 | 1.723 (0.093) | 1.047 | 2.951 |

## p-value size under the null DGP (rejection percent, MCSE in brackets)

| cfg | arm | p_percentile @5 | @10 | p_bc @5 | @10 | p_bca @5 | @10 |
|---|---|---|---|---|---|---|---|
| AN | bayes | 13.0 (7.0) n=23 | 13.0 (7.0) n=23 | 13.0 (7.0) n=23 | 13.0 (7.0) n=23 | - | - |
| AN | std | 12.5 (6.8) n=24 | 12.5 (6.8) n=24 | 12.5 (6.8) n=24 | 12.5 (6.8) n=24 | - | - |
| AT | bts | 30.0 (14.5) n=10 | 30.0 (14.5) n=10 | 30.0 (14.5) n=10 | 30.0 (14.5) n=10 | - | - |

## Joint pre-trend test (parallel pre-trends hold in both DGPs, so this is size)

| cfg | arm | family | n | rej @5 | rej @10 |
|---|---|---|---|---|---|
| AN | bayes | fixedmetric | 23 | 4.3 (4.3) n=23 | 4.3 (4.3) n=23 |
| AN | std | fixedmetric | 24 | 0.0 (0.0) n=24 | 0.0 (0.0) n=24 |
| AT | bts | studentized | 10 | 0.0 (0.0) n=10 | 0.0 (0.0) n=10 |

## Cells more than 2 MC standard errors from nominal

- BI/plainit main coverage@95 = 82.8 (n=29, MCSE 7.0)

*(The flag list above suppresses cells with n < 15. Two cells with n = 10 are past
2 MC standard errors and are named in the conclusions: `AT`'s non-displayed
coefficient-draw percentile family, at −3.63 MCSE on coverage and +3.63 on size.)*

---

# Conclusions, one surface at a time

Design A's rows are the τ = 0 run, which by the verified shift identity is endpoint-for-
endpoint the τ = 2 design the ruling specified.

Distances below are in MC standard errors computed under correct calibration,
√(p₀(1−p₀)/n). Read them as "how far, given how little data" — at n ≈ 15–28 a
1.7-MCSE gap is a hint, not a verdict.

## 1. The analytic clustered variance (the comparator) — *this is where the real signal is*

First, a fact that simplifies the design-B rows and was checked rather than assumed: on
design B the `att` and `att_it` estimators are **numerically the same thing**. On all 15
datasets both configs ran, the point estimates agree to 0.0 and the analytic standard
errors to a ratio of 1.000 (se 0.000). That is what the balanced panel forces — every
treated unit carries the same four post periods, so the observation weights and the cell
weights coincide realisation by realisation, not merely in expectation. The design-B
analytic result is therefore **one** measurement over all 29 datasets, not two.

| design | estimand | mean SE / SD | coverage @95 | distance |
|---|---|---|---|---|
| A, 24 equal clusters | `att` (windowed) | 1.024 (0.157) | 96.0% (n=25) | +0.23 |
| B, 12 unequal clusters | `att` ≡ `att_it` | 0.702 (0.104) | **82.8% (24 of 29)** | **−3.02** |

**Plain conclusion:** the analytic clustered variance is well calibrated on the balanced
design and **materially anti-conservative on the unequal-cluster design** — about 30% too
small, and its nominal-95 interval covers 83%. This is the only cell in the study past
2 MC standard errors on its own n, and it is the *analytic* path, not a bootstrap one.
Design B is deliberately adversarial (six treated clusters of sizes 1, 3, 5, 7, 9, 11
with the effect proportional to cluster size), but that is exactly the setting the help
already tells users to reach for the bootstrap in, and this quantifies why. (On the
first 15 datasets alone the same quantity read 86.7%; the extra 14 draws moved it down,
which is the ordinary behaviour of a rate at n = 15.)

## 2. `bootstrap(standard)` — percentile, bias-corrected, Wald

| design | family | coverage @95 | distance | coverage @90 | size @5 |
|---|---|---|---|---|---|
| A (n=24) | percentile | 87.5% | −1.69 | 87.5% | 12.5% (+1.69) |
| A (n=24) | bias-corrected | 87.5% | −1.69 | 87.5% | 12.5% (+1.69) |
| A (n=24) | Wald | 91.7% | −0.75 | 87.5% | — |
| B (n=15) | percentile | 93.3% | −0.30 | n/a | — |
| B (n=15) | bias-corrected | 93.3% | −0.30 | n/a | — |
| B (n=15) | Wald | 100.0% | +0.89 | 93.3% | — |

**Plain conclusion:** no surface is off by more than 2 MC standard errors. On design A
the percentile and bias-corrected families sit ~7 points below nominal at 95 (1.7 MCSE,
n = 24) while the Wald family from the same replications sits ~3 points below; the
scale is not the problem (mean SE / SD = 1.044, and the paired
se_boot/se_analytic = **1.039 ± 0.039**, i.e. the bootstrap standard error reproduces
the analytic one on the balanced design), so if the design-A gap is real it is a
tail-shape effect in the order statistics, not a variance-magnitude one. On design B the
bootstrap **repairs** the analytic shortfall: paired se_boot/se_analytic = **1.526 ±
0.102**, mean SE / SD moves from 0.702 (analytic) to 1.105, and coverage from 86.7% on
those same 15 datasets — 82.8% over all 29 — to 93.3%. That is the intended behaviour of
the feature, measured.

## 3. `bootstrap(bayes)`, `bootalpha(equal)`, with the 2026-09-06 per-cluster-multiplier fix

| design | coverage @95 (percentile) | mean SE / SD | paired se_bayes / se_analytic |
|---|---|---|---|
| A, balanced (n=23) | 87.0% (−1.77) | 0.950 (0.148) | 0.942 (0.042) |
| B, unequal (n=14) | 92.9% (−0.37) | 0.933 (0.199) | 1.276 (0.065) |

**Plain conclusion:** bayes tracks `standard` closely and is not off by 2 MCSE anywhere,
but it runs consistently **below** it: 0.942 vs 1.039 against the analytic on design A
(a 9% gap) and 1.276 vs 1.526 on design B (a 16% gap). The design-A gap is consistent in
size and sign with the one the 2026-08-30 report saw (0.877 vs 0.967 against the truth);
the two are different estimators on different draws, so read that as agreement in
direction, not a matched number. Design A is not expected to move under the fix anyway:
every cluster there holds two units, so the per-cluster multiplier and the old total-mass
split coincide up to scale. The new information is
design B, where they do not coincide: even with the fix, the Bayesian bootstrap recovers
only about half of the analytic shortfall that the pairs bootstrap recovers in full.
The bias-corrected family under bayes on design B covered 85.7% (−1.59, n = 14) —
the lowest bootstrap cell on that design, and the one most worth a top-up.

## 4. `bootdisplay(bca)` (design B, n = 15)

BCa completed all 12 delete-one-cluster jackknife fits on every dataset — 12 extra fits
on top of 199, about 6% more work, and `e(boot_jk_completed)` read 12 every time.
Coverage 86.7% at nominal 95 (−1.48), and
it is the **narrowest** bootstrap family on that design — mean width 1.344 against 1.430
(percentile), 1.351 (BC) and 1.590 (Wald). So BCa is doing what BCa is supposed to do,
tightening the interval relative to percentile, and on this design it is the only
bootstrap family that does not over-cover. At n = 15 that is a 1.5-MCSE hint, not a
finding; it is the second-highest-value top-up after design B's size arms.

## 5. `boottstat` percentile-t (design A, n = 10) — and one thing to look at

The displayed percentile-t interval covered **100% (10/10)** at nominal 95, with a mean
width of 1.754 against the analytic interval's 1.393 (26% wider). `AT` and `AN` were run
on the *same seeds* — the point estimates agree to 0.0 on all ten shared datasets — so
the comparison is paired: on those ten datasets `AN`'s percentile family covered **7 of
10** and its Wald family 8 of 10, against percentile-t's 10 of 10. The studentized joint
pre-trend test rejected 0 of 10.

**The thing to look at.** A `boottstat` run posts `e(boot_ci_percentile)`,
`e(boot_ci_bc)`, `e(boot_p_percentile)` and `e(boot_p_bc)` alongside the studentized
output, and those hold the **coefficient-draw** families — the same names, a different
statistic from the percentile-t interval the run displays. That every family is computed
whichever is displayed is documented general behaviour, so this is not a hidden surface;
the point is narrower and is about consumption. On these ten datasets the by-name
percentile family covered 70% and rejected the true null 30% of the time, while the
interval the user actually sees covered 100%. Someone who reads `e(boot_ci_bc)` by name
after a `boottstat` run — a plausible thing for a script to do — gets a materially
different interval from the reported one. n = 10, so treat the 70/30 figures as
directional; the *structural* point is exact and independent of n.

## 6. Joint pre-trend tests (parallel pre-trends hold, so this is size)

| family | arm | n | rej @5 | rej @10 |
|---|---|---|---|---|
| `fixedmetric` | `bootstrap(standard)` | 24 | 0.0% (−1.12) | 0.0% (−1.63) |
| `fixedmetric` | `bootstrap(bayes)` | 23 | 4.3% (−0.14) | 4.3% (−0.90) |
| `studentized` | `boottstat` | 10 | 0.0% (−0.73) | 0.0% (−1.05) |

`e(boot_jt_family)` read `fixedmetric` on every non-`boottstat` run and `studentized` on
every `boottstat` run, as designed. **Plain conclusion:** neither family shows any sign
of over-rejection; both look mildly conservative (0 of 24 at a nominal 10% where 2.4 were
expected), but nothing reaches 2 MC standard errors. This is the reassuring direction —
an over-sized pre-trend test would be the damaging failure — but the conservatism is
consistent across all three cells and deserves a larger n.
## 7. `att_it` and the two bootstrap legs (design B) — the clearest result

Because `att` and `att_it` coincide numerically here (§1), design B gives something
better than a coverage number for the re-identification leg: a clean isolation of the
**leg itself**. The same estimator, the same datasets, the same point estimates — the
only difference is that the `att` run replicates by cluster frequency reweighting while
the `att_it` run physically resamples and re-identifies duplicated clusters.

| | mean SE / SD | coverage @95 | distance |
|---|---|---|---|
| analytic (`att` ≡ `att_it`) | 0.702 (0.104) | 82.8% (n=29) | **−3.02** |
| `bootstrap(standard)`, frequency reweighting, percentile | 1.105 (0.228) | 93.3% (n=15) | −0.30 |
| `bootstrap(standard)`, physical re-identification, percentile | 1.159 (0.172) | **96.4% (n=28)** | **+0.35** |
| … same run, bias-corrected | — | 92.9% (n=28) | −0.52 |
| … same run, Wald | — | 100.0% (n=28) | +1.21 |

Paired within dataset against the analytic standard error: reweighting **1.526 ± 0.102**,
re-identification **1.571 ± 0.085** (both n = 15) — a 53–57% widening, six or more
standard errors from 1. Directly against each other on those same 15 datasets, the
re-identification leg's standard error is **1.046 ± 0.029** times the reweighting leg's.

**Plain conclusion, in two parts.** First, the strong one: on the unequal-cluster design
the analytic interval is badly anti-conservative (83% at nominal 95, the only 2-MCSE
cell in the study) and **the bootstrap fixes it** — the re-identification leg lands at
96.4%, within a third of an MC standard error of nominal, on the best-measured cell here
(n = 28). That is about 13 percentage points of coverage bought by the bootstrap.
Second, the careful one: the re-identification leg's extra cell-sampling variance
component shows up as roughly a **5% wider standard error** than reweighting on this
design, in the direction the help predicts ("reweighting … omits the cell-sampling
variance component and can understate those standard errors when treatment effects vary
across clusters" — and here they vary with cluster size by construction). But 1.046 sits
only 1.6 standard errors above 1, and the two legs ran on different `rseed`s, so part of
that spread is replication Monte-Carlo noise. **The direction is as designed; the
magnitude is not yet established.** Both legs covered 14 of the same 15 datasets, so
nothing here says a user is worse off with either — the cost of the leg is runtime, and
what it buys on this design is small but correctly signed.

# Surfaces off by more than 2 MC standard errors

1. **The analytic clustered interval on the unequal-cluster design** (`att`, which is
   numerically `att_it` here) — coverage 82.8% at nominal 95, n = 29, −3.02 MCSE. Not a
   bootstrap defect: the bootstrap arm on the same datasets is at 96.4%.
2. **The coefficient-draw percentile/BC families that `boottstat` posts under the usual
   `e()` names** — coverage 70% and size 30% **at n = 10 only**, against a displayed
   percentile-t interval that covered 100% on the same datasets. This is below the n ≥ 15
   flag gate and the rates are directional, not established; it is named here because the
   structural half of it — reading those matrices by name after a `boottstat` run returns
   a different family from the reported one — does not depend on n.

Nothing else in the study reaches 2 MC standard errors. The largest remaining gaps, in
order, are design A's percentile and bias-corrected coverage (−1.69, −1.77 at n ≈ 24)
and design B's bayes bias-corrected and BCa coverage (−1.59, −1.48 at n ≈ 15).

# What to run next, in priority order

The driver is resumable: `./run_bootse2.sh s3 <CFG>:<first>:<last>` appends and
`analyze_bootse2.py` globs, so every number above sharpens rather than being replaced.

1. `AN:19:28` then `AN:36:200`, and `BA:16:200` — the two flagship configs, at n = 25 and
   15. (`AN` holds draws 1–18 and 29–35: the time-slicer cut it at 18 and the top-up
   restarted at 29. The gap is only a gap in the seed sequence, not in validity.) Taking
   both to n = 150 would settle design A's −1.7 MCSE percentile gap and design B's BCa.
   ≈ 11 h.
2. `BN:1:150` — design B's size arms, the only requested surface with **no** measurement
   at all, and the sharpest test of the p-values under unequal clusters. ≈ 8 h.
3. `AT:11:150` — `boottstat`, currently n = 10, and the config that raised the
   posted-family question in §5. ≈ 8 h.
4. `BJ`, `BT` — design B's joint test and `boottstat`. ≈ 12 h.

# Files

Drivers and analyzer: `test/suite/mc/mc_bootse2.do`, `mc_bootse2_dgp.doh`,
`run_bootse2.sh`, `sched_bootse2.sh`, `analyze_bootse2.py`.
Probes: `test/suite/mc/smoke2/smoke2.do` (structure and timing), `smoke2/verify2.do`
(the three design checks).
Data: `test/suite/mc/out2/bootstrap/mcb2_{AN,AT,BA,BI}.csv` (193 arm rows),
`out2/bootstrap/tables.md` (regenerated tables), `out2/bootstrap/logs/`.
