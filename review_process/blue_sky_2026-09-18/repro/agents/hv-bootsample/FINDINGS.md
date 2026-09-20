# hv-bootsample FINDINGS

## VERDICT ON THE SEEDED CLAIM: **TRUE**, and the consequence surface is larger
## than the claim stated. Mechanism as described. Measured, not inferred.

### HV-BOOT-01 (high) — the bootstrap engine's cluster universe is the raw-`if`
### sample, so few-cluster SAFETY GUARDS are silently suppressed and bootdf() is wrong
Repro: `p1_guards.do` (log `p1_guards.log`). Oracle = physically `drop if missing(y)` and
re-run identically; wooldid's own markout makes the OBSERVED estimate bit-identical, so every
difference below is purely the sample mismatch.

D1 — 20 states, 8 treated, 4 treated states have y missing everywhere (16 live clusters,
4 live treated clusters):
| | full data | dead rows dropped |
|---|---|---|
| e(N) | 320 | 320  (identical sample) |
| point est | .65513021 | .65513021 (identical) |
| **e(boot_df)** | **19** | **15** |
| "only K treated cluster(s)" caution | **NOT PRINTED** | printed ("only 4") |

D2 — 12 states, 6 treated, 4 treated states dead (8 live clusters, 2 live treated):
| | full data | dead rows dropped |
|---|---|---|
| e(N) | 160 | 160 |
| **e(boot_df)** | **11** | **7** |
| "only N clusters are available" (<10) caution | **NOT PRINTED** | printed ("only 8") |
| "only K treated cluster(s)" (<6) caution | **NOT PRINTED** | printed ("only 2") |

D3 — 6 states, 5 dead, 1 live cluster: the full-data run PASSES the screen whose message reads
"Error: bootstrap() requires at least 2 clusters **in the estimation sample**"; the dropped-row
run correctly refuses r(459) with that message. The message names a sample it does not count.

Mechanism (verified in code): `wooldid_bootctl` ado:6459-6463 builds `smp2 = ((ifr)==1) &
!missing(clus)` and all three guards count on it; `wooldid_bootrun` ado:8908-8923 builds the
same unscreened `smp`, and `wdb_C` (ado:8956) → `wdb_r_nclus` (ado:9537) → bootdf(clusters)
df = G-1 (ado:9570) → `e(boot_df)` (ado:10233).

Why it matters: the <10-cluster and <6-treated-cluster cautions exist *precisely* to warn about
the setting the few-cluster literature says the pairs cluster bootstrap gets badly wrong. They
are the last line of defence and they go silent exactly when dead clusters pad the count.
The help's own text for `bootstrap()` promises "A caution fires when fewer than 6 treated
clusters are available for resampling (beside the existing caution at fewer than 10 clusters
in total)" — it does not.

Contradicted record: ado:666's comment says the six forwarded copies of the `if` expression mean
"the estimation sample, the permutation design and the bootstrap frame agree". They agree on the
exact-one rule (BS-37, not contradicted) but NOT on markout. ado:9367 shows the author knows the
resample frame holds rows e(sample) drops ("missing y, an event-study window, singletons"), so
the guard/df consequences read as an oversight rather than a ruling.

---

### HV-BOOT-02 (high) — `bootalpha(size)`: the Dirichlet concentration is read off the
### UNSCREENED frame, so the bootstrap SE is understated by a further factor sqrt(f)
Repro: `p2_size_degen.do` section B (log `p2_size_degen.log` lines 228-407).
14 states x 2 units x 12 years; y missing at random for ~half the cells (an outcome with
partial survey coverage). NO cluster is dead — every cluster keeps roughly half its rows.

| | full data | dead rows dropped (oracle) |
|---|---|---|
| e(N) | 164 | 164 (identical sample) |
| point estimate | .64896714 | .64896714 (identical) |
| **bootstrap SE** | **.0309309** | **.04522715** |
ratio = **0.6839**. Live fraction f = 164/336 = .48810; predicted ratio sqrt(f) = .6986.
MC noise on the ratio at reps(400) is about 5%, so the prediction is confirmed.
**The reported SE is 32% too small.**

Mechanism: `__wdbn` (ado:8946) is the cluster's row count in the map frame, which is built from
`smp` = raw-if ∧ nonmissing cluster. Under `bootalpha(size)` that count is BOTH the Gamma shape
and the divisor (help sthlp:548: "the Gamma shape n_c, its observation count"). Rows outside
e(sample) inflate n_c, so the Dirichlet concentration is too high by 1/f and the cluster-weight
CV^2 shrinks by f — variance x f, SE x sqrt(f).
This is NOT the understatement sthlp:548 already documents ("tying the Dirichlet's concentration
to the observation count rather than the cluster count shrinks the cluster-weight dispersion by
roughly the square root of the mean cluster size"). That one is about n_c vs G and is the
documented reason size is not a sound default. The factor here is an ADDITIONAL, undocumented
shrinkage that appears only because the count is taken off the wrong sample, and it scales with
how unbalanced the panel is — which the user cannot see anywhere in the output.

### HV-BOOT-03 (medium) — `bootsave()`'s observed row is on the wrong scale, contradicting the help
Repro: `p2_size_degen.do` section C.  Same run as HV-BOOT-02.
sthlp:557 promises: under `bootalpha(size)` "the _rep = 0 row holds 1/N, N the number of
observations in the **estimation sample**".
Observed: `_rep==0` holds w = .00297619 = **1/336**. e(N) = **164**. The oracle file holds
1/164 correctly. So the "observed run's weights on the same scale" row is off by 1/f — the file
cannot be used to recentre the size-flavoured draws as documented.
(The companion identity the same paragraph states, "it is the rows, sum over c of n_c times w_c,
that sum to 1", is likewise an identity only over the unscreened frame.)

### HV-BOOT-04 (high) — the phantom cluster needs NO missing data: a singleton unit does it
Repro: `p2_size_degen.do` section D (log lines 700-810).
10 states x 2 units x 8 years, plus one extra state holding a single unit observed in a single
period, with **every variable present**. reghdfe drops that observation as a singleton.

| | with the singleton row | without it |
|---|---|---|
| e(N) | 160 | 160 (the row never entered either way) |
| **e(boot_df)** | **10** | **9** |

So the trigger is not "missing data" at all. Anything the estimator itself drops — a singleton
FE group, a row outside the event-study window (ado:9367 names both) — mints a bootstrap cluster
that contributes nothing. Single-observation units are ubiquitous in unbalanced administrative
panels, which moves HV-BOOT-01's realism from "plausible" to "routine".

### CLEAN NEGATIVE — `bootstrap(bayes)` under the DEFAULT `bootalpha(equal)` is immune to
### partial within-cluster missingness
`p2_size_degen.do` section B3: the same 50%-missing panel, full vs dropped, reps(400),
SE = .14522768 in BOTH runs — **bitwise identical**, ratio exactly 1. Under `equal` the
multiplier is per-cluster and the estimator is scale-invariant in the weights, so the dead rows
carry no mass and the live rows' normalized weights are unchanged. This is the predicted result
(normalized sub-vector of a Dirichlet is a Dirichlet), and it usefully bounds the blast radius:
the default bayes flavour loses only the df and the guards, not the SE.

### NEGATIVE (partial) — one downstream guard DOES count correctly
`p2_size_degen.do` section A: 6 clusters in the frame, 1 in the estimation sample. The bootstrap
screen passes it and prints the false caution "only 6 clusters are available for the bootstrap",
but a separate cluster-robust guard then catches it: "Error: cluster(state) leaves 1 cluster(s)
in the estimation sample; cluster-robust inference needs at least 2." So the fully degenerate
SE = 0 case is not reachable at G_live == 1. At G_live == 2 that guard passes — see p4.

---

### HV-PERM-01 (high) — the permutation DESIGN universe is the same unscreened sample:
### onsets are assigned to clusters that contribute nothing, and the design is misreported
Repro: `p3_perm_assert.do` section A (log `p3_perm_assert.log`). `permsave()` names one
`onset_<id>` column per cluster the design assigns to, so it is direct evidence.
20 states, 8 treated, 4 treated states report no outcome (16 live clusters, 4 live treated).

| | full data | dead rows dropped (oracle) |
|---|---|---|
| e(N) | 320 | 320 (identical) |
| point estimate | .65513021 | .65513021 (identical) |
| **permsave() assignment columns** | **20** (`onset_1` … `onset_20`) | **16** (`onset_5` … `onset_20`) |
| **reported possible assignments** | **634,888,800** | **43,680** |
| "only K of the resampling clusters are ever treated" caution | **NOT PRINTED** | printed ("only 4") |
| att p (200 reps) | 0.0100 | 0.0100 (both at the add-one floor) |

Mechanism: `wooldid_permrun` ado:7479-7488 `smp = ((ifr)==1) & !missing(clus)`, and
`wooldid_permctl` ado:5957-5959 builds the byte-identical `insmp`, which drives the
never-treated/extended-missing screen, the cluster-nesting screen, the permonsetcohorts
bijection screen, the "at least 2 clusters and 2 onset groups" freedom screen (ado:6153-6165)
and the treated-cluster caution (ado:6066).

HONEST CAVEAT on severity: unlike the bootstrap, a defensible reading exists for permute —
randomization inference conditions on the assignment mechanism, which physically operated on
the clusters with no outcome data too. What is NOT defensible either way is the reporting: the
user is told the design admits 634,888,800 assignments when the analysis supports 43,680, and
the caution that tells them how few treated resampling clusters their inference rests on is
suppressed by clusters that contribute no data. I did not measure size, so I do not claim the
p-value is invalid; I claim the design description and the guard are wrong.

### HV-VAL-01 (medium) — wooldid REFUSES runs because of rows it would itself have dropped
Repro: `p3_perm_assert.do` section D. Realistic admin coding: 20 unusable records carry the
sentinel id `unit == -9` and no outcome (`y == .`).

| | rc | e(N) |
|---|---|---|
| full data, no `if` | **411** "Error: unit must consist of positive, non-zero integers, but does not" | — |
| the same 20 rows physically dropped | 0 | 240 |
| full data with `if junk == 0` supplied by the user | 0 | 240 |

The analysis sample is IDENTICAL in all three: `ifvar` zeroes every `y >= .` row (ado:686), so
those 20 rows could never have entered. `wooldid_validate` builds `local ifassert if (`if') == 1`
(ado:4956) off the RAW user `if` — EMPTY when no `if` is given — so the asserts at ado:4960,
4986, 5247, 5291, 5320, 5341 and 5359 run over the whole dataset. The message names no remedy:
it does not say that restricting with `if` would make the run legal, and a user reading it will
conclude their id variable is broken rather than that 20 unusable rows are in the file.
Sentinel codes (-9, -99) on records with no outcome are everywhere in administrative and survey
extracts, so this fires without anything pathological in the analysis data.

---

### HV-BOOT-01 addendum — wooldid PRINTS BOTH COUNTS, contradicting itself on one screen
The estimation header reports the correct live cluster count (it comes from the estimation
command on e(sample)); the bootstrap footer reports the reference distribution from the
unscreened count. From `p1_guards.log`:
```
line 136:  N clusters = 16; N = 320 (80 obs dropped from initial sample)
line 178:  ... referred to the t(19) distribution
line 351:  N clusters = 8;  N = 160 (80 obs dropped from initial sample)
line 389:  ... referred to the t(11) distribution
```
and from `p4_df_bca.log` (12 clusters in the frame, 2 in the estimation sample):
```
N clusters = 2; N = 60 (300 obs dropped from initial sample)
... referred to the t(11) distribution      <- t(1) is correct here
```
wooldid already has the right number on the screen. That is why this reads as an oversight
rather than a convention: no design could want the same output to assert both.

### HV-BOOT-05 (high) — worst case measured: t(11) where t(1) is correct
Repro: `p4_df_bca.do` section A. 12 states, only state 6 (treated) and state 12 (never treated)
report an outcome.

| | full data | oracle |
|---|---|---|
| e(N) | 60 | 60 |
| header "N clusters" | 2 | 2 |
| **e(boot_df)** | **11** | **1** |
| "only 2 clusters available" caution | **NOT PRINTED** | printed |
| "only 1 treated cluster(s)" caution | **NOT PRINTED** | printed |

t(11) 97.5th percentile 2.201 vs t(1) 12.706: any bootdisplay(wald) interval, and anything a
downstream `test`/`lincom`/esttab/coefplot reads out of e(V) with e(df_r), is **5.8x too narrow**.
Both runs report se ~ 1e-15 and p = 0.000 (a pre-existing property of a two-cluster Bayesian
bootstrap, not caused by this defect) — but on the full data that result arrives with every
few-cluster warning silenced.

---

### HV-PERM-02 (high) — `permexact` reports "exact shares with no Monte-Carlo error" over a
### reference set it pruned, and the p moves 4x with a preprocessing step that cannot
### change the estimation sample
Repro: `p5_perm_refusals.do` section D (log `p5_perm_refusals.log` lines 660-800).
7 units, unit 1 treated 2003, unit 2 treated 2004, units 3-7 never treated; units 1 and 7
report no outcome at all. `cluster(unit) permute(shuffle) permexact`.

| | full data | the same rows physically dropped |
|---|---|---|
| e(N) | 60 | 60 (identical sample) |
| point estimate | .51739116 | .51739116 (identical) |
| reported reference set | **"40 of 42 replications completed; 42 possible assignments"** | "5 of 5 replications completed; 5 possible assignments" |
| **exact att p** | **0.1000** | **0.4000** |

THE LOAD-BEARING CLAIM. The full-data footer says the p-values are "exact shares with no
Monte-Carlo error" and reports "42 possible assignments" — and in the same block prints
"Elements not estimable in every enumerated assignment (exact shares over a reduced reference
set) … main_att". Two of the forty-two enumerated assignments put both onsets on the two
clusters that carry no outcome, so they cannot be estimated and were pruned. The reported p is
therefore conditioned on estimability, and that conditioning exists ONLY because the enumeration
ranged over clusters that contribute nothing. Under the 5-assignment design there is nothing to
prune. So "exact, no Monte-Carlo error" is false as printed, and the footer names a reference
set (42) that is not the one the p was computed over (40), nor the one the analysis supports (5).

THE DEFENSIBILITY CAVEAT, same as HV-PERM-01. 0.1000 and 0.4000 are each an exact share of a
DIFFERENT randomization model (7 clusters randomized vs 5 clusters analysable), and a
randomization-inference purist can argue the 7-cluster model is the physically correct one.
What survives that argument: (a) which model wooldid uses is decided by whether the analyst
deleted unusable rows before calling it — a step that provably cannot change e(sample) — and
nothing in the output says which one ran; (b) the pruning in (a) above breaks the exactness
claim under EITHER model.

### HV-PERM-03 (high) — `permute()` HARD-REFUSES (r459) on rows wooldid itself drops,
### while the ordinary run on the very same file succeeds
Repro: `p5_perm_refusals.do` sections A, B, C.
All three screens live in `wooldid_permctl` and evaluate on `insmp` (ado:5957-5959), the
unscreened concept. Every message says "in-sample observation(s)".

| probe | file | plain `wooldid` | `+ permute()` |
|---|---|---|---|
| A: ONE stale record for a unit in state 7 carrying a different onset, y missing | 400 rows | rc 0, e(N)=399 | **rc 459** "requires clusters to be subsets of treatment-onset groups … 1 in-sample observation(s)" |
| A oracle: that one row deleted | 399 rows | rc 0, e(N)=399 | rc 0 |
| B: the same row carrying `.a` in the onset instead | 400 rows | — | **rc 459** "requires never-treated observations coded as plain missing" |
| B oracle | 399 rows | — | rc 0 |
| C: `permonsetcohorts`, ONE dead row giving a cohort a stray onset | 300 rows | rc 0, e(N)=299 | **rc 459** "permonsetcohorts asserts that coh IS the treatment-onset grouping … 1 in-sample observation(s)" |
| C oracle: that row deleted | 299 rows | rc 0, e(N)=299 | rc 0 |

C hits `permonsetcohorts`, which is ONE DAY OLD (T1, 2026-09-16) and whose assessment suite has
not run — see known_bydesign §7. A single unusable record anywhere in the file is enough to make
the whole engine unavailable, and the message points the user at their cohort definition rather
than at the row.

---

### HV-BOOT-06 (high, NOT SEEDED — found by following the code's own comment)
### The T2/T3 nesting-regime decision is made on `smp`, so ONE dropped row flips the
### bootstrap regime and collapses att_i's composition variance
`wooldid_bootrun` ado:8975-8986 says, of the `if smp` on the nesting predicate:
> "the `if `smp'' on BOTH the gegen and the assert is load-bearing, not tidiness. `smp' is
> exactly the row set the leg keeps before it draws … Evaluated on the raw frame instead, a
> cluster that nests i within the sample but straddles two i's only through rows the sample
> excludes would be misclassified … and a strictly COARSER cluster … would silently take the
> keep-i rule and collapse att_i's composition variance"

The comment is right about what `bsample` sees — the leg does `qui keep if `smp'` (ado:9345) and
draws on exactly those rows. What is wrong is that those rows are the wrong sample for the
PREDICATE. The nested `wooldid` inside each replication then drops more (ado:9367 says so in as
many words: "the replication drops rows the frame still holds (missing y, an event-study window,
singletons)"), so the units the equal-cell weighting actually sees are e(sample)'s, not `smp`'s.
A unit that straddles two clusters only through rows no replication will ever use is therefore
misclassified as non-nesting, and the trap the comment describes is live. The remedy is to decide
the regime on the OBSERVED run's e(sample), not on the resample frame.

Repro: `p6_nestregime.do` (log `p6_nestregime.log`). 16 states x 3 units x 10 years,
`cluster(state)` (coarser than `i`=unit), estimand `att_i`, `bootstrap(standard)` reps(300).
ONE stray record codes unit 19 to state 8 for a single year with **y missing** — a crosswalk /
boundary-change artefact on a row with no outcome. Inside e(sample) unit 19 is wholly in
state 7, so the cluster nests i and `wdb_nest` should be 1.

| | full data | that one row deleted |
|---|---|---|
| e(N) | 480 | 480 (identical) |
| att_i estimate | .21866955 | .21866955 (identical) |
| **att_i bootstrap se** | **.10710045** | **.11880241** |
| t | 2.04 | 1.84 |
| **p** | **0.040** | **0.093** |
| 95% CI | [.0045449, .41846671] | [-.01521439, .43350729] |

**The significance verdict flips.** One row that contributes nothing to the estimate moves
att_i from "not significant, CI covers zero" to "significant at 5%, CI excludes zero".
This is one-day-old code (T2/T3, 2026-09-16) whose assessment suite was deferred and has not run
(known_bydesign §7). Seed robustness in `p7_nest_seeds.do`.

### CLEAN NEGATIVE — `bootalpha(equal)` with WHOLLY DEAD clusters: SE unbiased
`p4_df_bca.do` section B. 18 states, 5 never-treated states report no outcome; reps(600).
SE full .11124138 vs oracle .10943693, ratio **1.0165** against MC noise of 2.89% — 0.6 sigma.
Confirms the theory: the normalized sub-vector of a Dirichlet is a Dirichlet, and wooldid's
estimand is scale-invariant in the weights, so wasted mass is harmless. Under the default bayes
flavour the loss is confined to df, guards and the reported design.

### NEGATIVE — the event-study window is NOT a row-dropping trigger here
`p4_df_bca.do` section D: `esprelength(2) espostlength(2)` left e(N) at 360, unchanged from the
unwindowed run. So of the three triggers ado:9367 names (missing y, an event-study window,
singletons), I confirmed missing-y and singletons and could NOT reproduce the window one.

### UNCERTAIN — `bootdisplay(bca)`: the jackknife runs one fit per FRAME cluster
`p4_df_bca.do` section C. 18 states, 5 dead. The footer says "The BCa acceleration comes from a
delete-one-cluster jackknife (**18 fits**)" when only 13 clusters can influence anything; the
loop is `forval g = 1/wdb_C` (ado:9493). Deleting a dead cluster reproduces the full estimate, so
five of the eighteen jackknife points are exact duplicates of the full fit. Because the standard
acceleration is built from deviations about the MEAN of the leave-one-out estimates, those five
duplicates shift the mean and therefore change `a`. Measured: full att p 0.007, CI
[.0764, .6446]; oracle p 0.021, CI [.0728, .5550]. Those runs also have different draw streams,
so the interval difference is confounded and I do not claim the size of the effect — only that
the acceleration is computed over the wrong set, which is structural.


---

## CAVEAT ADDED ON REVIEW — HV-BOOT-02 and the bayesboot parity claim
`bayesboot` is NOT installed on this machine, so I could not check how it builds its own n_c.
sthlp:548 sells `bootalpha(size)` as reproducing "what bayesboot does when passed
id(clustervar)". If bayesboot also generates weights over every row it is handed and lets the
estimation command drop the missing ones afterwards, then wooldid MATCHES it and the right
reading of HV-BOOT-02 is "the parity is preserved, and the understatement on an unbalanced panel
is worse than the help documents, by an undocumented sqrt(live fraction)". If bayesboot marks out
first, wooldid diverges from the command it claims parity with as well. Either way the measured
number stands: the reported SE is 32% smaller than a bootstrap of the estimation sample gives,
and HV-BOOT-03 (the bootsave `_rep=0` row) is a flat contradiction of sthlp:557 under both
readings. I mark HV-BOOT-02 `probable`, not `certain`, for this reason.

---

## HV-BOOT-06 CONFIRMED SYSTEMATIC — `p7_nest_seeds.do` / `.log`
Same design, four independent bootstrap seeds, reps(400). e(N) = 480 on BOTH sides every time.

| rseed | se, full data | se, stray row deleted | ratio |
|---|---|---|---|
| 101 | .10183318 | .12516709 | **.8135** |
| 202 | .09709323 | .11168878 | **.8693** |
| 303 | .09677665 | .10892175 | **.8884** |
| 404 | .10779581 | .11932307 | **.9033** |
| (p6, reps 300) | .10710045 | .11880241 | .9015 |

Five for five below 1, mean ratio .875. Both sides of each pair run the same seed over the same
16 clusters, so the sampling error is largely common and these ratios are near-deterministic
rather than five independent draws.

**The discriminating control settles the mechanism.** Give that same stray row a REAL outcome
(`replace y = 1.5 if stray`): the unit then genuinely straddles states 7 and 8 inside e(sample),
so the keep-i regime is the CORRECT one. Result: e(N) = 481, se = **.10190961** — indistinguishable
from the full-data run's .10183318 and nowhere near the oracle's .12516709. So the full-data run
is demonstrably executing the keep-i regime, which is the wrong one for its own estimation sample.
That rules out "different draw stream" and "MC noise" as explanations.
