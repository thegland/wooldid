# perm-invariance — FINDINGS
Agent: metamorphic testing of permute() invariance at fixed rseed.
Build: /private/tmp/.../bluesky/build (sha c1d8979).  All probes generate their own data.

---

## PI-1 (FOUND BY ME, not seeded) — a cluster with no usable outcome data stays in the permutation design, so rows the estimator discards change the p-value

**Surface:** `permute()` design population (`wooldid_permrun`, ado:7481-7490; `wooldid_permctl`,
ado:5959-5964) — the four-sample bug axis.

**Claim.** The permutation design sample is `((if)==1) & !missing(cluster)`. It does not require a
row to be usable for estimation. A cluster all of whose rows have a missing outcome therefore
enters the map frame, is handed an onset in the shuffle, and enlarges the reference set — while
contributing nothing to any fit. Two datasets with a **bit-identical estimation sample and
bit-identical point estimates** consequently report different permutation p-values and different
confidence intervals, with no caution. Excluding the same cluster by a missing CLUSTER value, or by
the user's `if`, DOES remove it from the design, so the three exclusion routes disagree.

**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/repro_phantomcluster.do`
Also probeH.do variants H0/H1ymiss/H2ydrop/H3cmiss/H4partial/H4bdrop and probeE.do's J0/J1/J2.

**Observed** (40 counties in 10 states, state 5 a 2006 adopter, 30 replications, rseed fixed):
```
                         y missing on state 5   state 5 rows dropped
 e(N) / e(N_clust)           360 / 9                360 / 9
 att point estimate          1.0659403              1.0659403      (bit-identical)
 e(perm_npossible)          16800                   5040
 p(pretreatment att)         0.645161               0.451613
 p(ES_att pre2)              0.580645               0.645161
 95% interval for att      [0.6364, 2.0193]       [0.2934, 1.9483]
```
Setting state 5's CLUSTER value to missing reproduces the DROPPED column exactly (probeH:
`dump_H3cmiss` vs `dump_H2ydrop` differ only in `e(N_dropped)`). Partial missingness inside a
cluster that survives is neutral (`H4partial` vs `H4bdrop`: identical).

**Sharper instance (probeE2 J0/J1/J2).** Make the phantom the *sole* carrier of an onset level
(state 7, the only 2008 adopter, outcome missing everywhere). The estimation sample then contains
no 2008 cohort at all, but the onset 2008 stays in the shuffle pool: **in 29 of 30 replications a
data-carrying cluster received the 2008 onset**, so the reference distribution is built largely
from designs with a cohort structure the observed data cannot produce. `e(perm_npossible)` is
16800 vs 1680 -- a factor of ten. Nothing failed: `e(perm_failed)=0`, `e(perm_blank)=0`,
identical column sets, so the mechanism is not broken replications, it is a different design.
p(pretreatment att) 0.710 vs 0.516; att interval [0.632, 1.921] vs [0.503, 2.005]; point estimates
bit-identical. `if state != 7` reproduces the DROPPED column exactly (J1 vs J2 equal to 1e-10).

**Expected.** Rows that cannot enter any fit should not change an inference surface. The
randomization distribution wooldid reports is supposed to be the distribution of the *realized*
statistic under reassignment. The count in `e(perm_npossible)` is the correct count for the design
wooldid actually uses -- the defect is not that the number is wrong, it is that the design silently
differs from the estimable one, and that the three routes to removing a cluster disagree: a missing
cluster value removes it, the user's `if` removes it, a missing outcome does not. wooldid knows the
rows left (`e(N_dropped)=40` on the kept-rows run) and says nothing. The help's only statement
about the design population (`sthlp:60`, BS-37) is that the user's `if` applies identically to the
estimation sample and the permute design; nothing says missing outcomes are treated differently
from a missing cluster id.

**Scope.** `wooldid_bootrun` builds its cluster map from the *same* predicate
(`ado:8912-8923`: `((ifr)==1)` then `replace smp = 0 if missing(clusv)`), so `bootstrap()` shares
the design population. That makes this a single shared choice rather than a permute-only slip, and
argues for a caution + documentation ruling rather than a numeric fix. (The repro also prints a
bootstrap `e(V)[1,1]` pair, 0.01910091 vs 0.01501891 — that is NOT a measurement of the effect:
at 30 replications the Monte-Carlo error of a bootstrap variance is of the same order, and the two
runs draw different resample sequences. The shared predicate is the scope evidence; the numbers
are only consistent with it.)

**Realistic:** yes. Units with no outcome data at all are routine in administrative and survey
panels (a state that does not report, a county below a disclosure threshold). The analyst's
assumption — "the estimator drops them, so they cannot matter" — is exactly the one that fails.

**Severity:** medium. Under the sharp null the reported test is still an exact randomization test
*for the design wooldid used*; what is defective is (a) the silent sensitivity of the p-value and
interval to rows the estimator discards, (b) the inconsistency between the three exclusion routes,
and (c) the absence of any caution when a whole cluster leaves the estimation sample but stays in
the randomization pool. It is medium rather than high because `bootstrap()` uses the identical
predicate, so this is one shared definition of the design population rather than a permute-only
slip. The argument for re-ranking, if the author wants it, is that an edit a user believes is a
no-op moves a reported p from 0.516 to 0.323 with nothing in the output to notice.
**Confidence:** certain on the measurement; probable on the ruling.

**Exact repro output** (from `repro_phantomcluster.log`):
```
estimation sample   A: N=360 G=9   B: N=360 G=9
att point estimate  A: +1.82acfb38ea209X-001   B: +1.82acfb38ea209X-001   identical = 1
perm_npossible      A: 16800   B: 5040   C: 5040   D: 5040
p(pretreatment att) A: 0.516129   B: 0.322581
p(ES_att pre1)      A: 0.580645   B: 0.451613
95% interval (att)  A: [0.4452, 1.4635]     B: [0.2666, 1.4431]
bootstrap e(V)[1,1] A: 0.01910091            B: 0.01501891   <- MC noise at B=30; not evidence
```
(A = state 5's outcome set to missing, rows kept; B = state 5's rows dropped;
 C = state 5's CLUSTER value set to missing; D = `if state != 5`.)

---

## NEGATIVE RESULTS — what I attacked and what came back clean

Every item below was measured on data I generated, with the full e() surface dumped in
exact hex (`%21x`) and compared element by element. "IDENTICAL" means rtol 0 unless noted.

### Tier 1 — the permsave() assignment matrix (probeA, probeB, probeD, probeE)
Invariant to: running the identical command twice in one session; running it in a **fresh
session** (probeA's `base` vs probeB's `B0`, rtol 0); a pseudo-random permutation of the
dataset's row order; a full reversal of the row order; reordering the option list;
relabelling `i` order-preservingly, non-contiguously (`unit*7+3`) and order-reversingly
(`1000-unit`); relabelling the cluster ids order-preservingly (`state*13+5`); duplicating
every row; every transformation of `y` and of the weights; shifting `(t, ttre)`;
`permstrata(a constant)`; changing the estimand list; changing the event-study window.
It moves ONLY under an order-REVERSING cluster relabel, which is the one case where it
legitimately must (the generation order is the ascending (stratum, cluster) order).
The **permreps prefix property** holds: the first 40 assignment rows of a `permreps(80)`
run are bit-identical to a `permreps(40)` run at the same seed.

### Tier 1 correctness, not just stability (probeD)
For replications 1, 7, 23 and 40 I read the assignment out of `ps_D0.dta`'s `onset_*`
columns, wrote it back onto the panel as a new onset variable, and re-estimated with the
engine's own replication option string. The stored replication coefficients came back
exact or 1 ulp off. The permsave assignment block is correctly labelled by cluster.

### Tier 3 — reported p-values, intervals, bands (probeA-probeG)
- Bit-identical under every tier-1-invariant transformation above.
- Under `permexact` the null distribution is a multiset invariant, so p must not move even
  when the assignment matrix does. Measured: reversing the cluster labels changes ONLY the
  ROW ORDER of `e(permdist)`/`e(permdist_b)`; every p, interval, band and kernel entry is
  bit-identical. Same on an **interleaved-strata** design (strata {1,3,5,7} and {2,4,6,8},
  so ascending cluster order crosses stratum order — the NFPERM-02 geometry): 36 = 6x6
  assignments, and reversing the stratum labels or the cluster labels leaves everything
  bit-identical. The shipped NFPERM-02 fix holds from the outside.
- `y -> a + b*y`: with b>0 every p and t-ratio bit-identical and sd/lb/ub/band scale by
  exactly b (rel 1e-16); with b<0 every reported p bit-identical and the band mirrors
  (`bandlo(-y) = -bandhi(y)`), i.e. the equal-tailed convention is sign-symmetric. Checked
  under the default, `permabsp`, `permside`, `permstepdown` and `permtstat`.
  `permside(left)` on y equals `permside(right)` on -y. `e(perm_p_lo)`/`e(perm_p_hi)`
  correctly swap. Step-down adjusted p's identical. (Sanity: `permside(left)` and
  `permabsp` DO move the p relative to the default, so these are not vacuous passes.)
- `permtstat` across nine orders of magnitude of outcome scale (y, 1000y, 1e-6 y): every
  p bit-identical, t differing by 2e-14. No absolute tolerance leaks into the studentized
  path at realistic scales.
- Duplicating every row: every coefficient, p, interval and band identical; only N,
  `escomp_*`, ll, rss, rmse and r2 move.
- `if` vs physically dropping the same rows: identical (BS-37's claim holds for the
  permute design, and again in probeE2's J1/J2).

### Independent oracles (python, computed from the dumped matrices)
- The reported p reproduces `2*min((#ge+1)/(n+1), (#le+1)/(n+1))` exactly on all 14
  coefficients of a sampled run, and the PLAIN share `2*min(#ge,#le)/n` on a `permexact`
  run — i.e. the observed assignment sits inside its own enumerated reference set and is
  not double counted.
- The interval is `[b_obs - qhi, b_obs - qlo]`, the documented shift interval.
- `permute(probabilities)`: a table with eight deterministic rows aimed at onsets the
  clusters do NOT have, written ascending, written descending, and with the cluster ids
  relabelled `100-state` — in all three every cluster drew exactly the onset its own row
  specifies across all 24 replications. The permprobs/onset pairing is value-keyed.

### RNG
- The documented replay route (`set rngstate e(perm_rngstate)`, then repeat the call) is
  bitwise exact on a run with no `rseed()`.
- A user `set seed` before a no-`rseed()` run is honoured and reproducible.
- With `rseed()`, the post-run stream is determined by `rseed()` alone (a subsequent
  `runiform()` is identical whatever the user's own pre-run seed was) — the standard
  `set seed`-inside-a-command behaviour the help describes.

### permonsetcohorts (shipped 2026-09-16, one day old)
Reversing the `i` codes, relabelling the clusters and shuffling the rows all leave every
number equal to 1e-10 with `e(perm_irebuilt)=1`. The documented "no-op when `cluster()`
equals `i`" claim is exactly true: the two runs differ ONLY by the presence of the scalar
`e(perm_irebuilt)=0`.

### subgroups() label swap (probeG)
Relabelling the subgroup values 1<->2 swaps the per-subgroup tables EXACTLY (`G0`'s
`e(wdidsgresults_1)` equals `G1`'s `e(wdidsgresults_2)` bit for bit, and conversely) and
changes nothing else; the pairwise contrast flips sign, mirrors its interval and keeps its
p bit-identical. An order-preserving relabel to {5,9} renames the keys and leaves the
contents bit-identical. No surface crossed the two subgroups. (Note the dual convention —
`e(b)` names them positionally, `main_SG1_*`/`main_SG2_*`, while the `e()` tables are keyed
by VALUE, `_1`/`_2` or `_5`/`_9`. Both map consistently through sorted order.)

### Other
- Never-treated coded as an out-of-window future onset (2020 on a panel ending 2009) gives
  results bit-identical to plain missing, permutation design included.
- `simultaneousbands` under `permute()` is refused with an informative r(198) naming
  `permstepdown` as the permutation multiplicity tool. Not a defect.
- Rescaling `aweights` by 1000 leaves every reported number equal within 1e-10; the
  replication accumulator moves 4.6e-10 at one cell, which is reghdfe's `tol(1e-9)`
  convergence floor. `e(escomp_refW)`/`e(escomp_trtW)` report the raw weight sum and so
  scale by 1000 — correct reporting, noted only so it is not mistaken for drift.
