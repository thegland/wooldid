# SKEPTIC report — refuting PTP-1 .. PTP-5 (pretrendpower)

Everything below was rebuilt from the claims' own terms on MY OWN data generator
(16 states, 4 onset cohorts of 4, onsets t=6/8/10 + never-treated, T=12, state random
intercepts, quadratic+cosine common year path, AR(1) rho=0.5) and MY OWN seed streams.
I read the prober's do-files but ran none of them.

Files: `s1_orient.do/.log`, `s2_ideal.py/.log`, `s3_mc_null.do/.log/.csv`,
`s4_mc_power.do/.log`, `s5_ex7_ati.do/.log`, `s6_nis_from_null.py/.log`,
`s7_hotelling.do/.log`.

---

## HEADLINE

The prober's NUMBERS replicate. The prober's MECHANISM does not.

Every symptom PTP-1 attributes to "inverting a 15-df cluster-robust S" is reproduced,
to within a point or two, by a textbook model with **no panel, no clustering and no
DiD at all**: `b ~ N(0, Sigma)` with `S ~ Wishart(Sigma, df_r)/df_r` drawn
independently of `b`. That is the *best* case for any covariance estimator.

| statistic (k=3 pre periods, df_r=15) | ideal model | my wooldid MC | prober |
|---|---|---|---|
| size of the joint pre-trend test | 0.0784 exact / 0.0786 sim | 0.1120, 0.1280 (2 streams, 250 each) | 0.1050 |
| GLS slope 95% CI coverage | 0.9110 | 0.9120 | 0.9050 |
| reported slope se / MC sd | 0.846 | 0.775 | 0.79 |
| median mdt_hat / mdt_true | 0.9070 | 0.8919, 0.8667 | ~0.83 |
| CV of the reported mdt_wald(80) | 0.198 | 0.227 | 0.23 |

The cause is that `wooldid_ptp` applies **known-covariance formulas to an estimated
covariance**; the damage is a deterministic function of `(k, df_r)` that is present at
EVERY df_r and has nothing to do with how the covariance was clustered. The prober's
"discriminating control" at G=80 is not exact either — the ideal model puts it at size
0.0546 and coverage 0.9432; 220 replications simply could not see that.

---
## SKP-NEW-1 — `pretrendpower()` uses KNOWN-covariance references throughout, so its joint test, its slope CI and its `mdt_wald` column are off by a deterministic function of `(k, df_r)` — nothing to do with clustering
severity **high**, confidence **certain**, realistic **yes**
(this is PTP-1's defect correctly diagnosed and extended; the two are one item, not two)
repro: `s7_hotelling.do` (also `s2_ideal.py`; the panel legs are `s3_mc_null.do`, `s4_mc_power.do`)

`wooldid_ptp` forms `W = b'S^-1 b` on the ESTIMATED pre-period covariance `S = e(V_es_*)`
and then refers it to `F(keff, df_r)` (`ado:10516`), takes `invF(keff,df_r,.95)` as the
cutoff that `mdt_wald` powers (`ado:10559`), and builds the slope CI from
`se = 1/sqrt(m'S^-1 m)` with `invttail(df_r,.025)` (`ado:10512-10514`). All three are
exact only when `S` is the TRUE covariance. With `S` estimated at `df_r` degrees of
freedom the exact reference for `b'S^-1 b` is Hotelling's: `W*(df_r-k+1)/(df_r*k) ~
F(k, df_r-k+1)`.

Exact size of wooldid's own cutoff, a test the run labels 5% (`s7_hotelling.do` PART A,
closed form; PART B confirms each by 100,000 Mata draws of `b ~ N(0,I)`,
`S ~ Wishart(I,df_r)/df_r`):

| k \ df_r | 15 | 23 | 39 | 79 | 279 |
|---|---|---|---|---|---|
| 1 | 0.0500 | 0.0500 | 0.0500 | 0.0500 | 0.0500 |
| 2 | 0.0611 | 0.0570 | 0.0540 | 0.0519 | 0.0505 |
| 3 | **0.0784** | 0.0673 | 0.0597 | 0.0546 | 0.0513 |
| 4 | **0.1034** | 0.0813 | 0.0670 | 0.0579 | 0.0521 |

and the companion damage at k=3 (100,000 ideal draws):

| df_r | size | slope CI coverage | reported se / MC sd | median mdt_hat / mdt_true |
|---|---|---|---|---|
| 15 | 0.0770 | 0.9108 | 0.8448 | 0.9083 |
| 79 | 0.0544 | 0.9429 | 0.9706 | 0.9822 |
| 279 | 0.0498 | 0.9488 | 0.9935 | 0.9955 |
| 15, k=4 | 0.1033 | 0.8831 | 0.7776 | 0.8686 |

`mdt_hat/mdt_true` has the exact law `sqrt(chi2(df_r-k+1)/df_r)`: at k=3, df_r=15 its mean
is 0.9132, its CV 0.198 and its 5-95% range [0.627, 1.221]. `sthlp:789` says "the minimum
detectable trend itself depends only on the design's precision, never on the estimated
coefficients", which reads as a design constant; it is an estimate with a 2:1 spread here.

A 16-state DiD with `i` = onset cohort, `fe(st)`, `cluster(st)` and `esprelength(3)` lands
exactly on the k=3, df_r=15 row (`s7_hotelling.do` PART C: `e(df_r)`=15, `e(roth_keff_att)`=3,
`e(N_clust)`=16). Ask for `esprelength(4)` on the same panel and the "5%" joint pre-trend
test has true size 0.1034 and the slope CI covers 0.883 — with no change in the number of
clusters at all.

Two qualifications I want on the record. (i) `F(q, df_r)` is Stata's own convention for
`test`/`testparm` after any clustered fit, so wooldid did not invent the liberal reference
for the joint test — my hand-rolled comparator (a plain pooled event study, `reghdfe ...,
absorb(st t) cluster(st)`, `testparm` on 3 pre-period dummies, same 250 datasets) rejects
at 0.0680, consistent with the 0.0784 ideal-model value. (ii) What is wooldid's own is that
`pretrendpower()` turns that same known-Sigma arithmetic into a REPORTED DESIGN QUANTITY —
`mdt_wald`, `bias_wald`, `roth_slope_ll/ul`, `roth_nissize_*` — where the error is
systematic and one-directional rather than one test's slightly-off p-value.

Fix shape: refer `W` to `F(keff, df_r-keff+1)` after the Hotelling scaling, use the matching
noncentral reference for `mdt_wald`, and say in the help that the MDT is an estimate.

---
## SKP-NEW-2 — under `bootstrap()`, `pretrendpower()`'s reference distribution is keyed only on `bootdf()` and never reflects `reps()`; under the default `bootdf(normal)` it is `chi2(k)` with no finite-sample correction at all
severity **low**, confidence **uncertain**, realistic **yes**
repro: `s8_boot.do`

Same 16-state panel, `bootstrap(standard)` at `reps(39)` and `reps(199)`:

| run | `e(df_r)` | `e(boot_df)` | `e(boot_N_V)` | joint pre-trend line |
|---|---|---|---|---|
| analytic | 15 | . | . | `F(3, 15)` |
| boot, default `bootdf` | 15 | 15 | 39 | `F(3, 15)` |
| boot, `bootdf(normal)` | . | . | 39 | `chi2(3)` |
| boot, default `bootdf` | 15 | 15 | 197 | `F(3, 15)` |
| boot, `bootdf(normal)` | . | . | 197 | `chi2(3)` |

`e(roth_W_att)` and `e(roth_slope_se_att)` are bit-identical between the two `bootdf`
settings at a given `reps()`, so only the reference moves — and it moves on `bootdf()`,
never on how many replications the covariance was built from. `bootdf(normal)` is the
DEFAULT whenever `cluster()` is absent (`sthlp`, Bootstrap inference entry), so a plain
`bootstrap(standard) reps(39)` run with `robust` gets the `chi2(3)` branch — the same
no-correction branch PTP-3 flags for poisson — on a covariance estimated from 39 draws.

HONEST LIMIT: I did not establish what the correct reference IS here. Treating the
bootstrap covariance as Wishart with `reps-1` df is the wrong idealization (resampling
more does not add information about Sigma beyond the G clusters), so I cannot convert this
into a size number the way SKP-NEW-1 does. Filed as a structural observation, uncertain.

---

# VERDICTS

## PTP-1 — **PARTIAL** (severity **high**, unchanged — see the note below)

**What survives.** The headline replicates, and on my data it is FIRMER than the prober
established it. `s4_mc_power.do`, own seeds, 250 null + 400 power replications: the median
reported `mdt_wald(80)` is 0.54610; a trend of exactly that size is rejected by wooldid's
OWN nominal cutoff in 294/400 = **0.7350** (se 0.0221, 95% CI [0.690, 0.777]), against a
documented 0.80. The prober's own CI [.7348,.8099] included 0.80; mine excludes it.
The design's true 80% point (wooldid's own lambda on the MC covariance of the three pre
coefficients) is 0.63007, so the reported column is 13% optimistic; at 0.63007 the nominal
test has power 0.8880, above target, because that test's size is 0.128.
Also replicated: no caution fires (0 "Caution" lines, `s1_orient.log` section A); size
0.1120 and 0.1280 on two independent 250-rep streams; slope CI coverage 0.9120;
reported se / MC sd 0.775; CV of the reported MDT 0.227.

**What is wrong: the mechanism, and with it the scope and the remedy.**
1. It is not a clustering phenomenon. `s7_hotelling.do` PART B reproduces the whole
   symptom set with no panel, no cluster and no DiD: `b ~ N(0,I)`, `S ~ Wishart(I,15)/15`
   independent — size 0.0770, slope CI coverage 0.9108 (mine on the panel: 0.9120),
   reported se / MC sd 0.8448, median mdt ratio 0.9083 (mine: 0.8919, 0.8667).
   Clustering contributes the difference between 0.077 and 0.112, not the effect.
2. "the identical geometry at G = 80 is exact" is false. The ideal model puts G=80
   (df_r=79, k=3) at size 0.0546 and coverage 0.9429 — their 220 replications could not
   resolve that, so their control does not discriminate.
3. The defect scales with k, which they never varied. At k=4 and the SAME 16 clusters the
   exact size is 0.1034 and the slope CI covers 0.883. A few-cluster caution — their
   implied remedy — would not address that, and would fire on designs (k=1, any df_r)
   where the test is exact.
4. Their `mdt_nis` leg is misattributed. `mdt_nis` is optimistic because of PTP-2's own
   normal-quadrature-at-a-t-critical-value, not because of `S`-inversion: at k=1, where
   the NIS test and the joint Wald test are the SAME test, wooldid prints `mdt_wald`
   1.09139 and `mdt_nis` 1.08237 on my panel (`s1_orient.log` section B).
5. "the MDT and the bias column it drives are ~17% optimistic" is right in size but the
   split matters: roughly 9 points of it is the exact sampling law
   `mdt_hat/mdt_true ~ sqrt(chi2(df_r-k+1)/df_r)` (median 0.9070 at k=3, df_r=15), and the
   remainder is the extra CRVE noise of a cohort-level contrast identified off 4 states
   per cohort.

**Correction, in one sentence.** `pretrendpower()` is optimistic and its tests are
oversized because it applies known-covariance formulas to an estimated covariance; the
damage is an exactly computable function of `(k, df_r)` present at every cluster count,
not a few-cluster cluster-robust failure that 80 clusters would cure.

**On severity: I am NOT downgrading this.** PARTIAL is about the prober's characterization,
not the defect. The defect is real, it is BROADER than claimed (every `df_r`, worse as `k`
grows, not curable by more clusters), and my own replication pins the power shortfall harder
than theirs did — 0.7350 with CI [0.690, 0.777], which excludes the documented 0.80, where
their CI did not. "Textbook estimated-covariance theory" names the cause; it does not make
the reported number less wrong. Severity stays **high**; what changes is that a few-cluster
caution is the wrong remedy and a Hotelling scaling is the right one.

## PTP-2 — **CONFIRMED** (severity **medium**, unchanged)

Reproduced exactly on my own panel (`s1_orient.log` section B):
`e(roth_nissize_att)` = 0.033052124846 and `2*normal(-invttail(15,.025))` = 0.033052124846,
difference 0.000000000000; `2*ttail(15,invttail(15,.025))` = 0.050000000000 exactly.
`ado:10369` takes a t critical value, `ado:10542` feeds it to `mvnormalcvqp()`.

Three corrections/extensions the prober did not have:
- **The right reference is a multivariate t.** On my 250 null replications at k=3,
  df_r=15 (`s6_nis_from_null.py`), the realised NIS rejection rate is **0.1000**
  CP95 [0.0658, 0.1441]; the same rectangle under `mvt(15)` averages **0.1094** (inside
  that interval); what wooldid reports averages **0.0769** (a 30% understatement, at the
  interval's edge). So the fix is the reference distribution, not the critical value.
- **The `mdt_nis` COLUMN inherits it, not just the size scalar.** At k=1 the two
  detectable-trend columns must be identical by construction — the joint F(1,df) test and
  "any coefficient individually significant" are the same test — and wooldid prints
  1.09139 vs 1.08237, `mdt_nis` being the optimistic one. Pure arithmetic, no simulation,
  and it is the cleanest single demonstration of the bug.
- **It is confined to the `ptpdfok` branch.** On the no-df branch `ptpniscrit` is
  `invnormal` and the quadrature is then self-consistent: my poisson run reports 0.1331
  for a genuinely normal-referenced test. So the poisson user gets the correct NIS size
  and the gaussian user the wrong one.

## PTP-3 — **CONFIRMED**, severity medium -> **low**

Mechanically exact on my own panel (`s1_orient.log` section C): poisson posts
`e(df_r) = .` while `e(N_clust) = 16`, the run prints `chi2(3)`, and `mdt_wald` is
0.202212 against an F-branch 0.230673 — ratio 0.877, the prober's 12.3% to three digits.

But the branch is INHERITED, not chosen. On the same data `ppmlhdfe ycnt i.t D,
absorb(st) cluster(st)` posts `e(df_r) = .`, and so does Stata's own
`poisson ycnt i.t D i.st, cluster(st)`. The entire poisson family reports z and chi2;
wooldid's `chi2(3)` line is internally consistent with the rest of that run and is
visible on screen. What is left is one sentence: `sthlp:789`'s unconditional "an F test
with the regression's residual degrees of freedom", which is false on that path. That is
a documentation-vs-code divergence, not a wrong number, so **low**.
Counterpoint to their framing: on that same branch the NIS size becomes CORRECT
(see PTP-2), so the chi2 branch is not uniformly the worse one.

## PTP-4 — **CONFIRMED**, severity **low**

Reproduced on my own 60-unit panel with `i` = the unit id, which is the natural reading of
`wooldid y i t ttre` (`s5_ex7_ati.log` PART1): help example 7 run verbatim prints
"...the joint test's covariance is degenerate and its p-value should not be read as a
calibrated test" and then computes the diagnostics anyway (`e(roth_ok)`=1, `e(df_r)`=59,
`e(N_clust)`=60, `e(roth_prep_att)`=0.5513).

Two corrections. (i) With a cohort-level `i` (4 values) it is worse than the prober says —
the pre-block is rank-deficient and `pretrendpower()` is skipped entirely with
"no event-study coefficient vector with a valid variance-covariance matrix was available"
(`s1_orient.log` section E) — but that is my 4-cluster artifact, not the example's fault,
and I do not file it. (ii) It borders on already-known: `ado:1436`'s own comment records
that "wooldid's own help twice recommends the clustering that trips this condition", and
PS-64 was decided BECAUSE of that. The genuinely new part is only that the worked example
itself was never revisited after the ruling.

## PTP-5 — **PARTIAL**, severity medium -> **low**

**The magnitude does not replicate.** My own at-i panel (60 units, `i` = unit,
`cluster(i)`, AR(1) rho 0.5, T=12, 150 null replications, `s5_ex7_ati.do` PART2) gives
joint pre-trend size **0.1000**, slope CI coverage 0.8933, reported se / MC sd 0.8121 —
not 0.28. So 0.28 is a property of their particular panel, not of the at-i geometry.

**And what remains is a ruled item.** `recon/known_bydesign.md` §1: "Joint pre-trend /
`adversarial` tests are oversized at or above `i` — 11-36% rejection at nominal 5% (G=40)
... Computed anyway with a named caution, in BOTH `jointtests` and `pretrendpower()`
(PS-64/BS-9 as amended, 2026-09-07)." My 0.1000 at G=60 sits at the bottom of that band.
The prober concedes the caution fires on this design. Of my measured 0.100, the ideal
model of SKP-NEW-1 already accounts for 0.0562 at df_r=59, so the at-i degeneracy
contributes less than half.

The one piece that is neither ruled nor rediscovery: the caution's text names the joint
test's p-value only, while `mdt_wald`, `bias_wald` and `roth_slope_ll/ul` are built on the
same covariance and are printed without qualification. That is a low-severity wording gap,
and it is the same gap SKP-NEW-1 covers more generally.
