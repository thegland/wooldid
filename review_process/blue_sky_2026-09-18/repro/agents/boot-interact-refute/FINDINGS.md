# boot-interact-refute — SKEPTIC verdicts on BOOTX-1..5, plus own findings

Build: `.../scratchpad/bluesky/build/wooldid.ado` ("version 1.4 2026-08-19"), asserted by
`which wooldid` in r00_setup.log:54. Every number below is from MY OWN generator and MY OWN
do-files in this directory; I read the prober's scripts but ran none of them.

## Panel used by most probes ("design A")
30 states x 2 counties = 60 counties, 16 periods, 10 treated states with onsets 7..12,
y = 2 + .05*t + unit effect + treatment effect + N(0,1). **No pre-trend in the DGP.**
i = county, cluster(state), G = 30, esfixedbaseperiod esprelength(6) espostlength(4).

---

## BOOTX-1 — PARTIAL (severity high -> MEDIUM)

**Mechanism: CONFIRMED exactly.** `r01_mech.do` hand-forms W = b'S^-1b on the pre block of
e(b_es_att)/e(V_es_att) and shows `Ftail(keff, e(df_r), W/keff)` equals the printed /
e(roth_prep_att) value to every digit at reps 20/50/100/300 and under bootdf(clusters),
bootdf(normal) (chi2) and bootdf(3). e(df_r) = 29 = G-1 at every reps(). `r11_boottstat.do`
shows a hand `test` over the six pre coefficients returns the SAME number
(analytic .000676436559 == e(roth_prep_att); boot reps(60) .4708198565 == e(roth_prep_att)).
The claimed DIRECTION is also right: the finite-draw noise in S inflates W and the excess
grows without bound as n falls toward k.

**Three corrections, none of which kills it.**

1. *Magnitude at replication counts a user actually runs.* The prober never measured size.
   `r03_oracle.py` idealises b ~ N(0,I) with S ~ Wishart(I, n-1)/(n-1) — exactly the divisor
   and centring `wooldid_bootves` uses (ado:14622-14625) — and asks what the printed p does.
   The oracle isolates the FINITE-DRAW excess only; it cannot adjudicate whether F(k, G-1) is
   the right reference at large n, because F(k, G-1) is the ordinary small-cluster convention
   for a cluster-robust Wald and the oracle assumes exact normality, under which no cluster
   correction is needed. So the honest baseline is the same reference's own Sigma-KNOWN limit
   (r03_limits.out), not 5%:

   | k  | n (draws) | df_r | oracle size (Sigma-known idealisation) | same reference, S known | excess | Hotelling |
   |----|-----------|------|------|------|------|------|
   | 6  | 20        | 29   | .175 | .0237 | 7.4x | .052 |
   | 6  | 50        | 29   | .063 | .0237 | 2.7x | .050 |
   | 6  | 100       | 29   | .040 | .0237 | 1.7x | .050 |
   | 6  | 200       | 29   | .031 | .0237 | 1.3x | .050 |
   | 6  | 999       | 29   | .025 | .0237 | 1.05x| .050 |
   | 11 | 50        | 29   | .109 | .0149 | 7.3x | .050 |
   | 11 | 60        | 39   | .101 | .0209 | 4.8x | .049 |
   | 8  | 50        | 29   | .078 | .0196 | 4.0x | .051 |
   | 6  | 50        | 199  | .096 | .0452 | 2.1x | .049 |

   Hotelling's (n-k)/((n-1)k)*W ~ F(k, n-k) is .049-.052 at every cell, which validates the
   oracle. Read in absolute terms the printed test over-rejects against nominal 5% only for
   n/k below roughly 8 — .109 at (k=11, n=50) and .175 at (k=6, n=20) are real over-rejections,
   while (k=6, n=200) is .031, i.e. BELOW nominal. So the practical exposure is a long pre
   window with a modest reps(), not reps() in general.
2. *More clusters is worse at fixed reps, in absolute terms.* At (k=6, n=50) the size rises
   from .063 at G=30 to .096 at G=200 (the excess factor falls, but the baseline rises faster).
   The prober's G=40 design sits near the mild end.
3. *The doc leg as stated does not stand.* sthlp:532's "Neither is the Wald test on e(V) that a
   hand test reproduces; that one is ... not printed" sits inside the **jointtests** paragraph
   ("Neither [of those two bootstrap forms] is ..."), and the same paragraph says of
   pretrendpower "The method is unchanged - only the covariance it consumes", which is exactly
   what the code does. The sentence's factual clause is nonetheless contradicted whenever
   pretrendpower() is also asked for — see BX-N1, which is the version of this leg that stands.

**Also measured, against the prober:** on well-behaved data at realistic reps the ptp Wald and
wooldid's own bootstrap RANK joint test agree (r01: reps 20 .585/.714; 50 .778/.784;
100 .739/.713; 300 .621/.580). The four-orders-of-magnitude disagreement is a low-reps event.

repro: `r01_mech.do` (mechanism), `r03_oracle.py` + `r03_limits.out` (size),
`r11_boottstat.do` (hand-test identity)

---

## BOOTX-2 — PARTIAL (severity high -> MEDIUM); the headline mechanism is REFUTED

**The boundary arithmetic is confirmed and trivial:** the only screen is
`diag0cnt(invsym(S))` (ado:10505), and a sample covariance of n draws is numerically singular
only while n <= k, so the guard is silent from n = k+1 on.

**The realistic leg reproduces but means the opposite of what the prober says.** On design A the
ANALYTIC run reports keff 5 of 6 and prints the rank-deficiency caution; every bootstrap arm
(reps 50, 200, four seeds each) reports keff 6 of 6 with no caution (r02_rank.log). The prober
reads this as the bootstrap hiding a real deficiency. `r07_nulldir.do` decides it:

- the analytic pre block's smallest eigenvalue is -2.17e-17 and its eigenvector is
  v ∝ 2·pre4 − pre7 (only two nonzero entries);
- v'S_analytic v = -2.17e-17 (exactly zero);
- **v'S_bootstrap v = 1.10e-01, which is 71% of that block's average bootstrap variance**, and
  the SD of v'b_r across the 199 listwise-complete replication draws in e(bootdist) is 0.331.

So the direction the analytic covariance calls zero-variance has ordinary variance in the data.
`r06_whyrank.do` locates it: the analytic pre block is rank 5 of 6 under cluster(state),
cluster(county) AND `robust`, and **full rank 6 of 6 under iid** (min eigenvalue 5.98e-2). The
exact deficiency is an artifact of the analytic sandwich covariance, not a property of the
estimator; the bootstrap's full rank is CORRECT, and it is the ANALYTIC caution that is the
false alarm. "The guard goes silent exactly where S^-1 is worst" is not what happens.

**The prober's numbers are reps(7) artifacts.** mdt_wald's systematic understatement is
sqrt((n-k-2)/(n-1)) — 7.4% at (k=6,n=50), 13% at (k=11,n=50) — not 3.8x.

**What survives, and is worth the author's time:** the printed slope, its se, its CI and
mdt_wald carry a Monte-Carlo component the printed se does not contain and nothing records.
Four seeds, same data, reps(50): slope .0505 / .0195 / .0011 / .0110 (across-seed SD ~ .021)
against printed se .029-.043 — the seed is ~55% of the printed se; printed joint p over the
same four seeds .271 / .714 / .399 / .679. At reps(200): slope .0274/.0303/.0202/.0284
(SD ~ .0044), p .490/.583/.550/.542. Oracle coverage of the nominal-95% slope CI: 92.8% at
(k=6,n=50), 88.7% at (11,50), 84.3% at (6,20). See BX-N2.

**Doc leg is wrong:** sthlp:532 promises only that the guard skips "if it is not invertible",
which is precisely what the code does. No draw-count screen is promised.

repro: `r02_rank.do`, `r06_whyrank.do`, `r07_nulldir.do`

---

## BOOTX-3 — PARTIAL (severity medium -> LOW)

Arithmetic confirmed on my data (`r05_tt_band.do` part C, 40 states / 16 periods,
esprelength(5) espostlength(5), colsof(e(b)) = 12): reps(8) -> rank(e(V)) = 7 = completed-1,
min eigenvalue -7.17e-17, and a hand `test` over five pre coefficients runs at df 5 with
drop = 0 and no message. reps(20) and reps(60) -> rank 12 of 12, min eigenvalue 9.8e-4 / 1.8e-3.
The analytic arm has rank 11 of 12 and Stata DOES drop one (df 4, drop 1).

Why it is low, not medium:
- the prober's own calibration shows the arithmetic is Stata's `bootstrap` prefix, not wooldid's;
- **the doc leg is weak**: sthlp:1340 promises a hand `test`/`lincom` reproduces the
  *bootdisplay(wald) p-values*, which are one-dimensional per-coefficient t/z, and it does; it
  never promises a k-dimensional test;
- the reachable version needs completed reps below the dimension of the block actually tested.
  A 72-column e(V) at rank 49 still gives a full-rank subblock for any 10-20 coefficient hand
  test, and e(V_es_att) with 16 horizons is full rank at reps(50).

Their strongest leg is real but is a ruling, not a defect: ado:6235 says verbatim "reps() is
mandatory (author ruling: no floor - low-rep runs are legitimate)" while ado:5890 enforces
permreps() >= 20.

repro: `r05_tt_band.do`

---

## BOOTX-4 — CONFIRMED (severity LOW, as filed)

My first thirteen own-built shapes never reached the precondition (`r08_tt2.do`: timetrends +
cluster on i at 40 and 24 states x 8/10/12 periods, with and without an ES window;
`r09_attrit.do` part A: 26 states x 18 periods with every treated unit on a distinct onset).
So the all-zero-covariance geometry is rarer than sthlp:1059's "typically" suggests. I then
read (did not run) their p15 for the geometry — 40 states, 12 periods, onsets 5/7/9 with ten
states each and ten never-treated, esprelength(3) espostlength(3) — and rebuilt it with MY
generator and MY seed (`r13b_disc.do` part B). It reproduces exactly:

- analytic: `Caution: estimation returned an all-zero variance-covariance matrix`,
  e(wdid_vcvstate) = **absent**, `confirm matrix e(V)` -> rc 111, and the printed table has
  every se / t / p / CI **missing** (att 2.6784309 with a blank se);
- bootstrap(standard) reps(40) rseed(4242): rc 0, e(boot_completed) = 40, and the printed table
  is **fully populated with bootstrap quantities** (att 2.6784309, se .86374121, t_wald 3.10,
  p 0.000, CI .95021064 to 4.2682074; every ES row likewise);
- yet e(wdid_vcvstate) is still **absent**, e(boot_v_posted) = **0**, `confirm matrix e(V)` ->
  rc 111, e(V_analytic) -> rc 111, e(V_boot) -> rc 0, and a hand `test ES_att_post0 = 0` fails
  **rc 321**.

So sthlp:1059's flat "every standard error, test statistic, p-value and interval in the run is
missing" is false under bootstrap(), exactly as filed. The severity call is right: the runtime
caution names e(V_boot) and the post-estimation failure is loud, so this is a documentation
divergence plus an inconvenience, not a wrong number.

repro: `r13b_disc.do` (part B; `r08_tt2.do` and `r09_attrit.do` record the thirteen shapes that
do NOT reach it)

---

## BOOTX-5 — PARTIAL (severity LOW, confirmed with a correction)

`wooldid_bootqs` (ado:~16770) is x = n*p; whole x -> mean of order statistics x and x+1;
otherwise order statistic ceil(x). At level 95 that IS the maximum for every n <= 19. But being
the maximum is not the defect: the expected coverage of the r-th of n order statistics is
r/(n+1), so

| n   | rule               | E[coverage] |
|-----|--------------------|-------------|
| 8   | v[8] = MAX         | .8889 |
| 12  | v[12] = MAX        | .9231 |
| 19  | v[19] = MAX        | **.9500** |
| 20  | avg(v[19],v[20])   | .9286 |
| 40  | avg(v[38],v[39])   | .9390 |
| 100 | avg(v[95],v[96])   | .9455 |
| 400 | avg(v[380],v[381]) | .9489 |

n = 19 is the maximum AND exactly right. Measured on design B (40 states, 16 periods,
esprelength(5) espostlength(5), simultaneousbands, one seed): analytic-simulated 2.6038;
bootstrap reps 8 -> 2.4543, 20 -> 2.3978, 40 -> 2.3535, 100 -> 2.6211, 400 -> 2.6646. So the
low-reps values do sit low, but the reps(8) value came out LARGER than reps(40) — the narrowing
is an expectation, not something visible in a run. And the claim "no caution keyed on the draw
count; only the generic 2/reps caution mentions reps at all" is half true: that generic caution
DID fire at reps(8) and reps(20) on my runs and not at 40+.

repro: `r05_tt_band.do` part B

---

# NEW FINDINGS (found while refuting; all mine unless marked)

## BX-N1 [medium/certain] — two differently-referenced "joint pre-trend test" p-values print in the same bootstrap run, and under `boottstat` the only one printed is the one the help says is never printed

`jointtests` prints `p-value for joint test that all pre-treatment att horizon estimates are 0:`
— the bootstrap RANK test on the fixed bootstrap metric, wooldid's own valid construction.
Twenty lines later `pretrendpower()` prints `joint pre-trend test: F(k, G-1) = ..., p = ...`
— a normal-theory Wald on the same bootstrap covariance. Same coefficients, same null, no
cross-reference, no statement that they are different tests.

Measured (r00_setup.log, design A at 4 counties/18 periods, reps(10), NO pre-trend in the DGP):
rank test **0.0909**, ptp Wald **0.0000** (F(6,29) = 8.5244). Measured again on 40 states /
22 periods / 6 treated / esprelength(9) at reps(200) (r10b_attrit.log): rank test **.00995**,
ptp Wald **.000814** — a factor of 12 at a replication count nobody would call small.

The sharp version is `boottstat` (r11_boottstat.log, design A, reps(60), rseed 2468):
- the run prints `Caution: the bootstrapped joint_horiz joint test is unavailable (no usable
  replication statistic - rank-deficient replication covariances, or a horizon missing from the
  replications).` and `e(joint_horiz)` is missing, `e(boot_jt_family)` = studentized. On this
  panel that is STRUCTURAL, not bad luck: every replication inherits the singleton-treated-cell
  degeneracy of BX-N3, so the studentized joint pre-trend test is unavailable at every reps();
- the run's own boottstat caution says "... and the joint pre-trend tests, the across-subgroup
  equality tests and the adversarial line tests are bootstrapped as test statistics";
- and the pretrendpower block nonetheless prints `joint pre-trend test: F(6, 29) = .957282,
  p = 0.4708`, which a hand `test ES_att_pre7 ... ES_att_pre2` reproduces to ten digits
  (.4708198565 both ways, df 6, drop 0) — i.e. **exactly the normal-theory Wald on e(V) that
  sthlp:532 says "is not printed"**, and the only joint pre-trend p the user sees.

Under `boottstat` the printed value is also bit-identical to the plain-bootstrap one
(.4708198565 in both), so the studentization the option exists to provide never reaches it.

repro: `r11_boottstat.do` (self-contained; prints the hand test beside e(roth_prep_att) for
analytic, bootstrap and boottstat arms)

## BX-N2 [medium/probable] — the Monte-Carlo component of the printed Roth slope / CI / mdt_wald / joint p is unrecorded and, at reps(50), is about half the printed standard error

The Roth GLS quantities are functions of S^-1, and under bootstrap() S is a draw. Nothing in the
printed line or in e() distinguishes sampling uncertainty from replication uncertainty, and
e(roth_slope_se_att) contains only the former. Four seeds on ONE dataset (design A, r02_rank.log):

| reps | slope by seed 101/202/303/404 | across-seed SD | printed se | printed joint p by seed |
|------|-------------------------------|----------------|------------|-------------------------|
| 50   | .0505 / .0195 / .0011 / .0110 | ~.021          | .029-.043  | .271 / .714 / .399 / .679 |
| 200  | .0274 / .0303 / .0202 / .0284 | ~.0044         | .033-.040  | .490 / .583 / .550 / .542 |

At reps(50) the seed moves the slope by about half the width of the standard error printed
beside it, and moves the printed joint p over a 0.44 range. Oracle coverage of the nominal-95%
slope CI (r03_oracle.out): 92.8% at (k=6, n=50), 88.7% at (k=11, n=50), 84.3% at (k=6, n=20).
The systematic part is E[S^-1] = Sigma^-1 (n-1)/(n-k-2), which understates the slope se and
mdt_wald by sqrt((n-k-2)/(n-1)) — 7.4% at (6,50), 13% at (11,50).

repro: `r02_rank.do` (four seeds x two reps counts), `r03_oracle.py`

## BX-N3 [medium/probable] — the singleton-treated-cell degeneracy reaches `robust`, silently rank-deficits the analytic event-study covariance, and cascades into pretrendpower() and the boottstat joint test with no caution about its cause

Found while deciding BOOTX-2. On design A — 60 counties in 30 states, i = **county**, so every
treated (i, relyr) cell is a single observation — the analytic pre block of e(V_es_att) is
exactly rank-deficient by one (r12b_corr.log):

- `robust`: Var(pre7) = .0849489, Var(pre4) = .0212372 (ratio exactly 4.000),
  corr(pre7, pre4) = **1 exactly**, and the whole of row 1 of the correlation matrix equals
  row 4. Smallest eigenvalue 2.90e-17.
- `cluster(state)` and `cluster(county)`: same structure; smallest eigenvalue -2.55e-17.
- `iid`: full rank, smallest eigenvalue 5.98e-2.
- The point estimates are NOT in a 2:1 ratio (b_pre7 = -.2489177, b_pre4 = -.0806095), so the
  covariance asserts a deterministic relation the coefficients do not have.
- The bootstrap contradicts it outright: for v proportional to 2*pre4 - pre7,
  v'S_analytic v = -2.17e-17 while v'S_bootstrap v = 1.10e-01 (71% of the block's mean
  bootstrap variance) and the SD of v'b_r over the 199 listwise-complete draws in e(bootdist)
  is 0.331 (r07_nulldir.log).

**Discriminator (r13b_disc.do part A).** Same data, i coarsened to the onset cohort with
`fe(county)` so the treated cells are no longer singletons: the pre block is **full rank** under
both `robust` (eigenvalues 3.17e-2 ... 5.53e-1, max off-diagonal |corr| 0.691) and
`cluster(state)` (1.22e-2 ... 3.20e-1, max |corr| 0.709). So the mechanism is the
**singleton-treated-cell degeneracy** — zero residuals in each treated cell — and it operates
under plain `robust`, where the clustering has nothing to do with it.

Why it is still reportable rather than a rediscovery: the project's record documents singleton
treated cells for `ccc_absorb()` and documents an all-zero covariance only on the
continuous-treatment / cluster-at-or-above-i geometry (sthlp:1076). Here the treatment is
binary, the VCE is `robust`, the failure is PARTIAL (rank 5 of 6, not all-zero), and nothing
names the cause. The consequence chain on the user's screen:

1. `pretrendpower()` prints its rank-deficiency caution and g-inverses the direction away, so
   the analytic slope, CI, joint test and detectable-trend columns are computed on 5 of 6 pre
   horizons — on a direction the bootstrap says carries ordinary evidence;
2. under `boottstat` each replication inherits the same degeneracy, so the run prints
   `Caution: the bootstrapped joint_horiz joint test is unavailable (no usable replication
   statistic - rank-deficient replication covariances, or a horizon missing from the
   replications)` and e(joint_horiz) is missing — leaving the ptp F(k, G-1) Wald as the ONLY
   printed joint pre-trend test (see BX-N1);
3. `wooldid y unit t onset` with i = the unit is the default way the command is called, so the
   singleton geometry is the common case, not an edge case.

Honest limit: singleton treated cells are necessary but not sufficient. A 40-state i = state
design with six treated states also has singleton cells and its analytic pre block is full rank
(keff 9 of 9, r10b_attrit.log), so something about design A's shape — two onsets, twenty treated
units sharing them, i finer than the onset — is doing additional work that I did not isolate.
That is why the confidence is `probable` and not `certain`.

repro: `r12b_corr.do` (the three covariance blocks), `r13b_disc.do` part A (the discriminator),
`r07_nulldir.do` (the bootstrap contradiction)

## Negative results worth recording
- **Listwise attrition behind e(V_es_att) does not fire on a plausibly hostile panel.** I
  expected the (n, k) ratio driving BOOTX-1/2 to be much worse than reps() suggests because
  `wooldid_bootves` drops any replication missing ANY column of the family. On 40 states with
  only 6 treated, 22 periods, esprelength(9) espostlength(6), reps(200): every one of the 15
  ES_att columns was available in all 200 replications; listwise complete = 200; e(boot_N_V) =
  200. No attrition at all. (`r10b_attrit.do`)
- **13 own-built `timetrends` x cluster(i) shapes produced a usable analytic covariance**
  (see BOOTX-4), including one where every treated unit has a distinct onset. The geometry that
  DOES reach the all-zero covariance is narrow: three onsets of ten states each in a 12-period
  panel with ten never-treated.
- **The bootstrap V_es construction matches Stata's own convention**: `wooldid_bootves`
  (ado:14622-14625) centres at the replicate mean and divides by (rows - 1). No divisor or
  centring defect.
