# FINDINGS — bootstrap(bayes), the Dirichlet engine
Agent: boot-bayes.  Method: reimplementation of the weight-generation process +
distributional testing.  All probes in this directory; every .do has its .log beside it.

## HEADLINE: the engine itself is CLEAN under reimplementation
I wrote an independent Mata implementation of the documented draw process and it
reproduces the multipliers wooldid actually applied **bitwise**, in both
`bootalpha(equal)` and `bootalpha(size)` (p02_reimpl.do, log lines 208-211/275-278:
`exact bitwise equal = 1`, `maxabsdiff = +0.0000000000000X-3ff`).  I also verified
bitwise that the Dirichlet multiplier MULTIPLIES the user weight rather than replacing
it (p04_wreplay.do).  Detailed negative results are at the bottom.

## FINDINGS
  F1  bootstrap(bayes) + safety(off) reports standard errors 4.4x apart for a
      BIT-IDENTICAL point estimate, and the cell-weighted standard error is exactly
      INVARIANT to cross-cluster treatment-effect heterogeneity.
      repro: p19_homogeneous.do (self-contained), p15_attvsattit.do (self-contained)
  F2  bootsave() under bootalpha(size) writes 1/(rows in the user-if sample) on the
      observed row, where the help promises 1/(estimation sample).
      repro: p13b_docsample_standalone.do (self-contained)
  F4  at cluster()==i the bayes percentile interval for att missed 7/30 at nominal
      95% where bootstrap(standard) missed 2, on the same panels at equal width;
      and every reported se there is below the true sampling sd, bayes att_i worst
      at 0.621.  repro: p22_coverage_selfcontained.do (self-contained)
  (not a finding) Stata 19's own rwgen/bayesboot IS an external oracle for this
      engine and wooldid matches it BITWISE. repro: p14b_oracle_standalone.do

## F1 (SEEDED - lead 5 asked what comes out under safety(off)) — at cluster() == i,
## `bootstrap(bayes)` + `safety(off)` reports standard errors 4.4x apart for a
## BIT-IDENTICAL point estimate (att vs att_it)

SCOPE, stated up front: every number below is measured at `cluster() == i`. The exact
invariance mechanism (a cluster-constant shift of a cluster's treated cells adds a
constant to every replicate of an equal-cell-weighted estimand) holds exactly only
there. With `cluster()` coarser or finer than `i` the Dirichlet reweights units WITHIN
a cluster, att_it's composition does move, and the magnitude is UNMEASURED - I did not
test it. `cluster() == i` is also the configuration in which the package's analytic
clustered SE is already known to undercover (recon sec 1), which is why the analytic
SE below is blind to the heterogeneity too.

repro: `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes/p15_attvsattit.do`
run  : `stata-mp -e do p15_attvsattit.do` (self-contained; generates its own panel)

On a fully balanced 24-cluster staggered panel every treated (i, relyr) cell holds
exactly one row, so the observation-weighted `att` and the equal-cell-weighted
`att_it` are the SAME FUNCTIONAL OF THE DATA. wooldid agrees: their point estimates
and their analytic standard errors are bit-identical.

    ANALYTIC  att    b = +1.047e22bacc981X+000   se = +1.96ee6cfe4673bX-004
    ANALYTIC  att_it b = +1.047e22bacc981X+000   se = +1.96ee6cfe4673bX-004
    identical point estimates: 1   identical ses: 1

    STANDARD (reps 400, rseed 31337)  se(att) = 0.34376133  se(att_it) = 0.34376133
                                      ratio att_it/att = 1.0000      <- correct
    BAYES  + safety(off) (same reps/seed)
                                      se(att) = 0.32877148  se(att_it) = 0.07458774
                                      ratio att_it/att = 0.2269      <- 4.4x apart

se(bayes)/se(standard) is 0.956 for att and 0.217 for att_it. So the Dirichlet engine
reproduces the pairs bootstrap for the observation-weighted estimand (as the help
claims) and collapses for the cell-weighted one - on numbers that are the same number.

Mechanism (verified, not conjectured): with cluster() == i the Dirichlet multiplier is
constant within a cohort, so it leaves that cohort's own cell effects essentially
unchanged; ALL of the between-cluster variability enters through the cell weights,
which `att_it`/`att_i`/`att_itime` discard by construction. The residual variance is
the regression-coefficient channel only.

Status vs what is already ruled: the DIRECTION is documented (sthlp:241, and the
refusal text at ado:6439 says the SEs "can be too small"), and the combination is
refused unless the user types `safety(off)` - which the refusal message itself tells
them to do. What is new here is (a) the MAGNITUDE, 4.4-5.6x on ordinary designs, i.e.
a nominal 95% interval with roughly the content of a 26% interval, and (b) the
demonstration that the two SEs disagree for a point estimate that is bit-identical,
which is checkable in one table by any user and is not what "can be too small" leads
a reader to expect.
Also measured on a second, unbalanced design (p09_cellwtd.do):
  att .2829/.3160 = 0.895 ; att_it .0590/.3160 = 0.187 ; att_i and att_itime .0596/.3327 = 0.179.

severity: medium-high   confidence: certain   realistic: yes
surface: bootstrap(bayes) x {att_it, att_i, att_itime} x safety(off)

## F2 (own) — `bootsave()` observed row under `bootalpha(size)` is 1/(rows in the
## user-if sample), not 1/(estimation sample), which is what the help promises

repro (SELF-CONTAINED): `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes/p13b_docsample_standalone.do`
(the original measurement is `p13_docsample.do`, which reuses `p01_panel.dta`.)

sthlp:557: "Under bootalpha(size) the column holds the per-observation multiplier
actually applied ... and the _rep = 0 row holds 1/N, N the number of observations in
the estimation sample."

Measured on p01_panel with 15% item nonresponse on the outcome (437 rows in the data,
375 estimated):

    ROWS_IN_DATA = 437
    e(N) estimation sample = 375
    bootsave _rep=0 row value = +1.2bef98e5a370fX-009 = .00228833   = 1/437
    1/e(N)                    = +1.5d867c3ece2a5X-009 = .00266667   = 1/375
    sum_c n_all_c * w_c (rep1) = 1.0000000000000000    (all data rows)
    sum_c n_y_c   * w_c (rep1) = 0.8595794591351535    (estimated rows)

ado:437 divides by `wdb_r_nrows', which wooldid_bootrun sets from the MAP FRAME - the
user-`if' sample intersected with nonmissing cluster() - which still holds every row
the estimator later drops. A user who follows the help and rescales the saved
multipliers by e(N) is off by the missingness rate (16.5% here).

The same map-frame count is the Dirichlet SHAPE alpha_c under bootalpha(size), so a
cluster with more item nonresponse carries more concentration than its contribution to
the estimate warrants. The MEAN is unaffected (E[surviving mass] stays proportional to
the cluster's estimated row count, so the replications remain centred), only the
dispersion, by about sqrt(n_estimated/n_rows) per cluster - ~7% here. Reported for
completeness, not as the substance of the finding.

severity: low   confidence: certain   realistic: yes (item nonresponse is ordinary)
surface: bootsave() x bootalpha(size) documentation

## (NOT A DEFECT - a method note, kept here because it retires the task premise)
## The claim "there is NO external oracle for bootstrap(bayes)" is FALSE, and the
## engine passes that oracle BITWISE

Stata 19 ships `rwgen` / `bayesboot` as official commands
(`/Applications/Stata/ado/base/r/rwgen.ado`, v1.6.5 21jan2025). At a common seed:

    bootalpha(size)   == `rwgen bayes S, reps(20) id(st)`                 max|diff| = EXACTLY 0
    bootalpha(equal)  == `rwgen bayes E, reps(20) id(st) priorpowers(-n_c)` * n_c,
                                                                        max rel diff 1.93e-16

repro (SELF-CONTAINED): `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes/p14b_oracle_standalone.do`
(the original measurement is `p14_oracle_bitwise.do`, which reuses p06's saved files.)

This is not a defect - it is a positive result the review should record, because it
means the Dirichlet engine has a cheap, exact regression test available to the project
that it is not currently using in this form, and it retires the "unverifiable engine"
framing. It also independently confirms the ado:8899 comment and sthlp:548.

severity: low (informational)   confidence: certain   realistic: n/a
surface: bootstrap(bayes) / bootalpha() verification

## NEGATIVE RESULTS — what I attacked and what came back CLEAN
(these are the bulk of the work; each is worth as much as a finding)

### N1. The draw generator, reimplemented — CLEAN (certain)
`p02_reimpl.do`. I wrote an independent Mata implementation of the process the code
and help document (per-replication outer loop, per-cluster-ascending inner loop,
`rgamma(1,1,alpha,1)` from the seeded stream, normalise by the replication's column
sum, divide by n_c under `size`) and compared it to the multipliers wooldid actually
applied, recovered via `bootsave()`. `exact bitwise equal = 1`, `maxabsdiff = 0`, in
BOTH `bootalpha(equal)` and `bootalpha(size)`. The declared RNG consumption order is
the real one. Anything wrong with index order, transposition, per-replication vs
per-cluster normalisation, or alpha selection would have shown up as a nonzero diff.

### N2. Distributional goodness-of-fit on the real multipliers — CLEAN (certain)
`p05_gof.do` / `p05b_analyze.do` / `gof.py`, R = 1000 real replications, C = 12
clusters with row counts 6..12.
  equal: rowsum(p) deviates from 1 by <= 4.4e-16; cluster means vs Dirichlet(1,...,1)
  max|z| = 1.84; per-cluster KS vs Beta(1,11) min p = 0.113 (Bonferroni 0.0042);
  mean off-diagonal correlation -0.09078 vs the exact theory -1/(C-1) = -0.09091;
  lag-1/2/3 serial correlation across replications <= 0.094 against 2/sqrt(R)=0.063;
  a 2000-shuffle exchangeability permutation test on the range of cluster means
  p = 0.494; pooled PIT mean .50111 var .08316 (theory .5 / .08333).
  size: rowsum(n_c w_c) - 1 <= 4.4e-16; KS vs Beta(n_c, N-n_c) min p = 0.0239
  (passes Bonferroni); corr(n_c, mean p_c) = 0.999 as the convention requires.
A wrong alpha, a shared Gamma draw across clusters, a reused stream position, or a
normalisation over the wrong margin would all have fired here.

### N3. Does the Dirichlet weight reach the estimator correctly? — CLEAN (certain)
`p04_wreplay.do`. On a design with one cluster carrying 500x the weight of the others
(so "multiply" and "replace" cannot coincide), I rebuilt `manw = bigpop * w_c` from the
bootsave file and re-ran the replication's own command by hand:
`wooldid y st year onset [aw=manw], att cluster(st) fastnose esrelativeto(-1)`.
Manual and `e(bootdist)` agree BITWISE for replications 1-3 in both alpha modes
(`manual_b=+1.1843451b7ddfcX+000  boot_b=+1.1843451b7ddfcX+000`). The negative control
(passing the multiplier as the weight instead of multiplying) differs in the first
significant digit. The multiplier multiplies; it does not replace.

### N4. `bootalpha()` validation — CLEAN (certain)
`p03_valid.do`. `equal` and `size` accepted. `EQUAL`, `Size`, `siz`, `eq`, `1`,
`equal2` are ALL refused r(198) with the full explanatory message - no case folding,
no Stata-style abbreviation, no silent default. `bootalpha()` with
`bootstrap(standard)` -> r(198) ("applies only to bootstrap(bayes)"); with no
`bootstrap()` at all -> the stray-companion r(198). `boottstat` + bayes -> r(198).
An empty `bootalpha()` takes the documented default `equal`.

### N5. Stata matrix ceiling on the C x R draw matrix — CLEAN (probable)
`p07_bigN.do`. The draw matrix is a Stata matrix and unclustered `bootstrap(bayes)`
sets C := number of observations, so I expected `c(max_matsize)` = 11000 to bite on an
ordinary micro panel. It does not: `st_matrix()` called from Mata builds 15000 x 10
without error (rc = 0 at 5000 / 11000 / 11001 / 15000). Not a defect. I did NOT test
the memory ceiling at, say, N = 500,000 x reps(999).

### N6. Reproducibility and the reps() prefix property — CLEAN (certain)
`p16_surfaces.do` R1/R2. Two runs with the same `rseed()` give bitwise identical
`e(bootdist)`; the first 50 rows of a `reps(120)` run are bitwise identical to a
`reps(50)` run at the same seed (so the outer loop really is per-replication).

### N7. F001 (the 2026-09-13 weighted-divisor fix) — CONFIRMED BUILT (certain)
`p16_surfaces.do` S4/S5 and `p08b_sescale.do`. The defect's signature was
|mean(draws) - estimate| / sd(draws) = 3.5 to 9.4 with every draw on one side. On a
weighted unbalanced design I measure -0.0190 (`bootalpha(size)`) and +0.0861
(`bootalpha(equal)`). The ruling is in the code, not just in a commit message.

### N8. Centring across the option surfaces — CLEAN (certain)
`p16_surfaces.do`. Event study (`esprelength(3) espostlength(3)`), `subgroups()` +
`jointtests` + pairwise contrasts: every column's (mean(draws) - b)/sd(draws) is within
+/-0.20, i.e. no displaced replication distribution anywhere. 26 columns checked.

### N9. "bayes silently keeps near-degenerate replications where standard fails" —
### my own hypothesis, NEGATIVE (certain)
`p11_fewcontrols.do`. 20 clusters, only 2 never-treated. `bootstrap(standard)` failed
61/500 replications (and said so); `bootstrap(bayes)` failed 0/500. I predicted bayes
would instead retain replications in which the control clusters drew near-zero
Dirichlet mass, producing heavy tails that silently inflate `e(V)`. It does not:
bayes draws have kurtosis 2.88 (thinner than normal), sd/IQR-sd = 0.957,
corr(|deviation|, control-cluster mass) = -0.057, and mean|deviation| is 0.090 for the
35 replications with control mass < 0.02 against 0.074 for the other 465. What IS true
and worth a line: bayes's interval there is 14% NARROWER than standard's (sd .0923 vs
.1075) and reports no failures, while standard discloses that 12.2% of its
replications were discarded, so the retained standard distribution is conditional on
an event the user is told about and the bayes one silently is not conditional at all.
Which is closer to correct is not settled by this probe.

### N10. The two documented quantitative claims about the engine — CLEAN (probable)
`p08b_sescale.do`, `p17_seratio_realistic.do`, reps 400-600 per cell.
  sthlp:532 "the Bayesian bootstrap ... shares [the pairs cluster bootstrap's] first
  two moments" => sd(bayes,equal)/sd(standard) should be ~1 (exactly sqrt(C/(C+1)) =
  0.984 at C = 30 for the weight distribution itself). Measured:
      0.9782 (C=30, unweighted, p04_panel)
      0.9258 (C=30, unweighted, p16_panel)
      0.9331 (C=30, [aw=pop], lognormal(0.4) population weights, p16_panel)
  All within Monte Carlo error (ratio-of-sd sd ~4% at these replication counts).
  sthlp:548 bootalpha(size) "shrinks the cluster-weight dispersion by roughly the
  square root of the mean cluster size" => sd(size)/sd(equal) ~ sqrt(C/N). Measured
  0.2406 and 0.2213 against sqrt(30/437) = 0.2620; 0.3003 against sqrt(30/360)=0.2887.
  Both documented statements hold.

### N11 (worth one line, UNCERTAIN, own). One dominant weight breaks the moment match.
On `p04_panel`, where one of 30 clusters deliberately carries 500x the aweight of the
rest, sd(bayes,equal)/sd(standard) falls to 0.6362 (p08b_sescale.do, [wtd] arm) while
the unweighted arm on the same data is 0.9782. The mechanism is not a bug: the pairs
bootstrap omits the dominant cluster entirely with probability (29/30)^30 = 0.36 and
the statistic is a ratio, so the two weight schemes stop agreeing beyond the first two
moments of the WEIGHTS. With ordinary lognormal population weights (max/min = 8.3 on p16_panel)
the ratio is back at 0.9331. I flag it only because a state panel weighted by population has
a California/Wyoming ratio near 70x, which is between the two designs I measured, and
sthlp:532's "shares its first two moments" does not prepare a reader for a 36% gap.
Not filed as a defect; filed as a measurement.

## F1 — ADDENDUM (the sharpest form of it): the bayes cell-weighted standard error is
## INVARIANT to the very variance component it is documented to omit

repro: `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes/p19_homogeneous.do`
(same design as p15, with the cluster-level effect shock switched from 1.2 to 0)

    [het=0  ] analytic se(att) = 0.099348   identical b: 1  identical se: 1
    [het=0  ] standard se(att) = 0.123593   se(att_it) = 0.123593   ratio 1.0000
    [het=0  ] bayes    se(att) = 0.113876   se(att_it) = 0.074588   ratio 0.6550
    [het=0  ] bayes/standard:  att = 0.9214   att_it = 0.6035
    [het=1.2] analytic se(att) = 0.099348   identical b: 1  identical se: 1
    [het=1.2] standard se(att) = 0.343761   se(att_it) = 0.343761   ratio 1.0000
    [het=1.2] bayes    se(att) = 0.328771   se(att_it) = 0.074588   ratio 0.2269
    [het=1.2] bayes/standard:  att = 0.9564   att_it = 0.2170
  (both blocks are printed by the SAME run of p19_homogeneous.do)

bootstrap(standard)'s att_it standard error rises 2.8x (0.1236 -> 0.3438) when the
treatment effect is given a cluster-level shock. The bayes one does not move at all:
0.074588 in both. The mechanism is exact, not approximate - with cluster() == i a
cluster-constant shift of the treated cells shifts every one of that cluster's cell
effects by the same constant, and `att_it` equal-weights a fixed number of cells per
cluster, so the shift is a constant added to every replicate and the replicate SPREAD
is unchanged. (The analytic clustered SE is blind to it too, for the same reason: the
cohort x relative-period dummies absorb the shock and the residuals are unchanged.)

Two consequences for the ruling as it stands:
 1. The refusal text conditions the understatement on effects that "vary across
    clusters". It is 0.60 of the correct value even with NO cross-cluster
    heterogeneity, so the conditional is too narrow.
 2. Where heterogeneity IS present the SE does not merely understate - it is exactly
    insensitive to it, so a user cannot reason "there is some of this variance in
    there, just not enough."

### N12. The `if` sample, the boundary cases and BCa — CLEAN (certain)
`p18_contreat_etc.do`.
  `if st <= 15 | st > 22` -> e(N_clust) = 23, the bootsave file carries exactly 23
  `w_` columns and the replication-1 multipliers sum to exactly 1 (+1.0000000000000X+000).
  The Dirichlet is built on the restricted sample, as it must be.
  2 clusters and 1 cluster both exit r(459) informatively: the 2-cluster case by the
  no-control-group screen ("you have no never-treated control group..."), the 1-cluster
  case by the engine's own screen ("bootstrap() requires at least 2 clusters in the
  estimation sample.").
  `bootdisplay(bca)` + bayes runs: acceleration -0.0107, z0 -0.0753, BCa interval
  [1.0399, 1.5389] against percentile [1.0586, 1.5713] - a sane, small correction.

### N13. NOT TESTED
  - `contreat()` + bayes: my probe DGP was collinear (`dose = post*dose0`) and BOTH the
    analytic and the bayes runs exited r(322) with the same informative message, so the
    combination is untested rather than clean.
  - the memory ceiling of the C x R draw matrix at, say, N = 500,000 unclustered.
  - `bootstrap(bayes)` + `poisson` beyond a smoke test (p03 S8 rc = 0).
  - `permonsetcohorts`-adjacent finer-than-i work (bayes takes no re-identification leg,
    so the T1/T2/T3 machinery does not apply to it).
  - **THE MAIN GAP: the cell-weighted understatement (F1) at cluster() != i.** Every
    measurement in F1 and F4 is at cluster() == i. With cluster() coarser than i (i =
    unit, cluster = state) or finer/cross-cutting, the Dirichlet reweights units WITHIN
    a cluster, so att_it's composition does move and the exact-invariance mechanism no
    longer applies in that form. The magnitude there is unmeasured. This is also the
    cell the recon note's section 7 calls out ("any finer-than-i arm that FAILS
    coverage"), so it should not be read as clean.

## F4 (SEEDED - lead 6 asked for a modest coverage MC) — at cluster() == i the
## bayes percentile interval for `att` missed 7 of 30 at nominal 95% where
## bootstrap(standard) missed 2, on the SAME 30 panels, with equal widths

repro: `p22_coverage_selfcontained.do` (self-contained; 30 simulated 30-cluster
panels, 79 replications each, true att = 1, homogeneous effects, cluster() == i).
Original run: `p12b_coverage.do` + `p20_cov_analyze.do` + `p21_covdetail.do`.

(reproduced exactly by the self-contained p22_coverage_selfcontained.do; the
analytic arm has no percentile family, shown as "-")

    arm       est     n  covered(pctl)  covered(wald)  mean se   true sd   se/sd
    analytic  att    30       -              26        .12554    .17784    .706
    analytic  att_i  30       -              25        .14698    .17942    .819
    bayes     att    30      23              28        .14392    .17784    .809
    bayes     att_i  30      24              24        .11138    .17942    .621
    standard  att    30      28              29        .15078    .17784    .848

Two things, both at cluster() == i:

(a) EVERY reported standard error here is below the true sampling sd of the
    estimate - the package's documented analytic understatement (recon sec 1) is
    not repaired by either bootstrap, and `bootstrap(bayes)` for the cell-weighted
    `att_i` is the WORST of the five at 0.621, worse than the analytic SE the user
    would be turning to the bootstrap to escape. That is F1 expressed as coverage.

(b) On the same 30 panels the bayes percentile interval for `att` covered 23/30
    where `bootstrap(standard)` covered 28/30, with mean percentile widths of
    0.6028 and 0.6146 - within 2% of each other - and identical point estimates.
    All 5 discordant panels go the same way (standard covers, bayes misses; never
    the reverse, p = 2^-5 = 0.031 on a sign test). The bayes interval is SHIFTED,
    not narrower: on draw 11, b = 1.2474, bayes gives [1.0038, 1.5676] and standard
    [0.9259, 1.4658]. bayes's misses are asymmetric (5 low-side, 2 high-side).

I am marking (b) UNCERTAIN: 30 panels give a +/-13pp band on a coverage
proportion, and I did not diagnose why the two percentile intervals sit at
different locations when both replication distributions are centred within 0.16
bootstrap-sd of the estimate. It is reported because it sits exactly in the
recon note's section 7 hunting ground ("any arm that FAILS coverage", "the effect
at a different cluster count"), it is the opposite direction from the open,
unruled `bootstrap(standard)` OVER-coverage result there, and 5/5 same-direction
discordances on paired data is not nothing.

severity: medium   confidence: (a) certain / (b) uncertain   realistic: yes
surface: bootstrap(bayes) coverage at cluster() == i
