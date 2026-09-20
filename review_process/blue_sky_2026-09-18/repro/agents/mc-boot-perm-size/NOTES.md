# MC size study of permute()/bootstrap() p-values under a sharp null — NOTES

WORKDIR: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-boot-perm-size

## State
- [x] Build asserted loadable via `adopath ++`. `which wooldid` -> pinned build.
- [x] t00_timing.do / .log: all five engines run on an 80-row state-year panel.
      TIMINGS (G=10, T=8, 80 rows, 1 estimand `att`):
        plain                              0.38 s
        permute(shuffle) permreps(79)     12.47 s
        bootstrap(standard) reps(99)      15.90 s
        bootstrap(standard) reps(99) boottstat 16.87 s
        bootstrap(bayes)  reps(99)        15.41 s
      => ~0.16 s per replication; per-rep overhead dominates, panel size barely matters.
- e() harvest names confirmed:
    perm: e(perm_p_att), e(perm_lb_att), e(perm_ub_att) (+ _att_i etc.)
    boot: matrices e(boot_p_percentile), e(boot_p_bc), e(boot_ci_percentile),
          e(boot_ci_bc), e(boot_ci_normal); colnames are `main_att` (underscore!)
          hidden full tables e(boot_wald_mainresults)/e(boot_pctl_mainresults)/
          e(boot_bc_mainresults); under boottstat e(boot_tpctl_mainresults)
    analytic joint controls: e(joint_horiz), e(adversarial_pval_att)  [need ES options]
- reps arithmetic: permreps(79) makes nominal 5% EXACTLY attainable for the
  equal-tailed add-one p: p = 2*(m+1)/80 <= .05 <=> m <= 1 -> 4/80 = .05 exactly.

## Design (committed)
Shape B "state-year panel": one row per (state, year). i = state = cluster.
This is the commonest applied DiD shape AND the region where the project's own
record says things are known bad (analytic SE undercoverage; bootstrap(standard)
over-covering att_i at cluster()==i) -> gives me two positive controls.
Sharp null: y = a_i + f_t + e_it, e AR(1) rho=.5 (Bertrand-Duflo-Mullainathan),
NO treatment effect anywhere.


## Calibration results (v01_validate.do / .log) — 2026-09-17
POSITIVE CONTROL 1 (analytic joint pre-trend, G=10, 60 draws):
  rej@5% = 0.433.  NF-1b's published number is 0.657 at G=10.  Same regime,
  same order of magnitude -> the DGP IS in the known-bad geometry.  (Mine is
  milder: fewer horizons and a milder within-state dependence.)
POSITIVE CONTROL 2 (analytic att / att_i, cluster(st)==i):
  coverage 0.917 both, rej@5% 0.083 -> the documented analytic UNDERstatement
  direction, milder magnitude than the 0.794-0.876 in review_process.
  NOTE: my v01 row-picking loop took the LAST 'att' row = `pretreatment:att`,
  not `main:att`.  Fixed in mc_main.do by matching on roweq=="main".
SANITY: e(perm_npossible) = 18900 = 10!/(2!2!2!4!). permute(shuffle) permutes
  ttre across i preserving cohort sizes.  Confirmed.
KEY STRUCTURAL FACT (useful): e(permdist_b) is BIT-IDENTICAL (max abs diff
  = +0.0000000000000X-3ff) across default / permabsp / permrawp at the same
  rseed.  The p-value convention does not perturb the RNG stream.
  Same run: default(equal-tailed add-one) .85 = 2*34/80; permabsp .7625 = 61/80;
  permrawp .759493670886 = 60/79.  All three formulas confirmed.
`adversarial` posts NO p without `esfixedbaseperiod` -- documented
  (sthlp:453, 758).  Not a defect.  Arm A9 adds esfixedbaseperiod.

## Timings per wooldid call, G=10, 80 rows, att+att_i
  plain                                  ~0.5 s
  permute(shuffle) permreps(79)          ~13 s
  permute ... permtstat                  ~17 s
  bootstrap(standard) reps(99)           ~20 s
  bootstrap(standard) reps(99) boottstat ~21 s
  bootstrap(bayes)  reps(99) (att only)  ~17 s

## Harness files
  dgp_snippet.txt  -- the sharp-null DGP program (mkdgp), pasted into every probe
  mc_main.do       -- the MC driver.  args: G d0 d1 "arms" hetv*100 tag
  runblock.sh      -- one block with the process-cap wait + deadline kill
  master1.sh       -- sequential chunk driver
  analyze.py       -- rejection rates + Wilson CIs from mcres_*.csv
  d01_specdiag.do  -- WHY does the ES-window arm report a different att?

## ACHIEVABLE NOMINAL LEVELS (must be quoted beside every rejection rate)
Rank-based p-values live on a grid, so "p <= .05" is not always a 5% rule.
  permute default (equal-tailed add-one), permreps(79):
      p = 2(m+1)/80 ;  p<=.05  <=>  m<=1  ->  exact level 4/80 = .0500
  permabsp, permreps(79):  p = (c+1)/80 ; p<=.05 <=> c<=3 -> 4/80 = .0500
  permrawp, permreps(79):  p = c/79     ; p<=.05 <=> c<=3 -> SAME region as
      permabsp at 5%.  The two conventions only separate at 1%, where permabsp
      CANNOT reject (floor 1/80 = .0125) and permrawp can (c=0 -> p=0), whose
      exact probability is 1/80 = .0125 > .01: the documented anti-conservatism
      of the no-add-one share.
  bootstrap percentile, reps(99): p = 2c/99 ; p<=.05 <=> c<=2 -> .0404
  boottstat tpctl, reps(99):      p = 2(1+m)/100 ; p<=.05 <=> m<=1 -> .0400
  bootstrap bc / wald: continuous, no grid effect.
So for percentile and tpctl the "p<=.05" rule is a 4% rule; the CI-excludes-0
criterion (recorded alongside) is a clean 95% rule and is the primary one for
those two families.

## reps() choice: 79 for BOTH engines (changed from 99 before any block ran)
permreps(79): default equal-tailed p<=.05 <=> m<=1 -> exact level 4/80 = .0500
reps(79)    : percentile p = 2c/79 <=.05 <=> c<=1 -> grid level 2*2/79 = .0506
              boottstat tpctl p = 2(1+m)/80 <=.05 <=> m<=1 -> 4/80 = .0500
At reps(99) the percentile rule would have been a 4.04% rule and tpctl a 4.00%
rule, which muddies a 5% size comparison. 79 aligns every rank rule on ~5%.

## d01 / d02 — a lead I chased and CLOSED as NOT a defect (2026-09-17)
d01_specdiag.do showed that on one noisy sharp-null draw the headline `att`
moved from +0.46203 (no event-study options) to -0.29540 with
esprelength(3) espostlength(3) and to -0.35225 with a window WIDER than every
realized horizon — a swing of ~0.8, larger than 2 reported SEs, with the sign
flipped.  Data mutation was ruled out (re-running the plain spec on the same
in-memory data and on a freshly reloaded copy gave a BITWISE identical
+1.d91d4163813feX-002).
d02_estimand_truth.do settles it: on an EXACT DGP (no error term at all)
  const effect  -> att = att_it = att_i = att_itime = 1 under ALL THREE specs
  ramped effect -> att = att_it = 1.5555556 and att_i = att_itime = 1.5 under
                   ALL THREE specs, matching the design truth exactly
so requesting event-study horizons does NOT change the estimand; it changes the
reference-period normalization and hence the ESTIMATOR, and the two unbiased
estimators simply differ by sampling noise on a 10-cluster panel.  NOT A DEFECT.
Worth knowing: the reported se on main:att swung 0.1704 / 0.2167 / 0.5188 across
the three specs on the same data.

## PRE-REGISTRATION (written 2026-09-17 00:26, BEFORE any block finished)
What each arm is GUARANTEED to deliver under my DGP, and what would be a finding.
My DGP randomizes the onset assignment uniformly over the cohort-size-preserving
assignments and generates y with NO dependence on the assignment.  Under that,
randomization inference is EXACT for ANY statistic, with NO asymptotics and with
no assumption about variance, serial correlation or cluster count.  So every
permutation arm carries a theorem; the bootstrap arms carry only asymptotics.

| arm | guarantee | nominal at my reps | a finding iff |
|---|---|---|---|
| A1 permute default (equal-tailed add-one) | EXACT | .0500 (4/80)   | rej >> .05 -> DEFECT |
| A2 permabsp                               | EXACT | .0500 (4/80)   | rej >> .05 -> DEFECT |
| A3 permrawp (c/n, no add-one)             | anti-conservative by <=1/(n+1) | same 5% region as A2 | rej >> .065 |
| A4 permtstat                              | EXACT | .0500          | rej >> .05 -> DEFECT |
| A10 permute-rewired joint pre-trend test  | EXACT (one-sided add-one rank) | .0500 (4/80) | rej >> .05 -> DEFECT |
| EX permexact (168 enumerated assignments) | EXACT | .0476 (8/168)  | rej >> .0476 -> DEFECT |
| PO permonsetcohorts (i rebuilt per rep)   | EXACT | .0500          | rej >> .05 -> DEFECT |
| A5 pctl / bc / wald                       | asymptotic only | .0506 p-rule, .05 CI-rule | report the number |
| A6 boottstat                              | asymptotic (2nd order) | .0500 | report; "never measured" |
| A7 bootstrap(bayes)                       | asymptotic | .0506 | report; no external oracle exists |
| A11 bootstrap-rewired joint test          | asymptotic fixed-metric rank | .0500 | report |
| A0 analytic att / A9 analytic joint       | KNOWN BAD | - | positive controls only |

## HV vs HVF — a distinction the project's own F049 measurement appears to miss
F049 (ado:6053-6056) reports "7 of 30 null datasets rejected at 5% with
treated-cluster error sd 3 vs 1, against 6 of 30 for the raw coefficient".
If the sd inflation is attached to TREATMENT STATUS, the randomization sharp
null is violated by construction (treatment changes the variance of y), so the
rejection rate measures POWER against a variance effect as much as size.  If the
inflation is attached to a FIXED cluster attribute drawn independently of the
onset assignment, randomization inference stays EXACT and any over-rejection is
a defect.  I run BOTH:
  HV  : sd 3 on the ever-treated states (the F049 shape; the applied
        "treatment is correlated with volatility" case -> expect inflation,
        which is the documented MacKinnon-Webb phenomenon, not a wooldid bug)
  HVF : sd 3 on the even-numbered ORIGINAL state ids, independent of the onset
        draw -> the sharp null HOLDS and exactness is a theorem.

## HARD STOP AND DECISION RULE (set 2026-09-17 00:30, before any block finished)
COMPUTE STOP 04:15.  REPORT WRITTEN BY 04:45.  Whatever has not run by then is
listed as DROPPED, with the reason, in the report.
DECISION RULE, fixed in advance: an arm with an EXACTness theorem is called a
DEFECT only if the 95% Wilson LOWER bound on its rejection rate exceeds its
achievable nominal level.  Detection bounds:
   n=100 draws: Wilson lower bound clears .05 at 11/100 (.11); clears .0476 at 11/100
   n= 60 draws: needs  8/60 (.133)
   n= 40 draws: needs  6/40 (.150)
   n= 30 draws: needs  5/30 (.167)
So at these sample sizes I can detect GROSS mis-sizing (2-3x nominal), not a
1-2 percentage-point distortion.  That limit is stated with every result.
BLOCK ORDER (value-weighted): chunk1 G=10 core (100) -> v02 verifications ->
A4/A6 (40) -> A10/A11 joint (40) -> G=40 core (30) -> PO (60) -> EX (40) ->
HVF (30) -> HV (30) -> top-ups.

## INTERIM SIGNAL at 7 draws (00:31) -- directional only, NOT a finding yet
  A7 bootstrap(bayes)   pctl 2/7 bc 3/7 rejections at 5% (one at 1%)
  A5 bootstrap(standard) pctl 0/7 bc 0/7 wald 0/7   on the SAME draws
  A1 permute             0/7
  A0 analytic att        2/7   (known-bad control, reproducing)
  A9 analytic joint      4/7   (known-bad control, reproducing; medp .0246)
  A7 wald                0/7 but median p .1535 vs A5 wald .2182
Candidate headline: bootstrap(bayes) is far more liberal than bootstrap(standard)
on identical data at G=10.  d03_bayes_vs_standard.do is written and ready to
diagnose whether the Dirichlet draw distribution is narrower or shifted.

## HARNESS CORRECTION (00:33) — e(perm_lb_<est>) is the NULL BAND, not the CI
d04_perm_dual.do chased what looked like a gross p-vs-interval contradiction:
on draw 8 wooldid reported e(perm_p_att)=.025 with e(perm_lb_att)=-1.2956 and
e(perm_ub_att)=+1.1604, an "interval" containing 0 by a mile.  Reconstructed
from e(permdist_b): the comparator is -1.3258, all 79 draws lie in
[-1.3009, 1.2912], so m0=0 and p = 2(1+0)/80 = .025 is CORRECT, and
-1.2956/1.1604 are exactly the 2.5/97.5 quantiles of the DRAWS.
The PRINTED TABLE is right: lb_95ci = -2.4863, ub_95ci = -0.0302, the shift
interval b_obs - q_hi / b_obs - q_lo, which excludes 0 as the p demands.
sthlp:1353 documents this: e(perm_lb_)/e(perm_ub_) hold the null band and are
DEPRECATED "since lb/ub elsewhere in wooldid means a genuine interval bound".
=> NOT A DEFECT.  My harvest was wrong.  mc_main_v2.do and the patched
mc_alt.do now take the permutation interval from e(wdidmainresults) cols 5/6.
CONSEQUENCE FOR THE DATA: in mcres_g10core.csv the A1 rows for draws 1-25
carry the null band in lo/hi.  Their p column is correct.  From draw 26 the
lo/hi are the real shift interval.  Analysis must not mix them.

## KEY INTERIM NUMBERS at 12-13 draws (00:36) — still preliminary
Rejection of a TRUE SHARP NULL at nominal 5%, G=10, i = state = cluster():
  bootstrap(bayes)      bc   5/12 = .417   pctl 4/12 = .333   (3/12 even at 1%)
  bootstrap(standard)   bc/pctl/wald       2/12 = .167
  ANALYTIC (the thing bootstrap is offered as a remedy for)  3/13 = .231
  permute(shuffle)                          1/13 = .077
On the SAME 12 draws the point estimate b is BITWISE identical across engines,
and the reported dispersion is
  true MC sd of the att estimator across draws   0.6921
  bootstrap(standard) se   0.4723   = 0.682 x truth
  bootstrap(bayes)    se   0.3964   = 0.573 x truth
  analytic cluster se      0.3800   = 0.549 x truth
  mean 95% CI width: standard 1.957, bayes 1.640, analytic 1.719;
  a correctly sized interval needs 2.713.
So at G=10 NEITHER bootstrap repairs the documented analytic understatement, and
bootstrap(bayes) is barely wider than the analytic se (1.10x) while
bootstrap(standard) is 1.29x.  sthlp:255 names "bootstrap() or permute()" as THE
remedies for the at-or-above-i undercoverage.  permute() delivers; bootstrap()
does not, and bootstrap(bayes) is worse than the analytic se it replaces.

## NEW, UNSEEDED (16 draws, 00:38): the analytic SE understatement is
## ESTIMATOR-specific, not a blanket property of cluster() >= i
Same data, same estimand `att`, same cluster(st) == i, 16 draws:
  spec                                     sd(b)   mean se   se/sd   rej@5%
  wooldid y st t ttre, cluster(st) att     0.6431  0.4407    0.685   4/16 = .250
  ... + esprelength(3) espostlength(3)
        esfixedbaseperiod jointtests       0.6064  0.6300    1.039   0/16 = .000
  corr(b_noES, b_ES) = 0.834
0.685 sits inside the help's own measured range "roughly 0.66 to 0.75"
(sthlp:255) for the DEFAULT spec, and the understatement DISAPPEARS once event
study horizons with a fixed base period are requested.  The natural mechanism:
with the pooled cohort-varying reference (the default, PS-34) the att carries a
variance component the cluster-robust covariance cannot see at cluster()==i;
with esfixedbaseperiod each coefficient is a difference from one reference
period and the covariance sees it.  If that holds up at 100 draws it locates the
documented defect much more precisely than "at or above the i level".

## d03 draw 1 (reps=399, att only, both engines, same data, same rseed)
  b_obs 0.7696
  sd of the replication draws:  standard 0.4401   bayes 0.3664   RATIO 0.832
  centring: mean(std draws) - b_obs = -0.0149 ; mean(bayes) - b_obs = +0.0431
  percentile p: standard .0402   bayes .0150
  reported se  : standard 0.4401  bayes 0.3664
  95% pctl CI  : standard [0.0083, 1.6759]  bayes [0.1197, 1.6436]
Theory for Dirichlet(1,..,1) over G clusters vs the multinomial pairs bootstrap
puts the variance ratio at G/(G+1) = 0.909, i.e. an SD ratio of 0.953.
MEASURED 0.832 at G=10 -> variance ratio 0.69, not 0.91.  ado:6498 tells the user
the Bayesian bootstrap "shares the pairs cluster bootstrap's first two moments".

## The few-cluster caution never fires on my design (ado:6491-6499)
The bootstrap's own "badly sized with few clusters (MacKinnon and Webb)" caution
is keyed on the count of TREATED resampling clusters being < 6.  My design has
G = 10 clusters of which 6 are ever treated, so r(N) = 6 and the caution does
NOT print.  A user on a 10-state staggered rollout therefore gets a
bootstrap(bayes) 95% interval that covers ~67-73% of the time with no warning of
any kind, while sthlp:255 is pointing them at bootstrap() as the remedy for the
analytic understatement.  The guard counts treated clusters; the few-cluster
literature's concern is the TOTAL cluster count.

## CONTRADICTION of a section-7 entry (worth reporting carefully)
known_bydesign section 7 records: "bootstrap(standard) at cluster() == i
OVER-covers: att_i 149/150 = .9933 vs nominal .95 (c8c4e2c, 2026-09-16)".
I measure the OPPOSITE DIRECTION on a neighbouring design: att_i percentile
95% CI excludes the true zero 17.6% of the time (coverage .824), G=10,
cluster(st) == i.  The designs differ - c8c4e2c's P5 had ONE state per cohort
("every cohort is its own resampling cluster") with cohort heterogeneity and
reps(199); mine has TWO states per cohort, a pure sharp null and reps(79).
So the right reading is NOT "the earlier measurement is wrong" but "the
over-coverage is a property of that particular geometry, not of
bootstrap(standard) at cluster()==i in general": move one state into each
cohort and the same engine goes from .993 coverage to .824.

## WATCH ITEM (not a finding yet)
At 18 draws the permute ECDF sits about 1.5x the 45-degree line at EVERY alpha
(.06/.17/.28/.39/.56/.78 against .05/.10/.20/.30/.50/.70).  A uniform excess
across the whole CDF is a different signal from a 5%-only excess, and the
theorem says the CDF should be AT OR BELOW the line everywhere.  At 18 draws
this is inside MC error (SE at .10 is .07).  RE-CHECK AT 50 DRAWS and say so
explicitly either way; do not rationalize it later.

## FRAMING CORRECTION (advisor, 00:42)
ado:6498's "shares the pairs cluster bootstrap's first two moments" sits inside
`if r(N) < 6' on the TREATED cluster count, so on my design (6 treated of 10) it
NEVER PRINTS.  It is therefore NOT a documented-vs-actual finding of its own.
The correct doc/guard anchors for MCBP-1 are:
  (a) sthlp:255 names "bootstrap() or permute()" as THE remedies for the
      at-or-above-i understatement, and bootstrap() is not one here;
  (b) the few-cluster guard keys on the TREATED cluster count (< 6), so a
      G=10 panel with 6 treated and 4 control - squarely inside the
      few-clusters literature - is warned about nothing at all.
The 0.83 SD ratio is the MECHANISM for MCBP-1, not a separate finding.
Lead with the ranking on identical data with a bitwise-identical point estimate:
  permute exact  >  standard 0.83x truth  >  analytic 0.69x  ~  bayes 0.69x.

## d03 COMPLETE (3 draws at reps(399), att only, same data and rseed per draw)
 draw  b_obs    sd_std  sd_bay  ratio   mean_std-b  mean_bay-b   p_std   p_bay
   1   0.7696   0.4401  0.3664  0.832    -0.0149      0.0431    0.0402  0.0150
   2   0.0542   0.5011  0.3817  0.762     0.0208     -0.0193    0.8766  0.9073
   3   0.1863   0.6901  0.5820  0.843     0.0527      0.0161    0.7960  0.7519
MEAN SD RATIO 0.812.  Theory (Dirichlet(1,..,1) vs multinomial over G=10
clusters) gives a variance ratio G/(G+1) = 0.909, i.e. an SD ratio 0.953.
CENTRING is fine for both engines (|mean draw - b_obs| <= 0.05 against SDs of
0.37-0.69), so the whole difference is DISPERSION, not location.
The pairs bootstrap's extra spread is the discreteness: with cohorts of two
states a resample can drop a whole cohort, which the Dirichlet never does.
Also visible: the "x" marks are FAILED standard replications - 1, 2 and 2 of 399.
bootstrap(bayes) never failed a replication.

## repro validation (rt_engine.log, 2 draws)
repro_engine_size.do runs clean.  It also confirms the load-bearing fact for
MCBP-1: the POINT ESTIMATE is identical across engines -
   |b_analytic - b_bootstrap(standard)| = 0 exactly
   |b_analytic - b_bootstrap(bayes)|    = 0 exactly
   |b_analytic - b_permute|             = 1.67e-16 max (one ulp; the permute
                                          fast-mode comparator refit)
so every difference in the reported inference is dispersion, not location.
(Bug found and fixed while validating: str10 truncated "boot_standard" to
"boot_stand".  Now str16.)

## d06 launched 00:47 -- G=40 engine comparison (20 draws, att + att_i)
Purpose: is the bayes over-rejection a FEW-CLUSTER phenomenon or general?
Same four engines, same DGP shape, 40 states instead of 10.

## CORE RESULT AT 50 DRAWS (01:05) -- these supersede every interim number above
G=10, i = state = cluster(st), TRUE SHARP NULL, permreps(79)/reps(79).
Rejection of the true null at nominal 5% (Wilson 95% in brackets):
  permute(shuffle)            att   .061 [.021,.165]   att_i .082 [.032,.192]
  bootstrap(standard) pctl/bc/wald  .122 [.057,.242]   (same for att and att_i)
  bootstrap(bayes)    bc            .204 [.115,.336]
  bootstrap(bayes)    pctl          .163 [.085,.290]
  bootstrap(bayes)    wald          .122 [.057,.242]
  ANALYTIC cluster se         att   .180 [.098,.308]   att_i .200 [.112,.330]
  ANALYTIC, ES spec (A9)      att   .061 [.021,.165]   att_i .082 [.032,.192]
  ANALYTIC joint pre-trend test     .755 [.619,.854]   (NF-1b published .657)
  ANALYTIC adversarial line test    .286 [.178,.424]
Reported dispersion against the TRUE MC sd of the estimator (0.5291 over 50 draws):
  analytic          0.4067 = 0.769 x truth   mean 95% CI width 1.840
  bootstrap(bayes)  0.4130 = 0.781 x truth   mean 95% CI width 1.686
  bootstrap(standard)0.4888 = 0.924 x truth  mean 95% CI width 2.031
  analytic ES spec  0.5791 = 1.011 x ITS OWN truth (sd 0.5728)
  a correctly sized 95% interval needs width 2.074
WATCH ITEM RESOLVED: the permute ECDF at 50 draws is .06/.12/.18/.27/.41/.65
against .05/.10/.20/.30/.50/.70 -- at or below the 45-degree line from alpha=.2
up.  The apparent uniform excess at 18 draws was noise.  permute() is CLEAN.
HEADLINE MUST BE RESTATED: at 12-18 draws bayes looked like .33-.42; at 50 draws
it is .16-.20.  The defensible statement is that bootstrap(bayes)'s standard
error is statistically indistinguishable from the ANALYTIC one it replaces
(0.781 vs 0.769 of truth) while bootstrap(standard) is materially better (0.924),
and that only permute() attains nominal size.

## d05 -- WHICH LEVER fixes the analytic att SE?  ANSWER: esfixedbaseperiod.
25 draws, G=10, cluster(st)==i, sharp null, plain (analytic) runs:
  S0  att only                              sd(b) 0.5474  se 0.3904  se/sd 0.713  rej .240
  S1  + esprelength(3) espostlength(3)      sd(b) 0.8085  se 0.6930  se/sd 0.857  rej .080
  S2  + esfixedbaseperiod ALONE             sd(b) 0.5224  se 0.4869  se/sd 0.932  rej .080
  S3  + window AND esfixedbaseperiod        sd(b) 0.5224  se 0.4869  se/sd 0.932  rej .080
S2 and S3 are IDENTICAL to the printed digits: once esfixedbaseperiod is on, the
event-study window changes nothing about the main att.  So the lever is
esfixedbaseperiod, not the window.  Under the DEFAULT pooled cohort-varying
reference the att estimator is BOTH noisier (sd .5474 vs .5224) AND has a
standard error understated by ~29%; under a fixed base period it is tighter and
the standard error is nearly right.
(25 draws -> the sd estimates carry about +/-15%; treat 0.713 vs 0.932 as the
signal and the third digit as noise.)

## v02 verifications (log_v02.log) -- all CLEAN
* e(perm_p_<est>) equals the PRINTED table's p BITWISE on every draw, under the
  default convention AND under permtstat.  My p harvest is the reported number.
* e(perm_lb_<est>) does NOT equal the table's lb - the documented null-band alias.
* PO design accepted: permonsetcohorts with i=cohort and cluster(st) finer than i
  gives e(perm_irebuilt)=1, e(perm_npossible)=369600, 79/79 completed, 0 failed.
* EX design: e(perm_exact)=1, e(perm_npossible)=168, e(perm_completed)=168,
  e(perm_failed)=0 -- the enumeration is COMPLETE, so the achievable level really
  is 2*4/168 = .0476 and the p is a plain share (p_att = .92857143 = 2*78/168).
* PO: e(perm_comparator_b) == e(b) exactly (.25352478 / .32194493), so the
  observed statistic is computed by the same functional as the replications
  and the exactness theorem applies to the PO arm as stated.

## d06 -- THE SECOND CLUSTER COUNT, G=40 (20 draws, same DGP shape)
Rejection at nominal 5% of a TRUE sharp null, ALL FOUR ENGINES, att and att_i:
   analytic .05   permute .05   bootstrap(standard) .05   bootstrap(bayes) .05
   (1 of 20 in every cell; the CI-exclusion rate is the same 1 of 20)
Reported dispersion against the true MC sd:
   att   (true sd 0.2552): analytic 0.897x   bayes 1.010x   standard 1.067x
   att_i (true sd 0.2774): analytic 0.842x   standard 0.971x
Point estimates again bitwise identical across engines (max |diff| = 0 for the
bootstraps, 1.11e-16 for permute).
=> THE BOOTSTRAP ENGINES ARE FINE AT G=40.  Everything in the G=10 result is a
FEW-CLUSTER phenomenon.  The analytic standard error stays understated at both
cluster counts (0.84-0.90 at G=40, 0.77 at G=10), which is the documented defect;
what is new is that bootstrap(bayes) does NOT repair it at G=10 while
bootstrap(standard) partly does, and that both repair it at G=40.
CAVEAT: 20 draws, so each G=40 rejection rate is 1/20 with a Wilson 95% of
[0.009, 0.236] - this rules OUT a gross distortion, not a 2-3 point one.

## A4 / A6 AT 40 DRAWS (01:32) -- both seeded "never measured" arms are CLEAN
G=10, cluster(st)==i, true sharp null, 40 draws:
  A4 permtstat        att   .075 [.026,.199]   att_i .050 [.014,.165]
  A6 boottstat (tpctl) att  .025 [.004,.129]   att_i .025 [.004,.129]
Replication completion: permute 79/79 on every draw; boottstat 77-79 of 79
(up to 2 replications failed, which is the documented degenerate-resample rule).
THE STRIKING COMPARISON: on the SAME design the NON-studentized bootstrap
families reject .122 while boottstat - the same engine, studentized - rejects
.025.  So the percentile-t family IS the bootstrap family that survives few
clusters (as MacKinnon-Webb theory predicts), and the default
bootdisplay(percentile) is not.  That is a practical recommendation the help
does not make: at few clusters use boottstat or permute(), not the default
percentile or bias-corrected families, and not bootstrap(bayes).

## A10 / A11 AT 20 DRAWS (01:49) -- THE 2026-09 REWIRE OF THE JOINT TEST WORKS
Same DGP, same G=10, same true sharp null, jointtests + adversarial with
esprelength(3) espostlength(3) esfixedbaseperiod:
  ANALYTIC joint pre-trend test (A9, 49 draws)      .755 [.619,.854]
  ANALYTIC adversarial line test (A9, 49 draws)     .286 [.178,.424]
  PERMUTE-rewired joint test    (A10, 20 draws)     .000 [.000,.161]
  PERMUTE-rewired adversarial   (A10, 20 draws)     .000 [.000,.161]
  BOOTSTRAP-rewired joint test  (A11, 20 draws)     .050 [.009,.236]
  BOOTSTRAP adversarial         (A11)               ABSENT on all 20 draws
     (documented: sthlp:759 - under bootstrap() the adversarial p is produced
      only with boottstat; without it e(boot_adv_absent) is set)
So the September 2026 rewiring of the joint and adversarial tests to
permutation / fixed-metric rank statistics REMOVES a 75% false-positive rate
and lands at or below nominal.  As far as the project's record goes this had
never been measured - F8 in the recon notes established that the p's ARE
permutation/bootstrap quantities, not that they are correctly SIZED.
SECOND OBSERVATION: A11's coefficient rows reject .05-.10 where A5 (the same
bootstrap without the ES options) rejects .122.  A11 carries esfixedbaseperiod.
That is the same lever d05 isolated: esfixedbaseperiod fixes the att standard
error, and with it every engine's size on the att row.

## A10 / A11 FINAL, 40 DRAWS (02:15)
  A10 permute joint pre-trend    .050 [.014,.165]   (analytic on the same DGP: .755)
  A10 permute adversarial        .000 [.000,.088]   (analytic: .286)
  A10 permute att / att_i rows   .050 / .050
  A11 bootstrap joint pre-trend  .025 [.004,.129]
  A11 bootstrap adversarial      ABSENT on all 40 draws (documented)
  A11 bootstrap att rows: pctl .075  bc .075  wald .025 ; att_i pctl .050 bc .075 wald .050
Every rewired test is at or below nominal.  The rewire is a genuine repair of a
75-percentage-point false-positive problem, and it had not been size-tested.

## *** POSSIBLE THEOREM VIOLATION -- PO DESIGN, 18-19 DRAWS (02:27) ***
PO design: 12 states, 3 onsets x 3 states + 3 never, one row per state-year,
sharp null, permute(shuffle) permreps(79).  BOTH arms:
  POonc (i = cohort, cluster(st) finer than i, permonsetcohorts)  .158 [.055,.376]
  POctl (i = st = cluster(st), the ordinary route)                .167 [.058,.392]
  ncomp 79/79 on every draw, 0 failures.
Randomization inference is EXACT on this design, so anything materially above
.05 is a defect -- and the Wilson lower bound already clears .05, just.
BUT: 18 draws.  Do NOT call this until 60+.  Note the control arm shows it too,
so if it survives it is NOT about permonsetcohorts; it is about permute() on a
12-state design with 3 never-treated, where my G=10 design (4 never, cohorts of
2) gave .061 over 50 draws.
NEXT IF IT SURVIVES: vary one thing at a time - never-treated count, cohort
size, cluster count - to locate it.

## PO at 40 pooled draws (02:40): att .100 [.040,.231], att_i .125 [.055,.261]
Elevated (2-2.5x nominal) but the Wilson lower bound on att is now BELOW .05.
IMPORTANT REFRAME: rather than treating each permute arm separately, POOL them.
The theorem applies to every permute arm, so the strongest statement is the
pooled rejection rate across all independent permute draws.  Running tally of
att-row rejections at nominal 5%:
   A1  permute default, G=10          3/49
   A4  permtstat, G=10                3/40
   A10 permute + ES + esfixedbase     2/40
   PO  permonsetcohorts, 12 states     4/40
   POctl permute, 12 states            4/40   (same data as PO, NOT independent)
Pooling the INDEPENDENT ones (A1, A4, A10, POonc): 12/169 = .071, and under a
true .05 the chance of 12 or more is about .15 -- CONSISTENT WITH NOMINAL.
Do the pooled calculation properly at write-up time.

## PO AT 63 POOLED DRAWS (02:48) -- THE SIGNAL IS GONE.  NOT A DEFECT.
  POonc  att .0635 [.0250,.1522]   att_i .0794 [.0344,.1727]
  POctl  att .0635 [.0250,.1522]   att_i .0794 [.0344,.1727]
The 16.7% at 18 draws and 13.3% at 30 draws were small-sample noise; at 63 the
Wilson lower bounds are .025 and .034, comfortably below .05.  The pre-registered
decision rule (defect only if the Wilson lower bound clears nominal) is what kept
this from being written up as a theorem violation at 30 draws.  RECORD THAT.
permonsetcohorts (T1, shipped 2026-09-16, never size-tested) is CLEAN, and so is
the ordinary route on the same data.

## PO at 100 pooled draws (02:57) -- superseded by the 140-draw line below
  POonc (permonsetcohorts, i=cohort, cluster(st) FINER than i, e(perm_irebuilt)=1)
        att .0396 [.0155,.0974]   att_i .0600 [.0278,.1248]
  POctl (ordinary route, i=st=cluster on the same data)
        att .0400 [.0157,.0984]   att_i .0606 [.0281,.1260]
  79/79 replications completed on every one of the 100 draws, 0 failures.
The finer-than-i permutation machinery shipped 2026-09-16 delivers nominal size
on a design where exactness is a theorem.  This is the first size measurement of
that feature; the project deferred its assessment suite to 2026-09-17.

## *** PO FINAL, 140 POOLED DRAWS (03:21) -- CLEAN ***
  POonc (permonsetcohorts; i=cohort, cluster(st) finer than i, perm_irebuilt=1)
        att   .0571 [.0292,.1087]     att_i .0714 [.0393,.1265]
  POctl (ordinary route on the same data, i = st = cluster)
        att   .0571 [.0292,.1087]     att_i .0714 [.0393,.1265]
  79/79 replications on every one of the 140 draws; 0 failures; 369600 possible
  assignments.  Both Wilson lower bounds are below the .05 nominal, so by the
  pre-registered rule this is NOT a defect.

## *** EX FINAL, 80 POOLED DRAWS (04:01) -- CLEAN, in fact CONSERVATIVE ***
permexact on 8 states (2 onset 2004, 1 onset 2006, 5 never), M = 168 assignments,
enumerated COMPLETELY on every draw (e(perm_exact)=1, e(perm_completed)=168,
e(perm_failed)=0), achievable nominal 2*4/168 = .0476:
  EXdef (equal-tailed plain share) att .0125 [.0022,.0675]  att_i .0250 [.0069,.0866]
  EXabs (absolute-value share)     att .0250 [.0069,.0866]  att_i .0500 [.0196,.1216]
All at or below nominal.  The default convention is on the conservative side
(P(X<=1 | n=80, p=.0476) = .10, so not a significant shortfall, just safe).

## HVF FINAL, 30 DRAWS (04:13) -- variance heterogeneity on a FIXED cluster
## attribute, independent of the onset draw, so the SHARP NULL HOLDS and
## exactness is still a theorem for the permutation arms
G=10, treated-or-not irrelevant: the even-numbered ORIGINAL state ids carry
error sd 3, the odd ones sd 1.  Rejection at nominal 5%:
  permute default      att .000 [.000,.114]   att_i .033 [.006,.167]
  permtstat            att .033 [.006,.167]   att_i .067 [.019,.213]
  bootstrap(standard)  att .100 [.035,.256]   att_i .067 [.019,.213]
  bootstrap(bayes)     att .100 [.035,.256]
  ANALYTIC             att .167 [.073,.336]   att_i .133 [.053,.297]
So under cluster-level variance heterogeneity that is NOT caused by treatment,
the permutation engine stays exact (as the theorem requires) while the analytic
and bootstrap engines carry the same G=10 over-rejection as in the homogeneous
baseline.  This is the clean half of the MacKinnon-Webb question.

## RUN CLOSED 04:34.  Final row counts:
  mcres_g10core.csv 983  (50 draws x A0 A9 A1 A5 A7)
  mcres_g10t.csv    161  (40 draws x A4 A6)
  mcres_g10jt.csv   481  (40 draws x A10 A11)
  mcalt_po.csv 241 + mcalt_po2.csv 321 (140 pooled PO draws)
  mcalt_ex.csv 161 + mcalt_ex2.csv 161 ( 80 pooled EX draws)
  mcalt_hvf.csv 301 (30 draws)   mcalt_hv.csv 226 (22-23 draws)
  d06_g40_engines.log (20 draws at G=40, four engines)
  d03 (3 draws at reps=399), d05 (25 draws, ES lever), v02, d01/d02, d04
DROPPED for time, and why:
  * permabsp / permrawp as their own MC arms.  permrawp is an EXACT arithmetic
    transform of permabsp ((80*p_abs - 1)/79, verified on live runs) and the two
    share the 5% rejection region at permreps(79), so the arm adds a 1%-level
    distinction only.  The conventions were verified live in v01.
  * permside(), bootdisplay(bca), bootalpha(size), ties, bootdf(normal).
  * G=40 for permtstat / boottstat / the joint tests (only the four main engines
    were run at G=40, via d06).
  * A larger G=40 sample: 20 draws only, so [0.009,0.236] per cell.

## REPRO FILES, ALL SYNTAX-VALIDATED ON 2 DRAWS BEFORE DELIVERY
  repro_engine_size.do     six engines at G=10 (analytic, permute, permtstat,
                           bootstrap standard, boottstat, bootstrap bayes),
                           prints rejection, CI-exclusion, mean and MEDIAN 95%
                           interval width, se/true-sd ratios, and the
                           bitwise-identical point-estimate check.  ~60 min at
                           NDRAW 60.  (rt2.log = its 2-draw validation.)
  repro_joint_size.do      joint + adversarial under analytic / permute /
                           bootstrap.  ~35 min at NDRAW 40.  (rt3.log)
  repro_perm_exactness.do  EX (permexact, 168 enumerated) + PO
                           (permonsetcohorts) + POC control.  ~45 min at
                           NDRAW 40.  (rt4.log)
  d05_es_lever.do          which lever fixes the analytic att se.  ~1 min.
  d03_bayes_vs_standard.do the Dirichlet-vs-pairs replication sd.  ~10 min.
  d02_estimand_truth.do    the ES/estimand negative result on an exact DGP. ~1 min.
  d06_g40_engines.do       the G=40 arm.  ~15 min.

## 2026-09-17 04:41 -- d02 RERUN, closing the last open premise of MCBP-ESFBP

The advisor flagged that d02's original spec list (none / w33 / w99) never included
`esfixedbaseperiod`, which d05 had since identified as THE lever that fixes the
analytic att SE.  Under a sharp null one cannot distinguish "different estimator of
the same estimand" from "different estimand", so the claim had an untested premise.

Patched d02_estimand_truth.do to add two specs and reran on the NOISELESS DGP
(where the design truth is known in closed form):

  TAUMODE=const   true att = 1
    none      main:att = +1.0000000000002X+000
    w33       main:att = +1.0000000000004X+000
    w99       main:att = +1.0000000000004X+000
    efbp      main:att = +1.0000000000000X+000     <-- NEW
    w33efbp   main:att = +1.0000000000000X+000     <-- NEW

  TAUMODE=ramp    true att = 1.5555556
    none      main:att = +1.8e38e38e38e38X+000
    w33       main:att = +1.8e38e38e38e40X+000
    w99       main:att = +1.8e38e38e38e40X+000
    efbp      main:att = +1.8e38e38e38e34X+000     <-- NEW
    w33efbp   main:att = +1.8e38e38e38e34X+000     <-- NEW

VERDICT: `esfixedbaseperiod` does NOT change the att ESTIMAND.  All five specs hit
the design truth to ~1e-15.  So MCBP-ESFBP stands as written: same estimand, two
estimators, and only the esfixedbaseperiod one gets a correct standard error.
