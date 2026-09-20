# pretrendpower(): does the reported minimum detectable trend have the stated power?
# Method: SIMULATION AGAINST THE CLAIM.  Agent label: PTP.

Status: COMPLETE.  Everything below is measured, with the do-file and log named.

## 0. What the displayed numbers are (verified exactly, not a finding)
`p02_algebra.do` / `.log`: hand-recomputing the whole of wooldid_ptp from the posted
`e(b_es_att)`, `e(V_es_att)` and `e(df_r)` reproduces `e(roth_slope_att)`,
`e(roth_slope_se_att)`, `e(roth_W_att)`, `e(roth_prep_att)`, `e(roth_power_att)[.,mdt_wald]`
and `[.,bias_wald]` to 1-2 ulp.  So the algebra is:
  m_j   = relyr_j - esrelativeto
  c2    = m'S^-1 m,  slope = m'S^-1 b_pre / c2,  se = 1/sqrt(c2)
  W     = b_pre' S^-1 b_pre,  p = Ftail(keff, df_r, W/keff)
  mdt_wald(P) = sqrt( npnF(keff, df_r, invF(keff,df_r,.95), 1-P) / c2 )
  bias(P)     = mdt(P) * (hbar - esrelativeto)
and `hbar` is the estimand-weighted mean post relative period over the FULL post support
(verified by hand on p01: 352/192 = 1.8333).
That algebra is the RIGHT algebra *if* S is the true sampling covariance of b_pre and the
null distribution of W/keff really is F(keff, df_r).  Both are testable only by simulation.

(placeholder - findings appended below as they are measured)

## 1. NEGATIVE RESULT (clean): the mvnormalcvqp vech ordering
wooldid hands `vech(S)'` (Mata's vech is COLUMN-major) to `mvnormalcvqp(L,U,M,V,q)`.
Stata's help never states the ordering, so a row-major consumer would silently
mis-assemble S for k>=3 pre-periods.  Decided exactly from the p01 posted S:
  correct S with t-crit bounds -> 0.12259799   (scipy oracle)
  row-major misreading         -> 0.09631734
  wooldid prints                  0.1225744646839422
=> ordering CORRECT.  The 2.4e-5 gap is the 1000-point quadrature.  CLEAN.
Repro: p11_vech.py (and TAG_A11 in p06_struct.do, an in-Stata brute-force check).

## 2. NEGATIVE RESULT (clean): DESIGN A -- the exact geometry
Balanced staggered panel, 80 units (48 treated in 3 cohorts, 32 never treated),
10 periods, iid homoskedastic normal errors, i = unit, no cluster() (e(vce)="ols"),
esfixedbaseperiod, esprelength(3) -> k = 3 pre coefficients, df_r = 279.
300 replications, delta = 0  (p04_simA.do, A0_null.csv):
    joint pre-trend F test rejects at 5% : 11/300 = 0.0367  CP95 [0.0184,0.0647]
    NIS pre-test rejects                : 40/300 = 0.1333  CP95 [0.0970,0.1771]
        (wooldid REPORTS this size as 0.1226 -- inside the interval)
    per-coefficient |b/se| > t_crit     : 0.0567 / 0.0467 / 0.0500 (nominal .05)
    Roth slope 95% CI covers the truth  : 288/300 = 0.9600 CP95 [0.9312,0.9792]
    median reported mdt_wald(80)        : 0.239843   (sd across replications 0.0096)
So on the geometry where the F test is exact by construction, the machinery is
correctly sized and the slope interval has its nominal level.  Nothing to report.

### 2b. DESIGN A, the power arm -- the MDT is EXACT here
300 replications carrying a linear pre-trend of slope exactly
median(mdt_wald(80)) = 0.2398425 (p04_simA.do, A1_pow80.csv):
    joint pre-trend F test rejects at 5% : 240/300 = 0.8000  CP95 [0.7502,0.8438]
    CLAIMED POWER                        : 0.80
    NIS pre-test rejects                 : 267/300 = 0.8900
    Roth slope 95% CI covers the true 0.2398 : 288/300 = 0.9600
Also: mdt_wald is IDENTICAL rep-for-rep between the null and the power arm
(53/99 bit-identical, the rest differing at 1e-15).  That is the saturation
property of i = unit -- the treated cells are fitted exactly, so s^2 comes from
the never-treated rows alone and the planted trend cannot move it.
=> the algebra and my harness both check out.  Whatever shows up on the other
   geometries is a property of S, not of the formula.
### 2c. DESIGN A, the 50% arm
300 replications at slope = median(mdt_wald(50)) = 0.1743356 (A2_pow50.csv):
    joint pre-trend F test rejects at 5% : 150/300 = 0.5000  CP95 [0.4420,0.5580]
    CLAIMED POWER                        : 0.50
Two points on the power curve, both on the nose.  DESIGN A IS CLEAN.

## 3. NEGATIVE RESULTS from the deterministic structural batch (p06_struct.do,
##    p06b_struct2.do; logs alongside)
- A0  planted trend of slope .25, noise 1e-8: e(roth_slope_att) = .25 exactly.
      The premise the power claim rests on (a linear pre-trend of slope d puts
      d*m into the pre-period coefficients) holds.
- A1  esrelativeto(-2), where the reported pre relyrs STRADDLE the reference so
      m carries both signs: slope = .25 exactly.  Sign construction correct.
- A8  BIENNIAL panel (t = 2000,2002,...): b_es = (pre6 -1, pre4 -.5), i.e.
      .25*m with m = (-4,-2) in t units; slope = .25, hbar = 4.2.  m is built
      in t units, correctly.
- A2  anticipation(1): e(roth_hbar_att) = 2.2666667, which is the mean post
      relative period on the SHIFTED onset (the raw-onset value is 1.8333333).
      Consistent with the event-study labels.  Correct.
- A4  the joint pre-trend p in the Roth block (.12418582) is the same number
      jointtests prints (0.1242).  One statistic, two entry points: consistent.
- A5  mdt_wald/mdt_nis with a planted trend of 0 vs 3.0 agree to ~1e-6
      (design-only dependence, as sthlp:791 claims).
- A6  esprelength(9) -> k=7: mdt_nis and the NIS size are omitted with the
      documented note; mdt_wald still computed.  Correct.
- A11 in-Stata brute force: mvnormalcvqp with vech(S)' gives 0.87134557 against
      a 400,000-draw Monte Carlo of 0.87141250; the row-major ordering returns
      missing.  Ordering CORRECT.
- B12 bootstrap(standard) + cluster(st): e(df_r) survives as 15 and the Roth
      block uses F(3,15) - consistent with the rest of the run.  CLEAN.
- B13 customvce(cluster st): df_r = 15, diagnostics computed.  CLEAN.

================================================================================
## FINDING PTP-1  (severity HIGH, confidence CERTAIN, seeded: YES -- the seed
## asked for a measurement "at a clustering geometry where the test is sound AND
## at one where it is not"; the design, the size-adjusted power, the inversion
## mechanism and the G-dependence are mine)
## The minimum detectable trend does not have the stated power on a canonical
## few-cluster DiD -- and wooldid prints no caution on that geometry.
================================================================================
DESIGN E: 16 states x 10 years, state-level outcome, policy adopted by 12 of the
16 states in three waves, AR(1) errors within state.
  wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3)
                        espostlength(2) pretrendpower(50 80)
  i = onset-cohort code (4 codes) with unit (state) effects supplied through
  fe() -- which is exactly the design sthlp:1334 recommends for pretrendpower:
  "a coarser i (cohorts, or treatment-onset groups) with unit fixed effects
  supplied through fe() is the design that satisfies them" -- plus the
  clustering any applied paper would use.  The at-i caution (ado:1413-1444)
  fires only when i is NESTED inside the cluster variable; here a cohort spans
  four states, so the predicate does not catch this geometry and NO caution is
  printed.  (Note: cluster(st) is finer than i but coarser than the cohort-by-
  period CELL, so it is not the "strictly finer than the cohort-by-period cells"
  remedy the caution's text names either -- it simply falls outside the screen.)
  df_r = 15, k = 3 pre-period coefficients.

MEASURED (p08_simE.do; E0_null.csv 400 reps, E1_pow80.csv 300 reps, E3_nis80.csv
300 reps; analysis in analyze.py and the inline scripts recorded in NOTES.md):

  (a) the test mdt_wald is supposed to power is OVERSIZED  [the robust core]
        joint pre-trend F(3,15) rejects at 5%, pooled over 600 null replications
        from two independent seed streams (E0_null.csv 400 + the self-contained
        repro's own 200):  63/600 = 0.1050  CP95 [.0816,.1323]
        (E0 alone 0.1025; the repro alone 0.1100)
        The size-corrected 5% cutoff is 4.4483 against a nominal F(3,15)
        critical value of 3.2874, i.e. 1.35x (bootstrap 95% CI on the cutoff
        [3.85, 4.90]).

  (b) at the REPORTED mdt_wald(80) = 0.5463 the achieved power is NOT the 0.80
      an honest 5% test would give
        with the nominal cutoff       : pooled 387/500 = 0.7740 CP95 [.7348,.8099]
                                        (E1 alone 226/300 = 0.7533; the repro's
                                        own 200-rep arm 161/200 = 0.8050)
        with the size-corrected cutoff: 180/300 = 0.6000  CP95 [.5421,.6559]
        CLAIMED: 0.80 at the 5% level.
      HONEST READING: the nominal-cutoff figure (0.774) is NOT decisively below
      0.80 - it is delivered by a test that rejects 10.5% of the time under the
      null.  The trend wooldid calls "detectable with 80% power by a 5% test" is
      detectable with 60% power by an honest 5% test.

  (c) mdt_nis is worse, and its shortfall is unambiguous even at the NOMINAL
      cutoff - this is the single cleanest number in the report, because it
      needs no size correction and so carries only binomial uncertainty
        at the REPORTED mdt_nis(80) = 0.42375
        nominal cutoff      : 214/300 = 0.7133  CP95 [.6586,.7638]  <- excludes .80
        size-corrected      : 157/300 = 0.5233  CP95 [.4652,.5810]

  (d) the Roth slope interval does not have its nominal level
        95% CI covers the true slope: 362/400 = 0.9050  CP95 [.8719,.9319]
        mean reported slope se 0.1437 vs the MC sd of the slope 0.1822 (se/sd .79)

  (e) the reported MDT is itself a noisy statistic, not a design constant
        mdt_wald(80) across the 400 null replications: median 0.5463, sd 0.1251,
        i.e. a coefficient of variation of 23%.  sthlp:791 says the MDT "depends
        only on the design's precision, never on the estimated coefficients",
        which reads as a design constant; it is proportional to an estimated
        scale with 15 degrees of freedom.

WHY, structurally (this is the part that is new):
  The per-coefficient inference on this design is, if anything, slightly
  CONSERVATIVE.  The individual |b_pre/se| > t_crit(15) rates are
  0.039 / 0.039 / 0.047 against a nominal 0.05, and the mean reported variances
  are 1.05-1.09 times their Monte-Carlo counterparts.  That makes the inversion
  story sharper: the marginal standard errors are not too small, yet the
  quadratic form over-rejects at twice its nominal size.  What breaks is the QUADRATIC FORM: pretrendpower() is built
  entirely on b'S^-1 b and m'S^-1 b / m'S^-1 m, and it inverts a cluster-robust
  S estimated on 15 degrees of freedom.  Inverting a noisy S is what produces
  both the 10% size and the 21% understatement of the GLS slope's standard
  error.  Nothing in the run tells the user this: the at-i caution does not
  fire, and the printed "F(3, 15)" line asserts the reference distribution
  that the measurement says is wrong.
  => This is NOT a rediscovery of the documented "analytic SEs understated at
  cluster() at or above i" (known_bydesign S1).  Here cluster() is strictly
  finer than i, the marginal SEs are correct, and the defect lives in the
  inversion.

================================================================================
## FINDING PTP-2  (severity MEDIUM, confidence CERTAIN, seeded: NO)
## e(roth_nissize_*) -- documented as "the actual size of the any-individually-
## significant pre-test" -- evaluates a NORMAL probability at a t critical value,
## so with few clusters it reports the pre-test as far more conservative than it is.
================================================================================
ado:10366 sets  ptpniscrit = invttail(df_r, .025)  and ado:10525/10589 feed that
bound to mvnormalcvqp(), i.e. a MULTIVARIATE NORMAL probability.  The test being
described rejects when |b/se| > invttail(df_r,.025), and b/se is referred to
t(df_r) everywhere else in wooldid (it is the p-value the event-study table
prints).  The two are the same only as df_r -> infinity.

Cleanest case, k = 1 (p06b_struct2.do, TAG_B15; 16 states, cluster(st), df_r=15,
esprelength(1)):  with ONE pre-period coefficient the "any individually
significant" test IS a single two-sided 5% t test, so its actual size is 0.05
by construction -- no simulation needed.
    wooldid prints  "size of the any-pre-period-individually-significant
                     pre-test: 0.0331"
    e(roth_nissize_att) = 0.03305212 = 2*normal(-invttail(15,.025)) exactly
    the true size       = 0.05
    understatement      = 34%
Same run, the two detectable-trend columns disagree although at k=1 they power
the SAME test:  mdt_wald = 0.77712, mdt_nis = 0.77069 (ratio 0.9917).

Measured at k = 3 on the same design (E0_null.csv, 400 reps):
    reported e(roth_nissize_att) (median over 12 runs) = 0.0781
    realised rejection rate of that pre-test          = 0.1125 CP95 [.0833,.1476]
The reported value sits below the lower confidence limit.

Scale of the error by df_r (2*normal(-invttail(df,.025)) against 0.05):
    df_r =  15 -> 0.0331  (-34%)     df_r =  39 -> 0.0431  (-14%)
    df_r =  23 -> 0.0386  (-23%)     df_r = 279 -> 0.0490  (-2%)
Direction matters: the number is printed to warn the user that the NIS pre-test
OVER-rejects ("its size ... exceeds its nominal level when there are several
pre-periods", sthlp:790).  With few clusters it understates that size, and at
k <= 2 it reports a size BELOW 5% for a test that is exactly 5%.
Repro: repro_nissize.do (self-contained; Part 1 is pure arithmetic, Part 2 a
300-replication Monte Carlo of the realised rate).

### PTP-1, the size of the error, pinned by a second power arm
p13_simE2.do put a trend of 0.6456 on the same design (E4_d6456.csv, 300 reps):
    power with the nominal F(3,15) cutoff      : 270/300 = 0.9000 CP95 [.8603,.9315]
    power with the size-corrected cutoff       : 223/300 = 0.7433 CP95 [.6900,.7918]
Using the POOLED 600-draw size-corrected cutoff (4.4483):
    delta = 0.5463 -> size-corrected power 0.6000 CP95 [.5421,.6559]
    delta = 0.6456 -> size-corrected power 0.7800 CP95 [.7288,.8256]
Model-free linear interpolation through those two points puts 80% power at
delta = 0.6566; a noncentral-chi-squared extrapolation (lambda proportional to
delta^2) gives 0.6567 - the two agree to four figures, so the figure is not an
artefact of the extrapolation model.
    trend detectable at 80% by an HONEST 5% test  ~ 0.657
    wooldid reports                                 0.546
So the reported minimum detectable trend is about 17% too small, and the bias
column (mdt x (hbar - ref)) understates by the same 17%.

================================================================================
## FINDING PTP-3  (severity MEDIUM, confidence CERTAIN, seeded: NO)
## Under poisson, pretrendpower() silently takes the chi-squared shortcut its own
## code comment calls optimistic -- even when cluster() leaves 15 df and wooldid
## has posted e(N_clust).  sthlp:789's description of the test is false there.
================================================================================
wooldid_ptp branches on e(df_r) (ado:10357-10366 and 10545-10548): present -> F
and npnF; missing -> chi2 and npnchi2.  The comment at ado:10543 says "a
chi-squared shortcut is optimistic by ~25 percent at 10 clusters".
Measured (p12_poiscluster.do; 16 states x 10 periods x 12 individuals, i = onset
cohort, fe(st), cluster(st)):
   gaussian path : e(df_r) = 15 -> "joint pre-trend test: F(3, 15)",  mdt 0.14024
   poisson path  : e(df_r) = .  -> "joint pre-trend test: chi2(3)",   mdt 0.02555
                   e(N_clust) = 16 IS posted on that same run
   the F-branch MDT on the same information would be 0.02914
   => mdt_wald is 12.3% too small on the poisson path at G = 16
   => the chi2(3) 5% cutoff has true size 0.0902 when the statistic is F(3,15)
      (and the run printed p = 0.0112 for that joint pre-trend test)
sthlp:789 says mdt_wald "powers the joint pre-trend test, an F test with the
regression's residual degrees of freedom".  On the poisson path it is a
chi-squared test with no degrees-of-freedom correction at all.
Honest framing: the poisson main table also reports z rather than t, so the
chi-squared branch is internally consistent with the rest of that run; what is
specific to pretrendpower() is (i) the help sentence, and (ii) that the MDT
column absorbs the optimism invisibly, where at least the "chi2(3)" label on the
joint-test line is visible.

================================================================================
## FINDING PTP-4  (severity LOW-MEDIUM, confidence CERTAIN, seeded: NO)
## The help's ONLY worked example of pretrendpower() puts the user on the
## geometry the run then says is not a calibrated test.
================================================================================
sthlp:1546, example 7:
  wooldid y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod
                      anticipation(1) pretrendpower(50 80)
Run literally, with the cohort variable actually named i (p06_struct.do TAG_A3):
  "Caution: pretrendpower() was combined with analytic clustered SEs clustered
   at the level of the cohort variable i or above (i is nested within the
   clustering variable i). ... the joint test's covariance is degenerate and its
   p-value should not be read as a calibrated test."
The caution's own source comment (ado:1436) already notes "wooldid's own help
twice recommends the clustering that trips this condition"; the PS-64 ruling
changed the code to caution rather than refuse, but the example was never
changed.  A user copying the documented example gets a warning telling them the
numbers the example was written to demonstrate are not calibrated.

================================================================================
## FINDING PTP-5  (severity MEDIUM-HIGH, confidence CERTAIN, seeded: YES -- the
## seed asked for "a clustering geometry where the test is not sound"; the
## MEASUREMENT of what that does to the MDT is new)
## At cluster() == i the MDT appears to deliver its 80% only because the test
## rejects 28% of the time under the null.
================================================================================
DESIGN B (p09_simB.do): the same 80-unit staggered panel, i = unit, cluster(unit),
AR(1) errors within unit, df_r = 79.  wooldid DOES print its at-i caution here.
  null arm B0_null.csv, 250 reps:
     joint pre-trend F(3,79) rejects at 5% : 70/250 = 0.2800 CP95 [.2253,.3401]
     the empirical 95th percentile of the null statistic is 2.158x the nominal
     F critical value
     Roth slope 95% CI coverage            : 213/250 = 0.8520 CP95 [.8018,.8936]
     mean reported pre-period variances are 0.43-0.68 of their Monte-Carlo
     counterparts; per-coefficient |b/se|>t rates 0.116 / 0.187 / 0.160
     m'S^-1 m: mean-reported 367.8 vs Monte-Carlo 220.3, so the reported MDT is
     0.774x the covariance-implied one
  power arm B1_pow80.csv, 220 reps at the reported mdt_wald(80) = 0.16971:
     with the nominal cutoff      : 175/220 = 0.7955 CP95 [.7360,.8467]
     with a size-corrected cutoff : 101/220 = 0.4591 CP95 [.3919,.5274]
     Roth slope CI covers the true 0.16971: 189/220 = 0.8591
So on this geometry the MDT LOOKS calibrated (0.7955 against a claimed 0.80) and
is not: the apparent power is bought entirely with a 28% false-positive rate.
Self-contained repro: repro_designB.do (no sibling files; p09_simB.do, which
produced the numbers above, sources ptplib.do by absolute path).
An author reading only "achieved 0.80" would conclude the column is fine.  The
(size, power) pair is the only honest summary, and it is (0.28, 0.46).

================================================================================
## PTP-1, THE DISCRIMINATING CONTROL: the same geometry at G = 80 is EXACT
================================================================================
DESIGN C (p07_simC.do): identical construction to Design E -- i = onset-cohort
code, unit effects through fe(), cluster() on the unit, AR(1) errors -- but with
80 clusters instead of 16 (df_r = 79).
  null arm C0_null.csv, 220 reps:
     joint pre-trend F(3,79) rejects at 5% : 10/220 = 0.0455 CP95 [.0220,.0820]
     empirical 95th percentile of the null statistic / nominal crit = 0.969
     Roth slope 95% CI coverage            : 209/220 = 0.9500 CP95 [.9123,.9748]
     mean reported slope se 0.06975 vs MC sd 0.06962  (se/sd = 1.0018)
  power arm C1_pow80.csv, 220 reps at the reported mdt_wald(80) = 0.235774:
     nominal cutoff      : 173/220 = 0.7864 CP95 [.7262,.8386]
     size-corrected      : 178/220 = 0.8091 CP95 [.7508,.8588]
     CLAIMED: 0.80.
CONCLUSION: PTP-1 is a FEW-CLUSTER effect, not a structural error in how
e(V_es_*) is assembled.  Everything is exact at G = 80 and broken at G = 16.
ACTIONABLE SHAPE OF A FIX: the existing screen (ado:1413-1444) keys only on
NESTING (is i inside the cluster variable).  It has no cluster-COUNT leg, so a
16-cluster run passes it silently.  A caution keyed on e(N_clust) - or a
small-sample correction to the joint pre-trend test and to the noncentrality the
MDT inverts - would have caught Design E.  The project's own se_coverage work
was scoped "G = 40 only" and a G = 200 follow-up was ruled SKIP; the direction
that actually matters for pretrendpower() is DOWNWARD in G, and it was not run.

### PTP-2, the self-contained repro was executed (repro_nissize.log)
  PART1 df_r = 15, npre = 1
        e(roth_nissize_att)                 = 0.03305212
        2*normal(-invttail(15,.025))        = 0.03305212   (identical)
        the size of ONE two-sided 5% t test = 0.05
        understatement                      = 33.9%
  PART2 300 replications, realised rejection rate of that same pre-test = 0.0600
        CP95 [0.0359,0.0932] -- consistent with the nominal 0.05 and EXCLUDING
        the reported 0.0331.  (The 0.06 is binomial noise around 0.05, nothing
        more; 300 replications cannot pin a 5% rate closer than that.)

### PTP-1, a second, USER-FACING consequence: one fit, one MDT, and it is noisy
The first run of repro_power.do (kept as repro_power_singledraw.log) read the
MDT off ONE fit, exactly as a user would:
    STEP1 mdt_wald(80) = 0.36167   (the median over 400 fits is 0.5463)
    STEP2 300 replications at exactly that slope: rejection rate 0.4167
          (claimed 0.80); NIS 0.5833
That 0.36167 is a -1.4 sd draw of a statistic with a 23% coefficient of
variation, so the 0.42 conflates the systematic shortfall with one unlucky fit -
but "one fit" is what every user has.  The re-run of repro_power.do targets the
MEDIAN over 25 fits instead, isolating the systematic part.

### PTP-1, the self-contained repro was executed (repro_power.log)
  STEP1a one fit (what a user sees): mdt_wald(80) = 0.36167,  df_r = 15
  STEP1b median over 25 fits of the SAME design: 0.55939
         and the range across those 25 fits is 0.29589 to 0.92505
         -- a 3.1x spread in a number the help describes as depending "only on
         the design's precision"
  STEP2  300 replications carrying a trend of exactly 0.55939:
         joint pre-trend test rejects at 5% : 0.7733   CLAIMED 0.80
         NIS pre-test rejects               : 0.8900
  (matching E1's 0.7533 at 0.5463; the size-corrected figure, which needs the
  null arm, is 0.58 -- see PTP-1 above.)

### A caveat that applies to every "size-corrected" figure in this report
The size-corrected cutoffs are the empirical 95th percentile of the null
statistic over 220-400 replications, so they carry roughly +/- 0.03-0.04 of
extra uncertainty on the power figures beyond the Clopper-Pearson intervals
quoted.  The NOMINAL-cutoff figures are exact binomial and are already decisive
on their own: E3's mdt_nis power is 0.7133 CP95 [.6586,.7638], which excludes
0.80 without any size correction at all.

================================================================================
## PROBED HARD, CAME BACK CLEAN
================================================================================
- the GLS slope, its se, W, p, mdt_wald and bias_wald are exactly the algebra
  described (1-2 ulp) -- p02_algebra.log
- the m vector: default reference, esrelativeto(-2) straddling the reference,
  and a biennial (gapped) relative-period axis all give b_pre = delta*m exactly
- hbar under anticipation() is on the SHIFTED onset, matching the event-study
  labels; bias = mdt*(hbar - ref) verified by hand
- hbar is the estimand's own weighting: simple mean of post relyr over treated
  post observations for att, cell-equal-weighted for att_it
- the Roth joint pre-trend p equals the jointtests pre-trend p on the same run
- mdt does not depend on the estimated coefficients (agrees to 1e-6 between a
  planted trend of 0 and of 3.0)
- Mata's column-major vech() IS the ordering mvnormalcvqp() expects (checked
  three ways: scipy oracle on the posted S, a 400,000-draw in-Stata Monte Carlo,
  and Roth's own pretrends package)
- mdt_nis agrees with Roth's `pretrends power 0.80` to 0.4% at df_r = 279, and
  the residual is exactly the t-vs-normal critical value ratio
- k > 4 correctly omits mdt_nis with the documented note
- bootstrap(standard) + cluster() keeps e(df_r) = G-1 and uses F(3,G-1)
- customvce(cluster) reaches the Roth block with the right df
- subgroups(): only the pooled families are diagnosed -- and sthlp:795 says so
  explicitly.  Not a finding.
- DESIGN A (iid, no cluster, i = unit, df_r = 279): size 0.0367, power 0.8000 at
  mdt_wald(80) and 0.5000 at mdt_wald(50), slope CI coverage 0.9600.  EXACT.
- DESIGN C (i = cohort, fe(unit), cluster(unit), G = 80): size 0.0455, power
  0.7864 / 0.8091, slope CI coverage 0.9500.  EXACT.

================================================================================
## UNTESTED
================================================================================
- whether the at-i caution's own remedy claim ("bootstrap-based clustered SEs
  give a joint test that is valid under this design") holds: each Monte-Carlo
  replication would need a full bootstrap, ~19 hours at the rep counts here.
- whether pretrendpower() under bootstrap() delivers its stated power (same cost).
- the G-dependence between 16 and 80 clusters: only the two endpoints measured.
- contreat()/ame and ame_it families: the hbar for `ame' reuses the att-side
  observation weighting and I did not check that against the ame estimand.
- rank-deficient pre-blocks (keff < k) and their generalized-inverse MDT.
- poisson + cluster + a Monte Carlo of the chi-squared branch's realised power
  (the arithmetic is measured; the realised power is not).

### PTP-5, the self-contained repro was executed (repro_designB.do / .log)
  STEP1 150 null replications
        size of the joint pre-trend test at 5% = 0.2800   (nominal .05)
        nominal F(3,79) 5% critical value = 2.72026
        empirical 95th percentile of the null statistic = 4.90010
        median reported mdt_wald(80) = 0.16926
  STEP2 150 replications at a trend of exactly that value
        power, NOMINAL cutoff        = 0.8333   (CLAIMED 0.80)
        power, SIZE-CORRECTED cutoff = 0.5133
  i.e. the file reproduces, on its own, the whole point: the column looks
  calibrated and is not.

### PTP-1, the self-contained repro was executed (repro_designE.do / .log)
(NOTE: the file's final sentinel still prints DONE_REPRO_DESIGNB in the shipped
log because the string was fixed after that run; the numbers are from
repro_designE.log.)
  STEP1 200 null replications
        size of the joint pre-trend test at 5% = 0.1100   (nominal .05)
        nominal F(3,15) 5% critical value      = 3.28738
        empirical 95th percentile of the null  = 3.89831
        median reported mdt_wald(80)           = 0.54662
  STEP2 200 replications at exactly that trend
        power, NOMINAL cutoff        = 0.8050
        power, SIZE-CORRECTED cutoff = 0.7200 (with this arm's own 200-draw
                                      cutoff; with the pooled 600-draw cutoff
                                      the figure is 0.60 -- see the caveat)
This independent seed stream reproduces the size (0.11 vs 0.1025) and the MDT
(0.5466 vs 0.5463) and is the reason the size-corrected figures in this report
are quoted from the POOLED null.

================================================================================
## SELF-CONTAINED REPRO FILES (all EXECUTED; logs alongside)
================================================================================
PTP-1  repro_designE.do   -> repro_designE.log
PTP-2  repro_nissize.do   -> repro_nissize.log
PTP-3  p12_poiscluster.do -> p12_poiscluster.log
PTP-4  repro_helpex7.do   -> repro_helpex7.log
PTP-5  repro_designB.do   -> repro_designB.log
also   repro_power.do     -> repro_power.log (+ repro_power_singledraw.log)
Each begins with `adopath ++ ".../scratchpad/bluesky/build"`, generates its own
data, and needs no sibling file.  The full Monte Carlos behind the pooled
figures are p04_simA.do (A), p07_simC.do (C), p08_simE.do + p13_simE2.do (E) and
p09_simB.do (B), with per-replication CSVs and analyze.py.
