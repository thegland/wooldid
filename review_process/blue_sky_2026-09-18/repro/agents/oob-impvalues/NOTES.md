# oob-impvalues — independent verification of the impvalues / Ibragimov-Muller claim

## Status
- M1_mech.do  DONE  (mechanism + noise sweep + verbatim output)
- M2_mc.do    RUNNING (R=200 Monte Carlo, 4 geometries x {impvalues, analytic})

## My own DGP (deliberately different from the prober's)
14 i-levels (cg): 12 treated over 4 onsets (2004/2005/2006/2007), 2 never-treated.
5 units per i-level = 70 units, 10 periods (2001-2010), 700 obs.
reg = 14 regions, units NEST in regions, regions CROSS cohorts, composition varied.
super = ceil(cg/2) = 7 groups, COARSER than i.
Effects spread widely: tau_cg = 0.30 + 0.12*(cg-1) (0.30..1.62) + mild dynamics.
Errors: unit RANDOM TREND us_u*(t-2005.5), us ~ N(0,0.15)  [NOT absorbed by cohort FE]
        + iid N(0,0.60).
File: mkdgp.doi

## M1 RESULTS (decisive)

### Harness sanity
noiseless analytic att b = 0.87818181818182 == hand-computed cell-weighted mean
of truetau over treated post cells (0.87818181818182). Harness correct.

### Mechanism: the G per-cluster "estimates" are copies, proven algebraically
Predicted, IF the 60 per-unit reads are 5 identical copies each of the 12
cohort-level reads:  b(cluster(unit)) == b(cluster(cg)) exactly, and
se(unit)/se(cg) = sqrt((K-1)/(nK-1)) = sqrt(11/59) = 0.431787770.

seed | b(unit)            | b(cg)              | se(unit)   | se(cg)     | ratio
4101 | 0.98191706734293   | 0.98191706734293   | 0.06401313 | 0.14825138 | 0.4317878
4102 | 0.93494865150714   | 0.93494865150714   | 0.05770185 | 0.13363475 | 0.4317878
4103 | 1.25573917858688   | 1.25573917858688   | 0.08541508 | 0.19781728 | 0.4317878

b identical to 14 digits; ratio equals the prediction. The 60 "per-cluster
estimates" carry exactly 11 degrees of cross-cluster information, but
e(df_r)=59 is posted and the t is taken on 59 df.

### Same-data se comparison (seed 4101)
cluster   IM se      analytic se
unit      0.064013   0.251175
cg        0.148251   0.221410
reg       0.038131   0.238648
super     0.193763   0.038210   <- analytic COLLAPSES when cluster is coarser than i

### Noise sweep (seed 5200, iid sd 0.2/0.6/2.0)
IM se cluster(unit): 0.0677 / 0.0794 / 0.1408
AN se cluster(unit): 0.2178 / 0.2646 / 0.4736
IM se cluster(cg)  : 0.1569 / 0.1840 / 0.3261
IM does move with noise (the pooled cohort taus carry noise), but under
cluster(unit) it is a fixed 0.4318 multiple of the cluster(cg) value —
refining the cluster label below i divides the SE by a pure label-count factor.

### What the user is told under geometry (a) — VERBATIM, no caution
Wooldid Estimation for Outcome: y; Standard Errors: Ibragimov-Muller p-values over cluster(unit) (t on G-1 df)
N clusters = 70; N = 700 (0 obs dropped from initial sample)
             |   estimate          se         t         p     lb_95ci
         att |  1.1132916   .06890962     16.16     0.000   .97540379
NOTHING else prints. No caution. (N clusters = 70 is also the raw count; only
60 treated units contributed a read, df_r = 59.)

## Documentation / prior-art context (established while the MC ran)
- sthlp:275 (the impvalues entry, rewritten by ruling T4 2026-09-15) frames the
  cluster-geometry problem purely as an ESTIMAND question: "IM-style inference
  changes wooldid's target estimand ... cluster() at the level of i with att_i is
  the pairing that is coherent by construction ... IM-style inference will yield
  incorrect results for NON-ATT estimands if clusters cross-cut cohorts".
  It says nothing about whether the SE is a standard error at all, and it
  explicitly scopes its warning to non-att estimands.
- sthlp:67 and sthlp:255 give the package's headline guidance: i = treatment
  onset cohorts, and "cluster BELOW the level of i". sthlp:1097 recommends
  impvalues as a valid fallback. Nothing reconciles the two.
- review_process/feature_tickets_2026-09-15.md T4 was RULED AND CLOSED as
  documentation-only, scoped to non-att estimands. It quotes IM's own
  requirement ("independent, approximately normal, COMMON MEAN") but does not
  ask whether the per-cluster reads are independent estimates. They are not.
- ado:12856 and ado:1528 both assert "IM dispersion inference is valid under
  containment" / "the estimator class that remains valid under cell saturation".
  Containment/saturation = cluster COARSER-OR-EQUAL to the cell = exactly the
  regime where the construction IS defensible. The code's own remedy caution
  (wooldid_amescreen) recommends IM only in that safe regime.

## Code sites confirmed (line numbers as the prober cited them, pinned build)
ado:2042  qui margins [w=w_att], dydx(targetflag) subpop(...) nose over(`cluster')
ado:2047-2048  grabs the second half of r(b) as the G per-cluster reads
ado:2052-2053  b = simple mean of the reads
ado:2054-2059  ttest reads == 0; se = r(se); ci = b +/- invt(r(df_t)) * se
ado:2066-2077  degeneracy scrub: fires only when se < 1e-10 * abs(b)  (relative only)
ado:4101 + wooldid_imdf: e(df_r) <- min IM df = q-1 (q = contributing clusters)

There is ONE margins call and ONE fitted model. The G numbers are G different
weighted averages of the SAME coefficient vector. Nothing is re-estimated.

## M2 MONTE CARLO RESULTS  (R = 200, same draws for every arm, my own DGP/seeds)
target = the estimand each arm actually reports, taken from a NOISELESS run of
that same arm (IM reports the cluster-equally-weighted mean, analytic reports
the cell-weighted att; each judged against its own target).

geometry                 arm  target    mean b   sd(b)    mean se  se/sd   cov    Wilson 95%
(a) cluster(unit) < i    IM   0.960000  0.923475 0.254055 0.069425 0.273  0.380  [0.315,0.450]
(a) cluster(unit) < i    AN   0.878182  0.842218 0.247357 0.237784 0.961  0.925  [0.879,0.954]
(b) cluster(reg) x i     IM   0.882224  0.845997 0.247890 0.039775 0.161  0.250  [0.195,0.315]
(b) cluster(reg) x i     AN   0.878182  0.842218 0.247357 0.243296 0.984  0.960  [0.923,0.979]
(c) cluster(cg) = i      IM   0.960000  0.923475 0.254055 0.160785 0.633  0.810  [0.750,0.858]
(c) cluster(cg) = i      AN   0.878182  0.842218 0.247357 0.153031 0.619  0.720  [0.654,0.777]
(d) cluster(super) > i   IM   0.956610  0.919752 0.253332 0.203793 0.804  0.935  [0.892,0.961]
(d) cluster(super) > i   AN   0.878182  0.842218 0.247357 0.044327 0.179  0.335  [0.274,0.403]

0 failures out of 1600 wooldid calls.
POSITIVE CONTROLS PASS: analytic at (a) 0.925 and at (b) 0.960 -- the two
geometries where the analytic VCE is the right object. So the harness is sound
and the (a)/(b) IM numbers are the option's, not mine.
(c)/(d) analytic under-cover as the help itself documents for cluster at/above i
(0.72 and 0.335) - not a harness fault.

VERDICT SHAPE: exactly the predicted pattern. impvalues is FINE where the code
comments say it is valid (cluster containing the cell: (d) 0.935, and (c) 0.810
which beats the analytic 0.720 on the same draws) and CATASTROPHIC where the
help's headline guidance sends users: (a) 0.380 and (b) 0.250 at a nominal 0.95.

## The algebraic identity, verified on ALL 200 independent draws (python on mcres.dta)
max |b_IM(cluster=unit) - b_IM(cluster=cg)|        = 1.11e-15
max |se_IM(unit)/se_IM(cg) - sqrt(11/59)|          = 5.55e-16
se ratio range over 200 draws: [0.431787769588372, 0.431787769588373]
i.e. the ratio is the constant sqrt((K-1)/(nK-1)) to the last bit, every draw.
So refining cluster() from the 12 cohorts to the 60 units -- which changes NO
data, only labels -- multiplies the reported standard error by 0.4318 and the
reported df from 11 to 59.

## Mean |t| against the TRUTH (should be ~0.8 for a well-calibrated t)
IM cluster(unit) : 3.141   P(|t|>3) = 0.415
IM cluster(reg)  : 5.700   P(|t|>3) = 0.625
IM cluster(cg)   : 1.356   P(|t|>3) = 0.080
IM cluster(super): 1.096   P(|t|>3) = 0.030

## A check I ran that did NOT support the claim (recorded so nobody re-runs it)
corr(IM se, honest analytic se) on the same draw is ~0 in EVERY geometry:
unit -0.020, cg -0.020, reg -0.006, super -0.018. It is ~0 even at (d), where
coverage is 0.935. At q = 6..12 the IM SE is a very noisy estimator, so a low
correlation is normal and is NOT evidence of the defect. Do not use it.

## M3 RESULTS

### A. NEAR-HOMOGENEOUS effects (effect spread scaled by 0.15), R=120
geometry              arm  target    mean b   sd(b)    mean se  se/sd   cov    covSE
(a) cluster(unit)<i   IM   0.144000  0.115322 0.249652 0.046004 0.184  0.342  0.043
(a) cluster(unit)<i   AN   0.131727  0.108318 0.244530 0.233253 0.954  0.933  0.023  <- control OK
(c) cluster(cg)=i     IM   0.144000  0.115322 0.249652 0.106543 0.427  0.708  0.042
(c) cluster(cg)=i     AN   0.131727  0.108318 0.244530 0.139650 0.571  0.675  0.043
=> LESS effect heterogeneity makes impvalues WORSE (0.380 -> 0.342 at (a);
   0.810 -> 0.708 at (c)), because the IM spread is driven by heterogeneity.
   My main MC is therefore the FAVOURABLE case for impvalues.

### B. wooldid refuses to estimate what impvalues claims to have estimated
wooldid y cg t ttre, att cluster(unit)   on ONE unit's rows  ->  r(459)
wooldid on that unit + all never-treated ->  b = 1.282, se = 0.2497
The honest one-unit-vs-controls SE is 0.25. impvalues' whole-sample IM se at
cluster(unit) is 0.069 - four times SMALLER than the SE of a single one of the
60 "per-cluster estimates" it claims to be averaging.

### C. honest analogue at cluster(cg): separate per-cohort fits vs wooldid's reads
(each separate fit = that cohort + the whole never-treated pool)
seed  wooldid read sd  separate-fit sd   IM se (wooldid)  IM se if separate
4101     0.5136            0.4920            0.14825          0.14202
4102     0.4629            0.4506            0.13363          0.13007
4103     0.6853            0.7089            0.19782          0.20464
=> At cluster = i the pooled reads reproduce the dispersion of genuine separate
   fits to within a few percent. The construction is a GOOD approximation of IM
   in that regime. This is the boundary, and it supports the innocent reading
   for (c)/(d) while leaving (a)/(b) indefensible.

### D. degeneracy caution on NOISY data - VERBATIM (sd(y) = 1.86)
Design: 13 i-levels x 6 units, 6 regions each holding exactly one unit per
i-level, so every region's read is identical.
  Caution: for one or more Ibragimov-Muller estimates the per-cluster estimates were numerically identical, so the cross-cluster (IM) standard error is not identified; wooldid has set the standard error, test statistic, p-value, and 95% confidence interval to missing for the affected estimate(s) while retaining the point estimate. A collapsed IM standard error usually reflects deterministic or noiseless data - please verify your input if this was not expected.
  Caution: standard errors were not identified for 1 of 1 reported estimates, so no variance-covariance matrix e(V) was posted. ...
  att |  1.0649753           .         .         .           .        .
The SAME data through the analytic path:
  att |  1.0649753   .54285793      1.96     0.107  -.33048539   2.4604361
The data is not deterministic; the message misattributes a structural property
of the construction (identical composition -> identical re-reads) to the input.

## REPRO
Self-contained: repro_OOB_impvalues.do in this directory (blocks 1-3, ~3 min).

## REPRO VERIFIED (repro_OOB_impvalues.log)
PREDICTED se ratio = 0.431787769588
seed 4101: |b(unit)-b(cg)| = 4.4e-16, ratio = 0.431787769588, df posted 59 vs 11
seed 4102: |b(unit)-b(cg)| = 3.3e-16, ratio = 0.431787769588, df posted 59 vs 11
seed 4103: |b(unit)-b(cg)| = 8.9e-16, ratio = 0.431787769588, df posted 59 vs 11
one unit's rows -> rc = 459
R=40 block: true sd 0.2938, IM se 0.0714, analytic se 0.2292,
            impvalues coverage 0.350, analytic coverage 0.850 (R=40 is noisy;
            the R=200 figures in M2 are 0.380 and 0.925).

## FILES
mkdgp.doi                 my DGP
M1_mech.do/.log           mechanism + noise sweep + verbatim geometry (a) output
M2_mc.do/.log, mcres.dta  R=200 MC, 4 geometries x 2 arms, same draws
M3_extra.do/.log,
  mcres_hom.dta           homogeneous-effect MC, degeneracy caution, separate fits
repro_OOB_impvalues.do/.log  self-contained 3-block repro (~3 min)
wilson.py                 binomial CI helper
