# BS-12: the adversarial pre-trend test's null

**Topic.** BS-12: the adversarial pre-trend test's null hypothesis — two-parameter free line vs. the one-parameter line through the reference period

## Finding

Both readings in the bug sweep are half right, and the sweep's framing needs three corrections. (1)
The shipped test is a coherent, exactly level-valid test — under known Gaussian V it has size
exactly alpha for the null it labels, and a pure reference-period level shift is a NULL POINT of
that null, not an alternative it fails to detect (measured size .0508 against mu = 2*1 with 40,000
draws). So this is not a size failure. (2) But the one-parameter form is not merely "better against
level shifts": because col(rel - esrel) is nested inside col([1, rel]), the noncentrality satisfies
lambda_1 >= lambda_2 for EVERY alternative without exception (0 violations in 20,000 random smooth
alternatives at each of three covariance shapes; it is a theorem of nested minimisation, not a
simulation result). The one-parameter form's only cost is one extra degree of freedom, bounded at
1.2-6.4 percentage points of power, and it is paid in full only in the regime where the reference
period's own level is nearly unidentified (s_r >> s_h). On the actual e(V_es_att) of a realistic
600-cluster design the implied s_r/s_h is 0.476 — squarely in the regime where the one-parameter
form wins — and it detects a reference-period anomaly of 2.18 horizon standard errors with 80% power
while the shipped test's power against any such anomaly, however large, is exactly .05. (3) The
blast radius is three consumers, not four: there is NO fixed-metric bootstrap adversarial test.
Without boottstat, wooldid_bootjnt posts e(adversarial_pval_*) = . and e(boot_adv_absent) = 1. The
three live consumers are the analytic path, boottstat, and permute(), and on the BS-12 design all
three return a large p (.9252 analytic, .8293 boottstat, .9016 permute) on a path that is uniformly
-5.0 with |t| of 40-68 and a joint pre-trend p of 0.0000.

## Mechanism

Under esfixedbaseperiod the estimator reports b_h = theta_h - theta_r for h != r, so the reference
period's own level is subtracted out of every reported coefficient and is not itself estimated. If
the underlying period effects lie on a line, theta_h = a + c*h, then b_h = c*(h - r): a ONE-
parameter family whose line is pinned to pass through (r, 0). The shipped test instead minimises
over the two-parameter family E[b] = g0 + g1*h. Because the constant vector 1 lies in col([1, rel]),
the whitened residual projection annihilates it, and the statistic is exactly invariant to b :+ c. A
uniform level shift of every reported coefficient is precisely what a bad reference period produces
— if theta_r alone is displaced by -c (an anomalous or anticipatory reference year), then b_h = +c
for every h — so that pathology is a null point of the shipped composite null, and the test returns
its nominal size against it however large c is. The extra restriction the one-parameter form
imposes, g0 + g1*esrel = 0, is exactly the statement that the reference period lies on the same line
as the rest of the path. The two statistics stand in an exact nesting relationship. Write z =
V^{-1/2} b in the whitened eigenspace, X2 = V^{-1/2}[1, rel], X1 = V^{-1/2}(rel - r), so col(X1) is
a codimension-1 subspace of col(X2). Then W_1 = W_2 + u^2 where u is the component of z along the
one direction in col(X2) orthogonal to col(X1), q_1 = q_2 + 1, and the noncentralities obey lambda_1
= lambda_2 + delta^2 with delta^2 >= 0 always. This is why the one-parameter form can never have
less noncentrality against any alternative: its only cost is one degree of freedom in the reference
distribution. Whether that cost is repaid depends on how much information the design carries about
the extra direction, and that quantity is delta^2/c^2 = 1'V^-1 1 - (m'V^-1 1)^2/(m'V^-1 m). Under
the canonical fixed-base covariance V = s_h^2 I + s_r^2 J (independent period effects plus the
common -theta_hat_r term), 1'V^-1 1 = (n/s_h^2)/(1 + n s_r^2/s_h^2) -> 1/s_r^2 as s_r grows. So
level-shift power is capped by the precision of the reference period itself: a thin reference period
is one you cannot check, and there the two tests coincide in content and the extra df is dead weight
(2-4 points of power). A well-estimated reference period — the ordinary balanced case, and the case
of the real e(V_es_att) measured here at s_r/s_h = 0.476 — puts the one-parameter form ahead on
essentially every alternative. The defect therefore is not a size failure and not confined to level
shifts: it is that the documented null answers a different question from the one a fixed-base event
study poses, and the direction it discards is the one that carries the reference period's own
credibility. The diagnostic consequence is the interpretation of the joint/adversarial pair. Under
the one-parameter form, 'joint pre-trend rejects but adversarial does not' means one thing: a clean
linear pre-trend through the reference period. Under the shipped form the same pattern is also
consistent with a reference-period anomaly of arbitrary size, which is a completely different
finding — and on the BS-12 design that is exactly what happens (joint p = 0.0000, adversarial p =
0.93). The change propagates to boottstat and permute() without disturbing their validity arguments,
because both are rank tests of whatever statistic wooldid_advW computes. Under boottstat the
recentred replication statistic W_r = advW(b_r - b_obs, V_r, rel) is still an exact draw from a
distribution satisfying the narrower null (E[b_r - b_obs] = 0 lies in span(m) as well as in col(X));
only the code comment's justification, which cites invariance to col(X), goes stale. Under permute()
the reference distribution is exchangeability under the sharp null of no effect at any horizon,
which holds for any statistic, so the test remains exact and only its power moves. The one path that
does NOT exist is a fixed-metric bootstrap adversarial test: without boottstat the p is not computed
at all.

## Evidence

- NULLS, STATED AS RESTRICTIONS. wooldid.ado:11757 sets X = J(rows(rel),1,1), rel and 11773 sets q =
  k1 - rankM. With k estimated horizons the shipped null is: there exist g0, g1 with E[b_h] = g0 +
  g1*h for every estimated h — equivalently every second difference of E[b] vanishes, k-2
  independent restrictions. Because b_h = theta_h - theta_r under esfixedbaseperiod, that is 'theta
  is linear over the ESTIMATED horizons, with theta_r free'. The one-parameter null adds g0 +
  g1*esrel = 0, i.e. E[b_h] = c*(h - r), k-1 restrictions: 'theta is linear over the estimated
  horizons AND the reference period'. The single extra restriction is exactly the reference period
  lying on the same line.
- THE REFERENCE PERIOD IS ABSENT FROM THE TESTED VECTOR, CONFIRMED. wooldid.ado:1510-1548 builds
  relyrcolumn only from relative periods that carry an estimable coefficient in e(b) (the regexm at
  1516 excludes the base level), so under esfixedbaseperiod b/V never contain the reference row.
  wooldid_advstat plants the reference rows at 10593-10594 into `advF' — the DISPLAY-line frame only
  — and passes the raw f2bmat/f2vmat/relyrcolumn to mata at 10628. The header comment at
  wooldid.ado:11710-11713 ('planted zero-variance reference rows receive zero weight under V^+ and
  so are mathematically inert here... omitting them from the computation is identical to including
  them') is therefore describing rows that are not there, and its equivalence claim is false in any
  case.
- BS-12'S MATA REPRODUCTION, EXACT. b=(1\1\1\1), rel=(-3\-2\0\1), V=.01*I(4), esrel=-1: shipped
  wooldid_advW(alt=0) returns W = 3.15544e-29, q = 2, chi2 p = 1.00000000. The one-parameter twin
  returns W = 400, q = 3, p = 0. The closed form b'V^-1 b - (m'V^-1 b)^2/(m'V^-1 m) with m = rel -
  esrel is exactly 400. Planting the reference row at variance eps and running the SHIPPED statistic
  gives W = 80 / 384.615 / 399.840 / 399.9984 / 399.99998 / 400.000 at eps = 1e-2 ... 1e-12 with q =
  3 throughout — so the shipped alt==0 statistic is NOT the eps->0 limit of the corresponding non-
  degenerate problem, and the one-parameter statistic IS. Verified across 200 random 6x6 full-rank V
  and b: max relative gap between the one-parameter statistic and the closed form 4.2e-14.
- INVARIANCE. Shipped W is invariant to b :+ c (3.16e-29 -> 2.02e-27 at c=7, i.e. numerical zero
  both times) and to b :+ 3*m. The one-parameter W is invariant to b :+ 3*m (400 -> 400) but NOT to
  b :+ c (400 -> 25600). That is the whole difference.
- LEVEL VALIDITY, MEASURED. Known-V Gaussian MC, 40,000 draws, rel = -4,-3,-2,0,1,2, esrel = -1, V =
  I + J (s_h = s_r = 1), alpha = .05 chi2 reference. Rejection rates: mu = 0 -> shipped .0500, one-
  param .0493; mu = .5*(rel-r) -> shipped .0498, one-param .0497; mu = 2*ones (pure reference-period
  anomaly) -> shipped .0508, one-param .2549. Analytic noncentral-chi2 predictions .0500/.0500,
  .0500/.0500, .0500/.2517 match to Monte-Carlo error. So the shipped test is exact-size for its own
  null under known V, and the level shift is a null point of it — the .05 is size, not lost power in
  the usual sense.
- POWER, THE THEOREM AND THE PRICE. lambda_1 >= lambda_2 always (nested minimisation over col(X1)
  subset col(X2)); empirically 0 of 20,000 random smooth alternatives (theta_h = c1 h + c2 h^2 + c3
  h^3, with and without an independent shock to theta_r) violated it, at s_r = 0.3, 1 and 3. Mean
  lambda_1/lambda_2 = 1.908 / 1.199 / 1.025 at those three s_r. The pure df cost at IDENTICAL
  lambda: power(q) - power(q+1) = .057/.064/.043/.010 at q=2 for lambda 4/8/12/20, falling to
  .012/.022/.023/.013 at q=8 — bounded above by about 6.4 points and typically 2-4.
- POWER, NET, AT FIXED SHIPPED POWER. Random smooth alternatives rescaled so the shipped test has
  power exactly .50 (lambda_2 = 6.4195, q2 = 4, q1 = 5): mean one-param power .7729 / .5515 / .4737
  / .4628 at s_r = 0.3 / 1 / 3 / 10; mean difference +.273 / +.052 / -.026 / -.037; share of
  alternatives where the one-param form loses .030 / .107 / .991 / .999; worst single loss -.038 in
  every column. Same picture rescaled to shipped power .80. So the crossover is s_r ~ s_h: below it
  the one-parameter form wins on nearly everything, above it the two tests are near-identical in
  content and the extra df costs 2-4 points.
- WHERE THE REAL DESIGN SITS. On the actual e(V_es_att) of the BS-12 DGP (7 horizons, cluster(id),
  600 clusters, 20 clusters spanning every treated cohort-period cell): mean horizon se = 0.0981,
  mean off-diagonal correlation among the horizon coefficients = 0.434 (the common -theta_hat_r
  component), level-shift noncentrality per unit shift delta^2/c^2 = 281.27, implied s_r/s_h =
  0.476. The one-parameter test reaches power .80 against a reference-period anomaly of c = 0.2136
  outcome units = 2.18 horizon standard errors; the shipped test's power against that anomaly, and
  against one of any size, is exactly .05.
- END-TO-END, ALL THREE LIVE SURFACES. DGP: 30 cohorts x 20 individuals x 12 periods; onsets 5/7/9
  for cohorts 1-8/9-16/17-24, 25-30 never treated; no treatment effect at any horizon; +5 added only
  at each treated cohort's own reference period (t = ttre - 1). Every event-study coefficient comes
  back at -4.99 to -5.10 with |t| between 40.8 and 67.9; joint pre-trend p = 0.0000. Adversarial p:
  0.7477 analytic with cluster(id) (W = 2.6887, q = 5, F(5,599) = 0.5377) and 0.9252 analytic with
  cluster(i) — against 0.9432 in the bug sweep on its own draw. The one-parameter statistic on the
  SAME e(b_es_att)/e(V_es_att) gives W = 7122.59, q = 6, F = 1187.10, p = 0. W1 - W2 = 7119.90, i.e.
  the entire signal lives in the one direction the shipped design absorbs.
- THE OTHER TWO SURFACES MOVE THE SAME WAY. boottstat (bootstrap(standard) reps(60), cluster(i)):
  e(adversarial_pval_att) = 0.8292682927, e(boot_advW) present with 60 usable replications.
  permute(shuffle) permreps(60): e(adversarial_pval_att) = 0.9016393443 with e(perm_advWobs) =
  1.701855. Both are exact/valid for their own reference distributions (recentring at b_obs under
  boottstat, exchangeability under the sharp null for permute) — what moves is power, not validity.
- THERE IS NO FIXED-METRIC BOOTSTRAP ADVERSARIAL TEST — BUG SWEEP CORRECTION. bootstrap(standard)
  reps(40) WITHOUT boottstat on the same data returns e(adversarial_pval_att) = . ,
  e(boot_adv_absent) = 1, e(boot_jt_family) = fixedmetric (that macro covers the joint and equality
  tests only), and no e(boot_advW). The code is wooldid.ado:6153-6208: `local wdb_advok = ("`bts'"
  != "" & ...)` gates the whole block on boottstat. The fixed-metric bootstrap twin wooldid_bootjqf
  at 13757 exists for joint/equality tests only. BS-12's refs '13437, 13442' land inside
  wooldid_permq/permtail in the current file; the fixed-metric-looking code it describes is
  wooldid_permadv at 13813-13860, which is the PERMUTATION path. Help sthlp:485 and sthlp:1061 are
  correct; sthlp:499's clause 'Without boottstat those tests are the fixed-metric rank tests
  described under bootstrap()' over-includes the adversarial tests and contradicts both the code and
  its own run-time caution at ado:6197.
- WOOLDID_PTP — WHAT IT ACTUALLY DOES, BUG SWEEP NUANCE. wooldid.ado:8940-8942 builds m = relyr -
  ptpref and 8997-8998 fit slope = (m'S^-1 b)/(m'S^-1 m), se = 1/sqrt(m'S^-1 m): one parameter, no
  intercept, through the reference. But ado:9000-9001's roth_W = b'S^-1 b is the b = 0 test, not an
  on-a-line test. The exact bridge, verified numerically: the one-parameter adversarial statistic
  evaluated on the PRE block equals roth_W - slope^2*(m'S^-1 m) with df n_pre - 1 (measured: roth_W
  1.97848539, slope^2*c2 0.55998443, difference 1.41850096, one-param advW 1.41850096, q = 2). So
  ptp's LINE is the one-parameter form; the package does carry two linear-pre-trend objects under
  two normalisations, but the object that differs from the adversarial test is ptp's slope, not
  ptp's printed joint test.
- BLAST RADIUS — REPORTED QUANTITIES. e(adversarial_pval_att / _att_it / _ame / _ame_it) and every
  _s<SG> variant, posted at ado:3923 and 3947 (analytic), 6186-6188 (boottstat), 7417-7419
  (permute); the printed lines at ado:9378 and 9520 (both say 'all (non-reference) ... on a single
  line'); e(boot_advW) and the bootsave advW column (ado:12481); e(perm_advW)/e(perm_advWobs). NOT
  affected: e(adversarial_line_*) (the unweighted display OLS at ado:10585-10608, which already
  carries the planted reference rows), the altmethod alt==1 branch, and every non-adversarial
  statistic.
- BLAST RADIUS — TESTS AND BASELINES THAT FLIP. test/test_advguard.do T2b (lines 103-131) uses a gap
  panel leaving only two post-treatment coefficients and ASSERTS that no p is posted and that the
  'composite-null covariance is degenerate' caution prints; under the one-parameter form q goes from
  0 to 1 and a p appears (verified in Mata: rel = (-2\0), esrel = -1, V = I, b = (1\1) -> shipped q
  = 0 / missing, one-param W = 2, q = 1, p = 0.1573). test/test_advstat.do T1 (141-165) reconstructs
  the p by hand and T5 (264-290) asserts the timetrends refusal — both touch the rank arithmetic
  that changes from rank(V) - 2 to rank(V) - 1. Five baseline CSVs carry adversarial values:
  test/baseline/adv_gapa.csv, adv_gapb.csv, adv_altb.csv, adv_r2001.csv, esf_adv.csv. Six .do files
  under test/ read e(adversarial_pval*).
- TEXT THAT BECOMES WRONG UNDER READING (a). sthlp:591-595 ('some line, with slope and intercept
  unknown'; 'the rank of the event-study covariance minus the number of line parameters
  identified'), sthlp:1118 ('all (non-reference) relative-time period treatment effect estimates'),
  the printed labels at ado:9378/9520, the false inertness comment at ado:11710-11713, the 'no
  testable residual dimensions (e.g. only two real horizons under the primary composite null)'
  comment at ado:10633-10636, and the boottstat recentring comment at ado:12645-12650 (invariance
  shrinks from col(X) to span(m); recentring at b_obs stays exact because E[b_r - b_obs] = 0
  satisfies the narrower null too, but the stated justification must be rewritten). Separately, the
  measured size claims at sthlp:604-614 were calibrated on the two-parameter statistic and would
  need re-measurement.
- THE EXISTING NULL DGPs REMAIN VALID NULLS. test/suite/mc/mc_adv.do generates y = a + .05*t +
  SLOPE*(t - 6 + 1)*(ttre < .) + e with ttre = 6, so the planted line is SLOPE*(h - (-1)) = SLOPE*(h
  - esrel): a line THROUGH the reference period. Config A1 (SLOPE = 0) and A2 (SLOPE = .15) are
  therefore one-parameter nulls as well as two-parameter nulls, so the existing size-calibration
  harness can be reused for both statistics without redesigning the DGPs.

## Recommendation

The ruling is the author's, and the discriminating fact is a question about what the reference
period is for, not about validity. Both statistics are exact-size tests of their own nulls under
known V; neither is broken. The choice is whether the adversarial diagnostic should privilege the
declared reference period. The case for (a), the one-parameter form: it is what the estimator's own
normalisation implies; lambda_1 >= lambda_2 against every alternative without exception, so there is
no alternative class against which it is less informative in noncentrality; the design where the
author's own realistic DGP sits (s_r/s_h = 0.476) is well inside the region where it wins; and it
makes the joint/adversarial pair diagnostic rather than ambiguous, since under the shipped form
'joint rejects, adversarial accepts' cannot distinguish a clean linear pre-trend from a reference-
period anomaly of any magnitude. wooldid_ptp's slope already uses this normalisation, and honestdid-
style sensitivity work assumes it. The case for (b), keep the two-parameter form: it is a robustness
property, not an accident. Its statistic is nearly insensitive to which period was chosen as the
base, so it will not fire merely because the user picked an unlucky anchor; the one-parameter form
deliberately makes the answer depend on that choice. And in designs where the reference period is
thin (s_r >> s_h), the extra parameter buys almost no information while the extra degree of freedom
costs 2-4 points of power on ~99% of smooth alternatives. The counter is that esrelativeto() and
anticipation() exist precisely so the user can move the anchor deliberately, and that a reference
period too thin to check is a fact the user should learn rather than have papered over. One
prerequisite before either ruling is safe: the size study in proposed_simulation. The direction the
one-parameter form adds is the constant direction, which carries the common -theta_hat_r component
and is plausibly the worst-estimated direction of a cluster-robust V with few clusters. The help
already documents that these line tests are severely oversized when cells are spanned by fewer than
about six clusters (sthlp:604-614), and those numbers were measured on the two-parameter statistic.
If the one-parameter statistic is materially more oversized at small G, that is the single fact that
could make (b) preferable on its merits, and it is currently unknown. Independently of the ruling,
three corrections should be recorded now because they are factual, not judgemental: (i) there is no
fixed-metric bootstrap adversarial test — sthlp:499's 'those tests are the fixed-metric rank tests'
over-includes the adversarial tests, contradicting the code at ado:6153-6208 and its own caution at
ado:6197 and sthlp:1061, and the bug sweep's downstream refs 13437/13442 point at the permutation
path; (ii) the inertness comment at ado:11710-11713 is false and describes rows that never reach the
function, whatever the ruling; (iii) the sthlp:604-614 size numbers need re-measurement in any case,
because they were calibrated before the current statistic's guards. If (a) is ruled, the
implementation carries four consequences the sweep only partly names: the timetrends sub-case (two
normalised periods force a linear path to theta ≡ 0, so the null collapses to b = 0 — set W = z'z, q
= k1, and skip the svsolve rather than handing it an empty design); test/test_advguard.do T2b flips
from a refusal to a posted p (q 0 -> 1) and must be respecified rather than re-blessed;
test/test_advstat.do T1 and T5 need their reconstructions and rank arithmetic updated; and five
baseline CSVs plus every suite config carrying an adversarial p need a re-bless under the pristine-
twice doctrine. Text: sthlp:591-595, sthlp:1118, the labels at ado:9378/9520 (drop '(non-
reference)'), the comments at ado:11710-11713, 10633-10636 and 12645-12650.

## Open questions

- Does the extra degree of freedom interact with the documented small-cluster oversizing? The direction the one-parameter form adds is the constant direction, which carries the common -theta_hat_r component and is a plausible candidate for the worst-estimated direction of a cluster-robust event-study covariance with few clusters. The help's own measured oversizing (threefold at two clusters per cell, an order of magnitude at one) was calibrated on the two-parameter statistic. This is the one unresolved fact that could reverse the ruling, and the known-V Monte Carlo run here cannot see it — it is the first item of the proposed simulation.
- Under reading (a) the primary and altmethod tests would sit under different normalisations. wooldid_advW's alt==1 branch (ado:11783-11816) appends the planted reference row to bP/relP as one ordinary unweighted observation in the OLS forecast operator, with a zero row and column in VPP — i.e. it treats the reference as an observation, not as a hard constraint. That is a defensible forecast-test design, but after (a) the two adversarial options would answer questions under two different normalisations, mirroring exactly the ptp/adversarial split the sweep complains about. Out of BS-12's scope; flagged, not investigated.
- The display line e(adversarial_line_*) is an unweighted OLS fit over the estimated horizons PLUS the planted reference rows with a free intercept, so it is already neither of the two tested nulls. It would remain a two-parameter fit under (a) while the test became one-parameter. That may be acceptable (it is documented as display-oriented and the sweep says leave it) but it is a third normalisation in the same output block.
- Both level-validity statements here are for known Gaussian V. With V estimated, both statistics are F approximations and the relative finite-sample size distortion is unmeasured. The shipped statistic's exactness under mu = c*1 is a property of the projection, not of the F reference, so it survives; the question is whether the two tests' F approximations degrade at different rates.
- The one-parameter test's power against a reference-period anomaly is bounded above by the precision of the reference period itself (delta^2/c^2 <= 1'V^-1 1 -> 1/Var(theta_hat_r)). No simulation was run on a design where the reference period is genuinely thin — e.g. esrelativeto placing r on a cohort's panel edge, or a gap panel where only a subset of cohorts realise r. Whether the loss regime (s_r >> s_h) is reachable in practice or only in construction is unsettled.
- The boottstat arm produced one replication with e(boot_advW) = 3.58e25 (a degenerate replication covariance) out of 60. That is orthogonal to BS-12 and did not affect the reported p (a rank test), but it means the boottstat adversarial reference distribution can carry pathological draws, and the guard at ado:11762-11766 that returns missing on degeneracy evidently did not catch this one.

## Proposed simulation

GOAL. Settle the one question the known-V analysis cannot answer: whether the one-parameter
statistic is more oversized than the shipped one under cluster-robust covariances with few clusters
and thinly spanned cells. Power is already settled analytically (lambda_1 >= lambda_2 always; the
exact df cost is tabulated), so this is a SIZE study first and a power confirmation second. HARNESS.
Extend test/suite/mc/mc_adv.do rather than writing a new one. Its null DGPs already serve both
statistics: y = a + .05*t + SLOPE*(t - 6 + 1)*(ttre < .) + e with ttre = 6 plants a line through the
reference period (h - esrel with esrel = -1), so config A1 (SLOPE = 0) and A2 (SLOPE = .15) are
valid nulls for BOTH the two-parameter and the one-parameter null. Add one post-estimation Mata hook
per replication that reads e(b_es_att)/e(V_es_att)/e(df_r) once and writes BOTH statistics and BOTH
p-values to the CSV, so every cell is a paired comparison on identical data with zero extra
estimation cost. Columns per row: config, shard, rep, seed, rc, k (horizons), rankV, W2, q2, p2, W1,
q1, p1, pjoint, and the level-shift information delta2 = 1'V^-1 1 - (m'V^-1 1)^2/(m'V^-1 m) plus
mean(sqrt(diag(V))) so the s_r/s_h location of each cell is recorded rather than assumed. DESIGN
GRID (size arm, 2000 reps per cell, all cells run at both statistics from the same fit). - Cluster
count G in {10, 24, 40, 200}. G = 10 and 24 are the regimes the help already flags; 200 is the
asymptotic anchor. - Clusters spanning each treated cohort-by-period cell in {1, 2, 6, 20}. This is
the axis the help says dominates size. Implement by fixing G and varying how many units per cell, as
mc_adv.do's A5 config already does at G = 10. - Horizon count k in {3, 5, 8} (set by
esprelength/espostlength), so q2 in {1, 3, 6} and q1 in {2, 4, 7}. k = 3 is where the df cost is
largest and where q2 = 1 makes the shipped test nearly a t-test. - esrelativeto in {-1, -3}. Under
-3 the pre relative years straddle the reference and m carries both signs, which exercises a code
path the -1 case does not. - Reference-period thinness in {equal, one third}: drop two thirds of the
observations at each treated cohort's own reference period. This is the axis that moves s_r/s_h and
is the only way to reach the regime where the one-parameter form is predicted to lose. Record the
realised s_r/s_h from delta2 rather than assuming it. - Null shape in {flat (A1), planted line
through the reference (A2)}. Both are nulls for both statistics; A2 checks that the one-parameter
form is not accidentally sensitive to a true through-reference trend. - Error structure in {iid,
AR(1) rho = .7 (A6), cohort-by-period common shock}. The last is the case the help says no analytic
variance estimator can handle; include it to confirm both statistics fail together rather than
differentially. STATISTICS TO REPORT PER CELL: empirical rejection at nominal .05 and .10 for p2 and
p1 separately, their difference with a paired Monte-Carlo standard error, the share of replications
where the two tests disagree at .05, and the share where each returned missing (degenerate). PASS
CRITERION for reading (a): size(p1) - size(p2) within Monte-Carlo error in every cell with 6 or more
clusters per cell, and no cell where size(p1) exceeds size(p2) by more than 2 percentage points at G
>= 24. DEGENERACY ARM (500 reps, cheap). Configs where the shipped test currently refuses: exactly
two estimated horizons (the test/test_advguard.do T2b gap panel), and timetrends +
esfixedbaseperiod. Record how often the one-parameter form posts a p where the shipped one refuses,
and — for the two-horizon case — the empirical size of that newly available q = 1 test. This is a
size question, not just a coverage-of-cases question, because the newly testable direction is the
one whose variance the reference period dominates. POWER ARM (500 reps per cell, confirmatory only).
Same grid restricted to G in {24, 200} and cells spanned by 20 clusters. Alternatives, each
parameterised in horizon standard errors so the cells are comparable across designs: (i) reference-
period anomaly, theta_r displaced by -c for c in {0.5, 1, 2, 4} horizon SEs — the shipped test must
sit at nominal size in every one of these, which is the empirical confirmation of the projection
result; (ii) realistic quadratic path theta_h = a h^2 at three a values chosen to put the shipped
test near power .3, .5 and .8; (iii) step at zero, theta_h = s*1(h >= 0), i.e. a genuine treatment
effect with no pre-trend, at two s values — this checks how much of the adversarial test's rejection
is being driven by the post-treatment effect itself under each form; (iv) a kink, theta linear
before the reference and linear with a different slope after, which is the alternative practitioners
most often have in mind. BOOTSTRAP AND PERMUTE ARM (100 outer Monte-Carlo draws x 200 inner
replications, one design only: G = 24, cells spanned by 6 clusters, k = 5, esrelativeto(-1), flat
null). One boottstat cell and one permute(shuffle) cell, each run at both statistics. The purpose is
not to re-derive validity — recentring and exchangeability are both statistic-agnostic — but to
confirm empirically that the add-one p under each is calibrated at .05 for both forms, and to count
how often a replication returns a missing or pathological W (the 3.58e25 draw observed here suggests
the degeneracy guard is not tight under resampling). BUDGET. The size arm dominates: 4 G x 4
spanning x 3 k x 2 esrel x 2 thinness x 2 shape x 3 error = 1152 cells at 2000 reps is too large.
Cut to a fractional design: run the full 4 x 4 (G x spanning) at k = 5, esrelativeto(-1), equal
thinness, flat null, iid errors — 16 cells, the primary question — then vary one factor at a time
off that centre point (k, esrel, thinness, shape, error structure), giving about 16 + 30 = 46 cells
at 2000 reps. At the observed cost of a wooldid fit on these panel sizes that is a few hours sharded
across four Stata instances, which matches the existing mc_adv.do sharding convention. WHAT WOULD
SETTLE IT. Reading (a) is safe if the size arm shows no differential oversizing beyond Monte-Carlo
error in the well-spanned cells and no more than a couple of points in the thin ones. Reading (b)
becomes the better choice if the one-parameter statistic is materially more oversized at G = 10 or
at one-to-two clusters per cell — because then the extra direction is not merely uninformative there
but actively misleading, and the shipped form's blindness to it would be a feature rather than a bug
in exactly the designs where the help already tells users not to trust the test.
