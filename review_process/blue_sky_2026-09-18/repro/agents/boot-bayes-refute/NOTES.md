# NOTES — skeptic pass on boot-bayes findings bb-01 / bb-02 / bb-03
workdir: this directory. Build loaded by `adopath ++ .../bluesky/build` (asserted by `which wooldid`).

## r01_saturated.do  (my own data, their structural shape: i=unit=cluster, singleton treated cells)
28 states x 14 years, onsets 2006/2009/2012 (7/7/6 states), 8 never-treated, cluster(st), i=st.
het arm adds a CLUSTER-CONSTANT shock to the post rows only. reps(200) rseed(77777).

- att and att_it ARE bit-identical in b AND analytic se.                     -> bb-01 part 1 REPRODUCES
- bootstrap(standard) reports se(att) == se(att_it) bitwise.                 -> bb-01 part 2 REPRODUCES
- bootstrap(bayes)+safety(off): se(att_it)/se(att_it|standard) = 0.5003 (het=0), 0.1671 (het=1.5)
                                                                             -> bb-01 part 3 REPRODUCES (6.0x, not 4.4x)
- bayes se(att_it) BITWISE IDENTICAL across het arms: +1.c432c1b72effdX-005 both; relmove 0.0000e+00
                                                                             -> bb-01 part 4 REPRODUCES exactly
- **CORRECTION THE PROBER MISSED**: the ANALYTIC clustered se is ALSO bitwise identical across the
  het arms (+1.2fdb96bdb85e9X-004 in both), and for BOTH att and att_it. The default inference in
  this geometry is exactly as blind to the cross-cluster heterogeneity as bayes is. The blindness is
  a property of the SATURATED SINGLETON-CELL geometry (treated rows are fitted exactly, residual 0,
  so a post-only shock changes no score), not of the Dirichlet engine.

## r02_scope.do  (the scope the prober declared untested)
B: i = onset cohort (5 groups of 8 states), cluster(state) FINER than i, 24 treated (i,t) cells of 8 rows
   ratio bayes/standard on att_it = 0.9654 (het=0), 0.8523 (het=1.5); bayes se(att_it) MOVES +53.7%
   with the heterogeneity. -> NO collapse. bb-01 does not generalise to a coarser i.
C: i = state, cluster(region) COARSER than i, singleton cells
   ratio bayes/standard on att_it = 0.6199 (het=0)  -> the collapse is NOT specific to cluster()==i.
   Analytic se differs only in the last hex digit across het arms (essentially invariant there too).

## r03_bootsave.do  (bb-02)  CONFIRMED
480 rows, 15% item nonresponse on y, e(N)=408, G=30, bootalpha(size).
_rep=0 = +1.11111111110c4X-009 = 1/480 (MATCH_NROWS=1); 1/e(N)=0.00245098 (MATCH_eN=0).
sum_c n_c*w_c on rep 1 = 0.9999999999999999 with the MAP-FRAME n_c, 0.8488 with the
estimated-row n_c -> the saved value is the right comparator for the saved columns;
only sthlp:557's phrase "estimation sample" is wrong.  bootalpha(equal) leg: 1/G, correct.

## r05_designD.do
D (i=st=cluster, THREE rows per (st,year) cell): ratio bayes/std att_it 0.6102 (het=0),
  0.2350 (het=1.5); bayes se(att_it) AND analytic se(att_it) both bitwise invariant across
  het arms.  -> singleton cells are NOT required; "cluster nests i" is the predicate.
E (same design, only the idiosyncratic sd changed 0.6 -> 1.5, no cluster heterogeneity):
  bayes se(att_it) 0.05519998 -> 0.13799996, relmove exactly +1.5000e+00 = the noise ratio.
  -> the bayes cell-weighted se is EXACTLY proportional to the within-cell noise and
  EXACTLY blind to the cluster-level component.  It is not a frozen constant.

## r07_centring.do  (own hunt - NEGATIVE result)
26-state unbalanced panel, unbalanced entry, aw=pop with pop p10 0.37 / p90 8.54,
att+att_it+att_i, reps(300).  (rep mean - b)/bootstrap sd:
  standard     att -0.179  att_it -0.013  att_i -0.073
  bayes_equal  att -0.054  att_it +0.027  att_i +0.007
The 2026-09-06 (equal-split) and 2026-09-13 (size divisor) centring fixes HOLD on weighted,
unbalanced data.  The ado comment says the old bug displaced the draws by 9.36 bootstrap sd;
nothing near that survives.

## r08_gteffects.do / r08b  (own hunt)
Event-study and pretreatment rows under bayes vs standard on the het panel, att only,
NO safety gate: ratios 0.907-1.018 (att 0.955).  No ungated collapse on the displayed rows.
(r08's automatic "se column" detection picked column 2, `onset', because strpos("onset","se")
hits - r08b re-runs with the correct column 6 for the gteffects CELLS.)

## r04_tailshape.do  (bb-03 b, mechanism)
Panel 1: standard R=999 repmean-b = -0.0023 sd, skew -0.042, pct width 0.786;
         bayes    R=999 repmean-b = +0.0029 sd, skew +0.015, pct width 0.755 (ratio 0.960).
         At R=79 the same panel gives bayes skew +0.486 and width 0.697 (ratio 0.871).
-> BOTH replication distributions are centred and near-symmetric at 999.  The prober's
   "the bayes interval is shifted, not narrower" does not hold up: on my panels the bayes
   percentile interval is NARROWER, and more so at their reps(79).

## r06_[abcd]  coverage MC (att, 30 clusters, true att = 1, reps(79), nominal 95%)
interim at 26 paired panels: standard pct 24/26 (.923), bayes 22/26 (.846), discordance 2-0,
mean pct width standard 0.642 vs bayes 0.568 (ratio 0.886) - NOT "within 2%".
mean se/true sd: standard .844, bayes .768, analytic .797.

## r09_analytic_blind.do
Four arms, one panel, identical but for the treatment-effect heterogeneity.
Reported ANALYTIC se(att): arm0 = arm1(unit shock sd1.5) = arm2(unit shock sd4.0) BITWISE
(+1.2e2d6dd8f2531X-003); arm3 (per-CELL shock sd1.5) differs by 3 ulps.
bootstrap(standard) se(att): 0.152361 (arm0) -> 0.935422 (arm2), factor 6.14.
Point estimate moves 0.905 -> 1.463.

## r10_sesd_het.do  (80 draws per arm)
unit-effect sd | TRUE sd of att | mean reported analytic se | se/sd || cell: TRUE sd, se, ratio
   0.0         |    0.1515      |        0.1415             | 0.933 || 0.8671  0.2166  0.250
   0.5         |    0.1921      |        0.1415             | 0.736 || 0.9318  0.2166  0.232
   1.5         |    0.3956      |        0.1415             | 0.358 || 1.5269  0.2166  0.142
The reported se is the SAME NUMBER in all three arms.  review_process/se_coverage_2026-09-14
measured 0.66-0.75 on DGPs with tau constant within cohort -> that is the middle arm.

## r11_gtcell_se.do  (60 draws, homogeneous, SELF-CONTAINED)
att   TRUE sd 0.15138  mean reported se 0.14064  se/sd 0.929
cell  TRUE sd 0.74714  mean reported se 0.21431  se/sd 0.287
cell 95% interval coverage of its own true effect = 0.450 (nominal 0.95)
No caution of any kind printed on those runs (grep -ci caution = 0).

## r06 coverage MC final-ish (59-60 paired panels)
standard pct 55/59 (.932) width 0.6422 se/truesd .906
bayes    pct 53/59 (.898) width 0.5977 se/truesd .847
analytic se/truesd .849 ; discordance 3-1 (sign test p .31)
-> bb-03(b)'s 16.7pp gap and 5/5 discordance do NOT reproduce; the gap is ~3.4pp and the
   bayes interval is NARROWER (width ratio .925), not shifted.

## r12_gtcell_finer.do  (scope control for N1: i = cohort, cluster(state) FINER than i)
(pending)

## r12_gtcell_finer.do  RESULT (scope control for N1)
Identical data and identical cell, i = onset cohort instead of state, cluster(st) FINER than i:
  att   TRUE sd 0.15138  mean reported se 0.15400  se/sd 1.017
  cell  TRUE sd 0.31685  mean reported se 0.32053  se/sd 1.012
  cell 95% interval coverage of its own true effect = 0.950 (nominal 0.95)  EXACTLY NOMINAL.
(The log shows the block twice: the file was launched once directly and once by a waiting
watcher loop, so two Stata processes ran it; both runs agree, the second simply read a
120-row postfile.  Either block is the result.)

## r06 FINAL (60 paired panels)
standard pct 56/60 (.933) width 0.64209 se/truesd .907
bayes    pct 54/60 (.900) width 0.59762 se/truesd .848 ; analytic se/truesd .850
discordance 3-1 (sign test p .31)
