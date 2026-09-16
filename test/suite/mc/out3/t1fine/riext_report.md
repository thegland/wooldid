# Randomization-inference extensions (mc_ri_ext)

outdir `out3/t1fine`; 4 shard file(s); 400 rows after de-duplication on (phase, arm, tau, shard, rep); 0 malformed line(s) dropped; 0 duplicate row(s) dropped.

The headline is SIZE: on every tau = 0 phase each p column must be uniform on its lattice, so the rejection count at a nominal alpha sits in the central ~95% region of Binomial(n, ATTAINED alpha). The attained alpha is enumerated from the run's own permreps (analyze_mc.py's DEVIATION D1/D2 arithmetic): at B = 199 the equal-tailed default attains .05 and .10 exactly; at small B it may not reach .05 at all, and a claim it cannot express is NOT JUDGED rather than failed.

`npos` is the size of the assignment space: a small space puts a floor under the achievable p. `pconv` names the p-value convention in force. `namesok` = 0 means a surface was read by POSITION after a name lookup failed - the numbers may be the wrong column's.

Cells with fewer than MINN=30 reps are reported, never judged.

| check | verdict | detail |
|---|---|---|
| rc == 0 on every row | PASS | rc!=0 on 0/400; threshold 0.00 (driver ACCEPTANCE: rc = 0 on every row) |
| permreps in force | NOTE | B = 199 (argument); attained alpha 0.05->0.0500 (eq) / 0.0500 (right), 0.10->0.1000 (eq) / 0.1000 (right) |

## Run inventory

| phase | arm | what | rows (rc=0) | rc!=0 | tau | pconv | level | npos | ncomp (modal) | nfail>0 | namesok=0 | mean secs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E2 | right | permside(right): one-sided add-one p | 200 | 0 | 0 | right | 95 | 18480 | 199 | 0 | 0 | 46.38 |
| E6 | fineclus | cluster() strictly finer than i under permonsetcohorts (T1) | 200 | 0 | 0 | equaltailed | 95 | 18480 | 199 | 0 | 0 | 28.12 |

## 0. Integrity (per-row facts; judged at any n - these are not statistical)

| check | verdict | detail |
|---|---|---|
| E2 ncomp == B (199) | PASS | ncomp != 199 on 0/200 (a lost replication puts that rep's p on a coarser lattice than the one judged here) |
| E2 nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw; driver ACCEPTANCE: nfail = 0 on every row) |
| E2 ndraw == B (199) for main_att | PASS | ndraw != 199 on 0/200 (main_att's accumulator lost draws) |
| E2 namesok == 1 | PASS | namesok = 0 on 0/200 (a surface read by POSITION after a name lookup failed, so the recorded p may be a different column's) |
| E2 npos >= 3000 | PASS | assignment space [18480]; the PREREGISTERED floor is ~3000 (PREREG_NEWFEAT_2026-08-24.md: spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p has a floor that is not the estimator's fault |
| E2 pconv == 'right' | PASS | pconv values off 'right': none (the convention decides the lattice and therefore the attainable alpha) |
| E2 level == 95 | PASS | level values off 95: none |
| E2 tau == 0 (it is a SIZE arm) | PASS | tau values off 0: none (a nonzero tau here would make every 'size' number below a power number instead) |
| E2 every p in [0,1] | PASS | no p outside [0,1] |
| E2 row count == M (200) | PASS | 200 rc-0 rows of 200 planned |
| E6 ncomp == B (199) | PASS | ncomp != 199 on 0/200 (a lost replication puts that rep's p on a coarser lattice than the one judged here) |
| E6 nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw; driver ACCEPTANCE: nfail = 0 on every row) |
| E6 ndraw == B (199) for main_att | PASS | ndraw != 199 on 0/200 (main_att's accumulator lost draws) |
| E6 namesok == 1 | PASS | namesok = 0 on 0/200 (a surface read by POSITION after a name lookup failed, so the recorded p may be a different column's) |
| E6 npos >= 3000 | PASS | assignment space [18480]; the PREREGISTERED floor is ~3000 (PREREG_NEWFEAT_2026-08-24.md: spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p has a floor that is not the estimator's fault |
| E6 pconv == 'equaltailed' | PASS | pconv values off 'equaltailed': none (the convention decides the lattice and therefore the attainable alpha) |
| E6 level == 95 | PASS | level values off 95: none |
| E6 tau == 0 (it is a SIZE arm) | PASS | tau values off 0: none (a nonzero tau here would make every 'size' number below a power number instead) |
| E6 every p in [0,1] | PASS | no p outside [0,1] |
| E6 row count == M (200) | PASS | 200 rc-0 rows of 200 planned |

## 1. SIZE under tau = 0 - the headline

(att and, on E1, ame surfaces plus the joint p's. `attained` is P(p <= alpha) on this run's lattice; the count band is the central ~95% region of Binomial(n, attained). chi2 is the 10-bin uniformity statistic, df=9, crit 21.666 at p=.01 - a secondary flag, never a verdict.)

| phase | column | surface | conv | alpha | attained | rej (n, Wilson 95%) | count band | chi2 |
|---|---|---|---|---|---|---|---|---|
| E2 | p_main | main_att | right | 0.05 | 0.0500 | 14/200 = 0.0700 [0.0422,0.1141] | [3,17] | 8.6 |
| E2 | p_main | main_att | right | 0.10 | 0.1000 | 22/200 = 0.1100 [0.0738,0.1609] | [11,29] | 8.6 |
| E2 | p_pre | pretreatment_att | right | 0.05 | 0.0500 | 6/200 = 0.0300 [0.0138,0.0639] | [3,17] | 8.5 |
| E2 | p_pre | pretreatment_att | right | 0.10 | 0.1000 | 16/200 = 0.0800 [0.0498,0.1260] | [11,29] | 8.5 |
| E2 | p_post0 | ES_att_post0 | right | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 4.5 |
| E2 | p_post0 | ES_att_post0 | right | 0.10 | 0.1000 | 18/200 = 0.0900 [0.0577,0.1378] | [11,29] | 4.5 |
| E2 | p_post1 | ES_att_post1 | right | 0.05 | 0.0500 | 7/200 = 0.0350 [0.0171,0.0705] | [3,17] | 5.8 |
| E2 | p_post1 | ES_att_post1 | right | 0.10 | 0.1000 | 17/200 = 0.0850 [0.0537,0.1319] | [11,29] | 5.8 |
| E2 | jointp | e(joint_horiz) | right | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 5.0 |
| E2 | jointp | e(joint_horiz) | right | 0.10 | 0.1000 | 18/200 = 0.0900 [0.0577,0.1378] | [11,29] | 5.0 |
| E6 | p_main | main_att | equaltailed | 0.05 | 0.0500 | 14/200 = 0.0700 [0.0422,0.1141] | [3,17] | 10.4 |
| E6 | p_main | main_att | equaltailed | 0.10 | 0.1000 | 26/200 = 0.1300 [0.0903,0.1837] | [11,29] | 10.4 |
| E6 | p_pre | pretreatment_att | equaltailed | 0.05 | 0.0500 | 13/200 = 0.0650 [0.0384,0.1080] | [3,17] | 19.5 |
| E6 | p_pre | pretreatment_att | equaltailed | 0.10 | 0.1000 | 17/200 = 0.0850 [0.0537,0.1319] | [11,29] | 19.5 |
| E6 | p_post0 | ES_att_post0 | equaltailed | 0.05 | 0.0500 | 12/200 = 0.0600 [0.0347,0.1019] | [3,17] | 3.4 |
| E6 | p_post0 | ES_att_post0 | equaltailed | 0.10 | 0.1000 | 22/200 = 0.1100 [0.0738,0.1609] | [11,29] | 3.4 |
| E6 | p_post1 | ES_att_post1 | equaltailed | 0.05 | 0.0500 | 9/200 = 0.0450 [0.0239,0.0833] | [3,17] | 3.5 |
| E6 | p_post1 | ES_att_post1 | equaltailed | 0.10 | 0.1000 | 22/200 = 0.1100 [0.0738,0.1609] | [11,29] | 3.5 |
| E6 | jointp | e(joint_horiz) | right | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 5.0 |
| E6 | jointp | e(joint_horiz) | right | 0.10 | 0.1000 | 18/200 = 0.0900 [0.0577,0.1378] | [11,29] | 5.0 |

| check | verdict | detail |
|---|---|---|
| E2 p_main (main_att) size | PASS | alpha 0.05 (attained 0.0500): 14/200 = 0.0700 [0.0422,0.1141] band [3,17]; alpha 0.10 (attained 0.1000): 22/200 = 0.1100 [0.0738,0.1609] band [11,29]; chi2=8.6 |
| E2 p_pre (pretreatment_att) size | PASS | alpha 0.05 (attained 0.0500): 6/200 = 0.0300 [0.0138,0.0639] band [3,17]; alpha 0.10 (attained 0.1000): 16/200 = 0.0800 [0.0498,0.1260] band [11,29]; chi2=8.5 |
| E2 p_post0 (ES_att_post0) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 18/200 = 0.0900 [0.0577,0.1378] band [11,29]; chi2=4.5 |
| E2 p_post1 (ES_att_post1) size | PASS | alpha 0.05 (attained 0.0500): 7/200 = 0.0350 [0.0171,0.0705] band [3,17]; alpha 0.10 (attained 0.1000): 17/200 = 0.0850 [0.0537,0.1319] band [11,29]; chi2=5.8 |
| E2 jointp (e(joint_horiz)) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 18/200 = 0.0900 [0.0577,0.1378] band [11,29]; chi2=5.0 |
| E6 p_main (main_att) size | PASS | alpha 0.05 (attained 0.0500): 14/200 = 0.0700 [0.0422,0.1141] band [3,17]; alpha 0.10 (attained 0.1000): 26/200 = 0.1300 [0.0903,0.1837] band [11,29]; chi2=10.4 |
| E6 p_pre (pretreatment_att) size | PASS | alpha 0.05 (attained 0.0500): 13/200 = 0.0650 [0.0384,0.1080] band [3,17]; alpha 0.10 (attained 0.1000): 17/200 = 0.0850 [0.0537,0.1319] band [11,29]; chi2=19.5 |
| E6 p_post0 (ES_att_post0) size | PASS | alpha 0.05 (attained 0.0500): 12/200 = 0.0600 [0.0347,0.1019] band [3,17]; alpha 0.10 (attained 0.1000): 22/200 = 0.1100 [0.0738,0.1609] band [11,29]; chi2=3.4 |
| E6 p_post1 (ES_att_post1) size | PASS | alpha 0.05 (attained 0.0500): 9/200 = 0.0450 [0.0239,0.0833] band [3,17]; alpha 0.10 (attained 0.1000): 22/200 = 0.1100 [0.0738,0.1609] band [11,29]; chi2=3.5 |
| E6 jointp (e(joint_horiz)) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 18/200 = 0.0900 [0.0577,0.1378] band [11,29]; chi2=5.0 |

## 2. Lattice check (a p off the lattice means that column lost draws)

(analyze_mcri2.py's check: with all B draws complete an equal-tailed p is an exact multiple of 2/200 and a one-sided p of 1/200. Off-lattice on main_att is a FAILURE - the exactness argument does not hold for the headline column; off-lattice elsewhere is a FINDING, because a horizon can be unidentified in some reassignments.)

| phase | column | on lattice | n | off-lattice examples |
|---|---|---|---|---|
| E2 | p_main | 200 | 200 | - |
| E2 | p_pre | 200 | 200 | - |
| E2 | p_post0 | 200 | 200 | - |
| E2 | p_post1 | 200 | 200 | - |
| E2 | jointp | 200 | 200 | - |
| E6 | p_main | 200 | 200 | - |
| E6 | p_pre | 200 | 200 | - |
| E6 | p_post0 | 200 | 200 | - |
| E6 | p_post1 | 200 | 200 | - |
| E6 | jointp | 200 | 200 | - |

| check | verdict | detail |
|---|---|---|
| E2 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |
| E6 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |

## 3. E3 - COVERAGE of the planted tau by e(perm_ci) at level(90)

(The only level() != 95 arm in any MC. A blank interval is the off-centre case mc_ri2 documents: the null band excludes 0, so no shift interval exists - those reps are excluded from the coverage denominator and reported separately.)

No E3 rows: the level(90) coverage arm did not run.
| check | verdict | detail |
|---|---|---|
| E3 present | NOTE | no rows for the level(90) arm in this outdir |

## 4. Secondary (information only, never fails the component)

| phase | mean p_main | mean p_pre | mean p_post0 | mean p_post1 | mean sd_main | mean b_main | n |
|---|---|---|---|---|---|---|---|
| E2 | 0.5258 | 0.5123 | 0.5227 | 0.5303 | 0.4246 | -0.0441 | 200 |
| E6 | 0.4989 | 0.5495 | 0.4951 | 0.4915 | 0.4246 | -0.0441 | 200 |

Under a correct null the mean equal-tailed p is not .5 but the lattice mean (0.5050 at B = 199); the one-sided mean is 0.5025. Read the means against those, not against .5.


## Verdict summary

- checks FAILED: 0
- claims NOT JUDGED (alpha unattainable at B=199): 0
- checks not judged for want of reps (MINN=30): 0
