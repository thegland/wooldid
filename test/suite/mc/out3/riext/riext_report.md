# Randomization-inference extensions (mc_ri_ext)

outdir `out3/riext`; 8 shard file(s); 1000 rows after de-duplication on (phase, arm, tau, shard, rep); 0 malformed line(s) dropped; 0 duplicate row(s) dropped.

The headline is SIZE: on every tau = 0 phase each p column must be uniform on its lattice, so the rejection count at a nominal alpha sits in the central ~95% region of Binomial(n, ATTAINED alpha). The attained alpha is enumerated from the run's own permreps (analyze_mc.py's DEVIATION D1/D2 arithmetic): at B = 199 the equal-tailed default attains .05 and .10 exactly; at small B it may not reach .05 at all, and a claim it cannot express is NOT JUDGED rather than failed.

`npos` is the size of the assignment space: a small space puts a floor under the achievable p. `pconv` names the p-value convention in force. `namesok` = 0 means a surface was read by POSITION after a name lookup failed - the numbers may be the wrong column's.

Cells with fewer than MINN=30 reps are reported, never judged.

| check | verdict | detail |
|---|---|---|
| rc == 0 on every row | PASS | rc!=0 on 0/1000; threshold 0.00 (driver ACCEPTANCE: rc = 0 on every row) |
| permreps in force | NOTE | B = 199 (argument); attained alpha 0.05->0.0500 (eq) / 0.0500 (right), 0.10->0.1000 (eq) / 0.1000 (right) |

## Run inventory

| phase | arm | what | rows (rc=0) | rc!=0 | tau | pconv | level | npos | ncomp (modal) | nfail>0 | namesok=0 | mean secs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E1 | contreat | contreat(dose) + contreatcontrols(basic), sharp null | 200 | 0 | 0 | equaltailed | 95 | 18480 | 199 | 0 | 0 | 209.85 |
| E2 | right | permside(right): one-sided add-one p | 200 | 0 | 0 | right | 95 | 18480 | 199 | 0 | 0 | 53.58 |
| E3 | level90 | level(90): COVERAGE of the planted tau by e(perm_ci) | 200 | 0 | 0 | equaltailed | 90 | 18480 | 199 | 0 | 0 | 50.96 |
| E4 | strata | permstrata(str): exact within strata | 200 | 0 | 0 | equaltailed | 95 | 3136 | 199 | 0 | 0 | 54.45 |
| E5 | tstat0 | permtstat at tau = 0 (the size companion of mc_ri2 D4) | 200 | 0 | 0 | equaltailed | 95 | 18480 | 199 | 0 | 0 | 70.63 |

## 0. Integrity (per-row facts; judged at any n - these are not statistical)

| check | verdict | detail |
|---|---|---|
| E1 ncomp == B (199) | PASS | ncomp != 199 on 0/200 (a lost replication puts that rep's p on a coarser lattice than the one judged here) |
| E1 nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw; driver ACCEPTANCE: nfail = 0 on every row) |
| E1 ndraw == B (199) for main_att | PASS | ndraw != 199 on 0/200 (main_att's accumulator lost draws) |
| E1 namesok == 1 | PASS | namesok = 0 on 0/200 (a surface read by POSITION after a name lookup failed, so the recorded p may be a different column's) |
| E1 npos >= 3000 | PASS | assignment space [18480]; the PREREGISTERED floor is ~3000 (PREREG_NEWFEAT_2026-08-24.md: spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p has a floor that is not the estimator's fault |
| E1 pconv == 'equaltailed' | PASS | pconv values off 'equaltailed': none (the convention decides the lattice and therefore the attainable alpha) |
| E1 level == 95 | PASS | level values off 95: none |
| E1 tau == 0 (it is a SIZE arm) | PASS | tau values off 0: none (a nonzero tau here would make every 'size' number below a power number instead) |
| E1 every p in [0,1] | PASS | no p outside [0,1] |
| E1 row count == M (200) | PASS | 200 rc-0 rows of 200 planned |
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
| E3 ncomp == B (199) | PASS | ncomp != 199 on 0/200 (a lost replication puts that rep's p on a coarser lattice than the one judged here) |
| E3 nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw; driver ACCEPTANCE: nfail = 0 on every row) |
| E3 ndraw == B (199) for main_att | PASS | ndraw != 199 on 0/200 (main_att's accumulator lost draws) |
| E3 namesok == 1 | PASS | namesok = 0 on 0/200 (a surface read by POSITION after a name lookup failed, so the recorded p may be a different column's) |
| E3 npos >= 3000 | PASS | assignment space [18480]; the PREREGISTERED floor is ~3000 (PREREG_NEWFEAT_2026-08-24.md: spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p has a floor that is not the estimator's fault |
| E3 pconv == 'equaltailed' | PASS | pconv values off 'equaltailed': none (the convention decides the lattice and therefore the attainable alpha) |
| E3 level == 90 | PASS | level values off 90: none |
| E3 every p in [0,1] | PASS | no p outside [0,1] |
| E3 row count == M (200) | PASS | 200 rc-0 rows of 200 planned |
| E4 ncomp == B (199) | PASS | ncomp != 199 on 0/200 (a lost replication puts that rep's p on a coarser lattice than the one judged here) |
| E4 nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw; driver ACCEPTANCE: nfail = 0 on every row) |
| E4 ndraw == B (199) for main_att | PASS | ndraw != 199 on 0/200 (main_att's accumulator lost draws) |
| E4 namesok == 1 | PASS | namesok = 0 on 0/200 (a surface read by POSITION after a name lookup failed, so the recorded p may be a different column's) |
| E4 npos >= 3000 | PASS | assignment space [3136]; the PREREGISTERED floor is ~3000 (PREREG_NEWFEAT_2026-08-24.md: spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p has a floor that is not the estimator's fault |
| E4 pconv == 'equaltailed' | PASS | pconv values off 'equaltailed': none (the convention decides the lattice and therefore the attainable alpha) |
| E4 level == 95 | PASS | level values off 95: none |
| E4 tau == 0 (it is a SIZE arm) | PASS | tau values off 0: none (a nonzero tau here would make every 'size' number below a power number instead) |
| E4 every p in [0,1] | PASS | no p outside [0,1] |
| E4 row count == M (200) | PASS | 200 rc-0 rows of 200 planned |
| E5 ncomp == B (199) | PASS | ncomp != 199 on 0/200 (a lost replication puts that rep's p on a coarser lattice than the one judged here) |
| E5 nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw; driver ACCEPTANCE: nfail = 0 on every row) |
| E5 ndraw == B (199) for main_att | PASS | ndraw != 199 on 0/200 (main_att's accumulator lost draws) |
| E5 namesok == 1 | PASS | namesok = 0 on 0/200 (a surface read by POSITION after a name lookup failed, so the recorded p may be a different column's) |
| E5 npos >= 3000 | PASS | assignment space [18480]; the PREREGISTERED floor is ~3000 (PREREG_NEWFEAT_2026-08-24.md: spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p has a floor that is not the estimator's fault |
| E5 pconv == 'equaltailed' | PASS | pconv values off 'equaltailed': none (the convention decides the lattice and therefore the attainable alpha) |
| E5 level == 95 | PASS | level values off 95: none |
| E5 tau == 0 (it is a SIZE arm) | PASS | tau values off 0: none (a nonzero tau here would make every 'size' number below a power number instead) |
| E5 every p in [0,1] | PASS | no p outside [0,1] |
| E5 row count == M (200) | PASS | 200 rc-0 rows of 200 planned |

## 1. SIZE under tau = 0 - the headline

(att and, on E1, ame surfaces plus the joint p's. `attained` is P(p <= alpha) on this run's lattice; the count band is the central ~95% region of Binomial(n, attained). chi2 is the 10-bin uniformity statistic, df=9, crit 21.666 at p=.01 - a secondary flag, never a verdict.)

| phase | column | surface | conv | alpha | attained | rej (n, Wilson 95%) | count band | chi2 |
|---|---|---|---|---|---|---|---|---|
| E1 | p_main | main_att | equaltailed | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 9.4 |
| E1 | p_main | main_att | equaltailed | 0.10 | 0.1000 | 20/200 = 0.1000 [0.0657,0.1494] | [11,29] | 9.4 |
| E1 | p_pre | pretreatment_att | equaltailed | 0.05 | 0.0500 | 6/200 = 0.0300 [0.0138,0.0639] | [3,17] | 14.3 |
| E1 | p_pre | pretreatment_att | equaltailed | 0.10 | 0.1000 | 19/200 = 0.0950 [0.0617,0.1436] | [11,29] | 14.3 |
| E1 | p_post0 | ES_att_post0 | equaltailed | 0.05 | 0.0500 | 12/200 = 0.0600 [0.0347,0.1019] | [3,17] | 8.9 |
| E1 | p_post0 | ES_att_post0 | equaltailed | 0.10 | 0.1000 | 20/200 = 0.1000 [0.0657,0.1494] | [11,29] | 8.9 |
| E1 | p_post1 | ES_att_post1 | equaltailed | 0.05 | 0.0500 | 13/200 = 0.0650 [0.0384,0.1080] | [3,17] | 16.8 |
| E1 | p_post1 | ES_att_post1 | equaltailed | 0.10 | 0.1000 | 36/200 = 0.1800 [0.1329,0.2391] | [11,29] | 16.8 |
| E1 | p_ame | main_ame | equaltailed | 0.05 | 0.0500 | 15/200 = 0.0750 [0.0460,0.1200] | [3,17] | 11.1 |
| E1 | p_ame | main_ame | equaltailed | 0.10 | 0.1000 | 29/200 = 0.1450 [0.1029,0.2005] | [11,29] | 11.1 |
| E1 | p_amepost0 | ES_ame_post0 | equaltailed | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 13.1 |
| E1 | p_amepost0 | ES_ame_post0 | equaltailed | 0.10 | 0.1000 | 31/200 = 0.1550 [0.1114,0.2116] | [11,29] | 13.1 |
| E1 | p_amepost1 | ES_ame_post1 | equaltailed | 0.05 | 0.0500 | 7/200 = 0.0350 [0.0171,0.0705] | [3,17] | 12.6 |
| E1 | p_amepost1 | ES_ame_post1 | equaltailed | 0.10 | 0.1000 | 18/200 = 0.0900 [0.0577,0.1378] | [11,29] | 12.6 |
| E1 | jointp | e(joint_horiz) | right | 0.05 | 0.0500 | 9/200 = 0.0450 [0.0239,0.0833] | [3,17] | 12.8 |
| E1 | jointp | e(joint_horiz) | right | 0.10 | 0.1000 | 19/200 = 0.0950 [0.0617,0.1436] | [11,29] | 12.8 |
| E1 | jointCp | e(joint_Choriz) | right | 0.05 | 0.0500 | 9/200 = 0.0450 [0.0239,0.0833] | [3,17] | 20.1 |
| E1 | jointCp | e(joint_Choriz) | right | 0.10 | 0.1000 | 21/200 = 0.1050 [0.0697,0.1552] | [11,29] | 20.1 |
| E2 | p_main | main_att | right | 0.05 | 0.0500 | 12/200 = 0.0600 [0.0347,0.1019] | [3,17] | 9.0 |
| E2 | p_main | main_att | right | 0.10 | 0.1000 | 23/200 = 0.1150 [0.0779,0.1666] | [11,29] | 9.0 |
| E2 | p_pre | pretreatment_att | right | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 7.1 |
| E2 | p_pre | pretreatment_att | right | 0.10 | 0.1000 | 18/200 = 0.0900 [0.0577,0.1378] | [11,29] | 7.1 |
| E2 | p_post0 | ES_att_post0 | right | 0.05 | 0.0500 | 7/200 = 0.0350 [0.0171,0.0705] | [3,17] | 4.3 |
| E2 | p_post0 | ES_att_post0 | right | 0.10 | 0.1000 | 15/200 = 0.0750 [0.0460,0.1200] | [11,29] | 4.3 |
| E2 | p_post1 | ES_att_post1 | right | 0.05 | 0.0500 | 5/200 = 0.0250 [0.0107,0.0572] | [3,17] | 12.6 |
| E2 | p_post1 | ES_att_post1 | right | 0.10 | 0.1000 | 14/200 = 0.0700 [0.0422,0.1141] | [11,29] | 12.6 |
| E2 | jointp | e(joint_horiz) | right | 0.05 | 0.0500 | 12/200 = 0.0600 [0.0347,0.1019] | [3,17] | 6.9 |
| E2 | jointp | e(joint_horiz) | right | 0.10 | 0.1000 | 22/200 = 0.1100 [0.0738,0.1609] | [11,29] | 6.9 |
| E4 | p_main | main_att | equaltailed | 0.05 | 0.0500 | 9/200 = 0.0450 [0.0239,0.0833] | [3,17] | 13.2 |
| E4 | p_main | main_att | equaltailed | 0.10 | 0.1000 | 20/200 = 0.1000 [0.0657,0.1494] | [11,29] | 13.2 |
| E4 | p_pre | pretreatment_att | equaltailed | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 1.9 |
| E4 | p_pre | pretreatment_att | equaltailed | 0.10 | 0.1000 | 19/200 = 0.0950 [0.0617,0.1436] | [11,29] | 1.9 |
| E4 | p_post0 | ES_att_post0 | equaltailed | 0.05 | 0.0500 | 11/200 = 0.0550 [0.0310,0.0958] | [3,17] | 16.3 |
| E4 | p_post0 | ES_att_post0 | equaltailed | 0.10 | 0.1000 | 23/200 = 0.1150 [0.0779,0.1666] | [11,29] | 16.3 |
| E4 | p_post1 | ES_att_post1 | equaltailed | 0.05 | 0.0500 | 7/200 = 0.0350 [0.0171,0.0705] | [3,17] | 13.2 |
| E4 | p_post1 | ES_att_post1 | equaltailed | 0.10 | 0.1000 | 16/200 = 0.0800 [0.0498,0.1260] | [11,29] | 13.2 |
| E4 | jointp | e(joint_horiz) | right | 0.05 | 0.0500 | 11/200 = 0.0550 [0.0310,0.0958] | [3,17] | 6.6 |
| E4 | jointp | e(joint_horiz) | right | 0.10 | 0.1000 | 17/200 = 0.0850 [0.0537,0.1319] | [11,29] | 6.6 |
| E5 | p_main | main_att | equaltailed | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 8.8 |
| E5 | p_main | main_att | equaltailed | 0.10 | 0.1000 | 19/200 = 0.0950 [0.0617,0.1436] | [11,29] | 8.8 |
| E5 | p_pre | pretreatment_att | equaltailed | 0.05 | 0.0500 | 9/200 = 0.0450 [0.0239,0.0833] | [3,17] | 15.7 |
| E5 | p_pre | pretreatment_att | equaltailed | 0.10 | 0.1000 | 19/200 = 0.0950 [0.0617,0.1436] | [11,29] | 15.7 |
| E5 | p_post0 | ES_att_post0 | equaltailed | 0.05 | 0.0500 | 11/200 = 0.0550 [0.0310,0.0958] | [3,17] | 5.1 |
| E5 | p_post0 | ES_att_post0 | equaltailed | 0.10 | 0.1000 | 24/200 = 0.1200 [0.0820,0.1723] | [11,29] | 5.1 |
| E5 | p_post1 | ES_att_post1 | equaltailed | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 13.9 |
| E5 | p_post1 | ES_att_post1 | equaltailed | 0.10 | 0.1000 | 15/200 = 0.0750 [0.0460,0.1200] | [11,29] | 13.9 |
| E5 | jointp | e(joint_horiz) | right | 0.05 | 0.0500 | 12/200 = 0.0600 [0.0347,0.1019] | [3,17] | 6.9 |
| E5 | jointp | e(joint_horiz) | right | 0.10 | 0.1000 | 22/200 = 0.1100 [0.0738,0.1609] | [11,29] | 6.9 |

| check | verdict | detail |
|---|---|---|
| E1 p_main (main_att) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 20/200 = 0.1000 [0.0657,0.1494] band [11,29]; chi2=9.4 |
| E1 p_pre (pretreatment_att) size | PASS | alpha 0.05 (attained 0.0500): 6/200 = 0.0300 [0.0138,0.0639] band [3,17]; alpha 0.10 (attained 0.1000): 19/200 = 0.0950 [0.0617,0.1436] band [11,29]; chi2=14.3 |
| E1 p_post0 (ES_att_post0) size | PASS | alpha 0.05 (attained 0.0500): 12/200 = 0.0600 [0.0347,0.1019] band [3,17]; alpha 0.10 (attained 0.1000): 20/200 = 0.1000 [0.0657,0.1494] band [11,29]; chi2=8.9 |
| E1 p_post1 (ES_att_post1) size | FAIL | alpha 0.05 (attained 0.0500): 13/200 = 0.0650 [0.0384,0.1080] band [3,17]; alpha 0.10 (attained 0.1000): 36/200 = 0.1800 [0.1329,0.2391] band [11,29]; chi2=16.8 |
| E1 p_ame (main_ame) size | PASS | alpha 0.05 (attained 0.0500): 15/200 = 0.0750 [0.0460,0.1200] band [3,17]; alpha 0.10 (attained 0.1000): 29/200 = 0.1450 [0.1029,0.2005] band [11,29]; chi2=11.1 |
| E1 p_amepost0 (ES_ame_post0) size | FAIL | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 31/200 = 0.1550 [0.1114,0.2116] band [11,29]; chi2=13.1 |
| E1 p_amepost1 (ES_ame_post1) size | PASS | alpha 0.05 (attained 0.0500): 7/200 = 0.0350 [0.0171,0.0705] band [3,17]; alpha 0.10 (attained 0.1000): 18/200 = 0.0900 [0.0577,0.1378] band [11,29]; chi2=12.6 |
| E1 jointp (e(joint_horiz)) size | PASS | alpha 0.05 (attained 0.0500): 9/200 = 0.0450 [0.0239,0.0833] band [3,17]; alpha 0.10 (attained 0.1000): 19/200 = 0.0950 [0.0617,0.1436] band [11,29]; chi2=12.8 |
| E1 jointCp (e(joint_Choriz)) size | PASS | alpha 0.05 (attained 0.0500): 9/200 = 0.0450 [0.0239,0.0833] band [3,17]; alpha 0.10 (attained 0.1000): 21/200 = 0.1050 [0.0697,0.1552] band [11,29]; chi2=20.1 |
| E2 p_main (main_att) size | PASS | alpha 0.05 (attained 0.0500): 12/200 = 0.0600 [0.0347,0.1019] band [3,17]; alpha 0.10 (attained 0.1000): 23/200 = 0.1150 [0.0779,0.1666] band [11,29]; chi2=9.0 |
| E2 p_pre (pretreatment_att) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 18/200 = 0.0900 [0.0577,0.1378] band [11,29]; chi2=7.1 |
| E2 p_post0 (ES_att_post0) size | PASS | alpha 0.05 (attained 0.0500): 7/200 = 0.0350 [0.0171,0.0705] band [3,17]; alpha 0.10 (attained 0.1000): 15/200 = 0.0750 [0.0460,0.1200] band [11,29]; chi2=4.3 |
| E2 p_post1 (ES_att_post1) size | PASS | alpha 0.05 (attained 0.0500): 5/200 = 0.0250 [0.0107,0.0572] band [3,17]; alpha 0.10 (attained 0.1000): 14/200 = 0.0700 [0.0422,0.1141] band [11,29]; chi2=12.6 |
| E2 jointp (e(joint_horiz)) size | PASS | alpha 0.05 (attained 0.0500): 12/200 = 0.0600 [0.0347,0.1019] band [3,17]; alpha 0.10 (attained 0.1000): 22/200 = 0.1100 [0.0738,0.1609] band [11,29]; chi2=6.9 |
| E4 p_main (main_att) size | PASS | alpha 0.05 (attained 0.0500): 9/200 = 0.0450 [0.0239,0.0833] band [3,17]; alpha 0.10 (attained 0.1000): 20/200 = 0.1000 [0.0657,0.1494] band [11,29]; chi2=13.2 |
| E4 p_pre (pretreatment_att) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 19/200 = 0.0950 [0.0617,0.1436] band [11,29]; chi2=1.9 |
| E4 p_post0 (ES_att_post0) size | PASS | alpha 0.05 (attained 0.0500): 11/200 = 0.0550 [0.0310,0.0958] band [3,17]; alpha 0.10 (attained 0.1000): 23/200 = 0.1150 [0.0779,0.1666] band [11,29]; chi2=16.3 |
| E4 p_post1 (ES_att_post1) size | PASS | alpha 0.05 (attained 0.0500): 7/200 = 0.0350 [0.0171,0.0705] band [3,17]; alpha 0.10 (attained 0.1000): 16/200 = 0.0800 [0.0498,0.1260] band [11,29]; chi2=13.2 |
| E4 jointp (e(joint_horiz)) size | PASS | alpha 0.05 (attained 0.0500): 11/200 = 0.0550 [0.0310,0.0958] band [3,17]; alpha 0.10 (attained 0.1000): 17/200 = 0.0850 [0.0537,0.1319] band [11,29]; chi2=6.6 |
| E5 p_main (main_att) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 19/200 = 0.0950 [0.0617,0.1436] band [11,29]; chi2=8.8 |
| E5 p_pre (pretreatment_att) size | PASS | alpha 0.05 (attained 0.0500): 9/200 = 0.0450 [0.0239,0.0833] band [3,17]; alpha 0.10 (attained 0.1000): 19/200 = 0.0950 [0.0617,0.1436] band [11,29]; chi2=15.7 |
| E5 p_post0 (ES_att_post0) size | PASS | alpha 0.05 (attained 0.0500): 11/200 = 0.0550 [0.0310,0.0958] band [3,17]; alpha 0.10 (attained 0.1000): 24/200 = 0.1200 [0.0820,0.1723] band [11,29]; chi2=5.1 |
| E5 p_post1 (ES_att_post1) size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 15/200 = 0.0750 [0.0460,0.1200] band [11,29]; chi2=13.9 |
| E5 jointp (e(joint_horiz)) size | PASS | alpha 0.05 (attained 0.0500): 12/200 = 0.0600 [0.0347,0.1019] band [3,17]; alpha 0.10 (attained 0.1000): 22/200 = 0.1100 [0.0738,0.1609] band [11,29]; chi2=6.9 |

## 2. Lattice check (a p off the lattice means that column lost draws)

(analyze_mcri2.py's check: with all B draws complete an equal-tailed p is an exact multiple of 2/200 and a one-sided p of 1/200. Off-lattice on main_att is a FAILURE - the exactness argument does not hold for the headline column; off-lattice elsewhere is a FINDING, because a horizon can be unidentified in some reassignments.)

| phase | column | on lattice | n | off-lattice examples |
|---|---|---|---|---|
| E1 | p_main | 200 | 200 | - |
| E1 | p_pre | 200 | 200 | - |
| E1 | p_post0 | 200 | 200 | - |
| E1 | p_post1 | 200 | 200 | - |
| E1 | p_ame | 200 | 200 | - |
| E1 | p_amepost0 | 200 | 200 | - |
| E1 | p_amepost1 | 200 | 200 | - |
| E1 | jointp | 200 | 200 | - |
| E1 | jointCp | 200 | 200 | - |
| E2 | p_main | 200 | 200 | - |
| E2 | p_pre | 200 | 200 | - |
| E2 | p_post0 | 200 | 200 | - |
| E2 | p_post1 | 200 | 200 | - |
| E2 | jointp | 200 | 200 | - |
| E3 | p_main | 200 | 200 | - |
| E3 | p_pre | 200 | 200 | - |
| E3 | p_post0 | 200 | 200 | - |
| E3 | p_post1 | 200 | 200 | - |
| E3 | jointp | 200 | 200 | - |
| E4 | p_main | 200 | 200 | - |
| E4 | p_pre | 200 | 200 | - |
| E4 | p_post0 | 200 | 200 | - |
| E4 | p_post1 | 200 | 200 | - |
| E4 | jointp | 200 | 200 | - |
| E5 | p_main | 200 | 200 | - |
| E5 | p_pre | 200 | 200 | - |
| E5 | p_post0 | 200 | 200 | - |
| E5 | p_post1 | 200 | 200 | - |
| E5 | jointp | 200 | 200 | - |

| check | verdict | detail |
|---|---|---|
| E1 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |
| E2 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |
| E3 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |
| E4 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |
| E5 p_main on the B=199 lattice | PASS | off-lattice on 0/200 reps (that column lost replications, so its p is not exact and the size verdict above is conditional) |

## 3. E3 - COVERAGE of the planted tau by e(perm_ci) at level(90)

(The only level() != 95 arm in any MC. A blank interval is the off-centre case mc_ri2 documents: the null band excludes 0, so no shift interval exists - those reps are excluded from the coverage denominator and reported separately.)

| quantity | covered (n, Wilson 95%) | count band | mean width | blank |
|---|---|---|---|---|
| main_att, tau = 0.0 | 177/200 = 0.8850 [0.8334,0.9221] | [171,189] | 1.4179 | 0/200 |

| check | verdict | detail |
|---|---|---|
| E3 shift-interval coverage of tau at level(90) | PASS | tau = 0.0: 177/200 = 0.8850 [0.8334,0.9221] vs nominal 0.90, count band [171,189], MARGINAL corridor [0.85,0.94] |

E3's p columns are NOT scored: with a planted tau they are power, not size, and mc_ri2 s5 documents pre-treatment p drift under a planted effect. Their means are in the secondary table below.

## 4. Secondary (information only, never fails the component)

| phase | mean p_main | mean p_pre | mean p_post0 | mean p_post1 | mean sd_main | mean b_main | n |
|---|---|---|---|---|---|---|---|
| E1 | 0.4709 | 0.4673 | 0.4984 | 0.4561 | 0.1901 | -0.0045 | 200 |
| E2 | 0.5239 | 0.5162 | 0.5142 | 0.5438 | 0.4262 | -0.0397 | 200 |
| E3 | 0.4685 | 0.5253 | 0.5052 | 0.4708 | 0.4262 | -0.0397 | 200 |
| E4 | 0.5230 | 0.5076 | 0.4987 | 0.5235 | 0.3678 | -0.0249 | 200 |
| E5 | 0.4683 | 0.5249 | 0.4932 | 0.4659 | 1.4684 | -0.0397 | 200 |

Under a correct null the mean equal-tailed p is not .5 but the lattice mean (0.5050 at B = 199); the one-sided mean is 0.5025. Read the means against those, not against .5.


## Verdict summary

- checks FAILED: 2
  - E1 p_post1 (ES_att_post1) size: alpha 0.05 (attained 0.0500): 13/200 = 0.0650 [0.0384,0.1080] band [3,17]; alpha 0.10 (attained 0.1000): 36/200 = 0.1800 [0.1329,0.2391] band [11,29]; chi2=16.8
  - E1 p_amepost0 (ES_ame_post0) size: alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 31/200 = 0.1550 [0.1114,0.2116] band [11,29]; chi2=13.1
- claims NOT JUDGED (alpha unattainable at B=199): 0
- checks not judged for want of reps (MINN=30): 0
