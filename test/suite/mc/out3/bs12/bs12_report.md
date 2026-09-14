# BS-12 adversarial size ladder (mc_adv_bs12)

outdir `out3/bs12`; 8 shard file(s); 1200 rows after de-duplication on (config, arm, slope, shard, rep); 0 malformed line(s) dropped; 0 duplicate row(s) dropped.

Every adversarial size number on record was measured on the FREE-INTERCEPT statistic that BS-12 replaced on 2026-09-10 (MC_REPORT_2026-08-24, REPORT_ADVBOUND_2026-08-26, the help). This campaign is the re-measurement; until it is judged, those numbers describe a statistic that no longer ships.

Headline: the rejection rate of `p` at nominal .05 and .10 on the NULL arm (slope == 0), per config and per arm, against the central ~95% binomial region. The `pm` arm is exact under the sharp null and is the paired reference.

Cells with fewer than MINN=30 reps are reported, never judged.

| check | verdict | detail |
|---|---|---|
| rc == 0 on every row | PASS | rc!=0 on 0/1200; threshold 0.00 (driver ACCEPTANCE: rc = 0 on every row) |
| replication count | NOTE | B = 199 (modal ncomp); attained alpha on the right-tailed add-one lattice 0.05->0.0500, 0.10->0.1000 |

## Run inventory

| config | arm | what | slope | rows (rc=0) | rc!=0 | N (modal) | ncomp (modal) | nfail>0 | advabsent | mean secs |
|---|---|---|---|---|---|---|---|---|---|---|
| AB1 | an | analytic p, e(adversarial_pval_att) | 0 | 200 | 0 | 2160 | - | 0 | 0 | 3.17 |
| AB12 | an | analytic p, e(adversarial_pval_att) | 0 | 200 | 0 | 2160 | - | 0 | 0 | 3.05 |
| AB12 | bt | boottstat p (bootstrapped W, containment cell) | 0 | 200 | 0 | 2160 | 199 | 0 | 0 | 534.48 |
| AB13 | an | analytic p, e(adversarial_pval_att) | 0 | 200 | 0 | 2160 | - | 0 | 0 | 3.09 |
| AB1C | an | analytic p, e(adversarial_pval_att) | 0 | 200 | 0 | 2160 | - | 0 | 0 | 12.26 |
| PM | pm | permutation p (sharp null, cohort level) - the exact reference | 0 | 200 | 0 | 2160 | 199 | 0 | 0 | 433.68 |

## 0. Integrity (per-row facts; judged at any n - these are not statistical)

| check | verdict | detail |
|---|---|---|
| AB1/an advabsent == 0 | PASS | advabsent>0 on 0/200 (the adversarial p is missing at rc 0, or e(boot_adv_absent) is set: the statistic never formed on that rep, so it contributes no test - a rate computed without it is conditional) |
| AB1/an p in [0,1] | PASS | no p outside [0,1] |
| AB1/an pjoint in [0,1] | PASS | no e(joint_horiz) outside [0,1] |
| AB1/an row count | NOTE | 200 rc-0 null rows; no --M target given, so completeness of the campaign is not checked |
| AB12/an advabsent == 0 | PASS | advabsent>0 on 0/200 (the adversarial p is missing at rc 0, or e(boot_adv_absent) is set: the statistic never formed on that rep, so it contributes no test - a rate computed without it is conditional) |
| AB12/an p in [0,1] | PASS | no p outside [0,1] |
| AB12/an pjoint in [0,1] | PASS | no e(joint_horiz) outside [0,1] |
| AB12/an row count | NOTE | 200 rc-0 null rows; no --M target given, so completeness of the campaign is not checked |
| AB12/bt advabsent == 0 | PASS | advabsent>0 on 0/200 (the adversarial p is missing at rc 0, or e(boot_adv_absent) is set: the statistic never formed on that rep, so it contributes no test - a rate computed without it is conditional) |
| AB12/bt p in [0,1] | PASS | no p outside [0,1] |
| AB12/bt pjoint in [0,1] | PASS | no e(joint_horiz) outside [0,1] |
| AB12/bt ncomp == B (199) | PASS | ncomp != 199 on 0/200 (lost replications put that rep's p on a coarser lattice than the one judged here) |
| AB12/bt nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw) |
| AB12/bt lattice observed | FINDING | 200 p value(s) are consistent with: add-one k/(B+1) at B=199. This analyzer bands against the add-one form (right); if the only consistent label is 'plain k/B', the attained alpha is floor(alpha*B)/B and the count bands in section 1 are off by one rank - change ARM_CONV and re-run. Reported, never scored |
| AB12/bt advw_big | FINDING | reps with a replication W above 1e10: 0/200; max W seen 9.894e+05. On a containment rung this is REPORTED, not scored (driver ACCEPTANCE) - it is KF-2's escape route, where a single degenerate draw carries the whole distribution |
| AB12/bt row count | NOTE | 200 rc-0 null rows; no --M target given, so completeness of the campaign is not checked |
| AB13/an advabsent == 0 | PASS | advabsent>0 on 0/200 (the adversarial p is missing at rc 0, or e(boot_adv_absent) is set: the statistic never formed on that rep, so it contributes no test - a rate computed without it is conditional) |
| AB13/an p in [0,1] | PASS | no p outside [0,1] |
| AB13/an pjoint in [0,1] | PASS | no e(joint_horiz) outside [0,1] |
| AB13/an row count | NOTE | 200 rc-0 null rows; no --M target given, so completeness of the campaign is not checked |
| AB1C/an advabsent == 0 | PASS | advabsent>0 on 0/200 (the adversarial p is missing at rc 0, or e(boot_adv_absent) is set: the statistic never formed on that rep, so it contributes no test - a rate computed without it is conditional) |
| AB1C/an p in [0,1] | PASS | no p outside [0,1] |
| AB1C/an pjoint in [0,1] | PASS | no e(joint_horiz) outside [0,1] |
| AB1C/an row count | NOTE | 200 rc-0 null rows; no --M target given, so completeness of the campaign is not checked |
| PM/pm advabsent == 0 | PASS | advabsent>0 on 0/200 (the adversarial p is missing at rc 0, or e(boot_adv_absent) is set: the statistic never formed on that rep, so it contributes no test - a rate computed without it is conditional) |
| PM/pm p in [0,1] | PASS | no p outside [0,1] |
| PM/pm pjoint in [0,1] | PASS | no e(joint_horiz) outside [0,1] |
| PM/pm ncomp == B (199) | PASS | ncomp != 199 on 0/200 (lost replications put that rep's p on a coarser lattice than the one judged here) |
| PM/pm nfail == 0 | PASS | nfail>0 on 0/200 (a failed replication means the estimate never formed on that draw) |
| PM/pm lattice observed | FINDING | 200 p value(s) are consistent with: add-one k/(B+1) at B=199. This analyzer bands against the add-one form (right); if the only consistent label is 'plain k/B', the attained alpha is floor(alpha*B)/B and the count bands in section 1 are off by one rank - change ARM_CONV and re-run. Reported, never scored |
| PM/pm advw_big == 0 | PASS | reps with a replication W above 1e10: 0/200; max W seen 22.17 (driver ACCEPTANCE: advw_big = 0 on PM) |
| PM/pm row count | NOTE | 200 rc-0 null rows; no --M target given, so completeness of the campaign is not checked |

## 1. SIZE on the null arm (slope == 0) - the headline

(`attained` is P(p <= alpha) the arm's lattice can reach: alpha itself for the analytic p, floor(alpha(B+1))/(B+1) for the add-one permutation and bootstrap p. chi2 is the 10-bin uniformity statistic, df=9, crit 21.666 at p=.01 - secondary.)

| config | arm | conv | alpha | attained | rej (n, Wilson 95%) | count band | chi2 |
|---|---|---|---|---|---|---|---|
| AB1 | an | continuous | 0.05 | 0.0500 | 9/200 = 0.0450 [0.0239,0.0833] | [3,17] | 14.2 |
| AB1 | an | continuous | 0.10 | 0.1000 | 22/200 = 0.1100 [0.0738,0.1609] | [11,29] | 14.2 |
| AB12 | an | continuous | 0.05 | 0.0500 | 145/200 = 0.7250 [0.6593,0.7822] | [3,17] | 1039.1* |
| AB12 | an | continuous | 0.10 | 0.1000 | 156/200 = 0.7800 [0.7176,0.8318] | [11,29] | 1039.1* |
| AB12 | bt | right | 0.05 | 0.0500 | 0/200 = 0.0000 [0.0000,0.0188] | [3,17] | 63.2* |
| AB12 | bt | right | 0.10 | 0.1000 | 0/200 = 0.0000 [0.0000,0.0188] | [11,29] | 63.2* |
| AB13 | an | continuous | 0.05 | 0.0500 | 36/200 = 0.1800 [0.1329,0.2391] | [3,17] | 129.9* |
| AB13 | an | continuous | 0.10 | 0.1000 | 64/200 = 0.3200 [0.2593,0.3875] | [11,29] | 129.9* |
| AB1C | an | continuous | 0.05 | 0.0500 | 10/200 = 0.0500 [0.0274,0.0896] | [3,17] | 18.2 |
| AB1C | an | continuous | 0.10 | 0.1000 | 23/200 = 0.1150 [0.0779,0.1666] | [11,29] | 18.2 |
| PM | pm | right | 0.05 | 0.0500 | 7/200 = 0.0350 [0.0171,0.0705] | [3,17] | 4.0 |
| PM | pm | right | 0.10 | 0.1000 | 18/200 = 0.0900 [0.0577,0.1378] | [11,29] | 4.0 |

| check | verdict | detail |
|---|---|---|
| AB1/an size | PASS | alpha 0.05 (attained 0.0500): 9/200 = 0.0450 [0.0239,0.0833] band [3,17]; alpha 0.10 (attained 0.1000): 22/200 = 0.1100 [0.0738,0.1609] band [11,29]; chi2=14.2 |
| AB12/an size (report-only) | INFO | alpha 0.05 (attained 0.0500): 145/200 = 0.7250 [0.6593,0.7822] band [3,17]; alpha 0.10 (attained 0.1000): 156/200 = 0.7800 [0.7176,0.8318] band [11,29]; chi2=1039.1*. NOT SCORED: cluster(i) containment, ruled uncalibrated 2026-08-26; driver header: 'expected to FAIL' |
| AB12/bt size | FAIL | alpha 0.05 (attained 0.0500): 0/200 = 0.0000 [0.0000,0.0188] band [3,17]; alpha 0.10 (attained 0.1000): 0/200 = 0.0000 [0.0000,0.0188] band [11,29]; chi2=63.2* |
| AB12/bt uniformity (sec) | FLAG | 10-bin chi2 = 63.2 > 21.666 on n=200: p is not uniform even where the tail rate is in band - investigate, never silent |
| AB13/an size | FAIL | alpha 0.05 (attained 0.0500): 36/200 = 0.1800 [0.1329,0.2391] band [3,17]; alpha 0.10 (attained 0.1000): 64/200 = 0.3200 [0.2593,0.3875] band [11,29]; chi2=129.9* |
| AB13/an uniformity (sec) | FLAG | 10-bin chi2 = 129.9 > 21.666 on n=200: p is not uniform even where the tail rate is in band - investigate, never silent |
| AB1C/an size | PASS | alpha 0.05 (attained 0.0500): 10/200 = 0.0500 [0.0274,0.0896] band [3,17]; alpha 0.10 (attained 0.1000): 23/200 = 0.1150 [0.0779,0.1666] band [11,29]; chi2=18.2 |
| PM/pm size | PASS | alpha 0.05 (attained 0.0500): 7/200 = 0.0350 [0.0171,0.0705] band [3,17]; alpha 0.10 (attained 0.1000): 18/200 = 0.0900 [0.0577,0.1378] band [11,29]; chi2=4.0 |

## 2. Paired against the exact reference (same datasets, same reps)

(PM/pm is the permutation p under the sharp null: exact by construction. Every other cell is run on the SAME dataset per rep, so the comparison is paired and free of dataset noise. Discordance is where the two arms disagree at alpha = .05.)

| config | arm | n paired | arm rej@.05 | pm rej@.05 | arm only | pm only |
|---|---|---|---|---|---|---|
| AB1 | an | 200 | 9/200 = 0.0450 [0.0239,0.0833] | 7/200 = 0.0350 [0.0171,0.0705] | 7 | 5 |
| AB12 | an | 200 | 145/200 = 0.7250 [0.6593,0.7822] | 7/200 = 0.0350 [0.0171,0.0705] | 138 | 0 |
| AB12 | bt | 200 | 0/200 = 0.0000 [0.0000,0.0188] | 7/200 = 0.0350 [0.0171,0.0705] | 0 | 7 |
| AB13 | an | 200 | 36/200 = 0.1800 [0.1329,0.2391] | 7/200 = 0.0350 [0.0171,0.0705] | 31 | 2 |
| AB1C | an | 200 | 10/200 = 0.0500 [0.0274,0.0896] | 7/200 = 0.0350 [0.0171,0.0705] | 8 | 5 |

The pairing is the point of the ladder: an analytic arm that rejects where the exact permutation p does not is over-rejecting on THESE datasets, not on an incomparable draw.


## 3. slope != 0 - reported separately, NOT scored

No slope != 0 rows in this outdir.


## 4. Secondary (information only, never fails the component)

| config | arm | e(joint_horiz) rej@.05 | mean p | advw_n (modal) | max advw_max | tstud_n (modal) |
|---|---|---|---|---|---|---|
| AB1 | an | 5/200 = 0.0250 [0.0107,0.0572] | 0.4652 | - | nan | - |
| AB12 | an | 79/200 = 0.3950 [0.3298,0.4641] | 0.0806 | - | nan | - |
| AB12 | bt | 1/200 = 0.0050 [0.0009,0.0278] | 0.6212 | 199 | 9.894e+05 | 199 |
| AB13 | an | 16/200 = 0.0800 [0.0498,0.1260] | 0.3185 | - | nan | - |
| AB1C | an | 5/200 = 0.0250 [0.0107,0.0572] | 0.4612 | - | nan | - |
| PM | pm | 7/200 = 0.0350 [0.0171,0.0705] | 0.4828 | 199 | 22.17 | - |

`advw_n` is the count of usable replication W statistics behind that rep's p and `advw_max` their maximum; `tstud_n` is e(boot_tstud)'s n for main_att. A mean p far from the lattice mean, or an advw_max in the 1e10s, says the tail is one draw.


## Verdict summary

- checks FAILED: 2
  - AB12/bt size: alpha 0.05 (attained 0.0500): 0/200 = 0.0000 [0.0000,0.0188] band [3,17]; alpha 0.10 (attained 0.1000): 0/200 = 0.0000 [0.0000,0.0188] band [11,29]; chi2=63.2*
  - AB13/an size: alpha 0.05 (attained 0.0500): 36/200 = 0.1800 [0.1329,0.2391] band [3,17]; alpha 0.10 (attained 0.1000): 64/200 = 0.3200 [0.2593,0.3875] band [11,29]; chi2=129.9*
- claims NOT JUDGED (alpha unattainable at B=199): 0
- checks not judged for want of reps (MINN=30): 0
