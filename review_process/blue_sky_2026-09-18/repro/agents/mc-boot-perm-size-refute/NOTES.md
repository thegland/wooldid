# SKEPTIC NOTES — refuting the "MC size of permute()/bootstrap() p-values at few clusters" findings

WORKDIR: .../probes/mc-boot-perm-size-refute
Build asserted: `which wooldid` -> pinned build (see s00_smoke.log line ~5).

## Source-level results already established (no MC needed)

### MCBP-GUARD -> PARTIAL (mechanism refuted, severity ~none)
- ado ~6473: `if r(max) < 10 dis "Caution: only N clusters are available for the bootstrap;
  bootstrap inference is unreliable with very few clusters."`  <-- a TOTAL-cluster caution
  exists, immediately above the treated-cluster one.
- ado 6491-6499: the `< 6` TREATED-cluster caution, explicitly author-ruled F050 2026-09-06,
  whose own comment says the total-cluster caution "above cannot see ... few TREATED clusters".
- EMPIRICAL (s00_smoke.log): G=9/6 treated -> "Caution: only 9 clusters are available..." fires.
  G=10/6 treated -> nothing. G=12 -> nothing.
- So the finding's stated mechanism ("keys on the TREATED cluster count") and its EXPECTED
  ("should key on the total cluster count, or at least both counts") are factually wrong:
  the code keys on BOTH. What survives is only that the total threshold is `< 10` and the
  prober chose exactly G=10 — a boundary/threshold-choice quibble.

### MCBP-JOINTSIZE -> REFUTED (rediscovery + misquoted band + suppressed runtime caution)
- ado 1392-1399 (SI-A1, author-ruled 2026-08-25): the analytic joint test under clustering at
  or above i is documented with "**measured rejection rates 11-66% at nominal 5%**, persisting
  at large cohort counts and large cells". The finding quotes the project band as "11-36%".
- known_bydesign.md §1 states "65.7% at G=10 (NF-1b)". The finding's own Wilson interval
  [.619,.854] CONTAINS .657.
- ado 1435-1443 (PS-64, author ruling 2026-09-08): on exactly this geometry wooldid PRINTS
  "Caution: jointtests was combined with analytic clustered SEs clustered at the level of the
  cohort variable i or above ... the joint test's covariance is degenerate and its p-value
  should not be read as a calibrated test. wooldid's randomization inference (permute()) or
  bootstrap-based clustered SEs give a joint test that is valid under this design..."
  The prober ran `qui` and never saw it; their EXPECTED ("a user ... is reading a statistic
  that fires three times in four") is the caution's own sentence.

### MCBP-BAYES novelty -> refuted at source
- ado 6491-6499 (the very lines the finding cites) already say "the Bayesian bootstrap shares
  its [the pairs cluster bootstrap's] first two moments" and cite MacKinnon-Webb badly sized.
  The claimed novel content ("the Bayesian/Dirichlet variant does not escape it either") is
  a restatement of the comment.
- bootstrap(bayes) + att_i is REFUSED r(198) under default safety (verified, s00_smoke.log):
  "Error: bootstrap(bayes) with att_i reweights clusters but cannot resample treatment cells,
  so those bootstrap standard errors omit the cell-sampling variance component and can be too
  small ...". The prober's own do-file notes this at line 72.

### MCBP-S7COVER -> design mismatch is larger than claimed
- `git show c8c4e2c` (the §7 record): M=150, reps(199), arms P4 cluster(unit) FINER than
  i=state and P5 cluster(state)==i. So the §7 design has MULTIPLE UNITS PER STATE, and its DGP
  has LIVE COHORT HETEROGENEITY in treatment effects ("the DGP's cohort heterogeneity is live
  on both arms"). It is not a sharp null.
- The finding says its geometry "differs only in having two states per onset cohort instead of
  one". It differs on at least four axes: sub-units per state, true effects present vs a sharp
  null, reps 199 vs 79, M 150 vs 49.

## Timing (s00_smoke.log timers)
 - analytic wooldid on G=10,T=10: ~0.30 s / call
 - bootstrap(standard) reps(199): ~37 s / draw
 - boottstat reps(199): ~41 s / draw

## Runs
- s00_smoke.do / .log          DONE  (build assert, caution census, timings)
- r01_analytic.do / .log       RUNNING (300 draws x 5 analytic specs + 120 draws multi-unit)

## r01_analytic.log — 300 draws, G=10, cluster(st)==i, SHARP NULL, my own generator

| spec | sd(b) | mean se | se/sd | rej@5% | df_r | N | Kposted |
|---|---|---|---|---|---|---|---|
| S0 att alone                      | 0.4742 | 0.3385 | **0.714** | **0.170** | 9 | 100 | 1 |
| S1 + esprelength(3) espostlength(3)| 0.5895 | 0.5144 | 0.873 | 0.123 | 9 | 100 | 8 |
| S2 + esfixedbaseperiod            | 0.6356 | 0.5144 | **0.809** | **0.140** | 9 | 100 | 2 |
| S3 window + esfixedbaseperiod     | 0.6356 | 0.5144 | 0.809 | 0.140 | 9 | 100 | 8 |
| S5 S3 + jointtests + adversarial  | 0.6356 | 0.5144 | 0.809 | 0.140 | 9 | 100 | 8 |

corr(bS0,bS2)=0.662, corr(bS1,bS2)=0.292.  df_r is 9 = G-1 on every arm, so the
reference distribution is t(G-1) and nothing here is a df mis-specification.

MULTI-UNIT (5 sub-units per state-year, i=state, cluster(st)==i, 120 draws, N=500):
| M0 att alone         | 0.4493 | 0.2840 | **0.632** | **0.258** |
| M2 + esfixedbaseperiod | 0.6019 | 0.3354 | **0.557** | **0.275** |

=> MCBP-ESFBP REFUTED on its actionable content:
   (1) at 300 draws on a G=10 sharp null, esfixedbaseperiod gives se/sd 0.809 and
       rejection .140 — NOT the 0.93-1.01 / .061-.080 claimed from 25 draws;
   (2) on a multi-unit panel it moves se/sd the WRONG WAY (0.632 -> 0.557) and
       raises rejection (.258 -> .275).  It is not a lever.
   (3) sd(esfbp) > sd(default) on both shapes, which is what sthlp:836 predicts
       ("pooled reference periods offer the possibility of greater precision");
       the finding's sd(S2) < sd(S0) was noise at n=25.
   (4) The default's se/sd 0.714 sits inside the documented 0.66-0.75 (sthlp:251),
       which also already says "at that level the analytic covariance is the wrong
       object rather than a mis-scaled one" — i.e. "the default SE is simply wrong"
       is the help's own sentence, not a new finding.
   (5) sthlp:832-836 describes esfixedbaseperiod as switching how ALL treatment
       effect estimates are formed, so "an event-study display option" is a misread.

JOINT (same 300 draws, S5 arm):  e(joint_horiz) rejection @5% = **0.6967** (n=300)
  ECDF .01/.05/.10/.20/.30/.50 = .567/.697/.750/.827/.887/.947
  e(adversarial_pval_att) rejection = 0.2367  (inside the documented 11-36%)
=> matches the project's own "65.7% at G=10 (NF-1b)" and the ado's own
   "measured rejection rates 11-66% at nominal 5%".  The finding's ".755, 2-7x the
   documented 11-36% band" misquotes the band and its own CI contains .657.

## r02_cautions.log — what wooldid actually SAYS
- jointtests at cluster(st)==i prints the PS-64 caution VERBATIM (line 121):
  "...the joint test's covariance is degenerate and its p-value should not be read
   as a calibrated test. wooldid's randomization inference (permute()) or
   bootstrap-based clustered SEs give a joint test that is valid under this design..."
- permute() at cluster==i prints its own caution naming MacKinnon-Webb permtstat and
  the one-observation-per-i-period limitation.
- Few-cluster caution census with bootstrap(standard), 6 treated states throughout:
  G=8  -> "Caution: only 8 clusters are available for the bootstrap..."   FIRES
  G=9  -> "Caution: only 9 clusters are available for the bootstrap..."   FIRES
  G=10 -> nothing;  G=11 -> nothing;  G=12 -> nothing
  G=10 with 4 treated -> "Caution: only 4 treated cluster(s) are available..." FIRES
  So BOTH a total-cluster and a treated-cluster caution exist; the finding's claimed
  "band of designs (6-9 treated clusters, any number of controls)" that gets nothing
  does not exist — only total>=10 AND treated>=6 is silent.

## r05_esfbp_altdgp.log — esfixedbaseperiod as an "SE lever" on two more shapes, 300 draws each
| shape | spec | sd(b) | mean se | se/sd | rej@5% |
|---|---|---|---|---|---|
| C8  (G=10,T=8, onsets 2004/05/06 = their onset shape) | default       | 0.5368 | 0.3977 | 0.741 | .153 |
| C8  | + esfixedbaseperiod | 0.5969 | 0.5291 | 0.886 | .140 |
| C12 (G=12,T=12, 4 cohorts)                            | default       | 0.4669 | 0.3362 | 0.720 | .173 |
| C12 | + esfixedbaseperiod | 0.5221 | 0.5233 | **1.002** | **.090** |

FOUR SHAPES now, esfixedbaseperiod's effect on se/sd:
  G10 T10 1-row  : 0.714 -> 0.809   rej .170 -> .140
  G10 T10 5-sub  : 0.632 -> 0.557   rej .258 -> .275     (WRONG DIRECTION)
  G10 T8  (theirs): 0.741 -> 0.886  rej .153 -> .140
  G12 T12        : 0.720 -> 1.002   rej .173 -> .090
It is shape-dependent, from -0.08 to +0.28. NOT a lever.
AND NOTE: even at se/sd = 1.002 the rejection is .090, i.e. 1.8x nominal — so the
finding's inference "se/sd ~ 1 therefore the SE is correct therefore 6-8% rejection"
does not follow; at G=10-12 the SE is estimated and the level stays oversized.
