# FINDINGS — Monte Carlo size study of permute() and bootstrap() p-values under a sharp null
Agent: mc-boot-perm-size.  Method: MONTE CARLO size/level study.  Build: pinned c1d8979.
All numbers are measured.  `report_tables.txt` is the machine-generated source of truth.

## 1. Design and the reason the arms are not symmetric
Sharp null: y_it = a_i + f_t + e_it with e_it AR(1) rho=.5 within cluster
(Bertrand-Duflo-Mullainathan), and NO treatment effect at any horizon.  Onsets
2004/2005/2006 assigned AT RANDOM to 20% of the clusters each, 40% never treated,
one row per (state, year), T = 8.  i = state = cluster() — the commonest applied
DiD shape, and the regime the project's own record calls understated.
Because the onset assignment is drawn uniformly over the cohort-size-preserving
assignments and y does not depend on it, RANDOMIZATION INFERENCE IS EXACT here for
ANY statistic, with no asymptotics and no assumption about variance, serial
correlation or cluster count.  So every permute() arm carries a THEOREM and the
bootstrap arms carry only asymptotics.  That asymmetry was pre-registered
(NOTES.md), together with the achievable nominal level of every rank rule and the
decision rule: an exact arm is a defect only if the Wilson 95% LOWER bound on its
rejection rate clears its nominal level.
permreps(79) and reps(79) were chosen so every rank rule's achievable level lands
on 5%: at reps(99) "p <= .05" is a 4.04% rule for the percentile family and a
4.00% rule for boottstat.

## 2. Positive controls — the harness is in the right geometry
| control | project's published number | measured here |
|---|---|---|
| analytic joint pre-trend test oversized at G=10 | .657 (NF-1b) | **.755 [.619,.854]**, 49 draws |
| analytic att SE / reference at cluster()==i | .66-.75 (sthlp:255) | **.769** of the true MC sd, 50 draws |
| e(perm_npossible), 10 states, cohorts (2,2,2,4 never) | 18900 | **18900** |

## 3. THE CORE RESULT — G=10, 50 draws, identical data, bitwise-identical estimates
Rejection of a TRUE sharp null at nominal 5%, and the reported dispersion as a
multiple of the estimator's true Monte-Carlo sd (0.5291):

| engine / family | rej @5% | Wilson 95% | se / true sd | mean 95% CI width |
|---|---|---|---|---|
| permute(shuffle) att | **.061** | [.021,.165] | (no se) | — |
| permute(shuffle) att_i | .082 | [.032,.192] | (no se) | — |
| bootstrap(standard) pctl / bc / wald | .122 | [.057,.242] | 0.924 | 2.031 |
| bootstrap(bayes) percentile | .163 | [.085,.290] | 0.781 | 1.686 |
| bootstrap(bayes) bias-corrected | .204 | [.115,.336] | 0.781 | 1.686 |
| bootstrap(bayes) wald | .122 | [.057,.242] | 0.781 | — |
| ANALYTIC cluster se att | .180 | [.098,.308] | 0.769 | 1.840 |
| ANALYTIC cluster se att_i | .200 | [.112,.330] | — | — |
| boottstat (percentile-t) | **.025** | [.004,.129] | — | — |
| permtstat | **.075** | [.026,.199] | — | — |
A correctly sized 95% interval on this design needs width 2.074.
The POINT ESTIMATE is bitwise identical across analytic / standard / bayes and
differs from permute's by one ulp, so every difference above is dispersion.

At G=40 (20 draws, same DGP shape) EVERY engine rejects 1/20 = .05 and the
dispersion ratios are analytic 0.897, bayes 1.010, standard 1.067.  The G=10
failure is a FEW-CLUSTER phenomenon.

## 4. Findings
See the structured return.  In one line each:
* **bootstrap(bayes) is not a remedy at few clusters.**  Its standard error is
  statistically indistinguishable from the ANALYTIC one it replaces (0.781 vs
  0.769 of the true sd), its percentile/bc families reject .163/.204, and
  sthlp:255 names bootstrap() as the remedy for exactly that understatement.
  Mechanism (d03, reps=399, same data and seed): the Dirichlet replication sd is
  0.812x the pairs bootstrap's where the Dirichlet-vs-multinomial variance ratio
  G/(G+1) predicts 0.953; centring is fine for both.
* **bootstrap(standard) only halves the distortion** (.122, se 0.924x truth).
* **The few-cluster caution never fires here**: ado:6491-6499 keys on the count of
  TREATED clusters being < 6.  My design has 6 treated of 10, so the user sees no
  warning at all.
* **boottstat is the bootstrap family that survives**: .025 against .122 for the
  non-studentized families of the SAME engine on the SAME data.
* **The analytic understatement is ESTIMATOR-specific, not a blanket property of
  cluster() >= i** (d05, 25 draws): att alone se/sd 0.713 (rej .240); with
  esfixedbaseperiod alone 0.932 (rej .080); with the event-study window alone
  0.857 (rej .080); window + esfixedbaseperiod is IDENTICAL to esfixedbaseperiod
  alone.  The lever is esfixedbaseperiod.
* **The 2026-09 rewire of the joint and adversarial tests works** (40 draws):
  analytic joint .755 -> permutation .050, bootstrap .025; analytic adversarial
  .286 -> permutation .000.  Never size-tested before.
* **CONTRADICTION of a section-7 entry, with a caveat.**  The record says
  bootstrap(standard) at cluster()==i OVER-covers att_i (149/150 = .9933).  On a
  neighbouring design I measure .878 coverage (.122 rejection).  Their P5 had ONE
  state per cohort; mine has two.  The over-coverage is geometry-specific.

## 5. Probed hard, came back CLEAN (with the detection bound)
* **permute() attains its nominal level.**  POOLED over five arms in three
  designs where exactness is a theorem — permute default, permtstat,
  permute+esfixedbaseperiod, permonsetcohorts (12 states, cluster finer than i),
  and permexact (8 states, all 168 assignments enumerated) —
  **17 / 325 = .052 [.033, .082]** against a .05 nominal.  This is the arm where a
  high rate would have been an unambiguous bug, and it is not high.
* **permonsetcohorts (T1, shipped 2026-09-16, never size-tested): .057 [.029,.109]
  on 140 draws**, with e(perm_irebuilt)=1 and 79/79 replications on every draw.
  Its ordinary-route control on the same data gives the same number.
* **permexact: .018 [.003,.094]** on 56 draws with a COMPLETE 168-assignment
  enumeration every time (e(perm_exact)=1, e(perm_completed)=168, 0 failures).
* The Wald family's p equals 2*t_tail(G-1, |b/se_boot|) to 1.1e-16 under both
  bootstrap modes — bootdf(clusters)-by-default is exactly what the code does.
* esprelength()/espostlength() do not change the `att' ESTIMAND: on an exact DGP
  all four estimands return the design truth under every window (d02).
* The permutation p and the printed shift interval are exact duals at level(95)
  with permreps(79) — verified algebraically from the quantile rule and on live draws.
* e(perm_p_<est>) equals the printed table's p BITWISE, under the default
  convention and under permtstat.
* Replication completion under a clean sharp null: permute 79/79 on every one of
  ~400 draws; bootstrap(standard) 77-79/79; bootstrap(bayes) 79/79.
* e(perm_lb_)/e(perm_ub_) carry the documented, deprecated NULL BAND, not an
  interval (sthlp:1353).  My harness had that bug, not wooldid.

## 6. Two results added at the end of the run

### boottstat is correctly sized but its intervals are 2.7x too wide (40 draws)
sthlp:550 says in as many words: "whether the resulting intervals attain nominal
coverage has not been measured".  Measured here at G=10, cluster()==i:
  coverage .975 (rejection .025 [.004,.129]) -- VALID
  95% interval width: median 5.69, mean 7.98, p90 14.4, max 45.1
  against the 2.074 a correctly sized interval needs, and against 2.03 for the
  non-studentized percentile family on the same draws.
So boottstat buys its validity with intervals that are typically ~2.7 times wider
than necessary and occasionally 20 times wider.  That is the honest trade at few
clusters and the help does not state it.

### variance heterogeneity: the two versions answer different questions
  HVF (sd 3 on a FIXED cluster attribute, independent of the onset draw, so the
       sharp null HOLDS and exactness is still a theorem), 30 draws:
        permute .000/.033   permtstat .033/.067
        bootstrap(standard) .100/.067   bootstrap(bayes) .100
        ANALYTIC .167/.133
  HV  (sd 3 on the EVER-TREATED states - the F049 shape - where the sharp null is
       violated by construction), 22-23 draws:
        ANALYTIC .364 [.197,.571]   <- collapses
        permute .043   permtstat .043   bootstrap(standard) .000   bayes .046
The project's F049 note (ado:6053) records "7 of 30 null datasets rejected at 5%
with treated-cluster error sd 3 vs 1, against 6 of 30 for the raw coefficient",
i.e. about 20-23% for the permutation arms.  I measure 1/23 = 4.3% for both, and
the engine that actually fails on that design is the ANALYTIC one at 36%.  My 23
draws give [.008,.210], which still contains .20, so this is UNDERPOWERED to
contradict F049 - but the point estimates disagree sharply and the analytic arm,
which F049 does not mention, is the one that breaks.  Worth re-running at 200.

## LATE (2026-09-17 04:41): the esfixedbaseperiod finding's last premise, closed

d02 originally tested only `none / esprelength+espostlength / wider window`.  It
never tested `esfixedbaseperiod`, which d05 had identified as the actual lever.
Reran d02 on the noiseless DGP with `efbp` and `w33efbp` added: att = 1 and
1.5555556 (the design truths) to ~1e-15 under ALL FIVE specs.  So the option
changes the ESTIMATOR, not the ESTIMAND -- the SE difference is not a
"different target" artifact.  Log: d02_estimand_truth.log.
