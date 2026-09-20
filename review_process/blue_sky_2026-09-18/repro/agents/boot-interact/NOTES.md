# boot-interact — NOTES (live)

Target: bootstrap() x everything. Method: interaction matrix + covariance provenance tracing.
WORKDIR: .../probes/boot-interact

## Code-read conclusions (to be CONFIRMED by run, not investigated further)
- **Seeded headline lead (ptp reads the wrong covariance) is GUARDED in code.**
  - ado:425 bootwire swaps e(V_es_*) -> bootstrap, preserves analytic as hidden e(V_es_analytic_*).
  - ado:449 shell calls wooldid_ptp AFTER bootwire; ptp (ado:10427-10437) gates on
    `cap confirm matrix e(V_es_analytic_`ty')` and SKIPS with a caution if absent.
  - ado:10225-10246 bootwire reposts e(df_r) = bootstrap df (or . under bootdf(normal)),
    keeping e(df_r_analytic). So ptp's `dfr()` at 449 is the BOOTSTRAP df. Provenance OK.
- Bootstrap joint tests are a genuine statistic-bootstrap: W_r (ado:15238 wooldid_bootjwr) uses
  the REPLICATION's own analytic e(V_es_*); W_obs (ado:15311 wooldid_bootjwobs) uses the
  OBSERVED analytic twin. Same metric on both sides + scale-free eigen guard. Consistent.
- ptp has a rank guard: `diag0cnt(invsym(S))` on the k x k PRE block only.
- bootwire uses `ereturn repost V=` ON PURPOSE to bypass the PD screen (ado:10037-10048);
  its own comment says downstream consumers (test, lincom, honestdid, ptp) "must screen for a
  non-positive-finite diagonal first". wooldid's own consumers do. USER's `test` does not.

## PARKED (no run)
- `V_es_analytic_att_it_s<sg>` is 22 chars + subgroup value; >10-char subgroup values would
  overflow Stata's 32-char e() name limit in an UNGUARDED `ereturn hidden matrix` at ado:10086.
  Needs 11+ digit subgroup codes = the author-excluded giant-integer class. Not run.

## Run log
(started)

## CONFIRMED (p01b_rank.log, p03_conseq.log)
- **Rank boundary measured.** 40-state staggered panel, cluster(state), esfixedbaseperiod,
  esprelength(4) espostlength(6) -> colsof(e(b)) = 12.
  reps(8): e(V) RANK 7 of 12, min eig -6.5e-17, e(boot_N_V)=8, e(boot_completed)=8.
  reps(19): RANK 12 of 12, min eig 8.4e-5.  reps(21)/reps(40): full rank.
  => rank(e(V)) = min(completed - 1, k). Boundary is completed >= k+1, NOT reps>=20.
  NOTHING in e() or in any caution reports the rank. `test` over 3 coefficients runs at df 3,
  drop=0, no message.
- Longer ES (esprelength(6) espostlength(6), k=14, 12 ES coefs):
  reps(8): e(V) rank 7/14, e(V_es_att) rank 7/12.  reps(12): 11/14 and 11/12.
  simultaneousbands critical value: analytic 2.6504 (simulated), reps(8) 2.0712 (n=8),
  reps(12) 2.4471 (n=12). The band at reps(8) is 22% NARROWER than the analytic one and
  e(esband_n_att) is the only record; no caution about the draw count.
  honestdid runs rc 0 on the rank-7-of-12 bootstrap e(V_es_att).
- Provenance baseline CONFIRMS the seeded headline lead is GUARDED: after bootstrap(),
  e(wdid_vcvstate)="bootstrap", e(V_analytic) and e(V_es_analytic_att) both exist,
  mreldif(e(V_es_att), e(V_es_analytic_att)) is large (0.017-0.15) i.e. genuinely swapped,
  e(boot_ves_failed) empty, e(df_r)=e(boot_df)=39 with e(df_r_analytic)=39 kept.

## FINDING CANDIDATE A (p07_testrank.log) — hand `test` on the low-reps bootstrap e(V)
40-state panel, 20 periods, cluster(state), esfixedbaseperiod, esprelength(6) espostlength(6),
jointtests. k = colsof(e(b)) = 14, 6 pre-period ES coefficients.
Hand `test` over the SIX pre-period coefficients (the natural hand pre-trend test):
  analytic          p = .01003909   df 6  drop 0     [= e(joint_horiz), they agree]
  bootstrap reps(4) p = .85675851   df 3  drop 1     rank(e(V)) = 3 of 14; wooldid's own
                                                     joint test ABSENT (e(boot_jt_absent)=1)
  bootstrap reps(8) p = 2.703e-06   df 6  drop 0     rank(e(V)) = 7 of 14 ; e(joint_horiz)=.1111
  bootstrap reps(12)p = .00006988   df 6  drop 0     rank 11 of 14 ; e(joint_horiz)=.0769
  bootstrap reps(30)p = .27820513   df 6  drop 0     rank 14 of 14 ; e(joint_horiz)=.3226
Only caution printed at reps(8): the 2/reps p-value GRID caution. Nothing about rank or about
a Wald test on a few-draw covariance.
sthlp:1340 explicitly invites the hand test: "e(df_r) is aligned with the reference
distribution actually reported, so a test or lincom run by hand after a bootstrap fit
reproduces the bootdisplay(wald) p-values under every bootdf() setting".

## FINDING CANDIDATE B (p10_ptpboundary.log) — pretrendpower() under bootstrap(), the
## rank guard stops firing at completed = k + 1 and the numbers are arbitrary there
40 states x 20 periods, cluster(state), esfixedbaseperiod, esprelength(6) espostlength(4),
pretrendpower(80).  6 pre-period coefficients (k = 6).
  analytic : KEFF 6/6  slope -.06623  se .02758  CI (-.1220,-.0105)  joint pre p .01004
                                                  mdt_wald .1101  bias_wald .6203
  reps(7)  : KEFF 6/6  slope +.08134  se .00722  CI (+.0667,+.0960)  joint pre p 2.36e-12
                                                  mdt_wald .02885 bias_wald .1625
             ONLY caution printed: the 2/reps p-value GRID caution. No rank caution.
  reps(8)  : KEFF 6/6  slope +.10420  CI (+.0542,+.1542)  p 2.703e-06
  reps(10) : KEFF 6/6  slope +.11512  CI (+.0649,+.1654)  p 1.788e-06
  reps(15) : KEFF 6/6  slope -.04258  CI (-.1167,+.0315)  p .3654
  reps(25) : slope -.05956  p .2981      reps(50): slope -.04690  p .3299
  reps(3/4/6) (completed <= k): KEFF 2/3/5 -> the rank caution DOES fire.
=> The SIGN of the reported pre-trend slope flips, its CI excludes zero on the wrong side,
   the joint pre-trend p goes to 1e-12, and mdt_wald understates the detectable trend 3.8x.
   sthlp:532 claims "the existing rank guard applies to the bootstrap block as it does to
   the analytic one" - it does, and it is not the relevant guard.

## MECHANISM (unifies A and B) — the missing Hotelling correction
e(V_es_att) is a SAMPLE covariance from n = e(boot_completed) draws, but e(df_r) stays on
the analytic scale (bootdf(clusters) -> G-1 = 39).  wooldid_ptp (ado:10545) computes
Ftail(keff, dfr, W/keff).  The correct reference for W = b'S^{-1}b with S from n draws is
Hotelling: (n-k)/((n-1)k) * W ~ F(k, n-k).  Inflation factor (n-1)/(n-k): 6.0 at n=7,k=6;
3.5 at n=8; 1.56 at n=15; 1.13 at n=50.
PROOF THAT IT IS THE SAME STATISTIC: at reps(8) the hand `test` over the six pre
coefficients returns p = 2.703e-06 (p07_testrank.log) and e(roth_prep_att) = 2.703e-06
(p10_ptpboundary.log) - identical to every printed digit.
wooldid's OWN jointtests are immune (rank test on the draws): e(joint_horiz) = .1111 at
reps(8).  So one run can print two bootstrap pre-trend p-values that differ by 1e5.

## p09 (realistic zero-caution case)
esprelength(8) espostlength(8) subgroups(grp): k = 72, 64 ES coefficients.
reps(50): rank(e(V)) = 49 of 72, min eig -5.9e-16. ZERO caution lines printed in the run
block (2/50 = .04 < .05 so ado:6477 does not fire). Hand test over the 64 ES coefficients:
p = 9.948e-33 on df 49, drop = 1 (Stata silently dropped 15). Analytic: p = 4.681e-17 df 36.

## p05 (cluster finer than / cross-cutting i) — CLEAN
i = st (10 states), cluster = county (40, finer): analytic se .2220, bootstrap .2244, ratio
1.011, 200/200 reps complete. cluster = xcut (7 cross-cutting groups): analytic .2371,
bootstrap .2331, ratio 0.983, 200/200 complete, e(df_r)=e(boot_df)=6. cluster = st (== i):
analytic .0975, bootstrap .2212, ratio 2.27, 169/200 complete (31 replications drew no
never-treated cohort and were excluded). att_it takes the physical leg and returns the same
numbers as att here. boottstat completes on all three regimes.

## p04b (provenance) — CLEAN
Analytic e() has 56 entries, bootstrap 84; NOTHING is analytic-only (purely additive).
e(pwsgs) IS rewired (mreldif 2.67; subgroup-contrast se .055 -> .387). e(escomp) and
e(adversarial_line_att) are unchanged (mreldif 0) - they carry composition diagnostics and
line-fit point estimates, not inference. e(adversarial_pval_att) is correctly BLANKED under
bootstrap without boottstat (1.109e-07 -> missing) with e(boot_adv_absent)=1 and a caution.

## p14 (ptp x bootdf variants, reps(8), same seed) — the reference distribution is the defect
slope .10420368 and se .02469802 are IDENTICAL in every variant (same covariance).
Only roth_prep moves: bootdf(clusters)/default df_r=39 -> 2.703e-06 ; bootdf(normal) df_r=.
-> chi2(6) -> 4.023e-10 ; bootdf(3) -> F(6,3) -> .04797 ; boottstat -> 2.703e-06.
Hotelling at n=8,k=6 -> .3006.
## p09 C60: 11 pre coefficients, reps(60): KEFF 11/11 (bootstrap) vs 9/11 (analytic, which
DOES print the rank caution). wooldid p .4984 vs Hotelling .6431, inflation 1.204, and ZERO
caution lines in the run block.

## p13 CALIBRATION (decisive for severity split)
Stata's OWN prefix: `bootstrap, reps(8) cluster(state): regress y i.t` on 40x20 ->
rank(e(V)) = 7 of 21; `test` of six time dummies gives p = 2.883e-81, df 6, drop 0, no
warning. reps(15) .001774 ; reps(50) .026803 ; reps(400) .093996 ; analytic .162645.
=> BI-2 (hand `test` on a low-reps bootstrap V) is GENERIC STATA behaviour, not a wooldid
bug. Report it as documentation (sthlp:1340 invites the hand test) + the no-floor asymmetry.
BI-1 stands alone: Stata never PRINTS such a test; wooldid does, as pretrendpower()'s
headline joint pre-trend line.
## p14 V5 bootstrap(bayes) reps(8): PREP .005246, slope .0400 - same class.
## p14 V6 no cluster(): bootdf default normal -> chi2(6), PREP .06305.

## p06a INTERACTION MATRIX (reps(20), cluster(state), esfixedbaseperiod espre(3) espost(3))
metric: e(wdid_vcvstate), e(boot_method), and sqrt(e(V)[1,1]) vs sqrt(e(V_analytic)[1,1]).
CORRECT (vcvstate=bootstrap, analytic twin kept, SE moved):
  1 (bare) 1.229 | 2 att_i 1.398 | 3 att_it 1.229 | 4 att_itime 1.398 | 5 all four |
  6 subgroups(grp) | 7 jointtests | 8 adversarial | 9 simultaneousbands | 10 gteffects |
  20 controls(x1) 1.290 | 21 interactivecontrols(x1) 1.125 | 22 anticipation(1) 1.138 |
  23 esrelativeto(-2) 1.247 | 24 boottstat | 25 bootdisplay(bc) | 26 bootdisplay(bca)
REFUSED r(198), all matching the documented incompatibility list:
  11 impvalues | 12 fastnose | 13 unconditionalse | 14 customvce() | 15 permute() |
  16 subgroups()+suppresssgprimaryeffects | 17 oldsyntax
NOT TESTED: 18 poisson (my y has negative values -> r(411), my generator's fault)
ODD: 19 timetrends (with cluster(state) == i): runs rc 0 and the printed table carries
  BOOTSTRAP standard errors, t, p and CI - but `ereturn repost V =` fails rc 503 (the
  estimate was posted without a V because the analytic one was all zero), so e(V) is NOT the
  bootstrap covariance, e(wdid_vcvstate) reads "absent", e(V_analytic) does not exist, and
  every downstream reader (test, lincom, esttab, coefplot, honestdid) sees nothing. A caution
  names e(V_boot). sthlp:1059 describes that design without carving out bootstrap(): it says
  "every standard error, test statistic, p-value and interval in the run is missing", which is
  FALSE under bootstrap().

## p15 (BI-4, confirmed)
analytic: vcvstate absent, e(V) rc 111. bootstrap reps(60): vcvstate absent,
boot_v_posted 0, e(V) rc 111, e(V_boot) rc 0, e(V_analytic) rc 111. Printed table fully
populated with bootstrap se/t/p/CI. Hand `test`/`lincom` fail rc 321 - LOUD, not silent.
So BI-4 is a doc divergence (sthlp:1059 unqualified for bootstrap) + inconvenience, LOW.

## p12b COMPLETE (flagship repro log)
reps(200): wooldid roth_prep = hand test = .49410733 ; Hotelling .50134147 ; jointtests
.45273632 ; inflation 1.03. Convergence confirmed: the three agree once n >> k.

## FINAL STATE 2026-09-17 02:55
FINDINGS.md is complete and current: BI-1 (high), BI-2 (low-medium), BI-4 (low), BI-5 (low),
plus the provenance baseline (clean) and the full interaction matrix.
p09b_zerocaution was still running at close; its numbers are already in p09_calib.log
(A block rank 49/72, p 9.948e-33, 0 cautions; C60 .4984 vs .6431, 0 cautions).
