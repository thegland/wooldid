# hv-bootsample REFUTATION — running notes

Build asserted: `which wooldid` = .../scratchpad/bluesky/build/wooldid.ado (r1_guards.log:61).
All data generated in-probe. No repo data, no repo generator, no test suite.

## Code facts established by reading (pinned build)
- ado:4941 / ado:394 pass `ifr(`if')` / `ifr(`wdp_if')` where `wdp_if` == the RAW user
  `if` (c_local at ado:615). The markout (`ifvar`) at ado:686 is NEVER forwarded.
- ado:8908-8923 `smp` = ((ifr)==1) & !missing(clusv). No markout. wdb_C = map-frame
  cluster count; `c_local wdb_r_nclus `wdb_C'` at ado:9537; boot_df = nclus-1 at 6600/6738/9570.
- ado:6459-6478 `smp2` likewise; the <2 refusal, <10 caution and <6-treated caution all
  key off it.
- ado:5957-5959 permctl `insmp` likewise (comment says it exists to stop FALSE r459s on
  bare-flag ifs — i.e. the author fixed the if-rule half of this and not the markout half).
- ado:4950-4956 `ifassert` is empty with no if(): the comment says so explicitly
  ("leaving the asserts running over the full dataset exactly as before").
- ado:15786-15790: BCa acceleration a = d3/(6 d2^1.5) from the JACKKNIFE ONLY.
  e(boot_accel) is posted (ado:10174). => a is MC-free; full-vs-oracle is a clean test.
- Help facts that bear on severity:
  * sthlp:532 says the 2-cluster refusal "fires when the options are screened, BEFORE
    the estimation runs" — the parse-time timing IS disclosed.
  * sthlp bootdf: "G is the number of cluster() groups" — does NOT say estimation sample.
  * sthlp:548 already says bootalpha(size) "is not a sound default for panel data" and
    gives "badly understated" SEs, and exists "for exact comparison against bayesboot".
  * sthlp:557 DOES say the bootsave observed row is "1/N, N the number of observations in
    the estimation sample".
  * ado:11681 ("exact shares with no Monte-Carlo error") and ado:11728 ("Elements not
    estimable in every enumerated assignment (exact shares over a reduced reference set
    ...)") are in the SAME footer block, printed one after the other.

## Measured so far
- r1 CASE A (20 frame states, 4 treated states with no outcome):
  A_FULL  N=384 Nclust=16 bootdf=19 dfr=19   (header "N clusters = 16", footer "t(19)")
  A_ORAC  N=384 Nclust=16 bootdf=15 dfr=15
  "only 4 treated cluster(s)" caution fires on the ORACLE, not on the full data.
  Same under bootstrap(bayes): AB_FULL bootdf=19, AB_ORAC bootdf=15.
  => HV-BOOT-01 reproduces on my own data.
- r1 CASE B (singleton unit, NO missing data anywhere):
  B_WITH N=160 Nclust=10 bootdf=10 ; B_WOUT N=160 Nclust=10 bootdf=9
  => HV-BOOT-04 reproduces: one single-period unit mints a phantom bootstrap cluster
     with no missing data at all.

- r1 CASE C (12 frame states, 2 live: 1 treated + 1 never):
  C_FULL N=40 Nclust=2 bootdf=11, header "N clusters = 2", footer "t(11)", NO cautions,
         se 4.018e-16, t 1.4e+15, p 0.000.
  C_ORAC N=40 Nclust=2 bootdf=1,  footer "t(1)", BOTH cautions print, se blank.
  => HV-BOOT-05 reproduces, and is slightly WORSE than claimed: on the full data the
     table prints a spuriously significant p=0.000 off a 4e-16 SE with no caution, while
     the oracle blanks the row AND cautions.
- r1 SE control (NOT claimed by the prober, and it holds up):
  A_FULL att .58092551 se .05645176 vs A_ORAC .58092551 se .05882665 (ratio .960)
  AB_FULL se .04893261 vs AB_ORAC se .04597457 (ratio 1.064, opposite sign)
  => under bootstrap(standard) and bootstrap(bayes)+bootalpha(equal) the dead clusters
     do NOT systematically move the SE; the damage is the guards and e(boot_df)/e(df_r).
     Correctly scoped by the prober.

- r3 PERM A (20 frame states, 4 treated states with no outcome):
  A_FULL  N=320 npossible=317,444,400 permsave onset cols = 20 (onset_1..onset_20) p=.00664452
  A_ORAC  N=320 npossible=43,680       permsave onset cols = 16 (onset_5..onset_20) p=.00664452
  "only 4 of the resampling clusters are ever treated" fires on the ORACLE only.
  CAVEAT: both p's sit exactly on the add-one floor 2/301, so this run says NOTHING
  about whether the p itself moves.  r3b re-runs it at a weak effect.
- r3 PERM B (permexact, 7 clusters, 2 dead):
  B_FULL npossible=42 completed=40 exact=1 incomplete=main_att p=0.2000
  B_ORAC npossible=5  completed=5  exact=1 incomplete=(none)  p=0.4000
  DISCLOSURE: the SAME output prints, in order, (i) "Caution: 2 of 42 enumerated
  assignments failed to estimate, so the enumeration is incomplete. The p-values stay
  exact shares ... but over the 40 assignments that completed", (ii) "Caution: under
  permexact the following elements were not estimable in every enumerated assignment
  ... a reduced reference set whose size is in e(perm_n): main_att", and (iii) the
  footer line "Elements not estimable in every enumerated assignment (exact shares over
  a reduced reference set; their n is in e(perm_n)): main_att".
  => HV-PERM-02's load-bearing claim ("'exact shares with no Monte-Carlo error' is false
  as printed") is NOT sustained: the qualification is printed three times in the same
  output.  Also: conditioning the reference set on estimability, where estimability is a
  function of the assignment and the observed assignment is a member of the retained set,
  is an EXACT conditional randomization test, not a broken one.
- r3 PERM C1: permute() r(459) "clusters to be subsets of treatment-onset groups" on a
  single stray row; deleted -> rc 0, N=160.  (The plain-run controls in r3 were
  mis-specified -- dots(0) needs bootstrap()/permute() -- and are redone in r3b.)

- r4 VAL A (20 sentinel rows unit == -9, y missing):
  A1 no if        -> r(411) "unit must consist of positive, non-zero integers"
  A2 if junk == 1 -> rc 0, e(N)=240
  A3 safety(off)  -> rc 0, e(N)=240    <-- a remedy the prober did not mention
  A4 rows deleted -> rc 0, e(N)=240
  => HV-VAL-01 reproduces; correction: safety(off) DOES lift it, and the message names
     neither that nor if().

- r5 (bootalpha(size)): frame rows 336, live rows 160, e(N)=160 on both sides, b
  bit-identical (+1.603ab04740d61X-001).
  SZ_FULL se .0270485748 ; SZ_ORAC se .0400620885 ; ratio .6752
  predicted sqrt(160/336) = .6901.  => HV-BOOT-02 reproduces, ~32% understatement.

- r6 (BCa): footer says "delete-one-cluster jackknife (18 fits)" on the full data and
  "(13 fits)" on the oracle, with e(N)=260 and e(N_clust)=13 on BOTH.  bootdf 17 vs 12.

- r2b STEP 1 (DIRECT REGIME READ-OUT, the instrument the original probe lacked):
  ado:406-407 prints the A10 caution iff wdb_reid & wdb_nest==1.  Naming i in fe()
  therefore prints the regime bit.  r2b_nest.log:
    line 187  REGIME.oracle  (no stray row)                 -> caution PRINTS  => nesting
    lines 249-320 REGIME.full (stray row, y MISSING)        -> NO caution      => keep-i
    lines 323-394 REGIME.control (stray row, y PRESENT)     -> NO caution      => keep-i
  One row that the estimator drops flips the regime.  HV-BOOT-06's mechanism is proven
  directly, not inferred from an SE.
- r2 DRAW-STREAM: max|w_full - w_oracle| = 0 EXACTLY over 401 saved rows, so the two runs
  draw the identical bsample sequence and any SE gap is the regime, not Monte Carlo.

- r5 CONTROL that isolates HV-BOOT-02: bootalpha(equal) (the DEFAULT) gives
  se = +1.efeac239772daX-004 on BOTH the full and the oracle data -- BIT-IDENTICAL.
  So dead rows move nothing except under bootalpha(size).  In the same fixture
  standard = .1428631773, equal = .1210734927, size(full) = .0270485748,
  size(oracle) = .0400620885: bootalpha(size) is already ~5.3x below standard, which is
  the documented sqrt(mean cluster size) shrinkage (frame n_c = 24, sqrt(24) = 4.9).
- r5 HV-BOOT-03: bootsave _rep==0 w column = .0029761905 = 1/336 exactly (frame rows),
  e(N) = 160; oracle file holds .00625 = 1/160.  sthlp:557 says 1/N with N "the number of
  observations in the estimation sample".  Divergence confirmed.

- r8 HV-BOOT-07 (MC-FREE): e(boot_jk_completed) = 18 with e(N_clust) = 13; footer
  "delete-one-cluster jackknife (18 fits)".  e(boot_accel) is bit-IDENTICAL across two
  different rseed()s (+1.d4f20b77b80d3X-005), proving it is a pure function of the
  jackknife.  Oracle accel = +1.d3b9024548fa4X-005.  Decimal .05724432 vs .05709506:
  a 0.26% difference.  So the mechanism is real and the numerical damage is negligible;
  the prober's p 0.007-vs-0.021 gap was stream confounding, as they flagged.

- r7 MY OWN NEW TRIGGERS (none suggested by anyone):
  A. ALWAYS-TREATED cohorts (onset 1996, window opens 2000 -- dropped by DOCUMENTED
     design, sthlp:73), NO missing data anywhere:
     A_FULL e(N)=288 Nclust=12 bootdf=15 ; A_ORAC e(N)=288 Nclust=12 bootdf=11.
     permute: AP_FULL npossible = 10,810,800 ; AP_ORAC = 5,940 (1820x).
  B. ZERO survey weights on two whole states (w==0, dropped at ado:686):
     B_FULL e(N)=288 Nclust=12 bootdf=13 ; B_ORAC bootdf=11.
  C. EVENT-STUDY WINDOW: C_NOW and C_WIN both e(N)=480, bootdf=11 -- the window did NOT
     drop rows.  Corroborates the prober's negative result on that third trigger.
  D. impvalues: D_FULL and D_ORAC both e(N)=288, Nclust=12, df_r=5 -- IDENTICAL.
     impvalues is CLEAN.  Negative result.
  E. permstrata(): E_FULL npossible=560 vs E_ORAC npossible=12, e(N)=288 both.

- r3b PERM p-sensitivity (the question r3 could not answer because both p's sat on the
  add-one floor): seed 11, weak effect: full p = .46589018 sd = .11977383 ;
  oracle p = .50582363 sd = .15011282.  The permutation NULL SD is 20% smaller on the
  full data.  Mechanism: the 20-cluster design gives a replication 8*(16/20) = 6.4 live
  treated clusters on average against the observed run's 4, so the null is tighter.

- r9 Q1 (the refusal-leak question): 12 frame clusters, exactly 1 live.
  Q1a bootstrap on the full file -> r(459) "cluster(state) leaves 1 cluster(s) in the
      estimation sample; cluster-robust inference needs at least 2. Use a finer cluster
      variable, robust, or fastnose."  (the ESTIMATION-path check, ado:1521)
  Q1b dead rows deleted          -> r(459) "bootstrap() requires at least 2 clusters in
      the estimation sample."  (the parse-time bootstrap refusal)
  Q1c ordinary analytic run      -> the same estimation-path error.
  => the bootstrap refusal leak has NO adverse consequence; the run is stopped either way,
     and the message the leak lets through is the MORE informative of the two.
- r9 Q2 (does the wrong e(boot_df) move a number?): 20 frame clusters, 15 live, reps 400.
  Q2a full        bootdf 19, lincom se .0427675382, CI [.4964515397, .6754785120]
  Q2c full + bootdf(14) (SAME draws, SAME se) CI [.4942377793, .6776922724]
  half-widths .08951349 vs .09172725, ratio .9759; t(19).975/t(14).975 = 2.0930/2.1448
  = .9758.  So on a realistic 20-cluster panel the hand-run lincom/test interval - the
  one the run's own footer tells the user reproduces its Wald numbers, and the one
  esttab/coefplot/honestdid read - is 2.4% too narrow.  Case C's version of the same
  thing is 5.8x.
- r3b C (HV-PERM-03, now with WORKING plain-run controls):
  C0  ordinary run on the stray file            rc 0, e(N)=160
  C0o ordinary run, stray row deleted           rc 0, e(N)=160   (identical sample)
  C1  permute() on the stray file               r(459) cluster-onset nesting
  C1s same + safety(off)                        r(459)  (NOT lifted)
  C2ord ordinary run on the .a-onset file       rc 0, e(N)=160
  C2/C2s permute() on it (from r3)              r(459) / r(459)
  C3clean permonsetcohorts, clean fixture       rc 0, e(N)=216
  C3ord ordinary run on the stray-onset file    rc 0, e(N)=216
  C3  permonsetcohorts on it                    r(459)
  C3s same + safety(off)                        r(459)  (sthlp:520 says so)
  C3o stray row deleted                         rc 0, e(N)=216   (identical sample)
- r2b STEP 2 (att_i SEs, dataseed 101, both rseeds, draw streams bit-identical):
  rseed 5551212  FULL .0978610043  ORAC .1156219131  CTRL .0977106359   ratio .8464
  rseed 8675309  FULL .0954973057  ORAC .1138693879  CTRL .0953086558   ratio .8387
  DRAW maxabsdiff = 0.0000000000 over 401 rows in both.
  The CONTROL (stray row with a real outcome, where keep-i is the CORRECT regime) lands
  on the full-data SE to three decimals and nowhere near the oracle - so the gap is the
  regime, not the extra row.
