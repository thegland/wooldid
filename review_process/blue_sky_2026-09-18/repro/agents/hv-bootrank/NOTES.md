# hv-bootrank — running notes

## Assignment
Verify the seeded claim: `reps(8)` posts a rank-7 e(V) over 20 coefficients;
a hand `test()` returns p=1.5e-14 with no warning. Stated mechanism:
`ereturn repost V =` at ado:10040-10066 bypasses `ereturn post`'s
positive-definiteness screen.

## Code reading done BEFORE running (facts, with line numbers)
- ado:10044-10053 — the repost site. The comment ALREADY states the cost:
  "e(V) can leave this program INDEFINITE, so any downstream consumer that
  inverts it or takes the square root of its diagonal (test, lincom,
  honestdid, the pretrendpower path) must screen for a non-positive-finite
  diagonal first". So the repost choice is deliberate and documented in code.
- ado:14498-14576 `wooldid_bootfin` — the covariance is the **LISTWISE
  complete-case** covariance over non-excluded columns. `nlw` = number of
  replications complete across EVERY non-excluded column; posted as
  `e(boot_N_V)`. Per-column `se[k]` uses a DIFFERENT row set (column k alone
  available). A11 ruled that divergence documented.
- ado:14563-14567 — any column with `excl==1`, `cnt<2`, **or `V[k,k]==0`**
  gets a `-1e100` diagonal and a zeroed row/col. The code comment itself
  says "a column can show a standard error in the table and carry the
  sentinel in e(V)".
- ado:6235-6260 — reps() has NO floor by author ruling; ado:5890 permute()
  has a hard floor of 20. Confirmed, and NOT in question.
- ado:6475-6477 — the reps() cautions key on the p-value GRID only
  (`wdbgrid > .05`). Nothing keys on rank, on q-vs-reps, or on `nlw`.
- ado:11786 — the footer line is
  `dis as result "  e(V) holds the bootstrap covariance (complete-replication
  listwise, e(boot_N_V) replications; ...)"` — `e(boot_N_V)` is INSIDE the
  string literal, so the footer prints the NAME, not the value. Neighbouring
  lines (11783/11784) correctly break out of the string to print
  `e(boot_completed)`. So the realized listwise count is never displayed.
- sthlp:532 — the help DOES say "Neither is the Wald test on e(V) that a hand
  `test` reproduces; that one is a normal-theory test on the bootstrap
  covariance and is not printed." So a hand `test` being a Wald on the
  bootstrap V is documented. What is NOT documented anywhere I can find is
  that its reference distribution is wrong when the draw count is small.

## Analytic calibration (python, scipy) — pB_size.txt
A chi2(q) Wald test on a B-draw bootstrap covariance is Hotelling-shaped:
W*(B-q)/(q*(B-1)) ~ F(q, B-q). Actual size of a nominal 5% test:
  B=8   q=3  -> 0.254      B=50  q=6  -> 0.105
  B=8   q=6  -> 0.735      B=50  q=10 -> 0.177
  B=13  q=6  -> 0.394      B=100 q=6  -> 0.074
  B=17  q=6  -> 0.283      B=200 q=10 -> 0.073
  B=30  q=6  -> 0.156      B=1000 q=10-> 0.054
**B here is e(boot_N_V), not reps().** That is the crux: reps() is what the
user sets and what every caution names; boot_N_V is what the covariance
actually rests on, and it can be far smaller.

## Status
- pA_core.do  : RUNNING (core repro + rank-vs-reps + ereturn post screen test)
- gen_panel.do: written (own generator, no repo data)
- gen_admin.do: written (realistic 30-state staggered adoption panel)

---
## MEASUREMENTS (as they land)

### pA_core.log / pA_core_partial.log — seeded claim, 200-unit panel, 40 clusters
    K=20  PROPS="b V"  VCVSTATE=bootstrap  BOOT_N_V=8  BOOT_FAILED=0
    RANKV=7 of KV=20   NSENT=0 (no -1e100 anywhere)
    eigs: 7 positive (3.18e-2 .. 1.36e-4), 13 at |1e-18|, 6 negative, min -2.65e-18
    test ES_att_post0 post1 post2 -> F(3,39)=90.97  p=1.177e-17  no constraint dropped
  NOTE the reference is an F, not a chi2: wooldid posts e(df_r)=39=G-1.
  NOTE the recon's demonstration is on POST coefficients, which have large TRUE
  effects; the tiny p there is a CORRECT rejection and proves nothing about size.
  FOOTER, verbatim from the log (line 373):
    "e(V) holds the bootstrap covariance (complete-replication listwise,
     e(boot_N_V) replications; ...)"
  -> prints the TOKEN. Two lines above it, "8 of 8 replications completed" and
     "resolution 2/8 = 0.250" interpolate their values correctly.

### pA2_lean.log — rank vs reps, 40-unit panel, 20 clusters, same K=20 spec
    reps   boot_N_V   rank(e(V))   K
      4        4          3        20
      8        8          7        20
     15       15         10        20
     30       30         10        20
     60       60         10        20
  rank = min(rank of the DISTINCT column space, boot_N_V - 1).
  The ceiling of 10 is not a defect: on a balanced panel att == att_it and
  att_i == att_itime EXACTLY, so 20 requested columns span only 10 dimensions.
  The point stands: nothing about `ereturn repost' enters. Rank deficiency is
  the arithmetic of a sample covariance over boot_N_V draws.

### size_table.txt — analytic size of a hand `test' (Hotelling)
  W(B-q)/(q(B-1)) ~ F(q, B-q), but wooldid refers it to F(q, e(df_r)).
  Nominal 5%, df_r=29:  B=8 q=6 -> .70 | B=13 -> .33 | B=20 -> .17 |
  B=30 -> .10 | B=60 -> .054 | B=200 -> .031 | B=1000 -> .025
  So the hand test is grossly OVERSIZED below B~55 and mildly CONSERVATIVE above.
  Severity therefore hinges entirely on how far boot_N_V falls below reps().

### pD_post.log / pD2_sentinel.log — the `ereturn post` screen, Stata 19.0 MP
  RANK_A=2 of 3   ereturn post -> POSTA_RC=0        (rank-deficient PSD ACCEPTED)
  RANK_B=2 of 5   ereturn post -> POSTC_RC=0        (bootstrap-shaped, ACCEPTED)
      test c1 c2 c3 c4 c5 -> "Constraint 1/3/4 dropped", chi2(2), p=.9796
      -> Stata DOES signal when q exceeds the rank. The seeded claim's
         "no constraint dropped" is therefore NOT a defect at q=3 <= rank 7.
  SENTINEL: V[1,1] = -1e100
      ereturn post   -> rc 0, e(V)[1,1] comes back VERBATIM as -1.00e+100
      ereturn repost -> rc 0, e(V)[1,1] comes back VERBATIM as -1.00e+100
      test x1 on that column -> "Constraint 1 dropped", chi2(0)=., p=.
      lincom x1 -> coefficient printed, se/z/p/CI all blank
  RANK-1 3x3 via ereturn post -> rc 0, V[1,1]=.0025 preserved (NOT zeroed);
      test x1 x2 -> "Constraint 2 dropped", chi2(1), p=.8554
  => THREE load-bearing code comments rest on a premise that does not
     reproduce on Stata 19.0 MP:
       ado:10046-10049 "the -1e100 sentinel rows ... would otherwise trip it"
       ado:14556-14560 "It survives ... because repost does not re-run the
                        positive-definiteness screen that ereturn post applies"
       ado:16956-16961 "ereturn post can zero a nonsymmetric or highly
                        singular VCE wholesale"
     All three describe `ereturn post' as a barrier it is not (for the
     symmetric cases wooldid actually produces).

### pB1_listwise.log — SPEC 1, 30-state staggered adoption, long ES window, NO subgroups
  wooldid y state t ttre, att att_it esfixedbaseperiod esprelength(6)
      espostlength(11) cluster(state) jointtests pretrendpower(50 80)
      bootstrap(standard) reps(60) rseed(4242)
    S1_COMPLETED=60   S1_REPS=60   S1_N_V=51      <-- 85% of reps, SILENT
    S1_K=38           S1_DFR=29    S1_VPOSTED=1
    S1_RANK=19 of 38  S1_NSENT=0   S1_MINEIG=-6.09e-17
    S1_RANK_VESATT=17 of 17        (the ES family block is FULL rank)
    S1_NAVAIL_MIN=51 MAX=60 MEAN=59.26
    S1_NCOL_BELOW_HALF=0           <-- no column trips the 50% caution
    S1_BD_LISTWISE_ALL=51          (my own recomputation from e(bootdist): matches)
  READ: a 15% listwise haircut with NOTHING said. rank 19 of 38 is the
  att/att_it duplication ceiling, NOT a draw-count deficiency (51-1=50 >> 19).
  At B=51, q=6 the hand test's size is ~.06 — essentially fine. So the POOLED
  long-window spec is a mild case. The subgroup spec (pB2) is the real test.
  NOTE pB1 aborted r(3000) after this block on a Mata syntax slip of mine
  (`cs[.,2]'` inside `mata:`); the S1 diagnostics above are complete and the
  file is patched. Lost: S1TEST, the pretrendpower dump and SPEC 3.

### pB2_subgroups.log — SPEC 2, THE COLLAPSE. Same panel + regional subgroups
  wooldid y state t ttre, att att_it esfixedbaseperiod esprelength(6)
      espostlength(11) cluster(state) subgroups(reg)
      bootstrap(standard) reps(60) rseed(4242)
    S2_COMPLETED=60  S2_FAILED=0  S2_REPS=60   S2_N_V=14   <-- 23% of reps
    S2_K=184         S2_DFR=29    S2_VPOSTED=1
    S2_RANK=13 of 184            S2_NSENT=0   S2_MINEIG=-2.15e-16
    S2_RANK_VESATT=17 of 17      (ES family block still FULL rank)
    S2_NAVAIL_MIN=39 MAX=60 MEAN=55.52
    S2_NCOL_BELOW_HALF=0         <-- min availability 39/60 = 65%, so NOT ONE
                                     column trips the 50% caution
    S2_BD_LISTWISE_ALL=14        (independent recomputation: matches)
  CAUTIONS PRINTED BY THE WHOLE RUN: **ZERO** (awk NR>240 | grep -c Caution = 0)
  FOOTER printed, verbatim:
    "e(V) holds the bootstrap covariance (complete-replication listwise,
     e(boot_N_V) replications; ...)"
    "e(V) and e(df_r) carry the standard-deviation family referred to the t(29)
     distribution, so test and lincom run by hand reproduce the
     bootdisplay(wald) numbers ... and so do the tools that read e(V) without
     being asked: esttab, coefplot and honestdid export the standard-deviation
     family, not the columns above."
  => the run NAMES esttab/coefplot/honestdid as consumers of e(V) and never
     says e(V) rests on 14 of the 60 replications. A hand test of q=6 pre
     coefficients is referred to F(6,29) on a 13-df covariance: actual size
     ~.30 against a nominal .05.

### pE_sediverge.log — THE CONSEQUENCE. se_printed vs sqrt(diag e(V))
  Same panel, `att' only + subgroups(reg), long window, reps(60), rseed(4242).
  boot_completed=60  boot_N_V=14  K=92  rank(e(V))=13  navail 39..60
  columns below the 50% caution threshold: 0   CAUTIONS PRINTED: 0
  (se_ratio_table.txt has the full 19-row table)
    main:att        b=1.007355  printed se .079677   sqrt(diag eV) .040155  x0.504
    ES_att:post9    b=1.580801  printed se .127551   sqrt(diag eV) .057044  x0.447
    ES_att:post4    b=0.971471  printed se .069095   sqrt(diag eV) .041152  x0.596 (navail 60!)
    17 of 19 ratios BELOW 1;  min .447  median .634  max 1.032
    median VARIANCE ratio .403
  => e(V) - which the footer itself says esttab, coefplot, honestdid, test and
     lincom all read - carries standard errors about 0.63x the ones wooldid
     PRINTED, i.e. t-statistics ~1.6x too large and 95% intervals ~37% too
     narrow, silently.
  The shrinkage is NOT confined to the thin columns: post4/post5/post6 have
  navail = 60 and still shrink to .60/.56/.54. That points at SELECTION - the
  14 listwise-complete rows are the resamples that drew every thin cluster,
  which are systematically less dispersed - rather than at sampling noise.
  pF_selection.do tests exactly that.

### pF_selection.log — SELECTION PROVEN (fsel_table.txt has all 84 rows)
  Within one run, from e(bootdist): for each column compare the sd over its own
  available rows (what wooldid PRINTS) with the sd over the 14 listwise-complete
  rows (what e(V) HOLDS), and against 2000 RANDOM 14-row subsets of the same
  column.  If the listwise rows were a random subset, p_lower would be uniform.
    F_BDROWS=60  F_NLW=14  columns examined = 84
    p_lower <= .05 in 19 of 84 columns (expected ~4.2 by chance)
    only 2 columns significantly HIGH
    65 of 84 ratios below 1; overall median ratio 0.889
    POOLED (non-subgroup) columns: 17/19 below 1, median ratio 0.634,
      12 of 19 with p_lower <= .05
      main_att        sd_full .079677  sd_lw .040155  ratio .504  p_lower .0010
      ES_att_post5    .110175 / .061694  ratio .560  p_lower .0060
      ES_att_post8    .115269 / .056349  ratio .489  p_lower .0005
      ES_att_post9    .127551 / .057044  ratio .447  p_lower .0010
      ES_att_post10   .131248 / .063692  ratio .485  p_lower .0005
  => the listwise-complete replications are NOT missing-at-random. They are the
     resamples that happened to draw every thin cluster, which are the LESS
     DISPERSED ones, so e(V) is biased DOWNWARD. Not a reconciliation nuisance:
     a factor of two on the headline att coefficient.

### pD3_zeroing.log — BOTH halves of the PS-62 premise, Stata 19.0 MP
    (A) rank-1 of 3, SYMMETRIC : ereturn post rc 0, V preserved, ALLZERO=0
    (B) NONSYMMETRIC (symdiff .01): ereturn post rc 0, e(V) ALL ZERO (ALLZERO=1)
    (C) NONSYMMETRIC via repost  : rc 0, e(V) ALL ZERO TOO  (ALLZERO=1)
    (D) -1e100 sentinel, symmetric: post AND repost both return it verbatim
  => "nonsymmetric OR highly singular" : the nonsymmetric half is real but
     repost shares it; the highly-singular half is false for both. The repost
     choice buys nothing against the only real hazard. Symmetry discipline
     ((S+S')/2, makesymmetric) is what protects e(V).

### navail_pattern.txt — the missingness is GEOMETRIC, not random (mechanism CONFIRMED)
  pE run, 92 columns, 60 completed replications, 0 failed.
  43 columns have navail = 60; 49 lose replications, and they are exactly the
  horizons carried by one cluster:
    ES_att_SG2_post7..post10   navail 39   (region 2's deepest post = state 2 alone)
    ES_att_SG1_post8..post10   navail 41   (region 1 = state 1 alone)
    ES_att_SG3_post4..post9    navail 41   (region 3 = state 3 alone)
    ES_att_SG4_post4..post8    navail 44   (region 4 = state 4 alone)
    ES_att_post10 51, post9 56, post8 59   (the pooled deep horizons)
  reg = mod(state-1,4)+1 puts the four single-state early adopters (onsets
  2008/2009/2010/2011) in four DIFFERENT regions, so those four dropout events
  are independent. P(a given cluster absent from a resample) = (1-1/30)^30 =
  .362, so a single-state horizon should be present in about .638*60 = 38.3
  replications: MEASURED 39, 41, 41, 44. And the listwise intersection should be
  about .638^4 * 60 = 10; MEASURED 14. The arithmetic matches the geometry.
  -> the selection story is not a guess: the listwise-complete rows are the
  resamples that drew every thin cluster.

### pG_export.log — DEMONSTRATED on the same estimate, same seed
    wooldid printed  main:att  b 1.007355  se 0.079677  t 12.643
                               95% CI [0.825545, 1.164670]
    esttab exported  main_att  1.007355  (0.040155)          <- 50.4% of it
    lincom main_att            se .04015483  CI [.92523, 1.08948]
    esttab ES_att_post9   1.580801 (0.057044)  vs printed 0.127551
    esttab ES_att_post10  1.270671 (0.063692)  vs printed 0.131248
    test of 6 pre coefficients -> F(6,29)=1.777 p=.1389 on a 13-df covariance
  boot_completed=60, boot_N_V=14, df_r=29. No warning from any of them.
