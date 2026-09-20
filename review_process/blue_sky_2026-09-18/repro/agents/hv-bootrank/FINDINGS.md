# hv-bootrank -- FINDINGS (live; updated as measurements land)

Build asserted on every run: `which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado
No repo data, no repo generator, no suite file touched.
Own generators: gen_panel.do, gen_admin.do.

---

## F1 -- VERDICT on the seeded claim: TRUE IN ITS FACTS, WRONG IN ITS MECHANISM

Seeded claim: reps(8) posts a rank-7 e(V) over 20 coefficients; test() returns
p=1.5e-14 with no warning. Mechanism given: `ereturn repost V =` at ado:10040-10066
bypasses `ereturn post`'s positive-definiteness screen.

Reproduced from scratch (pA_core.do / pA_core_partial.log, own 200-unit 12-period
3-cohort panel, 40 clusters nesting units):

    K=20            PROPS=b V        VCVSTATE=bootstrap
    BOOT_N_V=8      BOOT_FAILED=0    BOOT_VPOSTED=1
    RANKV=7         KV=20            NSENT=0
    eigenvalues: 7 positive (3.18e-2 ... 1.36e-4), 13 at |1e-18| or below,
                 6 of them negative (min -2.65e-18)
    test ES_att_post0 ES_att_post1 ES_att_post2
         -> F(3, 39) = 90.97,  p = 1.177e-17, no constraint dropped, no warning

K=20 OK, rank 7 OK, e(properties)="b V" OK, negative eigenvalues at roundoff OK,
hugely significant hand test with no warning OK. My p (1.177e-17 on F(3,39))
differs from the recon's (1.458e-14, df=3) only because the data differ -- note it
is an F, not a chi2: wooldid posts e(df_r) = G-1 = 39, so Stata's test uses
F(q, df_r).

THE MECHANISM IS WRONG, AND IT MATTERS FOR THE FIX.

1. Rank deficiency is not caused by the repost. It is arithmetic: the covariance is
   the sample covariance of the completed replication draws, so its rank is at most
   (draws - 1). Measured by varying only reps() on a fixed K=20 spec (pA2_lean.do):
   reps(4) -> e(boot_N_V)=4, rank 3.  rank = min(K, e(boot_N_V) - 1).
2. `ereturn post` would not have caught it. A rank-deficient PSD matrix is exactly
   what regress itself posts when a regressor is dropped for collinearity.
3. What the repost does let through is the -1e100 sentinel wooldid_bootfin writes
   (ado:14563-14567) -- an INDEFINITE matrix, a different object. The code comment
   at ado:10044-10053 already says so and already names the exposed consumers.
   On the reps(8) run NSENT=0, so no sentinel was involved.

ALSO MATERIAL TO SEVERITY: the recon tests ES_att_post0..post2, which have large
TRUE effects in any DGP with a real treatment effect. A tiny p there is a correct
rejection, not a wrong number. The test that carries evidence is on PRE-period
coefficients under a DGP with no pre-trends -- see F2.

REALISTIC? Not on its own. reps(8) is not a reported specification, and Stata's own
bootstrap prefix at reps(8) posts the same rank-7 matrix with the same silence.
sthlp:532 already states a hand test is "a normal-theory test on the bootstrap
covariance and is not printed". As stated: LOW severity.

---

## F2 -- (my own lead) e(V) rests on e(boot_N_V), not reps(), and the realized
##       count is never displayed

wooldid_bootfin (ado:14536-14552) builds e(V) as the LISTWISE complete-case
covariance over all non-excluded columns, while the reported standard errors come
from a PER-COLUMN pass over a different row set (A11 ruled the se/V divergence
documented). The listwise row count is e(boot_N_V).

From the code alone:
- Every reps() caution (ado:6475-6477) keys on the P-VALUE GRID (2/reps). None keys
  on e(boot_N_V), on rank, or on the number of constraints a user might test.
- The one line that mentions the count, ado:11786, is
  dis as result "  e(V) holds the bootstrap covariance (complete-replication
  listwise, e(boot_N_V) replications; ...)"
  -- e(boot_N_V) sits INSIDE the string literal, so the footer prints the token,
  not the value. The two neighbouring lines (11783/11784) correctly break out of
  the string to print e(boot_completed). So the count is stored and documented
  (sthlp:1340) but never shown.
- A chi2/F Wald on a B-draw covariance is Hotelling-shaped:
  W(B-q)/(q(B-1)) ~ F(q, B-q). wooldid posts e(df_r) (cluster count - 1), so test
  compares against F(q, df_r). Nominal-5% actual size (size_table.txt), df_r=29:

      B=boot_N_V |  q=6   |  q=14
             8   |  0.70  |  singular
            13   |  0.33  |  0.85 (q=10)
            20   |  0.17  |   --
            30   |  0.10  |  0.40
            60   |  0.054 |  0.116
           200   |  0.031 |  0.027

  The same B governs pretrendpower()'s reported slope CI, joint pre-trend p and
  detectable trend under bootstrap() (sthlp:532: Roth's formulas are fed the
  bootstrap covariance), and the honestdid handoff through e(V_es_*). Those are
  wooldid's own reported numbers, not hand tests.

Measurement of how far below reps() boot_N_V falls on realistic input: pB_listwise.do.

---

## F3 -- (my own lead) THREE load-bearing code comments assert an `ereturn post`
##       behaviour that does not reproduce on Stata 19.0 MP

The repost site (ado:10044-10049) justifies deviating from the standard idiom:

    * repost, not post, ON PURPOSE: `ereturn repost V=' does NOT re-run
    * `ereturn post''s positive-definiteness screen - MEASURED, a 2x2 with
    * V[1,1] = -1e100 reposts rc 0 and comes back verbatim - and the -1e100
    * sentinel rows wooldid_bootfin manufactures for excluded or
    * under-available columns WOULD OTHERWISE TRIP IT.

wooldid_bootfin repeats it (ado:14556-14560) and wooldid_vspl repeats a
stronger version (ado:16956-16961: "ereturn post can zero a nonsymmetric or
highly singular VCE wholesale").

MEASURED (pD2_sentinel.do / .log, Stata 19.0 MP, no wooldid involved):

    V[1,1] = -1e100, off-diagonals zeroed, symmetric:
      ereturn post   -> rc 0,  e(V)[1,1] = -1.00e+100  (VERBATIM)
      ereturn repost -> rc 0,  e(V)[1,1] = -1.00e+100  (VERBATIM)
      test x1 -> "Constraint 1 dropped", chi2(0)=., p=.
      lincom x1 -> coefficient shown, se/z/p/CI blank

    exact rank-1 3x3 PSD V (V = a a'):
      ereturn post -> rc 0, V[1,1]=.0025 preserved, NOT zeroed
      test x1 x2 -> "Constraint 2 dropped", chi2(1), p=.8554

    rank-2-of-5 bootstrap-shaped covariance (3 draws, 5 columns):
      ereturn post -> rc 0
      test c1..c5 -> "Constraint 1/3/4 dropped", chi2(2), p=.9796

So on this Stata, `ereturn post` accepts the sentinel, stores it verbatim, and
does not zero a highly singular VCE. The three comments describe a barrier that
is not there. No user-visible number is wrong because of this -- the route
chosen works -- but the stated justification for the nonstandard idiom is
false, and it is exactly the kind of "MEASURED" claim a maintainer would rely
on when refactoring. (Caveat: the author may have measured an older Stata; the
comments carry no version.)

Two corollaries that also correct the record:
- ado:14558-14560 says a hand test or lincom on a sentinel column "builds a
  quadratic form on a -1e100 variance with no error raised". No error is
  raised, but Stata does not fabricate a number either: `test` drops the
  constraint and reports missing, `lincom` prints blanks. The observable
  failure is clean.
- Stata's `test` DOES announce "Constraint N dropped" whenever the number of
  constraints exceeds the rank of the submatrix. So the seeded claim's "with no
  constraint dropped" is not itself a defect: at q=3 against rank 7, nothing
  should be dropped. The unsignalled case is q <= rank with a SMALL number of
  replications behind the covariance -- which is F2.

---

## F4 -- (my own lead, THE MAIN FINDING) the listwise complete-case rule makes
##       e(V) systematically NARROWER than the standard errors wooldid printed,
##       by a factor of two on the headline coefficient, in total silence

`wooldid_bootfin` (ado:14536-14552) computes two different things from the same
replication matrix:
  - the PRINTED standard error, `se[k] = sqrt(variance(A[use,k]))`, over the
    rows where column k alone is available; and
  - `e(V)`, the LISTWISE complete-case covariance over rows complete across
    every non-excluded column (that count is `e(boot_N_V)`).
The code comment says A11 "ruled this a documented consequence rather than a
defect" because "a hand-run test or lincom does not EXACTLY reproduce the
table's Wald numbers". Measured, it is not an exactness nuisance.

MEASURED (pE_sediverge.do / .log, se_ratio_table.txt; own 30-state panel):

    wooldid y state t ttre, att esfixedbaseperiod esprelength(6)
        espostlength(11) cluster(state) subgroups(reg)
        bootstrap(standard) reps(60) rseed(4242)

    boot_completed = 60   boot_N_V = 14   K = 92   rank(e(V)) = 13
    e(boot_navail): min 39, max 60   columns below the 50% threshold: 0
    CAUTIONS PRINTED BY THE WHOLE RUN: ZERO

    coefficient      b          printed se   sqrt(diag eV)   ratio   navail
    main:att       1.007355      .079677        .040155      0.504     60
    ES_att:post4   0.971471      .069095        .041152      0.596     60
    ES_att:post8   1.634908      .115269        .056349      0.489     59
    ES_att:post9   1.580801      .127551        .057044      0.447     56
    ES_att:post10  1.270671      .131248        .063692      0.485     51
    ... 17 of 19 pooled ratios below 1; min .447, median .634, max 1.032
        median VARIANCE ratio .403

IT IS SELECTION, NOT NOISE (pF_selection.do / .log, fsel_table.txt). Inside the
same run, from e(bootdist), each column's listwise-14 sd is compared with 2000
RANDOM 14-row subsets of that column's own available rows:

    columns examined 84;  p_lower <= .05 in 19 of them (expected ~4.2)
    only 2 columns significantly HIGH;  65 of 84 ratios below 1
    pooled columns: 17/19 below 1, median ratio .634, 12 of 19 at p_lower <= .05
      main_att      p_lower = .0010      ES_att_post8  p_lower = .0005
      ES_att_post9  p_lower = .0010      ES_att_post10 p_lower = .0005

The listwise-complete replications are not missing-at-random, and the reason is
specific: those 14 rows are the resamples in which states 1, 2, 3 and 4 were ALL
drawn at least once. In the full 60, each of those clusters is drawn zero times
in about 36% of resamples -- and a resample holding zero copies of a cluster is
exactly an extreme resample for every coefficient that cluster drives.
Conditioning on "at least one copy" TRUNCATES the multiplicity distribution of
those four clusters and deletes precisely the most extreme draws. So `e(V)` is
biased DOWNWARD for the coefficients those clusters drive.

That is why the shrinkage is heterogeneous rather than uniform: it is worst on
`main_att` and the deep POOLED post horizons, which states 1-4 dominate (ratios
.45-.60), milder on the subgroup columns (median .917), and one column,
ES_att_SG1_pre7, goes the other way at 1.36. A "the surviving resamples are
calmer in general" story would predict uniform shrinkage; the cluster-truncation
story predicts exactly the pattern measured.

THE MECHANISM IS CONFIRMED BY THE MISSINGNESS GEOMETRY, not assumed
(navail_pattern.txt). Of the 92 columns, 43 are available in all 60
replications; the 49 that lose replications are exactly the horizons a single
cluster carries:
    ES_att_SG2_post7..post10  navail 39   (region 2's deep post = state 2 alone)
    ES_att_SG1_post8..post10  navail 41   (region 1 = state 1 alone)
    ES_att_SG3_post4..post9   navail 41   (region 3 = state 3 alone)
    ES_att_SG4_post4..post8   navail 44   (region 4 = state 4 alone)
    ES_att_post10 51, post9 56, post8 59  (the pooled deep horizons)
With G=30, P(a cluster absent from a resample) = (1-1/30)^30 = .362, so a
single-state horizon should appear in about .638*60 = 38 replications --
measured 39, 41, 41, 44 -- and the intersection of four INDEPENDENT such events
(the four earliest adopters sit in four different regions, by
reg = mod(state-1,4)+1) in about .638^4*60 = 10 -- measured 14.

DEMONSTRATED, not argued (pG_export.do / .log, same spec, same seed):
    wooldid printed   main:att  b 1.007355  se 0.079677  t 12.643
                                95% CI [0.825545, 1.164670]
    esttab exported   main_att  1.007355   (0.040155)
    lincom main_att             se .04015483  95% CI [.92522917, 1.0894809]
    esttab ES_att_post9  1.580801 (0.057044)   vs printed se 0.127551
    esttab ES_att_post10 1.270671 (0.063692)   vs printed se 0.131248
    test of the six pre-period coefficients -> F(6,29) = 1.777, p = .1389,
      i.e. a q=6 Wald referred to F(6,29) on a covariance carrying 13 df
A regression table built from this estimate with esttab therefore reports a
standard error for the headline att that is 50.4% of the one the wooldid table
showed, and a lincom interval less than half as wide. Both commands run without
a word of warning.

(The `test` line above is ONE draw under a true null and is not a size
demonstration -- it shows only that the reference is F(6,29) on a 13-df object.
The size claim rests on size_table.txt and on pF's permutation p-values.)

THE DOCUMENTED FALLBACK CARRIES THE SAME BIAS. The footer offers
`e(V_boot)`, but ado:10058-10064 builds it from the SAME retained master copy:
`cap matrix vnew = bV` / `ereturn matrix V_boot = vnew`, where `bV` is the
matrix just reposted into `e(V)`. So `e(V_boot)` is `e(V)`, 14 rows and all. A
user who reads the footer and reaches for the alternative gets the same matrix.

THE RULING THIS CONTRADICTS. A11 (quoted in the code at ado:14500-14506) says
the divergence arises "whenever any column loses a replication", that
"sqrt(diag(e(V_boot))) is not the printed se column and a hand-run test or
lincom does not EXACTLY reproduce the table's Wald numbers", and rules it "a
documented consequence rather than a defect". Measured: it arises whenever
SEVERAL INDEPENDENT columns each lose about 35% of the replications, and the
factor is two, not a rounding tail. Per the coordinator note, contradicting a
ruled entry by measurement is the finding.

WHO THIS HURTS. The run footer names the consumers itself:
  "e(V) and e(df_r) carry the standard-deviation family referred to the t(29)
   distribution, so test and lincom run by hand reproduce the bootdisplay(wald)
   numbers rather than the columns printed above - and so do the tools that read
   e(V) without being asked: esttab, coefplot and honestdid export the
   standard-deviation family, not the columns above."
That sentence frames the difference as a FAMILY difference (percentile vs sd).
On this run it is also a factor-of-two SCALE difference. A regression table
built with esttab from this estimate carries a standard error for att of .040
where wooldid printed .080; t goes from 12.6 to 25.1; the 95% interval goes from
[.826, 1.165] to about [.925, 1.089]. coefplot draws the narrow bands. honestdid
starts from the narrow covariance.

Nothing announces it:
  - the reps() cautions (ado:6475-6477) key on the p-value GRID (2/reps) only;
    at reps(60) the grid is .033 and none of them fires;
  - the low-availability caution (ado:10296) keys on PER-COLUMN counts at a 50%
    threshold; here the thinnest column is available in 39 of 60 (65%), so it
    does not fire either -- yet the listwise intersection is 14;
  - the only line that mentions the realized count, ado:11786, is
      dis as result "  e(V) holds the bootstrap covariance (complete-replication
      listwise, e(boot_N_V) replications; ...)"
    with `e(boot_N_V)` INSIDE the string literal, so it prints the TOKEN, not
    the value -- verified in pB2_subgroups.log and pA_core_partial.log. The two
    neighbouring lines (ado:11783/11784) and the permute footer (ado:11651) use
    the correct `e(...)` interpolation, so this is a slip, not a house style.

SECOND CONSEQUENCE, same root. A hand `test` is a Wald on e(V) referred to
F(q, e(df_r)) with df_r = G-1 = 29, while the covariance carries 13 degrees of
freedom. Hotelling: W(B-q)/(q(B-1)) ~ F(q, B-q), so at B=14 and q=6 the actual
size of a nominal-5% test is about .30 (size_table.txt). wooldid's OWN joint
tests are immune -- the fixedmetric rank test (ado:6810-6845) uses the same
metric for the observed and replication statistics and carries its own rank
guard -- and e(V_es_att) came back FULL rank (17 of 17) on both collapsed runs,
so pretrendpower()/honestdid sit on a healthier, family-specific listwise set.
The asymmetry is the point: wooldid screens and discloses rank for its own
tests, and hands the user an e(V) with neither.

REALISTIC? Yes. A 30-state staggered-adoption panel, an event-study window long
enough to show the early adopters' late horizons, regional subgroups, and
reps(60). The trigger is structural: any specification with many columns of
which a few are carried by one or two clusters. Spec 1 of pB1 (the same panel
WITHOUT subgroups) gave boot_N_V = 51 of 60, a mild 15% haircut - so the effect
scales with the column count, and subgroups() is what makes it bite.

---

## Status: all probes run. Files in this directory:
  gen_panel.do / gen_admin.do  own generators (no repo data)
  pA_core.do / pA_core_partial.log   seeded-claim reproduction
  pA2_lean.do / .log                 rank vs reps (4/8/15/30/60)
  pD_post.do / .log, pD2_sentinel.do / .log, pD3_zeroing.do  ereturn post screen
  pB1_listwise.do / .log             pooled long-window spec  (boot_N_V 51/60)
  pB2_subgroups.do / .log            + subgroups            (boot_N_V 14/60)
  pE_sediverge.do / .log             printed se vs sqrt(diag e(V))
  pF_selection.do / .log             the selection permutation test
  size_table.txt, se_ratio_table.txt, fsel_table.txt

---

## F3b -- the PS-62 record measured only `repost`, and the half it asserted
##        about `post` is the half that is wrong

`review_process/held_for_ruling_2026-09-09.md:99` states:
  "Stata's `ereturn post` screen does zero a nonsymmetric or highly singular
   VCE wholesale, which is what `wooldid_amescreen , mode(vzero)` exists to
   detect. It does not arise on this path."
and then measures only `repost` ("probe at ~/.claude/jobs/be9644d0/.../probe_ps62.do,
2026-09-10 ... `repost` was measured to accept verbatim ...").

`post` itself was never measured. I measured it, both halves of the disjunction
(pD2_sentinel.do, pD3_zeroing.do; Stata 19.0 MP, no wooldid involved):

    (A) HIGHLY SINGULAR but symmetric, exact rank 1 of 3:
        ereturn post -> rc 0, e(V) preserved (V11=.0025, V12=.0015, V33=.0049),
        ALLZERO = 0.   The matrix is NOT zeroed.
    (B) NONSYMMETRIC (one off-diagonal moved by .01):
        ereturn post -> rc 0 and e(V) comes back ALL ZERO. ALLZERO = 1.
    (C) the SAME nonsymmetric matrix through ereturn repost:
        rc 0 and e(V) comes back ALL ZERO TOO. ALLZERO = 1.
    (D) the -1e100 sentinel, symmetric: post and repost BOTH return it verbatim.

So of the two hazards the record names:
  - "highly singular" is FALSE for both post and repost;
  - "nonsymmetric" is TRUE -- and equally true of `repost`, so choosing repost
    buys no protection from it at all.

The RULING stands (keeping the sentinel diagonal is right, and the "fabricated
certainty" argument against zeroing is untouched). What does not stand is the
REASON: the hazard the record says the path avoids by using `repost` is a hazard
`repost` shares, and the hazards it says `repost` uniquely tolerates are ones
`post` tolerates too. What actually protects e(V) here is wooldid's symmetry
discipline (the `(S + S')/2` at ado:10523, the `makesymmetric` calls), not the
choice of repost over post.

