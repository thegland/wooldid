# fastnose — FINDINGS

SUMMARY
  FN-03  HIGH    fastnose reports subgroup ame event-study rows the full run DROPS, as
                 exact zeros — including a POST-treatment horizon — and e(numpre_ame_s#)
                 goes 0 -> 4.  sthlp:292 promises the row sets are identical.
                 repro_FN03.do (self-contained), repro_FN03.log.   seeded (lead 3).
  FN-01  MEDIUM  permanalyticse is a complete NO-OP under fastnose; the run prints a
                 header and a footer saying the analytic SEs are in the table and in
                 e(V), and posts e(perm_analyticse)=1.  e(V) does not exist.
                 repro_FN01_FN02.do / p7_permase.do / p11_permseed.do.   not seeded.
  FN-02  LOW     "SEs and related materials are not estimated in fast mode" is printed
                 directly above a full gteffects table of se / p / 95% CI (correct, and
                 byte-identical to the analytic run's).  Documentation gap.  not seeded.
  N-1..N-10      ten negative results, including the resolution of seeded lead 2 (the
                 2026-09-13 ruling WAS built) and of seeded lead 3's binary half (clean).

Agent: FASTNOSE.  Method: exact identity + e()-surface differencing, analytic vs fastnose.
Build: /private/tmp/.../scratchpad/bluesky/build/wooldid.ado (c1d8979), asserted by
`which wooldid` in every probe log.

---

## FN-01  `permanalyticse` is silently defeated by `fastnose`, and the run then PRINTS
##        that the analytic standard errors are in the table and in e(V)
STATUS: CONFIRMED.  severity medium.  realistic: yes (moderately).
repro : /private/tmp/.../probes/fastnose/p7_permase.do   (block "P3")
log   : p7_permase.log lines 1160-1310

wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(20) ///
        permanalyticse fastnose esprelength(2) espostlength(2)

ACCEPTED, rc 0.  wooldid_validate (ado:5374-5412) refuses fastnose with unconditionalse,
impvalues and pretrendpower(), and ado:5882 refuses it with permtstat/permtci -- but
NOTHING screens permanalyticse, whose entire purpose is to restore the analytic standard
errors the permute() fast-mode default drops.  With fastnose in force `marginsdof' is
`nose' (ado:1973) so no analytic standard error is computed at all, yet:

  header  (ado:8215, e(wdid_vcespec)):
    "Standard Errors: randomization inference (permute(shuffle), 20 replications);
     se_analytic = the analytic cluster(cohort) standard error, shown for reference only"
  table   : the se_analytic column is entirely "." on every row
  footer  (ado:11679):
    "The standard-error column is the ANALYTIC standard error (se_analytic,
     permanalyticse), shown for reference only; ... e(V) and e(V_es_*) hold the analytic
     covariance, so test or lincom run by hand on them are analytic, not permutation,
     quantities."
  measured: e(wdid_vcvstate) == "absent";  `confirm matrix e(V)' -> r(111);
            `test main_att = 0' -> r(321).

So two printed statements are false, and a requested option is dropped in silence.  The
same run WITHOUT fastnose (block P2) gives vcvstate "analytic-blocks", e(V) present and
`test' rc 0; the permute default WITHOUT permanalyticse (block P1) prints the correct
"analytic standard errors not computed (permanalyticse restores them)".
The comment at ado:8246 states the assumption that breaks: "Fast mode only: under
permanalyticse an e(V) exists and the tools already work from it."
FIX SHAPE: refuse `fastnose' + `permanalyticse' the way permtstat/permtci are refused,
or let permanalyticse win over fastnose.

---

## NEGATIVE RESULTS (see probed_but_clean in the structured return)

## FN-02  fast mode prints "SEs and related materials are not estimated in fast mode"
##        and then prints a full table of standard errors, p-values and 95% CIs
STATUS: CONFIRMED.  severity low-medium (doc-vs-actual / UX, not a wrong number).
realistic: yes -- `gteffects' + `fastnose' is a natural exploratory pairing.
repro : p7_permase.do blocks G1 (analytic) and G2 (fastnose)
log   : p7_permase.log lines 731-782 (G1) and 840-890 (G2)

wooldid y cohort t ttre , fe(unit) gteffects esprelength(2) espostlength(3) [fastnose]

Under fastnose the header says "Standard Errors: not computed (fast mode)" and
"SEs and related materials are not estimated in fast mode."; every se/t/p/lb/ub cell of
the Main Results table is ".".  Immediately below it, the gteffects table prints
  r18 | estimate 1.721636  se .2322108  p 6.05e-13  lb_95ci 1.265295  ub_95ci 2.177978
and every one of its 18 rows x 4 inference columns is BYTE-FOR-BYTE the analytic run's.
Mechanism (not a bug in the numbers): the gteffects cells are read off the STAGE-1
regression's own e(b)/e(V) (ado:1765-1773 `matrix `gtsed' = vecdiag(e(V))'), and fastnose
only changes what `margins' is asked for.  sthlp:286 says fast mode computes
"no standard errors, test statistics, p-values, or confidence intervals ... for any
estimand", and the gteffects option's own paragraph (sthlp:1134-1145) never mentions
fastnose.  Nothing in the run reconciles the two.
POSITIVE SIDE: this is also a clean identity result -- gteffects inference does not move
between the two modes.
FIX SHAPE: one clause in the fastnose paragraph (and/or in the gteffects paragraph)
saying the gteffects cells keep the stage-1 analytic inference in fast mode.

## FN-03  fastnose REPORTS EVENT-STUDY ROWS THE FULL RUN DROPS, as exact zeros, under
##        contreat() + subgroups() -- and the help promises the row sets are identical
STATUS: CONFIRMED.  severity HIGH.  realistic: YES (continuous-treatment DiD with a
subgroup split and an event-study window is a completely ordinary specification).
repro  : /private/tmp/.../probes/fastnose/repro_FN03.do  (self-contained, builds its
         own 120-unit panel; also p4_sgrows.do ARM C)
log    : p4_sgrows.log lines 1870-1905 (analytic) and 2100-2260 (fastnose)

wooldid y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) ///
        subgroups(sgvar) esprelength(4) espostlength(5) [fastnose]

  ANALYTIC   e(b_es_ame_s1) : 1 x 5   colnames  post0 post1 post2 post3 post4
             e(b_es_ame_s2) : 1 x 2   colnames  post0 post1
             e(wdidsgresults_1) rows: ... main:ame ... ES_ame:post0 ... post4
             (no pretreatment:ame row, no ES_ame:pre# rows)

  FASTNOSE   e(b_es_ame_s1) : 1 x 9   colnames  pre4 pre3 pre2 pre1 post0 ... post4
                              values      0    0    0    0  1.3478682 ... (post = identical)
             e(b_es_ame_s2) : 1 x 6   colnames  pre4 pre3 pre2 pre1 post0 post1
                              values      0    0    0    0  2.0549423  2.03283
             e(wdidsgresults_2) GAINS  pretreatment:ame  and  ES_ame:pre4..pre1
  and the PRINTED subgroup table gains
             pretreatment | ame |  0
             ES_ame       | pre4 | 0 ; pre3 | 0 ; pre2 | 0 ; pre1 | 0

Those zeros are the regression's STRUCTURAL zeros: with contreatcontrols(untreatedZero)
the dose is zero outside treated cohort-periods by construction, so margins marks those
subgroup-by-horizon margins NOT ESTIMABLE and the SI-7 screen drops them in the full run.

MECHANISM.  ado:1974-1981 (BS-26) forces `sgmarginsdof' to df(`df_r') under fastnose for
the BINARY subgroup-by-horizon call at ado:2912, exactly so the omission marker survives
for wooldid_sglbl to filter on.  Its CONTINUOUS twin at ado:3290 still threads
`marginsdof', which fastnose sets to `nose' (ado:1973).  The carve-out was never applied
to the ame side.  (The binary side is clean: p4 ARMs A and B show identical att row sets
in both modes, including a genuinely ragged subgroup.)

WHY IT MATTERS BEYOND THE ROW COUNT
 * the printed table reads as a PERFECT pre-trend: four pre-treatment ame effects at
   exactly 0.  In a mode sold for specification search that is the worst possible
   direction for a fabricated number to point.
 * e(b_es_ame_s#) and e(numpre_ame_s#) are the documented honestdid/pretrends inputs;
   fast mode hands those consumers a different event-time geometry than the full run.
 * it contradicts sthlp:292 verbatim: "the one exception is the subgroup-by-horizon
   margins call, which keeps its variance so that the estimability screen described
   under subgroups() still drops a non-estimable horizon (fast mode therefore reports
   the same event-study rows as the full run, at a small cost per subgroup)".

FIX SHAPE: thread `sgmarginsdof' (not `marginsdof') at ado:3290, and null rows 2-6 of
that call's f9tab afterwards the way the binary twin does at ado:2978-2986.

### FN-03 SCOPING (measured, p4_sgrows.log lines 1640-1700)
The divergence is SUBGROUP-SPECIFIC.  In the ANALYTIC run the OVERALL table already
prints the ame structural zeros --
      pretreatment | ame | 0 (se .)
      ES_ame       | pre4..pre1 | 0 (se .)
-- because the MAIN ame event-study call (ado:2651) has no SI-7 filter at all; only the
SUBGROUP calls do.  So what fastnose changes is the SUBGROUP row set: the full run drops
those subgroup rows, fast mode keeps them.  Two consequences stand regardless:
  (i)  sthlp:292's "fast mode therefore reports the same event-study rows as the full
       run" is FALSE as written;
  (ii) e(b_es_ame_s#) changes SHAPE between the two modes (1x5 -> 1x9, 1x2 -> 1x6) and
       e(numpre_ame_s#) changes with it, and those are the documented honestdid /
       pretrends inputs.
A separate observation, not filed as a defect because it is the analytic behaviour in
both modes: wooldid is internally inconsistent about those structural zeros -- the
overall table prints them, the subgroup tables drop them.

---

# NEGATIVE RESULTS (attacked hard, came back clean)

## N-1  Seeded lead 2 is NOT reproducible on the pinned build -- the 2026-09-13 ruling
##      WAS implemented, in the code direction
The recon note (known_bydesign.md sec.7, item 11) predicted a live doc-vs-code divergence:
"Check whether fastnose actually posts e(b_es_att)/e(numpre_att); if it does not, this is
a live doc-vs-code divergence."  MEASURED (p1_smoke.log, p6_contract.log):
  fastnose posts e(b_es_att) (rc 0), posts e(numpre_att) = 3 (equal to the analytic run's),
  withholds e(V_es_att) (rc 111), and PRINTS the per-family Caution sthlp promises:
  "Caution: no horizon of the att event study has an identified variance, so e(b_es_att)
   posts without e(V_es_att). ..."
sthlp:1325's list is therefore correct for fastnose as it stands, and the ruling was
built by changing the CODE (ado:1982-1990 `f2bonly'; posting block ado:4147-4150), not by
editing the help.  (FN-03 above is a DIFFERENT, still-live divergence in the same option.)

## N-2  Seeded lead 3, BINARY half: fastnose and the full run agree on the att row set
p4_sgrows.log ARMs A and B, two realistic ragged subgroups (a policy variant only the
2009 cohort adopted; a handful of treated units in one cohort), esprelength(4)
espostlength(5):
   ARM A  s1 = pre4..post4 , s2 = pre4..post1   identical in both modes
   ARM B  s1 = pre4..post4 , s2 = pre2..post4   identical in both modes
   e(wdidsgresults_*) rownames identical in both modes.
BS-26's `sgmarginsdof df(df_r)' carve-out (ado:1981) does exactly what its comment says
on the binary side.  Only the continuous twin was missed (FN-03).

## N-3  Point identity is inside the documented envelope everywhere I looked
Full %21x differencing of every e() member (wdump.do), analytic vs fastnose:
  p8 (4 arms: cluster(cohort)/cluster(region)/nofastcollapse, subgroups, ES(5,6), att+att_it)
     worst relative difference anywhere  8.64e-15  (e(pwsgs), a subgroup pairwise contrast)
     worst on e(b)/e(wdidmainresults)/e(b_es_*)   1.04e-15
  sthlp:295 allows ~1e-15 absolute, up to ~1e-13 relative on subgroup pairwise contrasts.
  WITHIN CONTRACT.
  NO cell is ever filled by fastnose that the analytic run leaves missing (the differ
  counts FN-has/AN-missing = 0 on every arm).
  Turning the cell collapse OFF on both sides (nofastcollapse) reproduces the SAME 16
  differences with the SAME magnitudes -> the drift is margins' nose path, not the collapse.

## N-4  L3 (the cell-collapse degeneracy fallback) could not be triggered at all
wooldid_fcdeg reads r(V); under `nose' K==0 and it returns immediately, so
e(wdid_fcguard)=="fallback" is unreachable in fast mode while the analytic run can take
it and re-evaluate everything on the full data.  I could not MAKE the analytic run take
it: e(wdid_fcguard) is missing on every arm of p8 (4/5/6 clusters against 11 horizons,
two subgroups).  Arithmetic says why: fcdeg's cut is 1e-9 of the largest margin variance,
while the clustered-at-i subgroup pathology wooldid_sgdeg targets sits at rho ~ .999, i.e.
(1-rho)/(1+rho) ~ 3e-4 -- five orders of magnitude above the cut.  The structural asymmetry
is real; its reachable consequence is unmeasured.

## N-5  The published fastnose e() contract holds exactly (p6_contract.log)
  e(properties) == "b" (analytic: "b V");  e(wdid_vcvstate) == "absent" (analytic:
  "analytic-blocks");  e(V) absent (r111);  test -> r(321); lincom -> r(321);
  _se[] -> r(111); estat vce -> r(322); _b[] still works (rc 0).
  est store / est restore round-trip: rc 0, properties still "b".  Replay: rc 0.
  simultaneousbands ACCEPTED with the documented caution and no e(esband_*).
  makeplots: rc 0.
  Refusals all fire r(198): unconditionalse, impvalues, jointtests, adversarial,
  pretrendpower(), bootstrap(standard).  permtstat+fastnose r(198) (p7 block P4).
  subgroups() with no ES window: pairwise table prints with "." inference in both modes.

## N-6  The d3 ("dose constant in the cell") geometry is refused identically in both modes
p3_contreat_d3.log: `contreat(dose) contreatcontrols(untreatedZero)' with a cell-constant
dose errors r(322) ("a treated cohort-period dummy coefficient was dropped for
multicollinearity") in BOTH modes with the identical message, at i=unit and at
i=cohort+fe(unit).  wooldid_ccscreen's main scan is not gated on fastnose, so it stops the
run before the R6 dose-slope screen (which IS skipped under fastnose, ado:11958) could
matter on that geometry.  p12_r6.do pushes harder with four contreatcontrols() flavours
and a dose defined at the cohort-period level with 30 units per cell.

## N-7  A fastnose run leaves NO state that changes the next run (p10_contam.log)
Order A1 (analytic) -> F1 (fastnose) -> A2 (analytic, same spec) -> F2, on
`att att_it att_i att_itime subgroups(sg2) esprelength(3) espostlength(3) gteffects'
(+ jointtests adversarial on the analytic arms).  Full %21x differencing:
   A1 vs A2 : 0 value differences, 0 missingness differences  (BIT IDENTICAL)
   F1 vs F2 : 0 value differences, 0 missingness differences  (BIT IDENTICAL)
So `marginsdof'/`suppressjointtests'/`fcguard' and the Mata/global state fastnose touches
do not leak across calls.

## N-8  Point identity on the heaviest ordinary spec (p10, A1 vs F1)
4 estimands + subgroups(sg2) + ES(3,3) + gteffects:
   worst relative difference 1.13e-15 (e(pwsgs)); 4.88e-16 on e(b)/e(b_es_*).
   Surface dropped by fastnose: e(V), every e(V_es_*) (incl. _s# twins), every
   e(adversarial_line_*), and the macros e(vce), e(wdid_adversarial), e(wdid_vblocks),
   e(wdid_vblockn).  Surface ADDED: e(wdid_fastnose), e(wdid_suppressjt).  Nothing else.
   e(gteffects) present in BOTH with ZERO differing cells -- byte-identical inference
   columns, which is FN-02's other half.

### FN-03 CONFIRMED by the self-contained repro (repro_FN03.log)
repro_FN03.do builds its own 120-unit / T=11 panel, no project data, no project generator.
`wooldid y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit)
 subgroups(sgvar) esprelength(4) espostlength(4) [fastnose]`
  ANALYTIC  e(b_es_ame_s1) colnames: post0 post1 post2 post3   e(numpre_ame_s1) = 0
            e(b_es_ame_s2) colnames: post0 post1               e(numpre_ame_s2) = 0
  FASTNOSE  e(b_es_ame_s1) colnames: pre4 pre3 pre2 pre1 post0 post1 post2 post3
                                                               e(numpre_ame_s1) = 4
            e(b_es_ame_s2) colnames: pre4 pre3 pre2 pre1 post0 post1
                                                               e(numpre_ame_s2) = 4
  and e(wdidsgresults_1/2) gain a `pretreatment:ame' row under fastnose.
CONTROL, the SAME spec WITHOUT subgroups():
  ANALYTIC  e(b_es_ame) colnames pre4..post3, e(numpre_ame) = 4
  FASTNOSE  e(b_es_ame) colnames pre4..post3, e(numpre_ame) = 4     -> IDENTICAL.
So `e(numpre_ame_s#)', the argument the help tells you to hand honestdid, is 0 in the
full run and 4 in fast mode on the identical specification and data.

### FN-03, THE SHARPEST FORM: a POST-treatment horizon fabricated at exactly 0
Same repro, third arm, `contreatcontrols(basic)' instead of untreatedZero
(repro_FN03.log lines 1365-1690):
  ANALYTIC  e(b_es_ame_s2) = 1 x 1   post0 = -.300694594
  FASTNOSE  e(b_es_ame_s2) = 1 x 6   pre4 pre3 pre2 pre1 post0 post1
                                     =  0    0    0    0  -.300694594   0
The extra `post1' column is a POST-TREATMENT horizon that the full run drops as not
estimable and fast mode reports as an average marginal effect of exactly ZERO.  This is
no longer only about the structural pre-period zeros: fast mode publishes a
post-treatment treatment-effect estimate that the full run refuses to publish, on the
same data and the same specification, in the mode the help recommends for specification
search.  (s1 in the same arm: analytic post0..post3, fastnose pre4..pre1 = 0 prepended.)

### FN-01 QUANTIFIED with the RNG pinned (p11_permseed.log, p11dump.txt)
Same data, `set seed 777777' immediately before each call, permreps(60):
  R1 = permute(shuffle)                              (the fast-mode default)
  R2 = permute(shuffle) permanalyticse
  R3 = permute(shuffle) permanalyticse fastnose      (the suspect)
Full %21x differencing of the whole e() surface:
  R1 vs R3 differ in EXACTLY THREE items and in NO number at all:
     e(cmdline)            (trivially)
     e(wdid_vcespec)       R1: "... analytic standard errors not computed
                                 (permanalyticse restores them)"
                           R3: "... se_analytic = the analytic cluster(cohort) standard
                                 error, shown for reference only"        <-- FALSE
     e(perm_analyticse)    absent in R1, = 1 in R3                       <-- FALSE
  Every matrix and every cell -- e(perm_p), e(perm_ci), e(perm_kernel), e(b),
  e(wdidmainresults) -- is BYTE-IDENTICAL between R1 and R3.
  R2 (permanalyticse working) really does add e(V), e(V_es_*) and a populated se column.
CONCLUSION: `permanalyticse' is a COMPLETE NO-OP when `fastnose' is also given.  The only
thing the option produces is a header and a footer asserting the analytic standard errors
are in the table and in e(V), plus e(perm_analyticse)=1, which ado:11679 branches on to
print "e(V) and e(V_es_*) hold the analytic covariance, so test or lincom run by hand on
them are analytic".  e(V) does not exist and `test' returns r(321).
Severity therefore stays MEDIUM: nothing is mis-COMPUTED; a requested option is dropped
in silence and the run prints two false statements about it (and posts a machine-readable
flag that a replay or a downstream tool would believe).

## N-9  The R6 / d3 dose-slope screen IS skipped under fastnose (ado:11958) but I could
##      not reach it: an EARLIER, ungated screen stops the run in both modes every time
p12_r6.log, five arms (contreatcontrols basic / byTime / byCohort / untreatedZero, at
i=cohort+fe(unit) and at i=cohort alone), dose defined at the cohort-period level with
30 units per cell so the cohort-period dummies are well identified:
   A1 basic          r(322) BOTH modes      A2 byTime     r(322) BOTH modes
   A3 byCohort       r(459) BOTH modes      A4 untreatedZero r(322) BOTH modes
   A5 basic, no fe() r(322) BOTH modes
r(322) is wooldid_ccscreen's dropped-treated-dummy check, which is NOT gated on fastnose;
r(459) is the byCohort point-identification screen, which ado:12835 documents as running
regardless of fastnose.  Combined with p3 (same outcome on the cell-constant-dose
geometry), the R6 skip looks shielded in practice.  NOT filed as a defect.  What would
have shown up: an arm with rc 0 where the analytic run printed the "ame results as
missing" caution and fast mode printed ame numbers.

## N-10  FN-01 does not change any number (p11_permseed.log)
With the RNG pinned, permute()+permanalyticse+fastnose is byte-identical to the plain
permute() fast-mode default in every stored matrix and cell.  See the FN-01 quantification
above.  So the permutation inference is sound; only the labelling is wrong.

---

# CALIBRATION / CORRECTIONS applied after review

* FN-03's MECHANISM (that `nose' strips margins' not-estimable marker, leaving
  wooldid_sglbl nothing to filter on) is INFERRED from the code -- ado:1974-1981's own
  comment states it, and ado:3290 threads `marginsdof' rather than `sgmarginsdof'.  The
  direct margins-level demonstration (p13_nosemarker.do) is pending; my first version of
  it produced no marker in EITHER mode and was rebuilt.  The DIVERGENCE itself is
  measured three times (p4 ARM C, repro_FN03 untreatedZero arm, repro_FN03 basic arm).
* FN-03 severity HIGH, not critical: every standard error is "." in fast mode, so no
  inference is drawn from the fabricated rows, and the help tells you to re-run without
  fastnose before reporting.  The sharpest form is the `basic' arm's post1 = 0.
* FN-02 severity LOW: the numbers are correct and a note under the table says where they
  come from.  It is a documentation gap, not a hazard.
* N-6 and N-9 are the SAME negative (the d3 / R6 class is unreachable behind
  wooldid_ccscreen) measured from two different constructions -- count them once.

## N-11  Shape sweep: identity holds on every realistic panel shape completed
p9_shapes.do runs ONE spec (`att att_it att_i att_itime esprelength(3) espostlength(3)',
i = unit) on named realistic shapes and differences the whole e() surface in %21x.
Worst relative analytic-vs-fastnose difference per shape, and the count of cells fastnose
FILLS that the analytic run leaves missing:
     canon        2.73e-15   extra-filled 0
     unbal        0          extra-filled 0   (entry/exit/gaps: EXACTLY identical)
     rotate       8.51e-16   extra-filled 0
     admin        2.15e-15   extra-filled 0   (700 never-treated, one 260-unit cohort,
                                               four 2-4 unit cohorts, 11206 rows)
     tinytreated  7.09e-16   extra-filled 0
     rcs          6.03e-16   extra-filled 0   (repeated cross sections)
     biennial     2.54e-15   extra-filled 0
     holeyear     7.96e-16   extra-filled 0
     fewG5        2.76e-16   extra-filled 0
     pw_survey    0          extra-filled 0   (pweights: EXACTLY identical)
NO shape produced a dimension, rowname or colname difference in any e() matrix, and no
shape produced an e() item present under fastnose other than e(wdid_fastnose).
The surface dropped is always exactly e(V), e(V_es_att), e(V_es_att_it) plus the macros
e(vce), e(wdid_vblocks), e(wdid_vblockn).

## N-12  Option-matrix sweep: identity holds on 17 of 20 configurations, the other 3
##       being identical REFUSALS in both modes
p2_identity.do, one 80-unit T=10 panel, 20 specifications each run analytic and fastnose
with the full e() surface differenced in %21x.  Worst relative difference per spec:
  S1 bare att                                      0
  S2 four estimands                                0
  S3 ES(3,4)                                       1.26e-15
  S4 ES + esfixedbaseperiod + esrelativeto(-2)     6.00e-16
  S5 subgroups(sg2)                                6.79e-16 (e(pwsgs))
  S6 subgroups + ES                                1.95e-15 (e(pwsgs))
  S7 subgroups + ES + four estimands               1.95e-15 (e(pwsgs))
  S8 i=cohort + fe(unit) + ES + four estimands     0
  S9 controls(x1 x2) + ES                          3.41e-15   <- the largest anywhere
  S10 [aw=w] + ES + four estimands                 0
  S11 if x2 > -0.5 + ES                            1.18e-16
  S13 gteffects + ES                               0
  S14 cluster(region) + ES                         0
  S15 robust + ES                                  0
  S16 contreat + fe(unit) + ES (no subgroups)      0     <- matches the FN-03 control
  S18 interactivecontrols(x2) + ES                 0
  S20 subgroups(sgtiny) + ES                       0
  S12 poisson ES        r(322) in BOTH modes, same message
  S17 makecfxplots      r(198) in BOTH modes, same message
  S19 coarsecohortcontrols r(198) in BOTH modes, same message
No spec produced a dimension / rowname / colname difference in any e() matrix, and no
spec produced a cell filled under fastnose that the analytic run leaves missing.

### FN-01 / FN-02 CONFIRMED by the self-contained repro (repro_FN01_FN02.log)
repro_FN01_FN02.do builds its own 120-unit / T=10 panel.
  FN-01 control  permute(shuffle) permreps(40) permanalyticse            -> rc 0,
        e(wdid_vcvstate)="analytic-blocks", e(V) present, `test' rc 0.
  FN-01 defect   ... permanalyticse fastnose                             -> rc 0,
        e(wdid_vcvstate)="absent", e(V) ABSENT (r111), `test' r(321),
        e(perm_analyticse)=1, and
        e(wdid_vcespec) = "randomization inference (permute(shuffle), 40 replications);
                           se_analytic = the analytic cluster(cohort) standard error,
                           shown for reference only".
  FN-02  both gteffects runs rc 0; the fastnose run prints "SEs and related materials are
        not estimated in fast mode" and then a full gteffects se/p/CI table.

### p13_nosemarker: the plain-margins mechanism test is INCONCLUSIVE
Two constructions (dropping the cell, and setting y missing in it) both left the level out
of margins' over() list entirely, so no not-estimable marker appeared in EITHER mode and
the test says nothing.  FN-03's mechanism therefore rests on ado:1974-1981's own comment
("the SI-7 estimability drop reads margins' omission marker, which nose suppresses") plus
the fact that ado:3290 threads `marginsdof' where ado:2912 threads `sgmarginsdof'.  The
DIVERGENCE is measured three times and does not depend on the mechanism being right.

## N-13  The "roughly 2-4x" speed claim IS substantiated in the regime the help names
p5c_speed.do, 140 units x T=13 x 4 cohorts, i = unit (so K, the number of cohort-by-period
parameters margins differentiates, is large -- the help says K is what drives the cost).
One analytic / fastnose pair per spec, machine carrying 7-8 other Stata jobs throughout:
   bare att + ES(4,4)                      16.099s -> 7.892s    ratio 2.04
   four estimands + ES(4,4)                32.570s -> 11.347s   ratio 2.87
   att+att_it + ES(4,4) + subgroups(sg2)   67.969s -> 38.508s   ratio 1.76
   poisson + ES(4,4)                       r(322) in BOTH modes -- not a valid timing
So sthlp:288's "roughly 2-4x on typical specifications" HOLDS (2.04, 2.87).
The subgroup arm at 1.76 sits just below the band, which is exactly what BS-26's carve-out
predicts: the subgroup-by-horizon margins call keeps its full delta-method variance in
fast mode.  The help calls that "a small cost per subgroup"; measured here it is about a
40 percent erosion of the speed-up (2.9 -> 1.8).  NOT filed as a defect -- one timing on
a loaded machine, and "roughly"/"small" are not falsifiable at that resolution.  Recorded
so the author can calibrate the sentence if he wants to.
An earlier run on a SMALL panel (p5_speed.log: 400 units, T=14, i=cohort+fe(unit),
runtimes 1-3s) gave 1.73 / 1.11 / 1.72 -- fixed overhead dominates there and the numbers
do not test the claim, which the help scopes to large panels.
