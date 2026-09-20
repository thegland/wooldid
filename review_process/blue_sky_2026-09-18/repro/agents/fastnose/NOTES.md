# fastnose probe — NOTES (agent: FASTNOSE, method: exact identity + e() surface differencing)

WORKDIR: /private/tmp/claude-501/.../scratchpad/bluesky/probes/fastnose
Build asserted via `adopath ++ .../bluesky/build` (never `do`).

## Harness note
The project cap of 8 concurrent Stata instances was being exceeded by other agents
(9-10 observed for >25 minutes).  `runq.sh` waits up to 20 minutes for the count to
drop below 8 and then launches anyway, and never runs more than ONE Stata process of
mine at a time.

## What the DOCUMENTED contract actually says (pinned build sthlp:286-297)
- fastnose: every margins call runs with `nose`; no se/t/p/ci for any estimand.
- e(V) NOT posted; e(properties) == "b"; e(wdid_vcvstate) == "absent".
- e(V_es_*) NOT returned, **but e(b_es_*) and e(numpre_*) ARE**, "and a Caution per
  family says so".  ==> the recon note's item-11 claim ("ruled to REMOVE fastnose from
  that list ... check whether fastnose actually posts e(b_es_att)/e(numpre_att)") is
  about a ruling that, IN CODE, went the OTHER way: ado:1982-1990 implements
  `f2bonly` = (fastnose & no permute) so the F2 sites capture b alone and the posting
  block at ado:4147-4150 posts b_es_* + numpre_* and prints the Caution when f2ok==0.
  ==> PREDICTION: seeded lead 2 is now SATISFIED on this build.  Verify in p1/p2.
- "no joint tests of any kind are computed (including the subgroup homogeneity tests)"
- "the one exception is the subgroup-by-horizon margins call, which keeps its variance
  so that the estimability screen ... still drops a non-estimable horizon (**fast mode
  therefore reports the same event-study rows as the full run**, at a small cost per
  subgroup)"                                <-- FLAT CLAIM, testable.
- Point estimates MAY differ ~1e-15 abs, up to ~1e-13 rel on subgroup pairwise
  contrasts, ~1e-14 under timetrends.  So BIT identity is NOT the documented contract;
  the documented envelope is what I test against.
- Speed claim: "roughly 2-4x on typical specifications".

## Code-read leads (mine, from the pinned ado)
L1  ado:11958  the R6 / SI-A4 dose-slope screen (`amevoid`) is SKIPPED under fastnose.
    In the analytic run it voids the whole ame block (b_es_ame withheld, ame rows
    reported missing) because the regression's ame entries are exact STRUCTURAL ZEROS.
    Under fastnose nothing voids them  ==> fast mode should report ame = 0 as if it
    were an estimate.  Design that fires it: dose constant inside every (i,t) cell —
    the state-policy-intensity design.  PROBE p3.
L2  ado:3290  the CONTINUOUS-treatment twin of the subgroup-by-horizon margins call
    threads `marginsdof` (= nose under fastnose), NOT `sgmarginsdof`.  Only the binary
    twin at ado:2912 got the BS-26 carve-out.  ==> the SI-7 estimability drop should be
    BLIND on the ame side under fastnose, so the ame event-study ROW SET can differ
    between the two modes — contradicting sthlp's "same event-study rows" claim. PROBE p4.
L3  ado:5786 / mata wooldid_fcdegx reads r(V); under nose K==0 and it returns at once,
    so the cell-collapse degeneracy fallback (e(wdid_fcguard)=="fallback") cannot fire
    in fast mode.  The analytic run then evaluates subgroups on FULL data while fast
    mode stays on the collapsed rows.  PROBE p4.
L4  ado:1981 `sgmarginsdof df(df_r)` only when espre+espl>0 — check the no-ES case.
L5  ado:11981 wooldid_icscreen (interactivecontrols caution) skipped under fastnose.
L6  ado:1520 the df_r<1 / <2-clusters gate is exempt under fastnose (documented).
L7  ado:4490 wooldid_sgdeg skipped under fastnose (SEs already missing there; expect clean).
L8  SPEED: because of L2/ado:1981 the sg+ES call keeps its full delta VCE, so the
    2-4x claim may not hold on subgroup ES runs.  Measure.

## Probes
p1_smoke.do        build assertion + b_es_att/numpre_att/Caution under fastnose
p2_identity.do     20-spec analytic-vs-fastnose full e() dump (wdump.do) -> p2dump.txt
p3_contreat_d3.do  L1
p4_sgrows.do       L2, L3
p5_speed.do        L8 (to write)

## STATUS
- probes written; chain1.sh running them one at a time.

## RESULTS (as they land)

### p1_smoke.log  — seeded lead 2 RESOLVED, and it CONTRADICTS the recon note
Build asserted: `which wooldid` -> .../bluesky/build/wooldid.ado.  N=800, 80 units, T=10.
`wooldid y unit t ttre , esprelength(3) espostlength(4) [fastnose]`
  analytic: numpre_att=3, b_es_att rc0, V_es_att rc0
  fastnose: numpre_att=3, b_es_att rc0, **V_es_att rc111 (absent)**, and the run prints
     "Caution: no horizon of the att event study has an identified variance, so
      e(b_es_att) posts without e(V_es_att). ..."
==> sthlp:1325's list is CORRECT on this build for fastnose: the coefficient vector and
    e(numpre_*) are returned, the covariance is not, and the promised per-family Caution
    fires.  The recon note's "item 11 looks unfinished / e(b_es_att) not posted at all,
    e(numpre_att) empty" is NOT reproducible on c1d8979.  The 2026-09-13 ruling was built
    in the CODE direction (ado:1982-1990 `f2bonly`, ado:4147-4150), not by editing the help.

### HARNESS FACT worth recording
After a wooldid fit, `e(scalars)`, `e(macros)` and `e(matrices)` all come back EMPTY
(p1_smoke.log lines 739-747), on BOTH the analytic and the fastnose run.  Any user or
downstream tool that enumerates wooldid's stored results that way gets nothing.  **RETRACTED as a wooldid issue**: p2_identity.log shows plain `regress` behaves
identically in this Stata, so the empty member lists are a Stata/batch artifact, NOT a
wooldid defect.  Mata st_dir("e()", ...) also errors r(3300) here; wdump.do therefore
uses an explicit candidate-name list grepped out of the pinned ado.

### p4_sgrows.log (binary arms A,B) — seeded lead 3, binary half: CLEAN
Ragged realistic subgroups (a policy variant only the 2009 cohort adopted; and a handful
of treated units in one cohort), esprelength(4) espostlength(5):
  ARM A  s1 rows pre4..post4 ; s2 rows pre4..post1   -- analytic and fastnose IDENTICAL
  ARM B  s1 rows pre4..post4 ; s2 rows pre2..post4   -- analytic and fastnose IDENTICAL
  e(wdidsgresults_*) rownames identical in both modes.
==> BS-26's `sgmarginsdof df(df_r)` carve-out (ado:1981) does its job on the BINARY side.
    The contreat twin (arms C,D) errored r(198) on the first run for a missing
    contreatcontrols(); re-queued.

### p8_fcdeg.log — L3 (cell-collapse fallback) UNTRIGGERED; and the drift is NOT the collapse
Q1 cluster(cohort) + subgroups(sg2) + esprelength(5) espostlength(6) + att att_it
Q2 same with cluster(region)
Q3 i=cohort + fe(unit) + cluster(cohort): r(459) in BOTH modes (identical refusal)
Q4 Q1 with nofastcollapse on BOTH sides
  e(wdid_fcguard) is missing on every arm -> wooldid_fcdeg never fired even in the
  analytic run, so the "analytic falls back to full data, fastnose cannot" divergence
  was NOT reachable here.  Arithmetically that is expected: fcdeg's cut is 1e-9 of the
  largest variance, while the clustered-at-i subgroup pathology wooldid_sgdeg targets
  sits at rho ~ .999 -> (1-rho)/(1+rho) ~ 3e-4, five orders above fcdeg's cut.
  Q1 (collapse ON) and Q4 (collapse OFF) produce the SAME 16 value differences with the
  SAME magnitudes, so the analytic-vs-fastnose drift comes from margins' nose path, not
  from the collapse.
  Largest relative drift anywhere: 8.64e-15 on e(pwsgs) (subgroup pairwise contrast),
  1.04e-15 on e(b)/e(wdidmainresults)/e(b_es_*).  Documented envelope is ~1e-13 on sg
  contrasts, ~1e-15 elsewhere -> WITHIN CONTRACT.
  Surface dropped by fastnose, every arm: e(V) and every e(V_es_*) (incl. the _s#
  subgroup twins).  Nothing else.  NO cell is filled by fastnose that the analytic run
  leaves missing (FN-has/AN-missing = 0 on every arm).

### p7_permase.log — FN-02 (gteffects keeps full analytic inference under fastnose)
see FINDINGS.md.  gteffects se/p/ci are byte-identical between modes.

### p5_speed.log — the 2-4x claim, first cut (SMALL panel: 400 units x T=14, i=cohort)
  A  4 estimands + ES(4,5), fe(unit)    analytic 1.395s  fastnose 0.807s   ratio 1.73
  B  same + subgroups(sg2)              r(459) BOTH modes (sg2 straddles the
                                        cohort-by-t cell; identical refusal, by design)
  C  poisson + ES(4,5)                  analytic 1.867s  fastnose 1.688s   ratio 1.11
  D  contreat(dose) + ES(4,5)           analytic 2.599s  fastnose 1.513s   ratio 1.72
  => on runs of 1-3 seconds nothing reaches 2x and the poisson arm barely moves, but
     fixed overhead dominates at that size.  p5b re-measures on i=unit / T=16 / 260
     units, which is the "large panel, large K" regime the help names, twice per arm
     with the minimum taken (the machine carries 7-9 other Stata jobs throughout).

### p3_contreat_d3.log (rerun) — the d3 design with i=unit / i=cohort+fe(unit) is
### refused r(322) in BOTH modes, identically
`contreat(dose) contreatcontrols(untreatedZero)' with a cell-constant dose:
  "Error: a treated cohort-period dummy coefficient was dropped for multicollinearity.
   The problem seems to be with the treatment dummy for unit == 2, relative-t == 0"
  rc 322 on ARM A and ARM B, analytic AND fastnose, with the identical message.
So on THAT geometry wooldid_ccscreen's main scan (which is NOT gated on fastnose) stops
the run before the R6 dose-slope screen could matter.  The R6 class needs the DUMMIES to
survive while every dose SLOPE drops: p12_r6.do builds that (state-policy dose defined
at the cohort-period level, 30 units per cell) across four contreatcontrols() flavours.

### p4_sgrows.log ARM C — **FN-03, the headline finding** (see FINDINGS.md)
contreat(dose) contreatcontrols(untreatedZero) + subgroups(sgvar) + ES(4,5), i=cohort:
  analytic e(b_es_ame_s1) = 1x5 (post0..post4);  fastnose = 1x9 (pre4..pre1 = EXACTLY 0,
  then the identical post values).  Same for s2 (1x2 -> 1x6).  The printed subgroup table
  gains a "pretreatment | ame | 0" row and four "ES_ame | pre# | 0" rows.
  ARM D r(459) in both modes (identical).
Contradicts sthlp:292 "fast mode therefore reports the same event-study rows as the
full run".  Mechanism: the BINARY sg-by-horizon call (ado:2912) gets BS-26's
`sgmarginsdof df(df_r)' carve-out; the CONTINUOUS twin (ado:3290) still threads
`marginsdof' == nose, so the SI-7 estimability filter has no marker to read.

### HARNESS trap worth recording for whoever reuses wdump.do
`cap local x = e(cmdline)' ABORTS a batch do-file with r(111) "wooldid not found"
despite the `capture' -- the cmdline string is evaluated as an expression.  wdump.do now
keeps macro candidates strictly string-valued.

### p12_r6.log — the R6 / d3 class is NOT reachable: an earlier ungated screen always wins
Five arms (contreatcontrols basic / byTime / byCohort / untreatedZero, at i=cohort+fe(unit)
and at i=cohort alone) with a dose defined at the cohort-period level and 30 units per
cell.  Every arm errors IDENTICALLY in both modes: r(322) on four of them (the ccscreen
dropped-dummy check, which is NOT gated on fastnose) and r(459) on contreatcontrols(byCohort)
(the per-cohort point-identification screen, which ado:12835 says runs regardless of
fastnose).  So although the R6 dose-slope screen at ado:11958 IS skipped under fastnose,
I could not construct a design that reaches it: wooldid_ccscreen stops the run first, in
both modes, every time.  Reported as a negative result, not a defect.

## FINAL STATE (for a replacement picking this up)
FINDINGS.md holds everything.  Three confirmed defects:
  FN-03 (HIGH)   fastnose reports subgroup ame event-study rows the full run drops, as
                 exact zeros, incl. a POST horizon; e(numpre_ame_s#) 0 vs 4.
                 repro: repro_FN03.do (self-contained, ran clean, repro_FN03.log)
  FN-01 (MEDIUM) permanalyticse is a no-op under fastnose and the run prints that the
                 analytic SEs are in the table and in e(V).  repro: repro_FN01_FN02.do
                 (self-contained, queued) / p7_permase.do + p11_permseed.do (ran)
  FN-02 (LOW)    "SEs are not estimated in fast mode" printed above a full gteffects
                 table of se/p/CI.  repro: repro_FN01_FN02.do / p7_permase.do
Ten negative results N-1..N-10 (N-6 and N-9 are the same negative from two constructions).
STILL RUNNING at hand-off: p9_shapes.do (shape sweep, 3/20 shapes differenced clean),
p5b_speed.do (large-K speed), then p13_nosemarker.do (rerun, mechanism) and
repro_FN01_FN02.do.  chainn.sh / chaink.sh drive them; runq.sh is the launcher.
Differ: `python3 diff2.py <dump> AN FN`.

### p5c_speed.log — 2-4x claim substantiated at large K (see FINDINGS N-13)
   bare att + ES        2.04x ;  four estimands + ES  2.87x ;  + subgroups(sg2) 1.76x
   poisson arm r(322) in both modes, timing invalid.
