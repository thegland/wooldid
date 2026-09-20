# SKEPTIC pass - anticipation(), replay, level()   (findings AR-1 .. AR-4)

Every result below was produced by a script I wrote from the claim's own terms, on data I
generated, without running the other agent's do-files. Pinned build asserted via `which wooldid`
in each log.

Self-contained repros (each generates its own panel):
- `repro_SK_LEVEL.do`       -> `repro_SK_LEVEL.log`
- `repro_SK_REPLAYGUARD.do` -> `repro_SK_REPLAYGUARD.log`
- `repro_SK_REPLAYDIFF.do`  -> `repro_SK_REPLAYDIFF.log` + `R{J,B,C,P}_{live,replay}.txt`

Supporting probes: s1_render.do, s2_setlevel.do, s3_maps.do (+ s3_full.csv), s4_level.do,
s5_replay3.do, s6_cautions.do, s7_more.do, s8_antic.do, s9_levelprop.do, s10_mix.do.

---

## AR-4 - CONFIRMED (high). Mechanism correct; SCOPE understated by ~2.5x.

All four sub-claims reproduce on my own panel (repro_SK_LEVEL.log):

  ## SETLEVEL 99.9/99.99/99.91/32.22 : Stata's own `set level' accepts all four, rc 0
  ## OPT level(99.9)    -> rc 198, nothing posted
  ## OPT level(99.99)   -> rc 198  (the help's own documented upper bound)
  ## OPT level(99.91)   -> rc 0, e(level)=99.91
  ## SESS 99.9  (no level() typed) -> rc 198, LEVEL0=[99.90000000000001], LASTRC=[]
  ## SESS 32.22 (no level() typed) -> rc 198  (well below the claimed 64.01 floor)
  ## SESS 99.91 (no level() typed) -> rc 198  (although 99.91 PASSES as an option)
  ## C1 (`set level 99.9' then level(95)) -> rc 198 with e(N)=720, e(cmd)=wooldid,
        e(level)=95, b_att=3.284012 - a fully successful, fully posted estimation
        returning an error
  ## C2 c(level) AFTER = 95, not the caller's 99.9
  ## C3 LEVEL0=[99.90000000000001], CLEANUP=[fullmodel], LASTRC=[] all left set
  ## D1/D2 the failing call leaves LASTRC empty; a bare `wooldid' then REPLAYS the
        pre-failure results, rc 0, e(N)=720
  ## D3/D4 control: an ordinary failure DOES arm the guard and replay IS refused

MECHANISM - one correction. The string "level() can have at most two digits after the decimal
point" does NOT exist anywhere in wooldid.ado (grep is empty). It is Stata's own `set level'
message; wooldid's contribution is only the macro rendering it feeds in. The prober's CLAIM text
attributes it correctly; the title reads as though wooldid validates it.

SCOPE - the prober's rule covers only half the exposure. There are TWO different failure sets,
because c(level) is not the same double as the parsed literal (s3_maps.log):

  HEX L=99.91  parsed=+1.8fa3d70a3d70aX+006   c(level)=+1.8fa3d70a3d70bX+006   (one ulp apart)

Full enumeration of all 9000 legal two-decimal levels, both maps measured with TYPED literals
(s3_full.csv for map A, s11_mapb.csv for map B; cross-tabulated in this directory):

  * Map A - the level() OPTION path (`syntax [, Level(real 0)]' -> macro -> set level):
    1152 / 9000 fail (12.8%), first failure L = 64.01. Exactly the prober's rule and count.
  * Map B - the SESSION path (`set level L' -> `global WOOLDID_LEVEL0 = c(level)' -> restore):
    1768 / 9000 fail (19.6%), first failure L = 32.22, and `set level' itself refuses NONE of
    the 9000. Map A is a strict subset of map B (both = 1152, A-only = 0, B-only = 616), so
    616 levels - 99.91, 99.94, 32.22, ... - are ACCEPTED as a level() option but poison the
    session restore.

Correct statement: a level() you TYPE is refused on 12.8% of legal two-decimal values; a session
level you SET breaks every subsequent wooldid call on 19.6% of them, including 616 values
wooldid itself accepts as an option.

REALISTIC? Yes for 99.9 / 99.99 (the documented bound). Every integer level is safe, so the
trigger is a non-integer `set level', or a `foreach L in 90 95 99 99.9' sweep.

---

## AR-3 - CONFIRMED (severity revised medium -> low)

repro_SK_REPLAYGUARD.log, reproduced exactly on my own data:

  ## A good wooldid: e(cmd)=[wooldid] e(N)=720 LASTRC=[]
  ## B another command owns e(): e(cmd)=[regress] LASTRC=[]
  Error: command line input should feature all of y i t ttre; too few variables specified.
  ## C rc=198 LASTRC=[198]   <- the replay ATTEMPT armed the marker
  ## D estimates restore GOOD: e(cmd)=[wooldid] e(N)=720  (a perfectly valid stored model)
  Error: the last wooldid call failed with r(198); the stored results predate that failed call.
  ## D rc=198
  ## E after a refit, replay of the stored model works again: rc=0

No repo ruling covers the "foreign e(cmd)" case (`grep -rn NFRPL review_process/' returns only
echoes of the ado's own comment). Downgraded to LOW: no number is wrong, both messages print,
and any successful wooldid call clears the marker (## E).

---

## AR-2 - CONFIRMED (low), and the divergence is TWO lines, not one

repro_SK_REPLAYDIFF.do, diff of RB_live.txt vs RB_replay.txt:

  4d3
  <
  20a20
  >

i.e. replay prints one FEWER blank line at the top of the output AND one EXTRA blank line before
the `Bootstrap Inference' block. Identical pattern under permute() (RP: 15d3, 30a19). Net line
count is unchanged, which is why a line-count check misses it. The analytic case (RJ) has no
blank-line divergence at all once the cautions are set aside - so this is engine-path-specific.
AR-2 reported only the second half.

---

## AR-1 - PARTIAL (severity revised medium -> low): behaviour confirmed, mechanism refuted

Confirmed on my own data that replay drops results-describing cautions (diffs of RJ, RP, RC):

  * RJ - the jointtests + cluster-at-i "its p-value should not be read as a calibrated test"
    caution (6 lines) is gone; the rest is byte-identical.
  * RP - the permreps(25) resolution-floor caution and the permute exchangeability caution
    (12 lines) are gone.
  * RC - the strongest case. Live prints 24 lines of caution including "the within-cell average
    marginal effects (ame) ... are NOT IDENTIFIED in this design" and "standard errors were not
    identified for 1 of 2 reported estimates". The replay prints:

             |   estimate          se         t         p     lb_95ci     ub_95ci    relyr
    main     |
         att |   3.284012   .17220992     19.07     0.000   2.9394207   3.6286032        .
         ame |          .           .         .         .           .           .        .

    with NO explanation whatsoever for the blank ame row.

BUT THE STATED MECHANISM IS WRONG. AR-1 argues this is "an omission rather than a policy" because
wooldid_display already reconstructs other results-describing cautions from e(). I dumped the full
e() surface of that contreat run (s7_more.log, @@Cmac / @@Csca):

  macros : wdid_replayv wdid_esrel wdid_espre wdid_suppressjt wdid_cetype wdid_fcol
           wdid_semielas wdid_estimands wdid_contreat wdid_vcespec vcetype vce clustvar
           wdid_vcvstate estimate_type cmd cmdline properties depvar
  scalars: df_r N Nclus ll ll_0 r2 r2adj r2_within r2_withinadj rmse rss N_dropped N_clust
           level wdid_fcncell

There is NO flag recording the ame-not-identified condition, and none recording the
jointtests-uncalibrated condition. Two consequences, in the right order:

  (i) What makes the drop POLICY rather than an omission is that sthlp:1458 documents it
      explicitly. That alone answers the prober's framing.
  (ii) Separately, the sgnopre / roth precedent the prober cites does not transfer:
      wooldid_display reconstructs those from e(wdid_sgnopre*) / e(roth_*), and there is no
      analogous state for these cautions. So the remedy is new e() postings, not a display edit
      - a larger change than "wooldid_display already does this elsewhere" implies.

Downgraded to LOW because sthlp:1458
explicitly documents that cautions are not repeated; what remains is that the sentence's
JUSTIFICATION ("which describe the estimation run, not the results") is inaccurate for at least
five of them.

---

## NEW (mine): stale $WOOLDID_LEVEL0 masks the true return code and skips the cleanup epilogue

repro_SK_LEVEL.log, block 5:

  ## E0 baseline bad-cluster rc=111 (expect 111) LASTRC=[111]
  ## E1 poisoning call rc=198 LEVEL0=[99.90000000000001]
  ## E2 SAME bad-cluster rc=198 (was 111) LASTRC=[111] LEVEL0=[99.90000000000001]

Once WOOLDID_LEVEL0 holds an unrestorable rendering, `macro drop WOOLDID_LEVEL0' (ado:70) is
never reached, so the global is sticky for every subsequent call that fails BEFORE the worker's
level block (ado:574-576) - i.e. every parse-stage failure. For each of those:

  1. the true return code is replaced by 198 (111 -> 198 above), so `if _rc == 111' handling in
     a caller's do-file breaks, and
  2. the whole rc-gated `nobreak' epilogue (ado:72-95) is skipped - the frame restore, the
     emptycells restore, the `est drop ___wooldidfullmodel' cleanup, and the WOOLDID_LASTRC
     update all do not run. LASTRC is left at a value from an older call (111 above).

Both user-visible error messages still print, so it is diagnosable - but the shell's own safety
net, whose whole purpose is B18/B22 cleanup, is disarmed. A call that REACHES the worker's level
block refreshes the global and recovers (s4_level.log @@C2a-@@C2d, all rc 0).

---

## NEGATIVE RESULTS (attacked hard, clean)

1. anticipation(k) is an exact relabelling. anticipation(2) vs `gen ttre2 = ttre - 2':
   bit-identical e(b) (hex equal, rel dev 0.0000e+00) across att, att_it, att_i, att_itime
   (s8_antic.log @@ACMP/@@AWORST); identical again with an event study (esprelength(2)
   espostlength(3), @@ESWORST 0.0), identical again under `esfixedbaseperiod esrelativeto(-1)'
   (s10_mix.log @@FWORST 0.0, same column names, same e(N)=720), and identical again with a
   coarse cohort-level i plus fe(uid) (@@CWORST 0.0). Sanity check that the test could fire:
   |antic2 - antic0| = 8.57e-01.
2. anticipation(-2) is refused informatively, naming the delayed-onset design and the remedy.
3. level() propagation is exact. gteffects: 260 rows, every interval half-width ratio between
   level(95) and level(80) equals invttail(389,.10)/invttail(389,.025) to a max relative
   deviation of 1.5e-15, 0 rows outside 1e-6 (s9_levelprop.log @@GWORST). Main and event-study
   rows scale by the same factor. Column naming follows the documented `p' convention:
   lb_92p5ci, lb_99p91ci, lb_10p01ci, with e(level) matching.
4. CS-2 holds. With `set level 90' and no level() typed, bootstrap(standard) returns
   lb_90ci/ub_90ci, e(level)=90, and genuinely narrower intervals than the same run at 95.
5. Replay with no dataset in memory works (`drop _all' then bare wooldid, rc 0, e(N)=720).
6. Replay of the analytic table is byte-identical once the caution block is set aside
   (RJ: 33 of 33 remaining lines match exactly).
