# state-leak probe -- working notes

METHOD: state differencing. `wdsnap <tag>` (in common.do) writes snap_<tag>.txt
(c-class incl. full c(rngstate) and c(sortrngstate), every variable's type/format/
value label/variable label, datasignature, a planted row-order fingerprint
__wdrowid, _dta chars, every global) plus dirs_<tag>.log (matrix/scalar/estimates/
frames dir, mata describe, label dir).  `wdline <tag>` prints a compact greppable
one-fact-per-line state line.  `wddumpe <tag>` dumps e(b)/e(V)/e(wdidmainresults)/
e(b_es_att)/e(V_es_att)/e(gteffects) in exact %21x hex.

HARNESS FACTS ESTABLISHED
- pinned build resolves (`which wooldid` -> .../bluesky/build/wooldid.ado, v1.4).
- c(sortrngstate) EXISTS on this Stata and is snapshotable.
- default c(emptycells) is "keep"; worker sets "drop" (ado:602) and relies on
  program scoping.
- Stata's own S_* globals contain bare backticks; the snapshot must use
  `macval()' or it dies r(198).  (fixed in common.do)

QUEUE: ./queue.sh runs one do-file at a time, waiting for a free slot under the
project cap of 8 Stata instances.

ARMS
 p01_plain      reference: plain run alone in a fresh session, full snapshot diff
 p05_cleanmata  cleanmatrices vs plain: e(V) splice, Mata kernels, session poisoning
 p02_bb         attribution control (reghdfe/gtools alone); back-to-back same spec;
                RNG bleed on a plain run; coeftabresults / emptycells scoping
 p03_fail       parse fail / estimation fail / poisson fail / post-post fail; globals;
                level restore; repeated failures in a capture loop
 p06_permrng    permute(): seeding, e(perm_rngstate) replay, rseed determinism, permsave
 p04_diffspec   ten different specs each followed by the reference spec on identical data
 p07_bootrng    bootstrap(): same battery
 p08_frames_if  non-default frame (success + failure), plots/files/graphs left behind,
                replay with no data, `if _n<=' vs a flag variable, `if runiform()<'

STATUS: see queue.out

## RESULTS SO FAR (all measured on the pinned build, 60 units x 12 years)

### p02_bb  -- CLEAN
- reghdfe + gtools alone consume NO randomness (attribution control).
- a PLAIN wooldid run consumes NO randomness: c(rngstate) byte-identical before/after.
- back-to-back identical spec: e(b)/e(V)/e(wdidmainresults) BIT-IDENTICAL across
  (i) the same data wooldid just handed back, (ii) regenerated data, and
  (iii) a DIFFERENT SESSION (p03's run after five failures) -- cross-session bit equality.
- `set coeftabresults off' survives a wooldid run (program scoping holds).
- `set emptycells keep' survives a wooldid run (ado:602's reliance on program scoping holds).
- variable count, types, formats, labels, row order, c(level), frame: all unchanged.

### p05_cleanmata -- CLEAN (CONFIRMS known_bydesign s0)
- cleanmatrices run vs plain run: ereturn dumps BIT-IDENTICAL, e(wdid_vcvstate)=
  analytic-blocks and e(wdid_vblockn)=15 in both. The PS-62 covariance splice still
  happens under cleanmatrices.
- a plain run AFTER a cleanmatrices run: identical again -> the session's Mata kernels
  are back. Also tested with `set adosize 10000' (so the ado stays resident): still
  identical. gteffects after cleanmatrices: e(gteffects) present.
- cleanmatrices does destroy the user's own matrices and Mata objects (rc 111 both) --
  documented at sthlp:1114.

### p03_fail -- CLEAN
- parse fail r198, estimability fail r459, ppmlhdfe-parse fail r411, separation r322,
  three failures in a row inside `capture': after every one the ONLY residue is
  $WOOLDID_LASTRC (documented, NFRPL-01/02). No WOOLDID_CLEANUP, no WOOLDID_LEVEL0,
  no WOOLDID_REPMODE.
- level() restored after a failure that happens AFTER `set level' (B_LEVEL_RESTORED=1).
- emptycells restored on every failure path.
- _ppmlhdfe_d is NOT left behind on any poisson failure path tested.
- a clean run after the failures reproduces the reference numbers BIT-IDENTICALLY.

### p08_frames_if -- CLEAN
- run from a non-default frame: current frame is still `work' after success AND after
  an r459 failure.
- replay with no dataset in memory works (sthlp:1459 claim holds).
- replay after the data changed underneath works; replay after a `regress' falls
  through to the worker and errors r198 (documented design).
- `if _n <= 360' and `if flagvar == 1' give BIT-IDENTICAL e() -- the if expression is
  evaluated on the user's row order (ado:657), not after an internal sort.
- makeplots/saveplots leave one named graph (att_fullmodel) in memory and write
  <name>_ES.dta + <name>_ES_att.gph; graphs are created with name(..., replace)
  (ado:2291 comment), so re-runs overwrite.

### p09_misc -- CLEAN except one message
- simultaneousbands: c(rngstate) BYTE-IDENTICAL before/after (ado:10888 claim holds),
  and e(esband_crit_att) is the same from two different caller seeds.
- e(sample) sane: e(N)=count if e(sample)=720 under both fastcollapse and nofastcollapse.
- level(): `set level 80' honoured, level(99) used for that call only, c(level) back to 80.
- FINDING CANDIDATE: outcome missing on every never-treated row -> r459 "you have no
  never-treated control group ... modify the sample window manually", with no mention
  of missingness. Partial missingness is handled well (rc 0, counted in e(N_dropped)).

### p06_permrng -- CLEAN
- e(perm_rngstate) == c(rngstate) immediately before the call (doc holds).
- two identical `set seed 777' + permute runs reproduce EXACTLY (p, interval, estimate).
- rseed(4242) from `set seed 111' and from `set seed 999' give IDENTICAL results ->
  rseed() isolates the draws properly.
- the engine advances the session stream (documented, sthlp:542).
- permsave() file carries _dta[wooldid_engine/seed/rngstate] correctly.
- the observed att under permute() differs from the plain run's att by 4 ulps
  (~9e-16), which is the documented fastnose/fast-mode path, far inside rtol 1e-10.
- HARNESS NOTE: `set rngstate "<state>"' (quoted) is r(198); the token must be unquoted.

### p07_bootrng -- CLEAN
- e(boot_rngstate) == c(rngstate) immediately before the call.
- two identical `set seed 777' + bootstrap runs: every printed number BIT-IDENTICAL.
- rseed(4242) from `set seed 111' vs `set seed 999': BIT-IDENTICAL -> proper isolation.
- after a bootstrap run (and after a bayes run with bootsave), a plain analytic run
  reproduces the clean-session reference BIT-IDENTICALLY; e(boot_method)/e(perm_method)
  are empty on it (no stale engine markers).
- no leftover globals, no _dta chars, row order intact, c(level)/c(emptycells) restored.
- NOTE: two bootstrap runs back to back with no reseeding had 39 and 38 of 40
  completed reps (different draws, as expected); the OBSERVED estimate was identical
  in every arm.

### p04_diffspec -- CLEAN
Ten different specifications (verbose, subgroups, poisson, nofastcollapse, impvalues,
fastnose, level(80), customvce, unconditionalse[r198 by design], contreat[r198 by
design]) each followed by the SAME reference spec on regenerated data: all ten
reference dumps BIT-IDENTICAL to the in-session baseline.
- nofastcollapse (the path with no preserve/restore) leaves row order intact too.
- level(80) run: e(level)=80, c(level) back to 95 afterwards.
- LEAD FOR p17: under `fastnose' with NO ES window, e(b_es_att) is absent and
  e(numpre_att) is missing -- but that is confounded (no event study at all).
  p17 re-tests with esprelength(3) espostlength(3).

### p15_edge (partial) -- CLEAN so far
- A1/A2: REGRESSION CHECK for the 2026-09-05 fix ("a verbose bootstrap run left the
  LAST replication's model under ___wooldidfullmodel"). After `att_i verbose
  bootstrap(standard) reps(20)' and after `verbose permute(shuffle) permreps(20)',
  ___wooldidfullmodel restores with e(N)=720 and e(cmd)=reghdfe -- the OBSERVED
  model, not a replication (the physical re-id leg's replications have a different
  N by construction). ___wooldidrepmodel is gone in both. The fix HOLDS.
- B1/B2: a genuine POST-POSTING failure (saveplots to an existing .gph with no
  `replace' -> r602 at graph save, after ereturn post and after `set level 80'):
  c(level) is back to 95, $WOOLDID_LASTRC=602, replay refused r198.
- B3: the same on the DEFERRED (bootstrap) plot path (ado:259-300): r602,
  c(level)=95, c(frame)=default, no ___wooldidfullmodel left, replay refused,
  and the next good run clears the marker.
- C1: bootstrap + level(80): e(level)=80, e(boot_level)=80, table columns named
  lb_80ci/ub_80ci, c(level) back to 95 afterwards.

## HARNESS INCIDENT (recorded so a verifier is not confused)
Two Stata processes briefly ran p15_edge.do at once (an orphaned child of a queue
script I had killed, plus a fresh one), both writing p15_edge.log. Both were killed
and p15 was re-run alone. Any p15_edge.log content timestamped before 03:03:45 is
from the discarded pair. I also learned not to `pkill -f queue...` (the pattern can
match other agents' scripts); everything since is killed by explicit PID.

### p17_fastnose_es -- CLEAN; SETTLES known_bydesign s7 item 11 IN THE NEGATIVE
All arms with esprelength(3) espostlength(3) so "no event study" cannot confound:
  reference : e(b_es_att)=1  e(V_es_att)=1  e(numpre_att)=3  vcvstate=analytic-blocks
  fastnose  : e(b_es_att)=1  e(V_es_att)=0  e(numpre_att)=3  vcvstate=absent
  timetrends: e(b_es_att)=1  e(V_es_att)=0  e(numpre_att)=3  vcvstate=absent
  impvalues : e(b_es_att)=0  e(V_es_att)=0  e(numpre_att)=.  vcvstate=analytic-diagonal
  permute   : e(b_es_att)=0  e(V_es_att)=0  e(numpre_att)=.  vcvstate=absent
  fastnose + att_it also posts e(b_es_att_it) and withholds e(V_es_att_it).
=> sthlp:1325 ("when NO horizon has an identified variance -- fastnose, typically
   timetrends ... -- the covariance is not returned but the coefficient vector still
   is") is TRUE AS WRITTEN for fastnose and for timetrends.  The doc-vs-code
   divergence s7 item 11 hypothesised does NOT exist on the pinned build.
   impvalues' absence is separately documented at sthlp:284 ("Every IM margins call
   runs with nose, so ... e(b_es_), e(V_es_) and e(numpre_) are not returned"), and
   permute()'s at sthlp:448, with sthlp:500 saying permanalyticse restores them.
   My earlier p04 X6 reading (b_es absent under fastnose) was the confound of running
   fastnose with no ES window at all; corrected here.

### p15b_safetyoff -- FINDING 1 SHARPENED
  E0 never-treated rows present in memory: 180
  E1 `wooldid y unit year ttre, cluster(state) att'            -> r459, message names
     the DESIGN ("no never-treated control group ... modify the sample window")
  E2 `... safety(off)'  (the FIRST remedy the message offers)  -> rc 0,
     N clusters = 9, N = 360 (360 obs dropped), att = .415533924404, se = .0610866,
     p = .000137 -- and NO caution of any kind about the absent control group
  E3 the same data with the outcome intact                     -> att = .795871565626
     (so the safety(off) number is 48% below the truth, reported at p = .0001)
  E4 dropping the outcome-missing rows by hand                 -> r459 (correct)
  D3 plain verbose run leaves NO globals (the predicted WOOLDID_CLEANUP=keep leftover
     from ado:312 + ado:4947 does NOT survive: the first nested replication's shell
     drops it at ado:65).  Confirmed absent after verbose+bootstrap and verbose+permute
     too (p15 A1/A2: WDL|A1_post|GLOBAL|none, WDL|A2_post|GLOBAL|none).

### p16_snapdiff -- THE FULL BEFORE/AFTER SNAPSHOT DIFF (the method's deliverable)
Seven arms, each snapshotted immediately before and after ONE wooldid call.
`diff snap_<arm>_pre.txt snap_<arm>_post.txt` in full:

 s1 plain             : c(matastrict) off -> on ; c(sortrngstate) advanced. NOTHING ELSE.
 s2 nofastcollapse    : IDENTICAL (matastrict already on; this path does not even
                        touch the sort stream)
 s3 verbose           : c(k) 12 -> 13 and one new variable
                        `_est____wooldidfullmodel` (byte, label "esample() from
                        estimates store"); datasignature changes accordingly;
                        c(sortrngstate) advanced. No globals, no other change.
 s4 failing run (459) : ONE line added -- [GLOBAL] WOOLDID_LASTRC = 459. Nothing else
                        in the entire snapshot: not one variable, format, label,
                        characteristic, setting or estimate.
 s5 permute           : c(seed)/c(rngstate) advanced (documented), c(sortrngstate)
                        advanced. NOTHING ELSE -- no globals at all.
 s6 bootstrap+plots   : c(seed)/c(rngstate) and c(sortrngstate) advanced, plus one
                        new global `T_gm_fix_span = 0'.
 s7 poisson (success) : c(sortrngstate) only. _ppmlhdfe_d is gone, c(k) unchanged --
                        the documented cleanup works.

ATTRIBUTION (so none of this is filed as wooldid's that is not):
 - `T_gm_fix_span` is set by STATA's own graph code: it appears verbatim in
   /Applications/Stata/ado/base/g/grid.class and nowhere in wooldid.ado. Any command
   that draws a twoway graph leaves it. NOT a wooldid defect.
 - `matastrict` does not appear anywhere in wooldid.ado. It appears in reghdfe.mata,
   reghdfe3.ado, ftools.mata, _gtools_internal.mata and hdfe.ado -- wooldid's hard
   dependencies. p20_attrib.do measures whether reghdfe / gtools / margins alone flip
   it, in sessions where wooldid never runs.
 - `sortrngstate`: any `sort` consumes it; p20 measures a bare `sort` too.

### p20_attrib -- ATTRIBUTION SETTLED
  00 fresh session            matastrict=off  sortrng=20260916XZA0 (the seed itself)
  02 after a bare `sort'      matastrict=off  sortrng UNCHANGED
  03 before reghdfe           matastrict=off  sortrng=413035572XZA (wdgen's own sorts)
  04 AFTER reghdfe alone      matastrict=ON   sortrng unchanged
  06 after gtools alone       matastrict=on   sortrng unchanged
  08 after margins alone      matastrict=on   sortrng unchanged
  09/10 around wooldid        matastrict=on both sides; sortrng 72497204 -> 1944706868
=> c(matastrict) off->on is REGHDFE's, not wooldid's: `matastrict' does not occur
   anywhere in wooldid.ado, and one bare `reghdfe y post, absorb(unit year)' flips it
   in a session where wooldid never ran.  (And `mata: zz = 1 + 1' still works with it
   on -- matastrict governs function bodies, not interactive assignment -- so the
   user-visible consequence is small either way.)  NOT FILED.
=> c(sortrngstate): wooldid's internal sorts advance it, as any sorting command's do.
   The user's actual ROW ORDER is restored (roworder_intact=1 on every arm) and
   c(sortlist) is unchanged.  Recorded, not filed.
=> `T_gm_fix_span' after a plots run is STATA's: it is in
   /Applications/Stata/ado/base/g/grid.class and in no wooldid file.  NOT FILED.

### p14_replineparity -- ABANDONED for cost, replaced by p14b
16 spec pairs x 60 reps was running at ~4 minutes per pair under the shared load.
Spec 1 (the control) landed before it was killed:
  BIAS|01_plain|main_att| b = .792024015737  bmean = .791096499033
i.e. the replication mean sits essentially on the observed estimate, which is what
parity looks like.  p14b keeps R=25 and the eight specs most at risk of being dropped
from the hand-built replication option string at ado:318-363 (plain, esfixedbaseperiod,
anticipation, oostreatedcontrols, customweightmultiplier, interactivecontrols,
semielasticity, poisson).  Reading that string in the source, the only conditional
re-emissions that could silently drop an estimand-shaping option are
`if "`semielasticity'" == "eydx"' (ado:350) and
`if "`contreatelasticitytype'" != ""' (ado:351); everything else is a plain
`if "`opt'" != ""'.  validate normalises the semielasticity FLAG to the margins
operator name (dydx/eydx, used at ado:2030 etc.), so the == "eydx" test is the right
one -- that is a code-reading result, and p14b is the measurement.

## RESOURCE NOTE
From ~03:10 the shared 8-instance Stata cap was saturated (and at times exceeded, 9)
by other agents, so the tail of my queue (p18 if/e(sample), p19 session settings,
p12 rngstate-replay handle, p10 estimates ceiling, p11 data environment, p13 prefixes,
p14b trimmed rep-line parity, and the two repro runs) sat waiting for a slot.
Whatever of those has a .log beside this file ran; whatever does not, did not.

## UNTESTED / NOT REACHED (honest list -- updated at the end of the run)
Struck from this list because they DID run: p18 (A arms) and p18b (D/A/C arms),
p19 (baseline + dp comma), p12 (permute half), p16, p17, p20, and both repro files.
Still not reached:
- p18: whether the replication path's re-evaluation of the RAW `if' string (ado:9057,
  and the same shape at 7564) breaks `if e(sample)' -- the mechanism is confirmed in
  the source, the measurement is the missing piece.
- p19: `set dp comma', `set varabbrev off', `set type double', non-default
  cformat/pformat/sformat, a different sortseed, and a 30-run accumulation loop.
- p12: the documented "replay by `set rngstate' to e(perm_rngstate)/e(boot_rngstate)"
  handle (p06/p07 covered the seed-determinism half of the claim, not this half).
- p10: Stata's 300-stored-estimates ceiling under wooldid's own `est sto'; the
  estimates store/restore and estimates save/use round trips and replay across them.
- p11: xtset/tsset survival, _dta and variable notes, a pre-existing _ppmlhdfe_d on a
  non-poisson run, repeated makeplots, wooldid inside a user program that has issued
  `preserve', string variables in the data.
- p13: Stata's own bootstrap/statsby prefixes around wooldid; `svy:' refusal hygiene.
- p14b: replication-line parity for esfixedbaseperiod / anticipation /
  oostreatedcontrols / customweightmultiplier / interactivecontrols / semielasticity /
  poisson (spec 1, the control, measured clean before the full version was killed).
- A post-posting failure on the CONTINUOUS-TREATMENT deferred plot path
  (wooldid_cfxplots, ado:263-272) -- only the event-study twin was exercised.

## HARNESS HYGIENE (second incident, recorded for the coordinator)
My earlier `nohup ./queueN.sh ...' launchers did not always die when I thought they
had, so two of them (qfinal.sh and queue2.sh) kept re-launching files from their own
stale lists -- duplicating p14 and p10 and holding slots that blocked my own queue and
everyone else's. All of mine were found by PPID and killed at 03:12; the only launcher
scripts left running at that point belonged to other agents. Anyone reading a
p14_replineparity.log or p10_estimates.log from before 03:12 should treat it as
possibly interleaved output from two processes. I also used `pkill -f queue...` early
on, whose pattern can match other agents' scripts; everything after that was killed by
explicit PID.

### p18/p18b -- the `if'-expression surface under the replication engines
MECHANISM (source): the replication option line carries the user's RAW `if' STRING
(`ifr()'), and ado:9057 re-applies it to every nested call --
  `if `"`ifr'"' != "" local wtif0 `wtif0' if (`ifr')'
-- inside a frame that has already been resampled.  The same shape is at ado:7564 for
permute().  So a non-deterministic or position-dependent `if' is evaluated afresh in
every replication, on different rows than the observed run saw.
MEASURED so far (p18_ifreps.log, arms A1/A2):
  `if keepflag == 1'  + bootstrap(standard) reps(60) rseed(11)
      -> N=360, est=.833757148598, se=.168992017756, completed=57, failed=3
  `if _n <= 360'      + the same
      -> N=360, est=.833757148598, se=.168992017756, completed=57, failed=3
  IDENTICAL, so `_n' is NOT a hazard on this path: ado:8911 materialises `smp' once and
  the replication frame is a `keep if smp' copy in the ORIGINAL row order, in which
  `_n <= 360' selects every surviving row anyway.
p18b re-runs the `if e(sample)' pair (plain and bootstrap), the permute twin of A, and
the two-different-engines-back-to-back arm, with the wdse type-mismatch fixed (reading
e(perm_n) as a scalar is r(109) -- it is not one).

### p19_setstate (partial)
  BASE_RC=0, BASE_COLNAMES = estimate se t p lb_95ci ub_95ci relyr
  `set dp comma' DOES take effect even though c(dp) still reads "period" (a Stata
  reporting quirk): wooldid's own caution then prints "the smallest attainable
  per-coefficient permutation p-value is 0,095" with a comma.  And the substantive
  result is CLEAN -- `diff ereturn_set_base.txt ereturn_set_dpcomma.txt' is IDENTICAL
  for a run with esprelength(3) espostlength(3) subgroups(sg) jointtests, and
  e(wdidmainresults)'s column names and e(wdid_sglist) are unchanged.  So wooldid is
  decimal-separator invariant on that surface.  The engine arms under dp comma did not
  finish: I killed p19's Stata myself via TaskStop while re-prioritising the queue.
  The remaining arms (varabbrev off, type double, non-default cformat/pformat/sformat,
  a different sortseed, a 30-run accumulation loop and a 10-run loop with plots and
  subgroups) were still running when the queue was re-prioritised.

### p18b -- `if e(sample)' IS CLEAN (the strongest negative in this lane)
`regress y post if state <= 8' first (480 rows), then wooldid restricted two ways:
  D1 `if e(sample)'  plain      N=480 est=.808211020085 se=.058119209647
  D2 `if esflag==1'  plain      N=480 est=.808211020085 se=.058119209647   IDENTICAL
  D3 `if e(sample)'  bootstrap  N=480 est=.808211020085 se=.154166131168
                                 lb=.503987359545 ub=1.03328673672
  D4 `if esflag==1'  bootstrap  N=480 est=.808211020085 se=.154166131168
                                 lb=.503987359545 ub=1.03328673672          IDENTICAL
So although ado:9057 does hand the RAW `if' string to every replication, the
replications run on a frame that has already been reduced to the `smp' rows computed
ONCE at ado:8911, and re-evaluating `e(sample)' there selects the same rows.  The
"estimate on the sample my first stage used" pattern -- the most realistic version of
this hazard -- produces bit-identical inference.

### p18b A arms -- `_n' under permute() is clean too
  A3 `if keepflag == 1' + permute(shuffle) permreps(40) rseed(11)
      N=360 est=.833757148598 p=.0487804878049 lb=.193289575685 ub=1.29766584703
  A4 `if _n <= 360'     + the same
      N=360 est=.833757148598 p=.0487804878049 lb=.193289575685 ub=1.29766584703
  IDENTICAL, as under bootstrap.

### p18b C arms -- engines back to back
  C1 permute(shuffle)     -> e(perm_method)=shuffle, e(boot_method) empty
  C2 bootstrap(standard)  -> e(perm_method) EMPTY, e(boot_method)=standard,
                             e(perm_p_lo) ABSENT (no stale permutation surface)
  C3 permute(shuffle) again -> perm_method=shuffle, boot_method empty, and
                               e(bootdist) ABSENT, e(V_boot) ABSENT (no stale
                               bootstrap surface survives into the permute run)
  C4 plain analytic         -> both engine markers empty, e(perm_p_lo) ABSENT,
                               e(V_boot) ABSENT
  So alternating engines back to back leaves no stale inference surface in e() in
  either direction (p07 ARM G measured the bootstrap->plain direction bit-identical
  to the clean-session reference as well).

### BOTH REPRO FILES RUN CLEAN (verified end to end)
repro_missing_control_outcome.log:
  LEG1 rc 0  e(N)=720  att=.251427121182  se=.100824299164
  LEG2 rc 459 (the misdiagnosing message)
  LEG3 rc 0  e(N)=360  e(N_dropped)=360  att=.182200591071  se=.0610866413393
             p=.0175291582039, no caution
  LEG4 rc 459
repro_verbose_leaves_variable.log:
  PLAIN_LEFT_A_VARIABLE=0   VERBOSE_LEFT_A_VARIABLE=1
Both are self-contained (they build their own panel; no dependence on common.do).

### p16 dirs-diff (matrix/scalar/estimates/frames/mata/label), arm by arm
 s1 plain      : + Stata scalars `__gtools_gfile_ghdfeabs' and `__gtools_gfile_gregvcov'
                 (both = 65); + Mata globals `__00000L', `__00000M' and `k'
 s2 nofastcol. : nothing at all
 s3 verbose    : + the stored estimate ___wooldidfullmodel (reghdfe, y, 272 params,
                 "HDFE Linear regression") and the frame goes 720 x 12 -> 720 x 13
 s4 failing    : nothing at all
 s5 permute    : + Mata globals `__000017', `__000018'
 s6 boot+plots : + `__00000N', `__000015', `__000016'; `__000018' removed
 s7 poisson    : `__00000N' removed
The `__gtools_' scalars are self-attributing.  p21_matanames.do attributes the Mata
`k' / `__00000x' globals (reghdfe alone vs gtools alone vs wooldid alone in three
separate sessions) and asks the consequence question directly: with `mata: k = 12345'
set beforehand, does a wooldid run clobber it?  It also checks whether wooldid's own
top-level Mata global `wooldid_Mtmp' (ado:6689, 8121, 8791, 8795 -- the only bare
Mata assignments in the file) survives a `gteffects bootstrap()' run.

### p12_rngreplay -- the documented replay handle WORKS
sthlp:542 promises "a run without rseed() can still be replayed exactly by set rngstate
to that value before repeating the call".  Measured for permute():
  P1 permute(shuffle) permreps(50), no rseed  -> e(perm_rngstate), 5011 chars
  P2 `set rngstate <that>' then the identical call
  -> every printed number BIT-IDENTICAL (estimate, p_perm, lb, ub).
  P3: the same with the state QUOTED is r(198) -- `set rngstate' takes an unquoted
  token.  The help's wording does not show quotes, so this is a user trap rather than
  a defect, but it is worth a line in the help.
  The BOOTSTRAP half of p12 (B1/B2) was cut off by the 10-minute per-file deadline
  right after `set rngstate' returned rc 0; B2's run did not finish, so only the
  permute half of the sthlp:542 claim is measured.  p06/p07 separately measured the
  other half of the reproducibility promise (rseed determinism) for BOTH engines.
