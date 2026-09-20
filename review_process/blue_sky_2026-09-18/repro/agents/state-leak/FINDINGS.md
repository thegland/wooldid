# state-leak agent -- findings
Target: session state (RNG, globals, the user's data, sort order, back-to-back runs).
Method: STATE DIFFERENCING -- snapshot everything observable, run wooldid, snapshot, diff.

## HARNESS (common.do)
- `wdsnap <tag>` -> `snap_<tag>.txt`: every c-class setting including the FULL 5011-char
  `c(rngstate)` and `c(sortrngstate)`; every variable's storage type, display format,
  value label and variable label; `datasignature`; a planted row-order fingerprint
  (`__wdrowid`, asserted `== _n`); `_dta` characteristics; every global. Plus
  `dirs_<tag>.log`: `matrix dir`, `scalar dir`, `estimates dir`, `frames dir`,
  `mata describe`, `label dir`.
- `wdline <tag>` -- one greppable fact per line (batch logs hard-wrap at 80 columns).
- `wddumpe <tag>` -- e(b), e(V), e(wdidmainresults), e(b_es_att), e(V_es_att),
  e(gteffects) in exact `%21x` hex, for byte-level cross-run and cross-SESSION diffs.
- `wdgen` -- 60 units x 12 years, three treated cohorts + never-treated, carrying value
  labels, variable labels, non-default formats, a subgroup and a weight.
- Harness gotcha worth recording: Stata's own `S_*` globals hold compound-quoted values;
  a snapshot must fetch and sanitise them through Mata or it dies r(198). And
  `set rngstate "<state>"` (quoted) is r(198) -- the token must be unquoted.

## FINDING 1 (medium) -- an outcome missing on every never-treated row is
## misdiagnosed as a design failure, and the remedy the message names returns a
## silently wrong number
repro: repro_missing_control_outcome.do   (also p09_misc.do ARM E, p15b_safetyoff.do)
CLAIM  When the outcome (or i, t, or the weight) is missing on every never-treated
       row, wooldid exits r(459) with "you have no never-treated control group and no
       not-yet-treated control group whose treatment falls outside the sample window /
       To ignore this error, specify safety(off) or modify the sample window manually".
       The refusal is correct; the DIAGNOSIS is not.  Nothing in the message mentions
       missing data, and neither remedy it offers addresses the cause -- the sample
       window is fine, and safety(off) proceeds on a silently different sample and
       estimand (360 of 720 rows, 9 of 12 clusters, identified off not-yet-treated
       comparisons only), reporting an att 48% away from the complete-data answer.
       The header does disclose "360 obs dropped", but no caution says WHY the rows
       went, and none says the never-treated control group is gone.
       PRECISION NOTE: the safety(off) number is not "wrong arithmetic" -- it is the
       estimate for whatever design survives.  The defect is that the user is steered
       to it by a message that mis-names the cause.
OBSERVED  180 never-treated rows in memory, outcome missing on all of them:
       plain          -> r459, message above
       safety(off)    -> rc 0, N clusters = 9, N = 360 (360 dropped),
                         att = .415533924404, se = .0610866413393, p = .000137
                         and NO caution about the absent control group
       outcome intact -> att = .795871565626, se = .100824299164
       dropping the outcome-missing rows by hand -> r459 (the honest answer)
       partial missingness is handled well: 30% of control rows missing -> rc 0,
       e(N)=667, e(N_dropped)=53; one whole calendar year missing -> rc 0, e(N)=660,
       e(N_dropped)=60.
EXPECTED  The refusal is correct -- the realized sample really has no controls -- but
       the message should name the cause the user can act on.  wooldid already knows
       the number: the markout at ado:683 zeroes ifvar for missing y long before the
       control-group screen at ado:952 counts `ttre == . & ifstatement', and the help
       promises missing-value drops surface "in the header's dropped-observation
       count" (sthlp:898, 927, 1411) -- a count only printed at ~ado:1497, after the
       regression, so on this path it is never shown.
SURFACE  safety() / the r459 estimability screen / missing-value handling

## FINDING 2 (low) -- `verbose` leaves a variable in the user's dataset
repro: repro_verbose_leaves_variable.do   (measured in p16_snapdiff.do arm s3)
CLAIM  A `verbose' run adds `_est____wooldidfullmodel' (byte, label "esample() from
       estimates store") to the user's data and does not remove it: c(k) 12 -> 13 and
       the datasignature changes.  A plain run adds nothing (ado:4539 drops the stored
       estimate, and with it the marker variable).  The help documents the retained
       ESTIMATE (sthlp:1132: "This option also stores the underlying regression model
       used by wooldid to ___wooldidfullmodel") but not the
       variable, so a user who saves, exports, or loops `foreach v of varlist _all'
       after a verbose run ships a column they did not create.
SURFACE  verbose

## EVERYTHING ELSE IN THIS LANE CAME BACK CLEAN
Arm-by-arm evidence, with the mechanism each negative rests on, is in NOTES.md.
Headline negatives:
 1. A plain wooldid run consumes NO randomness. c(rngstate) is byte-identical before and
    after (5011 chars compared in full, not a hash). Attribution control: reghdfe and
    gtools called directly also consume none, so the negative is about wooldid.
 2. Back-to-back identical specifications are BIT-IDENTICAL -- on the data wooldid just
    handed back, on regenerated data, and ACROSS SESSIONS (a fresh session's first run
    equals a run made after five failures in another session).
 3. Ten different specifications each followed by the same reference spec: all ten
    reference dumps bit-identical to the in-session baseline. No stale-state effect.
 4. Nothing survives a failed run except `$WOOLDID_LASTRC`, which is the documented
    replay marker. Tested across parse-stage (r198), estimability (r459), ppmlhdfe
    refusal (r411), Poisson separation (r322), a genuine POST-POSTING failure on the
    plain path AND on the deferred bootstrap plot path (r602), and three failures in a
    row inside `capture`. c(level) and c(emptycells) are restored every time; no frame,
    estimate, variable or Mata object is left behind.
 5. WOOLDID_REPMODE cannot leak: the setters wrap it in `nobreak`, the worker consumes it
    READ-ONCE at ado:538, and a leaked value would affect at most one run. Measured absent
    after every completed and every failed engine run.
 6. `cleanmatrices` does not poison the session. The PS-62 covariance splice still happens
    on the cleanmatrices run itself and on the next run, with `set adosize 1000` and with
    `set adosize 10000`; the ereturn dumps are bit-identical to the plain run.
 7. rseed() isolates the engine draws exactly (same results from two different caller
    seeds); `set seed` before an engine run reproduces exactly; e(perm_rngstate) /
    e(boot_rngstate) equal c(rngstate) immediately before the call.
 8. `simultaneousbands` leaves c(rngstate) byte-identical and its critical value does not
    depend on the caller's seed.
 9. sthlp:1325's "fastnose returns the coefficient vector but not the covariance" is TRUE
    as written; known_bydesign s7 item 11's hypothesised divergence does not exist.

## FINDING 3 (low, UNCERTAIN -- attribution not settled) -- a plain run leaves three
## globals in the user's MATA workspace, one of them called `k'
repro: repro_mata_k.do  (measured in p16_snapdiff.do arm s1, `diff dirs_s1_pre.log
       dirs_s1_post.log'; the fuller attribution battery is p21_matanames.do)
CLAIM  One plain `wooldid y unit year ttre, cluster(state) att' on a session whose Mata
       workspace was EMPTY leaves it holding
         8 real scalar __00000L
         8 real scalar __00000M
         8 real scalar k
       and the Stata scalar namespace holding `__gtools_gfile_ghdfeabs = 65' and
       `__gtools_gfile_gregvcov = 65'.  A user with their own Mata variable named `k'
       -- an utterly ordinary name -- would have it silently overwritten.
WHY UNCERTAIN  The `__gtools_' pair is obviously gtools'.  For the Mata trio the
       attribution is NOT settled: wooldid.ado's only top-level Mata assignments are to
       `wooldid_Mtmp' (ado:6689, 8121, 8791, 8795), which are namespaced and do not
       fire on a plain run, so `k' is very likely reghdfe's or gtools'.  The probes
       that would settle it (repro_mata_k.do, p21_matanames.do -- reghdfe alone vs
       gtools alone vs wooldid alone, plus the direct `mata: k = 12345' clobber test)
       were written but never got a Stata slot: the shared 8-instance cap was
       saturated, and at times exceeded, by other agents for the last stretch.
       CODE EVIDENCE POINTS AWAY FROM WOOLDID: every `k' in wooldid's own Mata source
       is a DECLARED function-local (e.g. ado:14737 `string colvector k' inside
       wooldid_bootcellkeym), and an undeclared name inside a Mata function is a local,
       not a global, so wooldid cannot be creating the Mata global `k' that way; its
       only top-level Mata assignments are the namespaced `wooldid_Mtmp'.  The
       `__gtools_' scalars appearing in the same diff make gtools/reghdfe the likely
       source.  If that is confirmed this is NOT a wooldid defect and should be struck.
       Run repro_mata_k.do -- it takes under a minute.
SURFACE  Mata / scalar namespace on every run

## MECHANISMS EACH NEGATIVE RESTS ON (so the author knows what was actually checked)
- **WOOLDID_REPMODE cannot leak**: every setter wraps it in `nobreak` (ado:7672, 7786,
  9065, 9250, 9380, 9450, 9500) and the worker consumes it READ-ONCE at ado:538-539, so
  even a leaked marker would affect at most one run. Measured absent (`GLOBAL|none`)
  after every completed AND every failed engine run in p03/p07/p15.
- **WOOLDID_CLEANUP does not survive an engine run** even though the shell's second
  `wooldid_validate` (ado:312) re-sets it to `keep` under verbose (ado:4947): the first
  nested replication's own shell drops it at ado:65. Measured after verbose+bootstrap,
  verbose+permute and plain verbose: `GLOBAL|none` in all three.
- **`set emptycells drop` (ado:602) and `set coeftabresults on` (ado:552)** rely on
  Stata's program scoping. Tested by EFFECT, not by name: `set emptycells keep` survives
  a run, and with `set coeftabresults off` a post-wooldid `regress` still leaves no
  `r(table)`.
- **level()**: `set level 80` in the session is honoured, `level(99)` applies to that
  call only, c(level) returns to 80, and on the engine paths level(80) produces
  `lb_80ci/ub_80ci` columns and `e(boot_level)=80` with an interval strictly inside the
  95% one. Restored even when the run dies AFTER `set level` (r459 and r602 cases).
- **simultaneousbands** saves and restores `c(rngstate)` around its per-family
  `set seed` loop (ado:10888-10963): measured byte-identical, and the critical value is
  the same from two different caller seeds.
- **`if` is materialised once** at ado:657 on the user's row order, so `if _n <= 360`
  and `if flagvar == 1` give bit-identical e(). (The replication path re-evaluates the
  raw string per replication at ado:9057 -- see p18.)

## MORE HEADLINE NEGATIVES (added after the later probes)
10. **The full before/after SNAPSHOT DIFF is the strongest result here.** On a failing
    run (r459) the entire session diff is ONE line: `[GLOBAL] WOOLDID_LASTRC = 459`.
    Not one variable, format, label, characteristic, setting, matrix, estimate, frame
    or Mata object changes. On a plain successful run the diff is two lines, and BOTH
    are attributed elsewhere: `c(matastrict)` off->on is reghdfe's (one bare
    `reghdfe y post, absorb(unit year)` flips it in a session where wooldid never ran;
    `matastrict` occurs nowhere in wooldid.ado), and `c(sortrngstate)` advances as any
    sorting command's does while the user's actual row order is restored intact.
    On a permute run the diff is the RNG state (documented) and nothing else -- no
    globals at all. On a bootstrap+plots run it adds only `T_gm_fix_span`, which is set
    by Stata's own /Applications/Stata/ado/base/g/grid.class.
11. **`if e(sample)` is safe under the engines.** ado:9057 does hand the raw `if`
    string to every replication, but the replications run on a frame already reduced to
    the `smp` rows computed once at ado:8911, so `wooldid ... if e(sample)` and
    `wooldid ... if esflag==1` give bit-identical point estimates, bootstrap SEs and
    intervals. Same for `if _n <= 360` under both bootstrap() and permute().
12. **`set dp comma` is harmless.** The ereturn dumps are bit-identical to the default
    session for a run with esprelength(3) espostlength(3) subgroups(sg) jointtests, and
    wooldid's own cautions correctly print "0,095".
13. **The 2026-09-05 verbose-bootstrap regression has not returned.** After
    `att_i verbose bootstrap(standard) reps(20)` and after `verbose permute(shuffle)`,
    `___wooldidfullmodel` restores as the OBSERVED model (e(N)=720, e(cmd)=reghdfe) and
    `___wooldidrepmodel` is gone.
