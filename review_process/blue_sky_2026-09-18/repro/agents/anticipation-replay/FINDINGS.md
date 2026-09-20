# FINDINGS -- anticipation() / replay / level()   (method: EXACT EQUIVALENCE TESTING)
Workdir: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay

===============================================================================
AR-4  (HIGHEST VALUE)  level() refuses 1152 of the 9000 legal two-decimal levels,
       refuses them even when the user never typed level(), and on that path
       skips wooldid's own exit epilogue -- leaving the session level changed and
       the replay-staleness guard UNARMED.
Status CONFIRMED.  Repro: repro_AR4.do , p14_guardgap.do , p15_restorefail.do ,
                          p10_lvlbound.do , p11_map.do , p12_sessionlvl.do
-------------------------------------------------------------------------------
Mechanism. `syntax , Level(real 0)` (and `global WOOLDID_LEVEL0 = c(level)`) render
the double with 16 significant digits, so `level(99.9)` reaches wooldid_worker
ado:576 `set level `level'' as 99.90000000000001 and Stata refuses:
   "level() can have at most two digits after the decimal point"   r(198)
The refused set is {L : %.16g(L) != %.15g(L)}.  That rule is PREDICTED by the mechanism
(computed in Python: 1152 of the 9000 legal two-decimal levels, 12.8%, all at or above
64.01) and VALIDATED in Stata on 18 values -- 6/6 predicted failures failed, 12/12
predicted passes passed -- plus the directly observed rendering $WOOLDID_LEVEL0 =
99.90000000000001.  Among CONVENTIONAL choices only 99.9 and 99.99 are hit; that pair is
the realism claim, not the 12.8%.  99.99 is the upper bound sthlp:540 itself names.
Measured: `set level L` accepts 99.9 / 99.99 / 98.9 / 79.9 / 99.98 / 92.5 (rc 0 each);
`wooldid ..., level(L)` refuses 99.9, 99.99, 98.9, 79.9, 99.93, 99.96 and accepts
99.91, 99.92, 99.94, 99.95, 99.97, 99.98, 99.5, 92.5, 10.01, 88.8, 44.4, 12.3.
The message is FALSE about the value typed (99.9 has ONE decimal).

It fires with NO level() option.  wooldid_validate ado:5507 resolves an unspecified
level() to c(level) and c_locals it, so the worker's `if `level' != 0' is then true
and `set level` runs on the session value.  With `set level 99.9` in a profile.do,
a plain `wooldid y i t ttre, att robust` dies r(198) naming an option the user never
used (p12: 99.9, 99.99, 98.9 all r(198); 99.98, 92.5 fine).

Worse, a SUCCESSFUL run is reported as a failure (p15).  Session level 99.9,
`wooldid y i t ttre, att robust espre(2) espost(2) level(95)`:
   the full table prints, e(N)=500, e(level)=95, e(cmd)=wooldid  -- the run worked
   rc = 198                                            -- a do-file aborts here
   c(level) AFTER  = 95, not 99.9                      -- caller's level silently changed
   $WOOLDID_LEVEL0 = 99.90000000000001  (left set)     -- should be dropped
   $WOOLDID_CLEANUP = fullmodel         (left set)     -- epilogue never ran
   $WOOLDID_LASTRC  = <empty>                          -- guard NOT armed
Because the restore at ado:68-71 sits BEFORE the `nobreak' epilogue and is not
capture'd, the failing `set level $WOOLDID_LEVEL0' aborts the shell and the epilogue is
skipped.  What that costs depends on where the worker stood: in p15 (worker SUCCEEDED)
the losses are the level restore itself, the two global drops, and a false rc; `est dir'
showed nothing left behind, so do not read this as a leaked ___wooldidfullmodel.  In p14
(worker FAILED) the loss that matters is WOOLDID_LASTRC never being armed.

Contradicted guarantees (coordinator note: contradiction is first class):
 * sthlp:540  "between 10 and 99.99"  -- 99.99 is refused.
 * ado:571-573 "the shell restores the caller's value on EVERY exit path
   (WOOLDID_LEVEL0)"  -- it does not; the caller is left at the wrong level.
 * sthlp:1460-61 "Either way replay refuses with an informative error until a
   wooldid call completes successfully."  MEASURED OTHERWISE (p14): good run at 95;
   `set level 99.9'; the next wooldid FAILS r(198); $WOOLDID_LASTRC stays empty; a
   bare `wooldid' then REPLAYS the pre-failure results, rc 0, e(N)=500, e(level)=95,
   with no refusal and no warning -- precisely the misleading redisplay the ruled
   NFRPL-01/02 marker exists to prevent.
Realistic: yes -- `set level 99.9' in a profile.do, or reaching for the help's own
99.99 bound.  Severity HIGH.  Confidence CERTAIN.
Fix shape: stop handing `set level' a macro rendered from a real -- format the value to
two decimals as a STRING first (e.g. `set level `: display %9.2f `level''').  NOTE: an
`=round(`level',.01)' assignment does NOT help: round() returns the same double and the
`=' assignment renders it the same way (that is exactly how `global WOOLDID_LEVEL0 =
c(level)' produced 99.90000000000001).  The restore at ado:69 should also be capture'd.

===============================================================================
AR-1  replay drops the cautions that disown the numbers it reprints
Status CONFIRMED.  Repro: repro_AR1.do ; p3 r01/r02/r08/r10 ; p5 h01/h02 ; p7 r03
-------------------------------------------------------------------------------
sthlp:1458 -- replay redisplays "byte-for-byte as originally shown" and "does not
repeat estimation-time Caution messages (WHICH DESCRIBE THE ESTIMATION RUN, NOT THE
RESULTS)".  The parenthetical is false.  Measured losses, each describing the very
numbers replay reprints:
 1. jointtests at cluster() >= i: "...the joint test's covariance is degenerate and
    its p-value SHOULD NOT BE READ AS A CALIBRATED TEST."                    (r01)
    -> this is the runtime half of the ruled remedy PS-64/BS-9, recorded in
       known_bydesign S1 as "Computed anyway WITH A NAMED CAUTION".  The ruling does
       not survive the documented redisplay path.
 2. permute(): "with permreps(25) the smallest attainable per-coefficient p-value is
    0.077 ... a p-value printed at that floor is THE RESOLUTION LIMIT RATHER THAN
    EVIDENCE."                                                                (r08)
 3. bootstrap(): "with reps(25) the coefficient p-values fall on a grid of 0.080 ...
    the only value below .05 this run can produce is exactly 0 - WHICH RECORDS THAT
    NO COMPLETED REPLICATION REACHED THE FAR SIDE OF THE NULL."         (r10, h01)
 4. contreat() with a cohort-constant dose: "the within-cell average marginal effects
    (ame) ... are NOT IDENTIFIED in this design - the regression's ame entries are
    structural zeros, not estimates - so wooldid is reporting the ame results as
    missing", and "standard errors were not identified for 6 of 12 reported
    estimates."  On replay the table shows the blanks with NO explanation.   (p7 r03)
    -> this is the runtime half of AMESE R1 (known_bydesign S1, shipped 7054878).
 5. permute() at cluster()==i: the exchangeability / permtstat caution.       (r08)
 6. adversarial without esfixedbaseperiod: "the fitted line is reported without a
    p-value."                                                                (r01)
 7. histogramcohorteffects + an engine (shell ado:299/511): "histogram outputs are
    drawn from the ANALYTIC covariance ...; the reported tables ... are
    bootstrap-based / permutation-based."                                (h01, h02)
Omission, not policy: wooldid_display DOES replay results-describing cautions when a
state macro was stored for it -- the subgroup withheld-pre-row trio (e(wdid_sgnopre)
/wdid_sgnoprewhy/wdid_sgnoprerow, ado:4431-4437) and the Roth rank-deficiency and
interrupted-pass cautions (ado:11447-11477).  21 Caution sites live in
wooldid_worker, 15 in wooldid_validate, 5 in the `wooldid' shell, 8 in
wooldid_display.
Position: live prints the block at the TOP, before the header, not beside the number.
Everything from the header down replays byte-identically (r01: live 92, replay 72
lines; the 20 missing lines are exactly the caution block).
Realistic: yes.  Severity MEDIUM.  Confidence CERTAIN.

===============================================================================
AR-3  the replay guard is self-poisoning, and its message misattributes the failure
Status CONFIRMED.  Repro: repro_AR3.do ; p5 block C02
-------------------------------------------------------------------------------
A user who has run wooldid, then some other estimation command, and types a bare
`wooldid' to look at the DiD table again gets
   "Error: command line input should feature all of y i t ttre; too few variables
    specified."   r(198)
-- an error about a varlist, not "there are no wooldid results to replay" -- and that
REPLAY ATTEMPT is then recorded as a failed wooldid call ($WOOLDID_LASTRC = 198).
From then on even `estimates restore' of a perfectly valid stored wooldid model is
refused with "the last wooldid call failed with r(198); the stored results predate
that failed call. Refit wooldid before replaying."  Nothing is stale; the only
documented escape is a full refit.
(The blocking itself is documented at sthlp:1461 and the marker is ruled NFRPL-01/02.
What is new is that a replay ATTEMPT arms the marker, and that the message blames an
estimation failure that never happened.)
Realistic: yes.  Severity LOW-MEDIUM (the blocking itself is documented at sthlp:1461
and the marker is ruled NFRPL-01/02; what is new is only that a replay ATTEMPT arms the
marker and that the message blames an estimation failure that never happened).
Confidence CERTAIN.

===============================================================================
AR-2  replay is not byte-for-byte: an extra blank line before the engine block
Status CONFIRMED.  Repro: repro_AR2.do ; p3 r08/r10 ; p5 h01
-------------------------------------------------------------------------------
ado:48-51 (replay branch) prints `dis as result " "' before part(perm)/part(boot);
the live shell (ado:297/507) does not.  Live has 2 blank lines above "Randomization
Inference (...)" / "Bootstrap Inference (...)", replay has 3.
Severity LOW (cosmetic).  Confidence CERTAIN.

===============================================================================
NEGATIVE RESULTS (see NOTES.md for the full table)
 * anticipation(): 26 of 26 estimable configurations BIT-IDENTICAL (%21x on every
   visible e() scalar, macro and matrix) to the documented hand-edit.  k=0 an exact
   no-op; k<0 and non-integer k refused informatively; relative period 0 sits exactly
   at ttre-k.  The only exception is the one sthlp:826 already carves out
   (subgroups(ttre) keeps the RECORDED onset labels).
 * level() on every OTHER surface: exact.  level(L) == `set level L' bit for bit;
   every interval satisfies lb = b - crit*se at the requested level on the main,
   subgroup, pairwise, gteffects, impvalues, poisson, pretrendpower and
   simultaneousbands surfaces; labels correct including non-integer (lb_92p5ci);
   label round-trip identical for all 9000 two-decimal levels; level() perturbs no
   point estimate, SE, statistic or p-value and does not disturb the RNG.
 * replay of the TABLES: byte-identical everywhere tested, including the Roth
   pretrendpower block (the seeded "different code path" concern -- REPLAY-EXACT at
   the default level AND at level(92.5) and level(80)), with no dataset in memory,
   after estimates store/restore, with two stored models, and twice in a row.

===============================================================================
UNTESTED
 * A true keyboard Break (batch cannot send one).  `capture noisily' catches rc 1
   identically to any error, so the forced-failure block in p1 is the structural
   equivalent, but it is not the same event.
 * e(cfx*) / ame interval identity at level(): p4 L6 and p7 f06b both died on MY option
   errors (missing cfxlattice(); r(907) maxvar), so the continuous-treatment cell
   intervals at a non-default level have NO measurement.
 * anticipation() successfully PULLING an out-of-window cohort into the window: c22/c23
   hit r(459) in BOTH arms (no never-treated at all) and c04's oostreatedcontrols(3) did
   not reach the boundary.  Equal failure, not a differential test.
 * anticipation() + permonsetcohorts (needs i = cohort data; not built).
 * level() under unconditionalse / customvce().
 * level() interaction with saveplots' written plot data.
