# interactions recon - NOTES

Task: build OPTION-INTERACTION RISK MATRIX -> /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/interactions.md

## Map of the guard hubs (pinned build line numbers)
- wooldid (shell)        19-530   : replay branch 25-46; permute orchestration 108-300; bootstrap orchestration 311-520
- wooldid_worker         531      : syntax line 553 (must be byte-identical to validate's 4572)
- wooldid_validate       4568-5521: ALL parse-time cross-option guards + fast-collapse gate (5050-5150) + fastnose gates (5383-5410) + option-domain screens (5450-5500)
- wooldid_permctl        5805-6181: permute() screens
- wooldid_bootctl        6183-6524: bootstrap() screens (+ stray-companion screen, runs even with no bootstrap())
- wooldid_bootrun        8887-9543: rep loop; nesting regime at 8987-8994; 4 legs (bayes / phys / count / re-id)
- wooldid_permrun        7465
- rep option strings: wdp_ro (perm) at shell:112-149 ; wdb_ro (boot) at shell:322-364

## Confirmed-by-reading leads (details in the output file)
1. bootctl: the reps() resolution caution block (6480-6490) is NESTED INSIDE `if "`clus'" != ""`.
   Unclustered bootstrap gets no resolution caution at all.
2. T2 removed the refusal of bootstrap(standard)+cell-weighted estimand+NO cluster().
   bootrun sets wdb_nest=0 (keep-i) there and does plain `bsample` (obs resampling).
   Cell set is then ~fixed => the cell-sampling variance is omitted -- which is EXACTLY
   the defect that bootstrap(bayes)+cell-weighted still ERRORS on (bootctl 6462-6470).
   Asymmetric guard; the unclustered standard path is silent.
3. wdp_ro (permute) does NOT forward gteffects; wdb_ro (boot) does. permwire gets gt() -> blanks.
4. wdp_ro forwards customvce(); bootstrap REFUSES customvce(). Asymmetry.
5. permctl nesting check (6167) skipped entirely under permonsetcohorts.
6. bootrun re-id leg: cluster(<clus>) -> cluster(__wdbns) subinstr is a literal string edit (9376),
   only when boottstat AND clus!="" . P7 assertion 9394-9420 checks it landed.

## STATUS 2026-09-16: DELIVERABLE WRITTEN
/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/interactions.md (552 lines)
All ado:NNNN line numbers in it were re-verified by grep against the pinned build.

## Stata run: recon_check.do / recon_check.log (7 sentinels, all rc as expected)
CONFIRMED by that run:
 - A5 clustered bootstrap reps(10) prints the resolution caution; UNCLUSTERED does not
   (brace nesting: `if "`clus'" != "" {` opens bootctl:6459, closes :6501; the reps caution
    block is :6473-6478, inside it).
 - A1 `att_it bootstrap(standard) reps(20)` with NO cluster() runs rc 0, 20/20 completed,
   e(boot_nclus)=. ; the bayes twin still errors 198 at bootctl:6439.
 - regtol(3) under unconditionalse: rc 0, no message (silent no-op; regtol only reaches
   reghdfe/ppmlhdfe at ado:1461-1476).
 - contreatpoly(3) with no contreat(): rc 0, no message (silent no-op).
 - robust alongside cluster(): rc 0, header reads cluster(clus) (robust silently dropped at
   ado:1449-1451).

## POST-REVIEW REVISIONS (final)
 - A7 (permute + sample-dropping options) DEMOTED to Tier C addendum C1: under the sharp null the
   drop rule is symmetric in the assignment, so exactness survives; only power moves. NOT a defect.
 - New Tier A A7: boottstat + cluster() strictly finer than i on the re-id keep-i regime, WITH the
   `global WOOLDID_BOOTNOSUB 1` negative control (must reproduce the 2-1/G ~ sqrt(2) inflation).
   This is the core T3 deliverable and previously had no ranked row.
 - `robust` + `cluster()` relabelled NOT a defect (clustered VCE is robust; standard Stata).
 - A5 relabelled INCONSISTENCY, not REGRESSION (I never verified where the caution sat pre-RF-2).
 - A2: i is the SECOND positional argument; fe(cohid) with i=cohid is a repeated absvar, so
   controls(i.cohid) / coarsecohortcontrols(i.cohid) are the cleaner leak configs.
 - WOOLDID_REPMODE is consumed at ado:538-539 (worker), not in the shell replay branch.
