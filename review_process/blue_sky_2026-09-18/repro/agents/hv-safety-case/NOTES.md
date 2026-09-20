# hv-safety-case -- verification notes

Build asserted in every log: `which wooldid` ->
`/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`, `*! version 1.4 2026-08-19`.

## Files
- `mkdata.do`                  shared generator used by p1 only
- `p1_case.do` / `.log`        casing matrix
- `p2_ctw.do` / `.log`         first contreatwithin pass (has an RNG confound in leg B -- superseded by p4)
- `p3_ctw2.do` / `.log`        second pass (same confound in leg B)
- `p4_same.do` / `.log`        same-panel version, confound removed
- `p5_conseq.do` / `.log`      consequence surface (cluster count, opaque errors)
- `repro_safety_case.do`       SELF-CONTAINED repro, finding 1
- `repro_ctw_floatcaution.do` / `.log`  SELF-CONTAINED repro, finding 2 (verified standalone)
- `repro_ctw_ipositive.do` / `.log`     SELF-CONTAINED repro, finding 3 (verified standalone)
- `wooldid_pre_e26b3dd.ado`    pre-CS-98 source pulled with `git show e26b3dd^:wooldidinstall/wooldid.ado`

## Claim 1 -- safety() validated 500 lines after first consulted, case-sensitively
TRUE as to code ordering, FALSE as to the implied consequence.
ado:4617 is the first `safety`-keyed site; ado:5473 is the case-sensitive
`inlist("`safetyinput'","on","off")`. Measured (p1_case.log):
  safety(off)   rc 0      safety(on)   rc 0      safety( off )  rc 0 (Stata trims)
  safety(OFF)   rc 198    safety(Off)  rc 198    safety(oFF)    rc 198
  safety(ON)    rc 198    safety(On)   rc 198     -- all "safety() must be on or off."
Nothing is ever silently accepted. Every site in the window (4617, 4623, 4831,
4863-4, 4959, 4973) is a message or a refusal, and `!= "off"` means an unrecognised
casing takes the CONSERVATIVE branch. So no check silently took the wrong branch.
The one observable bite: `safety(OFF) unconditionalse att_i` dies at 4617 with
"Set option safety(off) to disregard this error" -- the user is told to do what they
believe they did, and the casing is never mentioned. Severity low.

## Claim 2 -- contreatwithin + contreat() rewrites safety to off
TRUE, and it is a REGRESSION plus a live doc-vs-code divergence, not just a design note.

ado:4828-4829 (CS-98, author ruling 2026-09-08, shipped in e26b3dd 2026-09-09):
    local safetyuser `safety'
    if "`contreatwithin'" == "contreatwithin" & "`contreat'" != "" local safety off
Before e26b3dd the rewrite sat at old-line 4657, i.e. BELOW these two checks
(old lines 4488 and 4493). The move pulled them under the relaxation and they were
never converted to `safetyuser'/`safetyinput', although the CS-98 comment says
"The user's own value is kept as `safetyuser' for the checks that must see what was typed."

Sites in ado:4830-5180 that read the rewritten working value:
  4831        cfx-plot caution wording          <- the intended target of CS-98
  4863/4864   espostlength(-999) advice lines   <- newly suppressed, cosmetic
  4927        si() handed to wooldid_permctl    <- moot, permute()+contreatwithin is refused
  4959        "i must be positive, non-zero integers" (error 411)  <- NEWLY SKIPPED
  4973        float-above-2^24 silent-merge caution (BS-34)        <- NEWLY SKIPPED
Everything else `safety`-keyed in the worker (801/836/908/955/1713/2413) was already
relaxed before CS-98, because the c_local blanket always exported the rewritten value.

sthlp:1195 says the rewrite leaves alone "every data-integrity assert"; sthlp:1208
lists "the cohort identifier not being a positive integer" among the refusals without
flagging it as working-value; sthlp:73 promises the float caution unconditionally.
All three are false on this build.

Nothing at runtime says safety was relaxed: no display line, and no e() key
(`e(wdid_safety)`, `e(safety)`, `e(wdid_safetystate)` all absent -- p3_ctw2.log C legs).
`e(cmdline)` echoes the user's own line, which does not contain safety(off).
Explicit `safety(on)` does NOT restore the guards (p4 leg A4, p5 leg C1/LEG3).

## Measured consequences
Float id (60 units -> 31 distinct after float rounding, onset assigned from the rounded
id so the onset-constancy screen does not fire and mask it):
  no contreatwithin : caution prints
  contreatwithin    : silent; run completes on 31 merged clusters
  SE with cluster(i): .07147848 (float, merged) vs .0684769 (double, intended) = 4.4% apart
i-positive-integer screen:
  non-integer i : r411 informative  ->  r452 "__000005: factor variables may not contain
                  noninteger values" under contreatwithin
  negative i    : r411 informative  ->  r452 "__000005: ... may not contain negative values"
  i == 0        : r411 informative  ->  rc 0, run completes.  HONEST NEGATIVE: on the panel
                  tested, the estimate is BIT-IDENTICAL to the same panel with ids shifted
                  by +1 or +100 (+1.aa36e6ab488d3X+000 in all three), so the aliasing of
                  i==0 with the internal never-treated code `_i'==0 did NOT move the number
                  here.  The loss is the guard, not (on this design) an answer.

## Probed and clean
- contreatcontrols() is validated token-exact and case-sensitively at ado:4792, which is
  ABOVE every `strpos(...,"byTime")` / "byCohort" / "byTreated" / "byTtre" consumer
  (ado:1200-1256, 1384, 5171). No case-sensitivity hole there.
- `safety( off )` with padding is accepted; Stata trims. No hole.
- No casing of safety() is ever silently accepted.
