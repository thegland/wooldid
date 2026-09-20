# hv-safety-case -- FINDINGS

Pinned build asserted in every log (`which wooldid` -> .../bluesky/build/wooldid.ado).

---

## HV-SAFETY-CASE-1  (seeded; verdict TRUE-BUT-DIFFERENT, severity low)

**Title.** safety() is validated ~500 lines after it is first consulted, but nothing is
silently accepted; the cost is one misleading error message.

**Claim.** `safety(OFF)`/`safety(Off)`/`safety(On)` are always rejected r(198)
"safety() must be on or off." at ado:5473. No intervening `safety`-keyed site (4617,
4623, 4831, 4863-4, 4959, 4973) silently takes the wrong branch, because every one of
them tests `!= "off"` and an unrecognised casing therefore takes the CONSERVATIVE
branch. The only observable defect is that a `safety`-keyed refusal at ado:4617 fires
first and tells the user to "specify safety(off)" without mentioning that their casing
is the problem.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case/repro_safety_case.do`

**Observed.**
```
safety(off)  rc 0     safety(on)  rc 0     safety( off ) rc 0   (Stata trims padding)
safety(OFF)  rc 198   safety(Off) rc 198   safety(oFF)   rc 198
safety(ON)   rc 198   safety(On)  rc 198     all: "Error: safety() must be on or off."

wooldid y i t ttre, att_i unconditionalse safety(OFF)
  -> r(198) "Error: unconditional standard errors may be unreliable with non-att
             estimand types."   "Set option safety(off) to disregard this error."
     (casing never mentioned)
```

**Expected.** The domain screen should run before any site keyed on the value, so the
user is told the real problem. The recon claim's alternative ("or is accepted
somewhere") is REFUTED: it is never accepted.

**Surface.** `safety()`. **Realistic.** true (uppercase option values are a common typo,
and `safety(On)` -- asking for the default -- is rejected too). **Confidence.** certain.

---

## HV-CTW-FLOAT  (seeded lead, extended; verdict TRUE, mechanism corrected; REGRESSION)

**Title.** contreatwithin + contreat() silently suppresses the float-above-2^24
silent-merge caution, contradicting sthlp:73 and sthlp:1195. Introduced by CS-98
(e26b3dd, 2026-09-09).

**Claim.** ado:4829 rewrites the working `safety` to off, and ado:4973 -- the BS-34
caution warning that `i`/`t`/`ttre` stored as float above 16,777,216 was rounded before
wooldid saw it, silently merging cohorts -- gates on that rewritten value. So a user who
never typed `safety(off)`, and even one who typed `safety(on)` explicitly, loses the
caution. Before e26b3dd the rewrite sat BELOW this check (old ado:4657 vs old ado:4493)
and the caution fired, so this is a regression, not a standing design.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case/repro_ctw_floatcaution.do`
(verified standalone; log alongside). Source evidence:
`git show e26b3dd^:wooldidinstall/wooldid.ado` -> local copy `wooldid_pre_e26b3dd.ado`.

**Observed.** 60 units, `gen float i = 20000000 + unit` (an 8-digit administrative code
in the default `gen` storage type); onset assigned from the ROUNDED id so the
onset-constancy screen -- which correctly gates on the user's own value -- does not fire.
Float rounding leaves 31 distinct ids out of 60.
```
LEG1  contreat, no contreatwithin          -> caution PRINTS,  rc 0
LEG2  contreat + contreatwithin            -> NO caution,      rc 0, SE .07147848
LEG3  contreat + contreatwithin safety(on) -> NO caution,      rc 0
LEG4  double ids + contreatwithin          -> rc 0,            SE .0684769
```
The run silently proceeds on 31 merged clusters instead of 60, and the clustered SE is
4.4% away from the panel the user meant, with nothing said. No display line and no e()
key reports the relaxation (`e(wdid_safety)`, `e(safety)`, `e(wdid_safetystate)` all
absent); `e(cmdline)` echoes the user's own line, which contains no safety(off).

**Expected.** sthlp:73 promises the caution whenever `i`, `t` or `ttre` is float above
2^24, with no condition attached. That promise is flatly false under
`contreatwithin` + `contreat()`.

The published safety inventory also places this site in NEITHER family. sthlp:1195 splits
the sites in two: the own-value half is "every data-integrity assert", followed by an
ENUMERATING parenthetical -- "(onset constancy within units, weight sign, the lattice and
cfxlattice domains, onset-period commensurability, the customweightmultiplier domain)" --
that does not contain the float caution; the working-value half is enumerated at
sthlp:1209 as the collinearity family, the realized relative-time range asserts and the
control-availability count, which does not contain it either. Measured behaviour puts it
squarely in the working-value family.

**Fix.** ado:4827's own comment states the rule -- "the user's own value is kept as
`safetyuser' for the checks that must see what was typed". ado:4973 was not converted.
Changing its gate to `safetyinput` (as ado:718/772/1355/5246/5290/5315/5340 already do)
restores the pre-e26b3dd behaviour.

**Surface.** `contreatwithin` + `contreat()`; `safety()`; the BS-34 float screen.
**Realistic.** true -- `gen id = <8-digit code>` defaults to float, and continuous-dose
designs with `contreatwithin` are exactly where administrative panels land.
**Severity.** medium. **Confidence.** certain.

---

## HV-CTW-IPOS  (seeded lead, extended; verdict TRUE; REGRESSION)

**Title.** contreatwithin + contreat() silently skips the "i must consist of positive,
non-zero integers" refusal, turning a named r(411) into a raw r(452) on an internal
tempvar name. Same CS-98 mechanism.

**Claim.** ado:4959 gates on the rewritten working `safety`, so under
`contreatwithin` + `contreat()` the screen never runs. sthlp:1208 lists this refusal
under "refuses, and safety(off) proceeds" without flagging it as a working-value site,
and sthlp:1195 places data-integrity asserts on the untouched side of the rewrite. Before
e26b3dd the check sat above the rewrite (old ado:4488 vs old ado:4657) and did fire.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case/repro_ctw_ipositive.do`
(verified standalone; log alongside).

**Observed.** One treated unit's `i` made bad; every leg otherwise identical.
```
                      no contreatwithin                with contreatwithin
non-integer i   r411 "Error: i must consist of    r452 "__000005: factor variables
                positive, non-zero integers"       may not contain noninteger values"
negative i      r411 same informative message     r452 "__000005: ... may not contain
                                                   negative values"
i == 0          r411 same informative message     rc 0, completes with no message
```
Honest qualifier on the `i == 0` leg: on the panel tested the estimate is BIT-IDENTICAL
(+1.aa36e6ab488d3X+000) to the same panel with ids shifted by +1 and by +100, so the
aliasing of `i == 0` with wooldid's internal never-treated code (`_i` = i if treated, 0
otherwise, ado:895 and inside wooldid_reltime) did NOT move the number here. What is lost
is the guard, not -- on this design -- the answer. I did not find a design where it does
move the number.

**Expected.** The documented r(411) refusal, whose whole purpose is to convert Stata's
factor-variable errors into a message naming the user's own variable. As with the float
caution, the published inventory assigns this site to neither family: sthlp:1208 lists
"the cohort identifier not being a positive integer" among the refusals without marking
it a working-value site, sthlp:1195's own-value enumeration does not name it, and
sthlp:1209's working-value enumeration does not name it. Measured behaviour puts it in
the working-value family.

**Fix.** Same as HV-CTW-FLOAT: ado:4959's gate should read `safetyinput`, per the rule
ado:4827 states in the very commit that moved the rewrite above it.

**Surface.** `contreatwithin` + `contreat()`; `safety()`; the ado:4959 identifier screen.
**Realistic.** true for non-integer and negative identifiers (decimal-coded keys, sentinel
codes such as -9); lower for `i == 0`. **Severity.** medium (uninformative error on
realistic input; doc-vs-code divergence; regression). **Confidence.** certain.

---

## Also noticed (not filed separately)

- ado:4863/4864, the `espostlength(-999)` advice lines, were also pulled under the CS-98
  rewrite. Cosmetic: the disclosure line above them prints either way by the P6 ruling.
- ado:4927 hands `si(`safety')` -- the rewritten value -- to wooldid_permctl. Moot on this
  build because `permute()` with `contreatwithin` is refused (sthlp:472/1235), but it is
  the one remaining live path where the rewrite could reach a screen it was not meant to.
- `contreatcontrols()` is CLEAN on the same axis: ado:4792 validates its tokens
  case-sensitively and exactly, above every `strpos(...,"byTime"/"byCohort"/"byTreated"/
  "byTtre")` consumer at ado:1200-1256, 1384 and 5171.
- Not chased: on the p4 panel, `att` was bit-identical with and without `contreatwithin`
  (B0 vs B3, +1.aa36e6ab488d3X+000, N=720 both). Plausibly by design when never-treated
  units are present -- `contreatwithin`'s documented effect there is mostly to skip the
  control-group checks -- but I did not verify it.
