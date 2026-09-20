# safety() guard-by-guard adversarial probe -- FINDINGS

Agent: SAFETY.  Method: for each safety-gated site, build the input that trips it,
then turn the guard off and characterise what comes out.
All probes are self-contained do-files in this directory; each has its .log beside it.

---

## SAFETY-01 (HIGH / realistic / certain)  the si8 maxvar preflight FALSELY REFUSES
##            an ordinary wide panel -- contradicting its documented zero-false-positive design

**Repro:** `p07_si8_falsepos.do` (maxvar 2048) and `p08_si8_confirm.do` (Stata's DEFAULT maxvar 5000).

**Claim.** The variable-count preflight at `ado:1355-1377` is NOT a lower bound on the
regression's variable demand on the default engine path, so it refuses runs that
complete correctly.  `ado:1341-1348` states the opposite ("The preflight uses the EXACT
treated interaction-cell count ..., a lower bound on the regression's variable demand,
so an early refusal can only fire when the late r(907) was certain anyway"), and the
project record (retracted defect, 2026-09-16) states the design "accepts false negatives
to guarantee zero false positives".  Measured: a **false positive**.

**Observed (p08, DEFAULT maxvar = 5000).**  4,800 units x 8 years = 38,400 rows, 20 of the
units ever treated, `wooldid ync unit t ttre, att cluster(unit)`:

    Error: this specification requires at least 5078 variables (about 70 treated
    cohort-by-period interaction terms) but maxvar is currently 5000, so estimation
    would fail with r(907).
    r(907);

The identical call with `safety(off)` returns **rc 0**, e(N) = 38,400, and a valid ATT
table.  p07 reproduces it at `set maxvar 2048` (Stata/BE's fixed ceiling, which a BE user
cannot raise at all): demand quoted as 2178, refusal r(907), `safety(off)` rc 0, N = 15,200.

**Expected.**  The refusal says estimation *would fail with r(907)*.  It does not fail.

**Diagnosis (engine-conditional -- this is the part that decides the fix).**
`local si8demand = si8cells*(1 + ("`contreat'" != "")) + si8ilev + si8tlev + 200` (`ado:1369`),
where `si8ilev` is the number of distinct `i` values and `si8tlev` the number of distinct
`t` values.  Those two terms are a real variable demand ONLY on the `reg`/`poisson` path
(`ado:1482-1483`), which writes `i.`i' i.`t'` out as regressors.  On the DEFAULT path they
are ABSORBED -- `reghdfe ... , absorb(i.`i' i.`t' `fe' ...)`, `ado:1476` -- and never enter
the design matrix, so charging for them is pure overcharge, and the overcharge is the
panel's unit count: the very dimension the message tells the user is "comparatively cheap".

Pinned in one run (`p10_engine_split.do`, 2,300 units, `set maxvar 2048`, quoted demand
2,578):

    default safety                      -> r(907), refused
    cluster(unit) safety(off)           -> rc 0, e(N) = 18,400   [reghdfe: FALSE POSITIVE]
    customvce(cluster unit) safety(off) -> r(907)                [reg path: demand is real]

So the fix is not "drop the level counts" -- that would reintroduce a late r(907) for
`unconditionalse`/`customvce()` users.  It is to add them only on the branch that
materializes them.  The `si8cells + 200` part is a genuine lower bound on both branches.

**Why it matters / realism.**  Any panel with more than ~4,700 units and few treated ones
trips it at Stata's default `maxvar`: school districts (~13k), ZIP codes, firms, hospitals,
census tracts, individuals in a household panel.  On Stata/BE (maxvar fixed at 2048) the
threshold is ~1,800 units and the user has no way to raise it.  The remedy the message offers -- coarsen `i`, move units into `fe()` -- is harmless for
plain `att` (the same parameter is estimated), but it is NOT harmless for the cell-weighted
family: `att_i`, `att_it` and `att_itime` weight cells and units equally, so redefining `i`
redefines the estimand.  A user chasing this refusal into a coarser `i` can therefore end up
reporting a different parameter than the one they asked for, to escape a constraint that was
not real.

**Surface:** `safety()` / si8 maxvar preflight (ado:1341-1377).  **Seeded:** partly -- the
seed asked about a `lattice()` false NEGATIVE; the false POSITIVE is mine.

---

## SAFETY-02 (MEDIUM / realistic / certain)  contreatwithin silently disables two screens
##            the help says gate on the user's own safety() value

**Repro:** `p02b_ctw_map.do`.

**Claim.** `sthlp:1188-1200` says "Roughly half the sites gate on your own value and are
therefore untouched by that rewrite - **every data-integrity assert** ...", and `sthlp:913`
says contreatwithin disables "the estimability-relaxation checks; the data-integrity
screens, such as onset constancy within units, still run unless safety(off) is given
explicitly".  Two data-integrity sites gate on the WORKING value and are therefore
silently disabled by `contreatwithin` + `contreat()` (the rewrite at `ado:4829`):

  * `ado:4959` -- `assert i > 0 & int(i) == i`, r(411).  The help lists this one explicitly
    under "**refuses, and safety(off) proceeds**", i.e. as a user-value site.
  * `ado:4973` -- the float-above-2^24 caution on `i`, `t`, `ttre`.  This site is not in the
    help's working-value list (that list is "the whole post-estimation collinearity family,
    the realized relative-time range asserts, and the control-availability count").

**Observed.**  With `t`/`ttre` stored as float 8-digit calendar codes (`tdate = 20000101 +
(t-2000)*10000`): the binary run prints the caution twice
(`p02b_ctw_map.log:922,928`); adding `contreat(dose) contreatcontrols(basic) contreatwithin`
prints it **zero times** (same log, the H2 block), on the same variables.  The caution
exists precisely because the damage -- adjacent periods/cohorts merged by float rounding
before wooldid ever sees them -- is otherwise silent and produces a wrong answer.
For `ado:4959`: an `i` coded 0 is refused r(411) on the binary path and accepted with no
message under contreatwithin.

**Expected.**  Per the help, both are data-integrity screens and should survive the
contreatwithin rewrite.

**Not a wrong number, and I checked.**  With 0-indexed cohort codes the numbers are
IDENTICAL to the 1-indexed run (att -.95083572, ame -.92666667 both ways; and on the
binary path att 1.5800311 under `safety(off)` == the 1-indexed default answer).  `_i == 0`
does not collide with the never-treated sentinel, because never-treated rows always carry
`_t == 0`.  A non-integer `i` under contreatwithin dies r(452) -- an uninformative crash
where the named r(411) would have fired, but not a silent number.

**Surface:** `contreatwithin` / `safety()`.  **Seeded:** the contreatwithin rewrite was
seeded; which two sites it disables, and the doc contradiction, are mine.

---

## SAFETY-03 (MEDIUM / realistic / probable)  the commensurability bypass manufactures a
##            14-sigma pre-trend with no message

**Repro:** `p06_unlocks.do`, block (b).

**Framing (read this first).** The issue is not which rounding convention is right for a
mid-year effective date -- both are defensible.  It is that wooldid adopts one silently, and
the adopted one folds the true onset year into the PRE window, which is exactly the failure
`sthlp:73` names as the pre-fix behaviour ("silently mislabeled the onset period as
pre-treatment").  The visible symptom is a pre-trend that the DGP does not contain.

**Reachability.** The silent path needs a pre-window request.  I ran `esprelength(2)`;
with the default `esprelength(0)` the shifted `_t` stays fractional and the run is expected
to die r(452) instead -- not tested.

**Claim.** `safety(off)` on a fractional onset (`ttre = ttre + 0.5`, a mid-year effective
date on an annual panel) does not merely proceed: it silently re-labels the true onset
period as pre-treatment, and the run reports a large, highly significant pre-treatment
`att` on a DGP whose pre-trend is exactly zero.

**Observed.** On an otherwise clean 60-unit panel, `esprelength(2) espostlength(2)`:

    ttre (correct)      att = 1.1373096   pre-period att = -0.00211341 (t = -0.04, p = .965)
    ttre + 0.5, safety(off)  att = 1.096691    pre-period att =  0.68206619 (t = 14.39, p = .000)
    ttre + 1 (honest recode) att = 1.096691    pre-period att =  0.68206619

The bypassed run is bit-for-bit the `ttre + 1` answer: the onset year is absorbed into the
pre window.  Nothing is printed.  Default safety refuses with r(459).

**Expected.** Either a refusal (the default) or, at minimum, a statement of the convention
adopted.  `sthlp:73` itself names this failure mode as the pre-fix behaviour ("silently
mislabeled the onset period as pre-treatment") and the safetysites table lists the screen
under "refuses, and safety(off) proceeds" without saying what "proceeds" does here.
Note the refusal message itself never mentions `safety(off)` -- only the help does -- so
the user most likely to hit this is one who set `safety(off)` for the runtime reason
`sthlp:1167` recommends and later merged in a fractional onset.

**Surface:** `safety(off)` / onset-commensurability screen (ado:5340-5366).  **Seeded:** no.

---

## SAFETY-04 (LOW / realistic / certain)  the negative-weight screen refuses on rows the
##            run discards anyway

**Repro:** `p06_unlocks.do`, block (c), the `ymiss` / `wneg2` pair.

**Claim.** `ado:5246` asserts `w >= 0` on the RAW `if` domain (`ifassert`), not on the
estimation-aware domain.  A negative weight sitting only on rows with a missing outcome --
rows deterministically dropped by markout -- refuses the run r(411).  `safety(off)` on the
same data runs to completion.  The project already recognised this class and fixed it for
the subgroup screen (VAL-G15A, `ado:5253-5257`: "the raw-if-domain check that lived here
falsely rejected runs whose only negative values sat on deterministically-discarded rows
(e.g. missing outcomes)") and for `customweightmultiplier()` (`ado:718`, which uses
`if `ifstatement'`).  The weight screen itself was not moved.

**Observed.** `WC_RC = 411` (default safety, negative weights only on missing-`y` rows);
`WD_RC = 0` (same call, `safety(off)`).

**Severity note.** Negative sampling weights are rare in practice, which is why this is
low, not medium.

**Surface:** `safety()` / weight-sign screen.  **Seeded:** no.

---

## SAFETY-05 (LOW / certain)  the esfixedbaseperiod pre/post refusal degrades to a bare
##            Stata error under safety(off)

**Repro:** `p05_bypass_numbers.do`, block (B).

Default safety: `Error: cannot find sufficient pre- and post-treatment periods for the
treated units.` r(459).  Same call with `safety(off)`: `insufficient observations`, r(2001)
-- no wooldid concept named.  A crash, not a wrong number, so low.

---

# NEGATIVE RESULTS (probed hard, came back CLEAN)

* **safety() token handling is clean.**  `safety(ON|OFF|Off|oFF|0|1|yes|no)` all exit
  r(198) with `Error: safety() must be on or off.`  The seeded 500-line gap between the
  first consumer (ado:4618) and the case-sensitive validation (ado:5473) is benign: an
  uppercase token compares unequal to "off" at every gate, so those 500 lines treat it as
  safety-ON -- the conservative direction -- and the run then dies with the right message.
  `safety(off )` / `safety( off)` are accepted (Stata strips the spaces).  `p01_tokens.do`.

* **safety(off) moves no number on runs that succeed either way.**  Seven configurations
  (plain att; all four estimands at cluster(region) with a 2x2 event-study window;
  subgroups(); controls()+fe()+jointtests; impvalues; esfixedbaseperiod+adversarial+makeplots)
  run default and again with `safety(off)`: `mreldif(e(b))`, `mreldif(e(V))` and
  `mreldif(e(wdidmainresults))` are all exactly **0**.  `p04_offchanges.do`, CFG1-CFG4,
  CFG7, CFG8.  I also traced the one place the code branches on safety while mutating data
  (`ado:12086/12108`, the reference-period forcing in `wooldid_reltime`): the
  `safety(off)` branch and the caution branch both end with `_t` missing at `t == mint`,
  so they are number-identical by construction, not by luck.

* **The no-never-treated bypass gives the CORRECT number.**  Every unit treated inside the
  window (onsets 2003/2005/2007, exact DGP): default r(459); `safety(off)` reports
  att = 1.6666667 against a design truth on the realized sample of 1.6666667, with
  240 of 800 observations dropped and the drop disclosed in the header.  The refusal is
  conservative, and the bypass is sound on this shape.  `p05_bypass_numbers.do`, block (A).

* **A collinear co-timed policy control does NOT trip r(322).**  `controls(policy2)` where
  policy2 is exactly one treated unit's post indicator: reghdfe drops `policy2`, not the
  cell dummy, and the ATT equals the design truth (1.05) under both safety settings.
  `p05_bypass_numbers.do`, block (C).

* **A region-by-year FE containing one treated state does NOT trip r(322) either** -- the
  state's rows are removed as reghdfe singletons before the screen can see them, and the
  reported ATT (1.0714286) equals the design truth on the reduced sample.
  `p03_r322.do`.

* **`unconditionalse` + a cell-weighted estimand** (refused at ado:4617, bypassed by
  safety(off) with no caution at all) produced an SE indistinguishable from the honest
  ones on a clean 60-unit panel: att_i se = 0.0505284 unconditional vs 0.0506031
  conditional-robust vs 0.0491569 clustered.  One dataset is not a coverage study, so this
  is "did not reproduce here", not "the refusal is wrong".  `p06_unlocks.do`, block (a).

* **bootstrap(bayes) + att_i** and **subgroups() splitting a unit over time** both refuse
  under default safety and proceed under `safety(off)` exactly as documented.
  `p06_unlocks.do`, blocks (d), (e).

---

## Repro-file note for verifiers

`REPRO_si8_false_refusal.do` is the FULLY SELF-CONTAINED demonstration of SAFETY-01: it
generates its own 4,800 x 8 panel inline, uses no include, no repo data and no test suite,
and runs at Stata's own default `maxvar` (5000).  Blocks (1) and (2) are the finding:

    MAXVAR = 5000
    ROWS   = 38400
    (1) default safety -> "Error: this specification requires at least 5078 variables
        (about 70 treated cohort-by-period interaction terms) but maxvar is currently
        5000, so estimation would fail with r(907)."   RC_DEFAULT = 907
    (2) same call + safety(off) -> RC_SAFETYOFF = 0, N = 38400, att = 1.440716
        (DGP effect 1.5; the reported se 0.0029 is the known analytic understatement at
        cluster() = i, not part of this finding)

Block (3) of that file re-runs the same data through the `reg` path and is SLOW (a
4,800-level factor on 38,400 rows).  It is only confirmatory: the decisive engine split is
`p10_engine_split.do`, which settles it in seconds at 2,300 units and `set maxvar 2048`.
