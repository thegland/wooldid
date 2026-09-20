# SKEPTIC pass on "safety(): what the guards prevent and what safety(off) unlocks"

Every number below comes from a do-file in this directory that generates its own data and
asserts the pinned build with `which wooldid`.  I did not run the original prober's scripts
as evidence.

---
## SAFETY-01 -- CONFIRMED (and stronger than claimed).  severity high
**Repro:** `REPRO_SAFETY01_si8_falsepositive.do` (+ `r1a_si8_defaultmaxvar.do`,
`r1b_si8_control.do`).

The si8 preflight (ado:1355-1377) refuses, with r(907), runs that complete correctly on the
same `maxvar`.  Three things I established that the original finding did not:

1. **Bit-identity, not merely "it runs."**  4,800 units x 8 years, 20 ever-treated, at
   Stata's DEFAULT `maxvar` 5000.  Default safety: `"...requires at least 5107 variables
   ... so estimation would fail with r(907)"`, r(907).  `safety(off)`, same ceiling:
   rc 0, N=38400, att hex `+1.8c965d1a10262X+000`, se hex `+1.3ed047a4ebb78X-007`.
   The SAME call under DEFAULT safety with `set maxvar 20000`: **the same two hex values**
   (`ATT_BITIDENTICAL=YES`, `SE_BITIDENTICAL=YES`).  The refused run is not degraded, not
   approximate, and not lucky - it is the answer wooldid's own default path gives.
2. **The trigger is the NEVER-TREATED CONTROL COUNT.**  With the SAME 20 treated units:
   4,000 units -> rc 0; 4,900 units -> r(907) quoting demand 5207.  The 920 extra units are
   pure never-treated controls; they add precision, and on the default engine they add no
   regressor (reghdfe absorbs `i.i`, ado:1476).  A user is refused for making their control
   group bigger.
3. **The overcharge is engine-conditional**, at the default ceiling: 5,200 units,
   `maxvar` 5000 -> default safety r(907) quoting 5507; `safety(off) cluster(unit)` rc 0
   N=41600 (reghdfe, FALSE POSITIVE); `safety(off) customvce(cluster unit)` r(907) - on the
   reg path (ado:1482-1483 writes `i.i i.t` out as regressors) the demand IS real.  So the
   fix cannot be "drop `si8ilev`+`si8tlev`"; it has to branch on the engine.

This contradicts ado:1341-1348 ("a lower bound ... an early refusal can only fire when the
late r(907) was certain anyway"), sthlp:1163 ("if the count of variables the specification
requires exceeds the current maxvar"), and the project record of 2026-09-16 ("The design
accepts false negatives to guarantee zero false positives", opened and retracted the same
day).

Aggravating, and measured in `r8_calib.do` block E: **the preflight also refuses the remedy
the help itself prescribes for this specification.**  `i`=unit with `cluster(unit)` is the
case where sthlp:251 sends the user to `bootstrap()`/`permute()` for valid inference; the
preflight runs on rep 0 (`wdrepmode != 1`), so `wooldid y unit t ttre, att cluster(unit)
bootstrap(standard) reps(2)` on the 4,800-unit panel also exits r(907) quoting 5107.
Aggravating further: wooldid's own two message lines offer only "coarsen `i` / use `fe()` /
`safety(off)`" - the "Use set maxvar" hint that appears below them is Stata's r(907)
boilerplate, not the guard's - and coarsening `i` changes the estimand (`att_i` becomes a
mean over cohorts rather than over units), so the guard's own advice is to answer a
different question.

---
## SAFETY-02 -- CONFIRMED as a doc-vs-code divergence; no number moves.  severity medium
**Repro:** `REPRO_SAFETY02_ctw_screens.do`.

- Float 2^24 caution (ado:4973): with `t` and `ttre` held as float 8-digit calendar codes,
  the caution prints **twice** on the binary run, **twice** with `contreat()` alone, and
  **zero times** with `contreat() contreatcontrols(basic) contreatwithin` - with no
  `safety()` typed by the user.  The `contreat()`-without-`contreatwithin` arm is the
  control that pins `contreatwithin` as the cause.
- Positive-integer `i` assert (ado:4959): 0-indexed cohort codes are refused r(411) on the
  binary path and accepted silently under `contreatwithin`.  A non-integer `i` dies r(452).
- sthlp:1199 says "Roughly half the sites gate on your own value and are therefore untouched
  by that rewrite - **every data-integrity assert** (...)"; neither of these two is in that
  parenthetical, and neither survives.  The universal is false.

**Correction to the original finding: no number moves.**  I coded a *treated* cohort as 0
(the case where `_i == 0` could collide with the never-treated sentinel) and compared against
the 1-indexed recoding: `att`, `att_i`, `ame` and `ame_i` are **bit-identical**
(`ATT_MOVED=NO`, `AME_MOVED=NO`), because never-treated rows always carry `_t == 0` while a
treated cohort 0's rows carry `_t > 0`, so the cells never collide.  The consequential half
of this finding is the float caution alone; the `i`-integrality half costs only the named
error (r(411) -> a bare r(452)).

**Control that came back clean:** F13, the onset-constancy screen the `contreatwithin` help
entry names by example (sthlp:913), **does** survive the rewrite - r(459) with no `safety()`
typed.  It gates on `safetyinput` (ado:772) as documented.

---
## SAFETY-03 -- PARTIAL.  severity revised medium -> low
**Repro:** `REPRO_SAFETY03_fracOnset.do`.

What survives: under `safety(off)` a fractional onset proceeds, and the pre-treatment att row
picks up the onset-period effect - `pretreatment att` = 0.739 against -0.157 on the correct
integer onset, on a DGP with exactly zero pre-trend.  The numbers do match an honest `ttre+1`
recode bit-for-bit.

What is refuted:
- **"silently"**: the event-study rows are LABELLED with the fractional horizons.  The
  printed table shows `pre1.5 / pre.5 / post.5`; `e(b_es_att)`'s column names are
  `ES_att_pre1.5`, `ES_att_pre.5`, `ES_att_post.5`; and the `relyr` entries in
  `e(wdidmainresults)` are -1.5, -0.5, +0.5.  A half-integer relative-time grid on an annual
  panel is announced in three places.
- **"adopts the round-up convention"**: there is no rounding anywhere.  `_t = t - ttre`
  exactly (ado:12081), and the +0.7 arm proves it - its rows come out `pre1.7 / pre.7 /
  post.3`, not `pre2 / pre1 / post0`.  The coincidence with `ttre+1` is a pure relabelling
  identity that holds for ANY offset in (0,1), because adding a constant to `ttre` moves every
  row's horizon by the same constant and so partitions the rows identically.  No convention
  was chosen, silently or otherwise.  The pre/post split is exactly `t >= ttre`, which is what
  the help says the treatment indicator means.
- sthlp:73 documents this failure mode by name and names `safety(off)` as the bypass;
  sthlp:1207 lists it under "refuses, and safety(off) proceeds".

Their untested caveat CONFIRMED: default `esprelength(0)` dies r(452).  I add that offset
**+0.3 also dies r(452)** at `esprelength(2)` while +0.5 and +0.7 run - the shifted `_t` has
to land on exact binary integers, so the bypass is a floating-point lottery.

---
## SAFETY-04 -- CONFIRMED.  severity low
**Repro:** `REPRO_SAFETY04_05.do` (part 1).

Negative `aweight` placed only on rows whose outcome is missing (`NEG_WEIGHTS_ON_USABLE_ROWS
= 0`): default safety exits r(411) "weights must not be negative"; `safety(off)` runs
(N=475).  New evidence: the `safety(off)` answer is **bit-identical** to the default-safety
answer with `if !missing(ymiss)` typed (`D2_EQ_D3=YES`), which proves those rows really are
deterministically discarded and the workaround is one token.  The screen at ado:5246 tests
over `ifassert` (the raw `if` domain) while its two siblings were explicitly moved to the
estimation-aware domain for exactly this reason (VAL-G15A comment, ado:5253-5257; and
ado:718 for `customweightmultiplier()`).  Note this sits next to the already-ruled VAL-G06
("negative weights error r(411) **at parse**"); the novel part is only the domain
inconsistency with the two siblings, which is why low is right.

---
## SAFETY-05 -- PARTIAL (real, but not r(2001)).  severity low
**Repro:** `REPRO_SAFETY04_05.do` (part 2).

The named r(459) "cannot find sufficient pre- and post-treatment periods for the treated
units" does degrade to a bare Stata error under `safety(off)`.  But it is not one rc:
- with `esprelength(2) espostlength(2)` -> **r(198) "invalid syntax"**
- with no ES window options at all -> **r(461) "no observations in subpop() subpopulation"**
- the original report's r(2001) is a third variant.

r(198) "invalid syntax" is materially worse than r(2001): it tells the user their
syntactically valid command is mistyped and sends them hunting through their own do-file.

---
## NEW (mine) -- SKEP-01: at `i` = unit the analytic ATT SE degrades with the SIZE OF THE CONTROL POOL, reaching SE/SD 0.21 and 37% coverage
**Repro:** `r7_se_arms.do` (three arms) and `r8_calib.do` block D (the calibration arm).
60 replications each, same DGP, nominal 95%, truth 1.5.

| arm | spec | mean SE / SD(b) | coverage |
|-----|------|-----------------|----------|
| A | `i`=unit, `cluster(unit)`, 20 treated of **300** units | **0.212** | **0.367** |
| D | `i`=unit, `cluster(unit)`, 20 treated of **60** units | 0.508 | 0.733 |
| C | `i`=unit, `robust`, 20 of 300 | 0.646 | 0.783 |
| B | `i`=coarse cohort + `fe(unit)`, `cluster(unit)`, 20 of 300 (the help's intended practice) | 0.891 | 0.900 |

Arm B validates the harness.  The point estimate is identical across arms A/B/C on the same
draws (1.493863), so only the variance is at issue.

**Mechanism, which is what makes this a regime finding rather than a different fixture.**
With `i` = unit every treated (unit, period) cell is a singleton carrying its own dummy, so
its residual is exactly zero and the CRVE omits the treated-side variance entirely.  What is
left is the control-side (time-FE) variance, which shrinks as controls are added, while
SD(b) stays pinned at roughly 1/sqrt(treated observations).  So SE/SD is governed by the
CONTROL-TO-TREATED RATIO, not by the cluster count, and it goes to zero as the control pool
grows.  Arms D and A are that prediction measured: holding the 20 treated units fixed and
growing the control pool from 40 to 280 moves SE/SD from 0.508 to 0.212 and coverage from
0.733 to 0.367.  Adding never-treated controls - free precision, and something an applied
researcher does without a second thought - makes the reported standard error arbitrarily
optimistic.

sthlp:251 states the understatement as "roughly 0.66 to 0.75 of the reference on the
diagnostic fixtures"; the project record pairs that with coverage 0.794-0.876 at G=40.  Arm C
(robust) lands in that band; arm A is three times worse on SE/SD and a 95% interval covers
barely a third of the time.  Nothing is printed on the `att` rows - the granular-`i` caution
that does print is about RUNTIME, and the SI-A1 caution fires only with `jointtests`.  The
help's first remedy, "a finer `cluster()`", does not exist when `i` is already the unit; and
per SAFETY-01 the remaining remedy, `bootstrap()`, is itself refused r(907) on the wide panel.

The phenomenon is ruled (known_bydesign S1) and I am not claiming it as new.  The
contribution is the magnitude and its driver in the few-treated / many-controls regime, which
S7 names as open ("the effect at a different cluster count").  Confidence `probable`: the
gradient is unambiguous, but the comparison is against a differently-constructed fixture.

---
## NEW (mine) -- SKEP-02: the printed `relyr` column rounds fractional horizons together
**Repro:** `REPRO_SAFETY03_fracOnset.do` (leg L3; see `r3_fracOnset.log` for the table).

Under the `safety(off)` fractional-onset bypass, `e(wdidmainresults)` stores relyr
-1.5 / -0.5 / +0.5 correctly, but the printed table renders that column as `-2 / -0 / 0`:
two distinct horizons print as the same value, and one prints as "-0".  Cosmetic, and only
reachable behind a bypass the default refuses.

---
## CLEAN -- probed hard, found nothing
- **`safety(off)` does not suppress the disclosure of a sample change.**  On an unbalanced
  panel where the never-treated units exit after 2005, wooldid drops 2006-07 for want of a
  control.  Default safety prints the caution; `safety(off)` does not - but both print
  `N = 360 (40 obs dropped from initial sample)`, both post `e(N_dropped) = 40`, and the att
  is bit-identical.  sthlp's "It has never suppressed a disclosure about what was estimated"
  HOLDS.  (`r4_own.do` block E.)
- **No tempvar leak after a mid-run r(907).**  `safety(off) customvce()` at a tight ceiling
  dies inside estimation; the dataset varlist is unchanged afterwards and `est dir` is clean.
  (`r4_own.do` block F, `F_LEAKED=NO`.)
- **A treated cohort coded 0 does not collide with the never-treated sentinel.**  Bit-identical
  att/att_i/ame/ame_i against the 1-indexed recoding, both under `safety(off)` on the binary
  path and under `contreatwithin`.  (`r2`, `REPRO_SAFETY02`.)
- **The F13 safety(off) bypass moves no number on the default spec.**  An onset recorded only
  from the onset period onward gives an att bit-identical to the correctly filled data, because
  with the default pooled reference those rows land in the reference pool either way.
  (`r5_f13_and_esf.do` block G.)
