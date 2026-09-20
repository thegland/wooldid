# SKEPTIC NOTES — refuting hv-safety-case-1 / hv-ctw-float / hv-ctw-ipos

Build asserted in every log: `which wooldid` -> scratchpad/bluesky/build/wooldid.ado.
All data generated inline in each do-file. Project test suite NOT touched.

## Source facts I verified myself (grep/sed on the pinned build)
- `wooldid_validate` spans ado:4568-5522 and is called FIRST thing by the worker (ado:568).
  It c_locals `safety` AND `safetyinput` back to the caller (ado:5509-5511), so the worker,
  wooldid_reltime etc. all see the REWRITTEN working value.
- ado:4828-4829  `local safetyuser `safety'` / `if contreatwithin & contreat != "" local safety off`
- ado:4959 (i positive-integer refusal) and ado:4973 (BS-34 float caution) both gate on
  `"`safety'" != "off"` -> both sit BELOW the rewrite.
- ado:5473 the domain screen `!inlist("`safetyinput'","on","off")` is the LAST block of validate.
- Every other safety-keyed site inside 4568..5473 tests `!= "off"` (conservative) except ado:4831,
  which is display text only.
- Regression provenance CHECKED MYSELF: `git log -S'16777216' -- wooldidinstall/wooldid.ado`
  -> a651e6f, 2026-09-06 (BS-34). At `e26b3dd^` the i-check is at 4488, the float check at 4503,
  and `local safety off` at 4657 — BELOW both. e26b3dd is 2026-09-09. Regression framing HOLDS.

## r1_case.do / r1_case.log — hv-safety-case-1
safety(off) rc0; safety(on) rc0; safety() rc0 (runs as ON); safety("off") rc0 (syntax strips quotes).
safety(OFF)/Off/oFF/ON/On all rc198 "Error: safety() must be on or off." NEVER silently accepted.
BITE1: `att_i unconditionalse safety(OFF)` -> rc198 with "unconditional standard errors may be
unreliable with non-att estimand types. / Set option safety(off) to disregard this error."
  NEW DETAIL the prober missed: with CORRECT safety(off) the same line gives a DIFFERENT error
  ("unconditionalse requires robust or clustered SEs"), so the advice the user is given is not
  merely redundant, it is a dead end.
BITE2: `[aw=wt], att unconditionalse robust safety(OFF)` prints the aweight caution
("Set option safety(off) to suppress this caution") and then dies "safety() must be on or off."
=> CONFIRMED, low.

## r2_float.do / r2_float.log — hv-ctw-float
float id = 20000000+_n -> 31 distinct ids out of 60 (log line 76).
L1 contreat, no ctw            -> caution PRINTS (log 132), rc0
L2 + contreatwithin            -> NO caution, rc0
L3 + contreatwithin safety(on) -> NO caution, rc0
L5 binary, no ctw              -> caution PRINTS (log 356), rc0
e(safety)/e(wdid_safety) absent; e(cmdline) echoes the user's own line. => CONFIRMED, medium.

## r3_ipos.do / r3_ipos.log — hv-ctw-ipos
noctw: 0.5 / -3 / 0 all r(411) "Error: id must consist of positive, non-zero integers, but does not"
ctw  : 0.5 -> r(452) noninteger; -3 -> r(452) negative; 0 -> rc 0, NOTHING SAID
safety(off) binary: identical to the ctw column. => CONFIRMED, and the r(452) shapes are exactly
what the DOCUMENTED bypass already produces, so the novelty is only that ctw triggers it unasked.

## r4_izero.do / r4_izero.log — MY OWN EXTENSION (the prober said they found no such design)
Same panel, leg A: treated unit's id = 0; leg B: every id shifted +1 (identical partition).
Binary + safety(off), att att_i att_itime att_it esprelength(3) espostlength(2):
e(wdidmainresults) differs ONLY in pretreatment:att_i and pretreatment:att_itime:
  A -1.a60a5c2483380X-004   B -1.aa91b1d1c7a80X-004   (se, t, p differ too)
In the balanced leg B, pretreatment:att_i EQUALS pretreatment:att, as a balanced panel requires.
In leg A it does not. Leg A is the wrong one.
MECHANISM (mine): ado:1602 `replace posttreatment = 2 if t<ttre & ttre!=. & `_i'!=0 & `_t'!=0`.
`_i' is 0 for EVERY row of a cohort labelled 0 (ado:928-929), so that cohort's pre-window rows
are never promoted from 0 to 2; ado:1620 then builds w_per_i by(posttreatment i), pooling that
cohort's pre-window AND reference rows into one denominator, so its pre-period att_i weight sums
to < 1 and the cohort is silently DOWNWEIGHTED in every pre-treatment att_i/att_itime row.
ctw legs in r4 died r(459) for an unrelated design reason (byCohort baseline dose slope) -> r5.

## r5_izero_strong.do / .log — THE DECISIVE LEG
20 cohorts x 3 units x 8 periods, onset t=5 (cohorts 1-7) / t=6 (8-14), 15-20 never treated,
a distinctive pre-trend planted in cohort 1 only. Leg A labels cohort 1 as id 0; leg B adds
1 to every id. Same panel, same partition.
  LEG_A_guard_RC=411  LEG_B_guard_RC=0   (default safety refuses A -- the guard works)
  LEG_A_bin_RC=0      LEG_B_bin_RC=0     (safety(off))
  LEG_A_ctw_RC=0      LEG_B_ctw_RC=0     (contreatwithin+contreat, NO safety() typed)
A_ctw printed NO caution and NO error (grep of log 303-423 matches only the results header).

THE IDENTITY THAT PINS WHICH LEG IS WRONG: the panel is balanced, so pretreatment:att_i
must equal pretreatment:att.
  leg B  line 660 pretreatment:att     +1.68c9f81132080X-005   +1.bb8109559ea4fX-003
         line 662 pretreatment:att_i   +1.68c9f81132080X-005   +1.bb8109559ea4fX-003   EQUAL
  leg A  line 228 pretreatment:att     +1.68c9f81132080X-005   +1.bb8109559ea4fX-003
         line 230 pretreatment:att_i   +1.467a7e6eb4400X-006   +1.bb31e3fcd5953X-003   NOT EQUAL
0.04404 vs 0.01993 -- leg A's att_i pre-treatment estimate is 2.2x too small.
Same identity holds for ame vs ame_i in the ctw legs (B 872/874 equal, A 440/442 not).
Every other row (main:*, pretreatment:att, pretreatment:att_it, every ES_*) is bit-identical
between A and B: ES rows substitute att_it weights (ado:1657 comment), so they are clean.

## r6_extra.do / .log
(a) float vs double identifiers under contreatwithin, same options:
    float  main:att +1.126646badd3d4X+000 (1.0719)  main:ame +1.f73a67f9513a3X-001
    double main:att +1.36e4d075fbf0eX+000 (1.2145)  main:ame +1.9ed711fa4f992X-001
    13% / 21% apart with the BS-34 caution suppressed.
(b) treated units with no pre-period: esfixedbaseperiod no-ctw r(459) named;
    +contreatwithin r(198). Pooled path keeps its named r(459) both ways.

## r7_bare198.do / .log
  1 default                r(459) "cannot find sufficient pre- and post-treatment periods"
  2 +contreatwithin        r(198), output is literally "invalid syntax" and nothing else
  3 safety(off), no ctw    same bare r(198)  -> mechanism is the skipped ado:908 assert
  4 binary + safety(off)   same bare r(198)
  5 ctw + verbose          reghdfe RUNS, then "invalid syntax" in the margins stage
So contreatwithin merely reaches a general safety(off) cost unasked -- but it burns a full
reghdfe first and the user is told nothing.

## Build asserted: all 7 logs echo scratchpad/bluesky/build/wooldid.ado from `which wooldid`.
