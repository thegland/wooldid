# SKEPTIC VERDICTS — hv-safety-case-1 / hv-ctw-float / hv-ctw-ipos
All evidence below is from MY OWN do-files and data, in this directory. The prober's
do-files were read but never run.

---------------------------------------------------------------------------
## hv-safety-case-1 — CONFIRMED (severity low, confidence certain)
repro: r1_case.do / r1_case.log

An unrecognised casing of safety() is NEVER silently accepted. Measured:
  safety(off) rc0 | safety(on) rc0 | safety() rc0 (runs as ON) | safety("off") rc0
  safety(OFF)/Off/oFF/ON/On -> rc198 "Error: safety() must be on or off."
Source reason (mine): every safety-keyed site between ado:4568 and the domain screen at
ado:5473 tests `!= "off"`, and the one `== "off"` site (ado:4831) is display text only.
So the prober's own refutation of the recon claim is right.

The residual defect reproduces exactly as described:
  BITE1 `att_i unconditionalse safety(OFF)` -> r(198) "unconditional standard errors may
    be unreliable with non-att estimand types. / Set option safety(off) to disregard
    this error." -- the casing is never mentioned.
  BITE2 `[aw=wt] att unconditionalse robust safety(OFF)` -> prints "Set option safety(off)
    to suppress this caution", then dies "Error: safety() must be on or off."

CORRECTION I ADD (the prober understated it, though not enough to move severity):
the advice is not merely redundant, it is a DEAD END. Typing the correct `safety(off)` on
the BITE1 line produces a DIFFERENT error -- "Error: the option unconditionalse requires
specifying either robust or clustered standard errors." (r1_case.log, LEG_bite1ok). So a
user who follows the instruction verbatim lands on a second refusal and still never learns
that the casing was wrong.

---------------------------------------------------------------------------
## hv-ctw-float — CONFIRMED (severity medium, confidence certain)
repro: r2_float.do / r2_float.log ; number consequence in r6_extra.do / r6_extra.log

60 units, `gen float id = 20000000 + _n` -> 31 distinct ids out of 60 (r2_float.log:76).
  L1 contreat only                  -> caution PRINTS (log:132), rc0
  L2 + contreatwithin               -> NO caution, rc0
  L3 + contreatwithin + safety(on)  -> NO caution, rc0
  L5 binary, no contreatwithin      -> caution PRINTS (log:356), rc0
e(safety), e(wdid_safety) absent; e(cmdline) echoes the user's own line.

NUMBER CONSEQUENCE, measured by me (r6_extra.log:165/236, same panel, same options,
identifiers float vs double):
  float  main:att = +1.126646badd3d4X+000 (1.0719)   main:ame = +1.f73a67f9513a3X-001
  double main:att = +1.36e4d075fbf0eX+000 (1.2145)   main:ame = +1.9ed711fa4f992X-001
i.e. 13% on the ATT and 21% on the AME, with nothing said. The caution exists for exactly
this and is the only signal the user would ever get.

REGRESSION PROVENANCE, checked by me rather than taken on trust:
  git log -S'16777216' -- wooldidinstall/wooldid.ado  ->  a651e6f, 2026-09-06 (BS-34)
  git show e26b3dd^:wooldidinstall/wooldid.ado -> i-check 4488, float check 4503,
    `local safety off` at 4657, i.e. BELOW both.  e26b3dd is 2026-09-09.
So the caution DID work when it shipped and stopped working three days later. Confirmed.

---------------------------------------------------------------------------
## hv-ctw-ipos — CONFIRMED, SEVERITY RAISED medium -> HIGH (confidence certain)
repro: r3_ipos.do (the guard) and r5_izero_strong.do (the number)

Guard legs reproduce exactly (r3_ipos.log):
  no contreatwithin: i = 0.5 / -3 / 0 all r(411) "Error: id must consist of positive,
    non-zero integers, but does not"
  + contreatwithin : 0.5 -> r(452) "__000005: factor variables may not contain noninteger
    values"; -3 -> r(452) negative; 0 -> rc 0 with nothing said
  binary + explicit safety(off): byte-for-byte the same column as contreatwithin.

THE CORRECTION. The finding's "HONEST QUALIFIER" -- "I did not find a design where it
moves the number" -- is wrong. Such a design exists and is ordinary.

r5_izero_strong.do: 20 cohorts x 3 units x 8 periods, onset t=5 (cohorts 1-7) and t=6
(8-14), cohorts 15-20 never treated, a distinctive pre-trend in cohort 1 only.
Leg A labels cohort 1 as identifier 0; leg B adds 1 to every identifier. SAME PANEL,
SAME PARTITION, different labels, so every reported number must agree.

  binary + safety(off) (r5 log 221-302 vs 653-734):
    pretreatment:att_i   A +1.467a7e6eb4400X-006 (0.01993)   B +1.68c9f81132080X-005 (0.04404)
    pretreatment:att_itime  same pair.  se/t/p/CI all move with them.
    -> the pre-treatment att_i estimate is off by a factor of 2.2.
  contreatwithin + contreat(), NO safety() option typed at all (r5 log 425-578 vs 857-1010):
    pretreatment:att_i   A -1.e5e49fbd85000X-004   B -1.c475586cc1ce0X-004
    pretreatment:ame_i   A -1.8ff7eb0f66b71X-001   B -1.8e4c52036407dX-001
    (att_itime and ame_itime move identically; se/t/p/CI too)
  Same panel with NO relaxation: leg A r(411), leg B rc 0 -- so the guard is the only
  thing standing between an ordinary user and these numbers.
Every other row -- main:att, main:att_i, main:att_it, main:att_itime, pretreatment:att,
pretreatment:att_it and every ES_* row -- is bit-identical. The regression is untouched;
only the att_i / att_itime PRE-period weighting is wrong.

MECHANISM (mine, and it predicts exactly that row set):
  ado:928-929  `gen double _i = i if _t>0 & _t<. & ifstatement` ; `replace _i = 0 if _i==.`
               -> a cohort labelled 0 has _i == 0 on EVERY row.
  ado:1602     `replace posttreatment = 2 if t<ttre & ttre!=. & _i!=0 & _t!=0 & esamp==1`
               -> that cohort's pre-window rows are never promoted from 0 to 2.
  ado:1620     `gegen w_per_i = total(w*esamp), by(posttreatment i)`
               -> its pre-window rows share a denominator with its reference rows, so
                  w_att_i = w/w_per_i sums to LESS THAN ONE over the pre-window and the
                  cohort is silently downweighted in every pre-treatment att_i row.
  Post-treatment rows are posttreatment==1 for everyone, which is why main:* is clean.

So the consequence is not "a named r(411) becomes a raw r(452)"; that is the least of it.
The consequence is a SILENTLY WRONG pretreatment att_i / att_itime NUMBER on the i==0 shape, and the
pre-treatment row is the one people read to decide whether the design is credible.

---------------------------------------------------------------------------
## NEW (mine) — sk-izero-preweight: a cohort identified as 0 is silently downweighted
   in every pre-treatment att_i / att_itime / ame_i / ame_itime row
repro: r5_izero_strong.do (see numbers above). severity HIGH, confidence certain.
Broader than hv-ctw-ipos: it fires under the PLAIN DOCUMENTED safety(off) bypass as well
as under contreatwithin. sthlp:1208 lists "the cohort identifier not being a positive
integer" under "refuses, and safety(off) proceeds" -- "proceeds" is all it says. Nothing
warns that proceeding yields a label-dependent pre-treatment estimate. A guard whose only
documented cost is "you are on your own" should not be the sole barrier to a wrong number;
either the weight classifier should use the user's own `i' (or a treated flag) rather than
`_i != 0', or the refusal should be unconditional.

---------------------------------------------------------------------------
## NEW (mine) — sk-ctw-bare198: contreatwithin + contreat() + esfixedbaseperiod on a panel
   with no usable pre-period dies with a BARE r(198) and no message
repro: r7_bare198.do / r7_bare198.log. severity medium, confidence certain.
40 units, treated units have ttre == 1 (no pre-period) -- an extremely common user mistake.
  default safety           -> r(459) "Error: cannot find sufficient pre- and post-treatment
                              periods for the treated units." (ado:908)
  + contreatwithin         -> r(198), output is literally "invalid syntax", nothing else
  + safety(off) (no ctw)   -> same bare r(198)  [so the mechanism is the skipped ado:908
                              assert, and contreatwithin merely reaches it unasked]
  verbose shows the run gets all the way through reghdfe and dies in the margins stage.
The POOLED (non-esfixedbaseperiod) path keeps its named r(459) under contreatwithin, so
this is specific to the esfixedbaseperiod branch.

---------------------------------------------------------------------------
## ADDENDUM — which leg is WRONG, not merely different
The r5 panel is balanced (3 units per cohort, equal pre-window widths, unit weights), so
pretreatment:att_i must EQUAL pretreatment:att. Leg B satisfies that identity to the last
hex digit; leg A does not.
  leg B r5 log 660/662: pretreatment:att and pretreatment:att_i both +1.68c9f81132080X-005
                        (se both +1.bb8109559ea4fX-003)
  leg A r5 log 228/230: pretreatment:att +1.68c9f81132080X-005 but
                        pretreatment:att_i +1.467a7e6eb4400X-006 (se +1.bb31e3fcd5953X-003)
  0.04404 (correct) vs 0.01993 (reported) -- 2.2x too small.
Same identity for ame vs ame_i in the contreatwithin legs (B 872/874 equal, A 440/442 not).
Leg A_ctw printed no caution and no error of any kind.
Scope note: only the `pretreatment:` rows move. Every ES_* row is bit-identical because the
event-study horizons substitute att_it weights (ado:1657).
