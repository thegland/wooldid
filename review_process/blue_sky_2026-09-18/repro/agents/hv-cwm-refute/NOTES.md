# hv-cwm-refute — SKEPTIC notes (customweightmultiplier x att_it/att_i/att_itime)

STATE: verdicts settled. Every number below comes from MY OWN panels and MY OWN do-files;
the prober's repro was read but never run.

## Files (all in this directory)
- `r1_core.do/.log`               arms A (no cwm) / B (cwm=pop) / S (cwm=3) / D (cwm=popg), gteffects oracle
- `r2_pre_and_missing.do/.log`    ES arms E1-E4 (the secondary pre-row claim) + M0-M4 (missing cwm)
- `r3_note_honesty.do/.log`       impvalues / subgroups / contreat-ame x cwm
- `r4_impvalues_cwm.do/.log`      impvalues with a cwm that varies WITHIN cluster (kills my own R3 lead)
- `r6_engines_carry_cwm.do/.log`  bootstrap(standard) SEs with and without cwm
- `repro_cwm02_bayes_gate.do/.log`     SELF-CONTAINED repro for the hv-cwm-02 verdict
- `repro_cwm_missing_error.do/.log`    SELF-CONTAINED repro for my new finding
- `panel_r*.dta`, `oracle_r1.dta`      generated panels (no project data, no project suite)

## hv-cwm-01 -> PARTIAL (arithmetic right, "silent"/"undocumented" wrong, severity low)

ARITHMETIC CONFIRMED on my data (r1_core.log):
  ARM A no cwm:    att 2.2656454  att_it 2.1660737  att_i 2.0400197  att_itime 2.0290046
  ARM B cwm(pop):  att 3.2087182  att_it 3.0353185  att_i 2.9297245  att_itime 2.8995162
  w_att_it in ARM A: min=max=1.  In ARM B: 1..40, matching the raw-data oracle
  sum_cell(sw*pop)/sum_cell(sw) to maxreldif 3.7397e-16.  ARM D (cwm=popg, varying within
  cohort over time): oracle maxreldif 2.0947e-16.  So the cell weight is exactly
  sum_cell(w*v)/sum_cell(w), as claimed.

THREE PARTS OF THE CLAIM THAT DO NOT HOLD:

1. "No option combination recovers the documented estimand in the same run" — FALSE.
   In ARM B (cwm in force) the simple unweighted mean of e(gteffects)'s `estimate` column
   over relyr>=0 is
       +1.1541e74632e42X+001  =  2.1660737126
   which is BIT-IDENTICAL to ARM A's reported att_it (+1.1541e74632e42X+001).  One option,
   `gteffects`, already in the prober's own arms, returns the documented estimand exactly,
   in the same run.

2. "No message names the substitution" — FALSE on both channels.
   - Runtime: "Note: all Estimates Shown Apply pop Weights Post-Estimation" fires and says
     ALL estimates shown.  It is accurate, not misleading: the user asked for pop weights
     and all four estimands got them.
   - Help: sthlp:1107 says "the weights NORMALLY USED to produce the requested type of
     treatment effect (ATT, ATT_i, etc.) will be multiplied by this variable", naming ATT_i
     explicitly.  That sentence is the code verbatim: the normal att_it weight is w/W_it
     (ado:1644 without the multiplier) and with it the code writes w*v/W_it.  sthlp:1107 is
     an exact, not an approximate, description of the mechanism, and the entry goes on to
     say the point is to allow "certain bespoke types of treatment effects".

3. "att and att_it now differ by 0.07%, so the heterogeneity contrast reads insensitive" —
   an artifact of the prober's panel, not a property of the option.  On my panel the gap
   WIDENS under cwm: 4.6% (2.2656 vs 2.1661) with no cwm, 5.7% (3.2087 vs 3.0353) with it.

WHAT SURVIVES: sthlp:218 (att) carries the "multiplied by the customweightmultiplier()
variable if one is specified" caveat; sthlp:220-221/226-227/232-234 (att_it/att_i/att_itime)
state their estimand with no such caveat and no cross-reference.  That asymmetry is real and
worth a one-clause fix.  It is a missing cross-reference in three option entries, not a
silent substitution: severity LOW.

SCALE INVARIANCE (my own check, ARM S, cwm = a constant 3): all four estimands and every SE
are HEX-IDENTICAL to ARM A.  w_att_it = 3 in every cell.  So the multiplier is exactly a
no-op when it carries no information, as it must be.

## The secondary claim (att_i pre row vs main row) -> REFUTED as cwm-specific

r2_pre_and_missing.log, event study esprelength(3) espostlength(3):

  E1  no cwm, time-constant sample weight
        w_att_i  per-cohort POST/PRE shares:  .375/.375   .375/.375   .250/.250   reldif 0
        w_att    per-cohort POST/PRE shares:  .510/.312   .358/.403   .132/.285
                                              reldif .1506, .0324, .1185
  E4  no cwm, TIME-VARYING sample weight
        w_att    reldif .1854, .0118, .1568
        w_att_i  reldif 0, 0, 0
  E2  cwm(pop), time-constant multiplier:     w_att_i reldif 0, 0, 0
  E3  cwm(popg), time-varying multiplier:     w_att_i reldif .0248, .0012, .0275

So "the pretreatment row weights the cohorts differently from the main row" is ALREADY TRUE
of the DEFAULT `att` estimand with NO multiplier and even a constant sample weight — 15.1%
for the first cohort in E1, larger than the 11.6% the prober reports for att_i.  The pre and
post blocks cover different periods, so any weight not constant within a cohort-block gives
them different cross-cohort shares.  cwm does not create the phenomenon; a TIME-VARYING cwm
merely moves att_i into the class `att` was always in.  Not a defect.

## hv-cwm-02 -> CONFIRMED (independently)

repro_cwm02_bayes_gate.log, multiplier varying WITHIN the (i,t) cell:
  K1 [aw=sw*epop] (the SAMPLE-WEIGHT route the Dirichlet draw takes, ado:9449):
       w_att_it min 0.99999999999999989  max 1.00000000000000000   -> cancels
  K2 customweightmultiplier(epop) (the NUMERATOR-only route, ado:1644):
       w_att_it min 0.71891796227648885  max 51.94960583131855003  -> does not cancel
  K3 bootstrap(bayes) + att_it + cwm(pop): r(198), the gate fires with its named message
  K4 bootstrap(bayes) + att only + cwm(pop): r(0)
The gate's comment (ado:6437-6439) reasons about the Dirichlet draw, which travels as the
sample weight; cwm travels as a numerator.  No inconsistency.  The recon's dichotomy is false.

## MY OWN findings

NEW-1 (filed): a missing customweightmultiplier() removes the row from the ESTIMATION SCREEN
(ado:704), i.e. from the regression, not only from the margins weights — while sthlp:218 says
the multiplier "enters the margins weights only and is deliberately not folded into the
regression weight".  Measured (r2, M-arms):
   M0 cwm=pop, nothing missing:   e(N)=560, underlying regression e(N)=560
   M1 cwm missing on 148 rows:    e(N)=412, underlying regression e(N)=412; identical to M3
   M2 the same rows set to ZERO:  e(N)=560, regression 560, DIFFERENT estimates
   M3 the same rows deleted:      e(N)=412, estimates identical to M1
The drop count IS printed ("N = 412 (148 obs dropped from initial sample)"), which is why
this is low and not medium.  The sharp edge is the error case: when the multiplier is an
administrative variable that exists only for programme participants, EVERY control row is
marked out and the run stops with
   "Error: you have no never-treated control group and no not-yet-treated control group
    whose treatment falls outside the sample window."
on a dataset with NINE never-treated counties and a strictly positive multiplier
(min 240, max 640, zero negative rows).  Nothing names customweightmultiplier(); both
remedies the message offers are wrong, and the actual remedy (fill the multiplier in with 0
on control rows) is named nowhere.  repro: repro_cwm_missing_error.do, arms A-D.

## Probed hard, found NOTHING (negative results worth recording)

- impvalues x cwm.  r3 showed BIT-IDENTICAL estimates with and without cwm under impvalues,
  which looked like the note printing a falsehood.  r4 killed my own lead: in r3 the
  multiplier was constant within the cluster, where an equal-weight cross-cluster
  aggregation makes it cancel arithmetically.  With the multiplier varying WITHIN cluster
  (r4: 6 states x 4 counties, pop 50/20/4/1 inside each state, cluster(state)) cwm moves the
  impvalues numbers properly: att 2.6825278000 -> 3.0889815582, att_it 2.7164835512 ->
  3.1155623408.  Also moves the impvalues event-study rows (post0 3.5082549 -> 4.0039385).
  CLEAN.
- contreat() ame family x cwm.  cwm moves BOTH families, not just att:
  att 2.6903257834 -> 3.5593013358 and ame 0.3843827889 -> 0.4448615390,
  ame_it 0.5169106248 -> 0.8176591677.  No att-only asymmetry.  CLEAN.
- subgroups() x cwm.  The main rows are bit-identical to the no-subgroup cwm run
  (+1.c7972fa8ebdf2X+001), and the per-subgroup note fires on every subgroup block.  CLEAN.
- bootstrap(standard) x cwm: the replications DO carry the multiplier.  se(att_it)
  .13647570134077 without cwm vs .17354532636481 with it, same rseed(4242), same reps.
  The inference is computed for the estimand that is reported, not for the unweighted one.
  CLEAN.  (r6_engines_carry_cwm.log)
- Constant multiplier: exact hex no-op on all four estimands and all SEs.  CLEAN.

## NOT filed (measured, judged below the bar)
- safety(off) is advertised in the negative-cwm refusal as the bypass.  With a few negative
  rows it does bypass (r(0), arm F of repro_cwm_missing_error).  With many negative rows, on
  an earlier panel of mine, it died r(402) "negative weights encountered" with no wooldid
  context.  Negative multipliers are not realistic input and ado:5638-5644 already calls the
  signed case caveat emptor; recorded here, not filed.
- permute() x cwm was attempted in r6 and both arms errored r(198) (permute() requires
  cluster()); showse then printed the STALE e() of the previous run, so the two identical
  "P0"/"P1" blocks in r6_engines_carry_cwm.log are a harness artifact and prove nothing.
  Closed separately in r7.

## r7 (permute x cwm) — CLEAN
`cluster(county) permute(shuffle) permreps(80) rseed(4242)`, with and without cwm(pop):
p_perm is 2/81 = .024691358024691 in both arms (the add-one floor — no replication beyond
the observed statistic either way, so the p is uninformative here), but the permutation
SHIFT INTERVALS move: [1.6901723179182, 3.9089828591455] without cwm vs
[0.99426003996833, 3.4422436392532] with it for att, and likewise for att_it.  The
permutation distribution therefore does see the multiplier.  CLEAN.
