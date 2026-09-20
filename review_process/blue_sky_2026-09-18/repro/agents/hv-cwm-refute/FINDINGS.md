# hv-cwm-refute — FINDINGS (skeptic pass over hv-cwm-01 / hv-cwm-02)

Build: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build (adopath ++, `which wooldid` asserted in every log).
All data generated in-file. No project test suite, no project .dta, no repo writes.
The prober's do-file was read, never run: every number below is from my own panels.

---

## VERDICT hv-cwm-01 — PARTIAL (severity revised medium -> LOW)

CONFIRMED: the arithmetic. With `customweightmultiplier(v)` the estimand weight for a
treatment cell is `sum_cell(w*v)/sum_cell(w)` (ado:1644 numerator `w * v`, ado:1619
denominator `total(esamp*w)` with no `v`), so `att_it` is a v-weighted, not an unweighted,
mean of the cells.  My panel (14 counties, 10 years, 4 establishments/county-year,
pop 1..40): att_it 2.1660737 without cwm, 3.0353185 with `customweightmultiplier(pop)`;
`w_att_it` runs 1..40 and matches a raw-data oracle `sum(sw*pop)/sum(sw)` to 3.7e-16.
A time-varying multiplier matches its oracle to 2.1e-16.  A constant multiplier is an exact
hex no-op on all four estimands and all SEs.
  repro: r1_core.do (`ARM*`, `*_simplemean`, `*_ORACLE_maxreldif_cellw`)

REFUTED: three load-bearing parts of the claim.

1. "No option combination recovers the documented estimand in the same run."  It does.
   In the cwm arm, the simple unweighted mean of `e(gteffects)`'s `estimate` column over
   `relyr >= 0` is `+1.1541e74632e42X+001` = 2.1660737126, BIT-IDENTICAL to the no-cwm run's
   reported att_it.  `gteffects` — an option the prober himself ran — returns the documented
   estimand exactly, in the same run.

2. "No message names the substitution."  Two do.
   - Runtime: `Note: all Estimates Shown Apply pop Weights Post-Estimation` (ado:11369),
     which says ALL estimates shown.  It is true, not misleading.
   - Help: sthlp:1107 — "the weights normally used to produce the requested type of
     treatment effect (ATT, ATT_i, etc.) will be multiplied by this variable" — names ATT_i
     explicitly and describes ado:1630/1644 verbatim: the normal att_it weight is `w/W_it`,
     and with the multiplier the code writes `w*v/W_it`.  The entry's stated purpose,
     "certain bespoke types of treatment effects", is precisely the reported behaviour.
   The claim's "three entries that cannot all be true at once" is ordinary help structure:
   the estimand entries describe the default, the option entry describes the modification.
   The code says the same thing about itself: `wooldid_ptp` rebuilds the att_it-family
   weight and comments (ado ~10404) "cell totals use the raw weight, without any custom
   weight multiplier - THE DOCUMENTED ASYMMETRY".  The exact mechanism the prober flagged as
   an accidental, unnamed substitution is called documented design at its second
   implementation site, and that site reproduces it deliberately so the two agree.

3. "att and att_it now differ by 0.07%, so the heterogeneity contrast reads insensitive."
   An artifact of the prober's panel.  On mine the gap WIDENS under cwm: 4.6% (2.2656 vs
   2.1661) without, 5.7% (3.2087 vs 3.0353) with.

WHAT SURVIVES, and the correction: sthlp:218 (att) carries the multiplier caveat;
sthlp:220-221 / 226-227 / 232-234 (att_it / att_i / att_itime) state their estimand with no
caveat and no cross-reference to customweightmultiplier().  That asymmetry inside the help
is real and is worth one clause in each of the three entries.  It is a MISSING
CROSS-REFERENCE, not a silent substitution: severity LOW.

SECONDARY CLAIM (the att_i pretreatment row weights cohorts differently from the main row)
— REFUTED as cwm-specific.  Event study, esprelength(3) espostlength(3):
  no cwm, constant sample weight: `w_att` per-cohort POST/PRE share reldif .1506 / .0324 /
    .1185 — the DEFAULT estimand already does this, with no multiplier anywhere, and by MORE
    than the 11.6% the prober attributes to cwm.
  no cwm, time-varying sample weight: `w_att` reldif .1854 / .0118 / .1568.
  `w_att_i` reldif is exactly 0 in both of those arms and in the constant-cwm arm, and
    .0248 / .0012 / .0275 only under a TIME-VARYING cwm.
The pre and post blocks cover different periods, so any weight that is not constant within a
cohort-block gives them different cross-cohort shares.  cwm does not create the phenomenon.
  repro: r2_pre_and_missing.do (arms E1-E4)

## VERDICT hv-cwm-02 — CONFIRMED (independently reproduced), severity low/none

My own arms, multiplier varying WITHIN the (i,t) cell:
  `[aw=sw*epop]` (the route `bootstrap(bayes)` takes: `repw = basew*drawv`, ado:9449) ->
     `w_att_it` min 0.99999999999999989, max 1.00000000000000000 — cancels exactly.
  `customweightmultiplier(epop)` (the numerator-only route, ado:1644) ->
     `w_att_it` min 0.71891796227648885, max 51.94960583131855003 — does not cancel.
  `bootstrap(bayes)` + `att_it` + cwm -> r(198), the gate fires with its own message.
  `bootstrap(bayes)` + `att` only + cwm -> r(0).
The gate comment (ado:6432-6439) reasons about the Dirichlet draw, which travels as the
sample weight; cwm travels as a numerator.  No inconsistency.  The recon's dichotomy is false.
  repro: repro_cwm02_bayes_gate.do

---

## NEW FINDING (mine, not seeded)

### hvS-cwm-miss-01 [low/certain/realistic] — a missing customweightmultiplier() removes the
row from the REGRESSION, and when the multiplier is participants-only the resulting error
misstates the user's data

CLAIM. sthlp:218 says the multiplier "enters the margins weights only and is deliberately not
folded into the regression weight".  ado:704 marks a row with a missing multiplier OUT of the
estimation screen, so the multiplier changes which observations the underlying regression
uses.  When the multiplier is an administrative variable that exists only for programme
participants — a normal shape for such data — every never-treated row is marked out and
wooldid stops with "Error: you have no never-treated control group and no not-yet-treated
control group whose treatment falls outside the sample window" on a dataset that plainly has
nine never-treated counties.  No message names customweightmultiplier(); neither remedy the
error offers addresses the cause; the remedy that works is named nowhere.

SCOPE.  The ado:704 markout is shared by every variable in that one `foreach` -- cluster(),
the controls, ccc_absorb(), contreat(), lattice(), interactivecontrols() -- so a treated-only
controls() variable produces the identical r(459).  For a regressor that is unremarkable:
missing regressors drop observations everywhere in Stata.  customweightmultiplier() is the
one member of that list whose help entry promises it is post-estimation only (sthlp:218,
"enters the margins weights only and is deliberately not folded into the regression weight"),
which is what makes the sample change surprising here and nowhere else.

OBSERVED (repro_cwm_missing_error.do, strictly positive multiplier, min 240 max 640, zero
negative rows, 9 never-treated counties, 270 never-treated rows):
  arm A  no multiplier                        -> r(0), e(N) = 600
  arm B  customweightmultiplier(progsize)     -> r(459), the false "no never-treated control
                                                 group" error, preceded only by
                                                 "Caution: estimability checks found time
                                                 periods exist in sample without a control
                                                 unit, causing them to be dropped."
  arm C  ... safety(off) (remedy 1 of 2)      -> r(0), but the run then proceeds on 231 of
                                                 600 observations off not-yet-treated
                                                 controls, with the header's
                                                 "369 obs dropped from initial sample" as
                                                 the only trace and the cause still unnamed.
                                                 Remedy 2, "modify the sample window
                                                 manually", does nothing for this at all.
  arm D  multiplier filled in with 0 on the   -> r(0), e(N) = 600, att 2.2629781
         control rows (the real remedy)
And on scattered missingness (r2_pre_and_missing.do, arms M0-M3):
  M0 nothing missing:            e(N) 560, underlying regression e(N) 560
  M1 multiplier missing, 148 rows: e(N) 412, underlying regression e(N) 412, and the four
     estimands bit-equal to M3 (the same rows physically deleted)
  M2 the same rows set to ZERO:  e(N) 560, regression 560, different estimates
The M1 header does print "N = 412 (148 obs dropped from initial sample)", which is why this
is LOW rather than medium; in arm B no header prints at all, so there the user gets no clue.

EXPECTED.  Either sthlp:218's "enters the margins weights only" should say that a missing
multiplier also removes the observation from the regression (sthlp:1107's "missing values
drop the observation" is the nearest thing and does not say from what), or the r(459) should
name the screening variable that emptied the control group.  A user told their data has no
never-treated control group, when it has nine such counties, will go and rebuild their panel.

  repro: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-cwm-refute/repro_cwm_missing_error.do

---

## NEGATIVE RESULTS (probed hard, clean)

- `impvalues` x cwm.  r3 showed BIT-IDENTICAL estimates with and without cwm and looked like
  the runtime note printing a falsehood.  r4 killed my own lead: in r3 the multiplier was
  constant within the cluster, where an equal-weight cross-cluster aggregation cancels it
  arithmetically.  With the multiplier varying WITHIN cluster (6 states x 4 counties,
  pop 50/20/4/1 inside each state, cluster(state)) cwm moves the impvalues numbers properly:
  att 2.6825278000 -> 3.0889815582, att_it 2.7164835512 -> 3.1155623408, and the impvalues
  event-study rows move too (post0 3.5082549 -> 4.0039385).
- `contreat()` ame family x cwm.  Both families move, not just att: att 2.6903257834 ->
  3.5593013358, ame 0.3843827889 -> 0.4448615390, ame_it 0.5169106248 -> 0.8176591677.
- `subgroups()` x cwm.  Main rows bit-identical to the no-subgroup cwm run; the per-subgroup
  note fires on every block.
- `bootstrap(standard)` x cwm.  The replications carry the multiplier: se(att_it)
  .13647570134077 without vs .17354532636481 with, same rseed and reps.  The inference
  belongs to the estimand that is reported.
- `permute(shuffle)` x cwm.  p_perm sits at the add-one floor 2/81 in both arms and is
  uninformative, but the shift intervals move ([1.690, 3.909] -> [0.994, 3.442] for att), so
  the permutation distribution does see the multiplier.
- Constant multiplier: exact hex no-op on all four estimands and every SE.

## RECORDED, NOT FILED
`safety(off)` is advertised in the negative-multiplier refusal as the bypass.  With a few
negative rows it does bypass (r(0)).  With many negative rows, on an earlier panel of mine,
it died r(402) "negative weights encountered" with no wooldid context.  Negative multipliers
are not realistic input and ado:5638-5644 already calls the signed case caveat emptor.
