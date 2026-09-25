# Econometric review of the help and the guide, and the oldsyntax removal, 2026-09-25

**Status line, dated: opened 2026-09-25. The documentation corrections are built (16f5392); the oldsyntax
removal is built (aeec885); every code, message-wording and number-moving item is HELD for the author's
ruling (section 5), and the deletion candidates await his choice (section 6).** The author asked on
2026-09-25 for (1) the oldsyntax option removed as a simplification, without running the suite; (2) opus
agents to review the help and related materials for candidate features to delete, to be put to him by
multiple choice; (3) separate econometrically oriented agents to review the same materials for
econometric errors, with corrections implemented; and (4) when the message-review page is next spun up,
an editable artifact of the main help page.

## 1. oldsyntax removed (aeec885)

The option renamed the event-study row for relative period 0 from post0 to contemp and was refused by
both replication engines. An opus agent removed it in a worktree: both syntax lines, the validator's
normalisation, the four row-naming branches, wooldid_sglbl's old() argument, the permctl and bootctl
refusals, the Mata label parser and covariance splice, the help (syntax table, two refusal lists, the
entry), the guide (feature list, refusal lists, two engine entries, the safety(off) inventory, the entry;
chapter 14 regenerated), test_bootstrap_gates legs 5f and m4 (no renumbering), the ext_cov10_oldsyntax
config and the suite README. Seven old-versus-new comparison runs were byte-identical in e(b), e(V), row
names, scalars, macros and printed output; a call naming the option stops with r(198). The suite was not
run, at the author's direction; the two suite legs were edited unrun. Notes: oldsyntax_removal_notes.md
in this directory. Left as they stand: test/suite/compare_oldnew.py still accepts contemp in its
event-study row pattern (dead, harmless); a comment in coverage_check.py says 85 options (it counted 83
before and 82 now); the untracked external bundle still defines the removed config. Cards m148 and m211
were deleted from the message-review record.

## 2. Method

Round 1 (wf_8ab761b9-a0e, 26 opus agents): six reviewers, one econometric lens each (identification and
estimands; analytic and bootstrap standard errors; permutation inference; continuous treatment; event
study, pre-trends, adversarial test and bands; poisson), read both documents in full and checked every
substantive claim against the ado, sorting each finding into document-misdescribes-code,
document-econometric-error, code-econometric-error or the author's ruled design choice, and searching
review_process/ for rulings. They filed 49 findings (9 high, 17 medium, 23 low). Each of the 16 high and
medium findings that the cap admitted went to an independent opus refuter told to refute by default and
to check the rulings: all 16 were confirmed. Round 2 (wf_2ff3fd6d-ccd, 88 opus agents) refuted the other
33 (25 confirmed, 4 the author's ruled choices, 4 refuted), then for every confirmed documentation
finding an opus editor drafted exact anchored edits from the refuter's corrected fix and an opus
verifier checked each edit against the code, the register and the rest of both documents (62 accepted,
46 accepted with a fix, 0 rejected, 9 further edits the verifiers found missing). Three deletion lenses
(user value, maintenance burden, redundancy) and a synthesiser produced the candidate list. Refuters ran
Stata probes where reading could not settle a claim; their probe files are under the session scratch
and their findings are quoted in round1_confirmed.txt and the two results files here.

## 3. Findings ledger

| finding | severity | kind | verdict | site | outcome | the text at issue |
|---|---|---|---|---|---|---|
| dose:CT-1 | high | CODE_ECONOMETRIC_ERROR | CONFIRMED | wooldid.ado:13536 | held: code/number-moving | if `amelive' == 0 & `amedrop' > 0 { |
| dose:CT-2 | high | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid.sthlp:793 | 3 doc edits applied, 2 items held | pieces without it are dropped for collinearity, the reported average marginal effect then  |
| dose:CT-3 | high | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:817 | 2 doc edits applied, 1 items held | When specified, {it:wooldid} will still present ATT-type average treatment effects estimat |
| ident:F1a | high | CODE_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:847 | held: code/number-moving | reported treatment effects will not vary as a direct function of these controls; the effec |
| ident:F1b | high | CODE_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:889 | held: code/number-moving | An alternate version of coarsecohortcontrols() that admits only continuous variables and p |
| ident:F4 | high | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:674 | 4 doc edits applied | The difference in SG level estimates returned will be the triple diff estimate. |
| perm:PERM-1 | high | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:383 | 5 doc edits applied, 4 items held | it protects only when the cluster-robust standard error is informative about the treated c |
| se:SEB-01 | high | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.ado:4949 | 5 doc edits applied, 3 items held | A cluster() strictly finer than the cohort-by-period cells, bootstrap() or permute() gives |
| dose:CT-4 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:800 | 6 doc edits applied | these different g() options set the reference effect of Contreat_jt relative to which Cont |
| dose:CT-5 | medium | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid_agent_guide.md:92 | 1 doc edits applied | the joint pre-trend and adversarial tests are not calibrated, the `ame` standard errors un |
| es:ES-1 | medium | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid.sthlp:963 | 3 doc edits applied | {it:ES_att_post#} (the event-study horizons, numbered outward from the reference period) |
| es:ES-2 | medium | DOC_ECONOMETRIC_ERROR | RULED | wooldid.sthlp:566 | no change | the line is chosen to fit the estimates as well as possible, a null that is adversarial to |
| es:ES-3 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:427 | 7 doc edits applied, 2 items held | / `simultaneousbands` / sup-t band over a family's horizons / two horizons with variance / |
| ident:F2 | medium | DOC_ECONOMETRIC_OMISSION | CONFIRMED | wooldid.sthlp:612 | 6 doc edits applied | If the user specifies the option {it:esfixedbaseperiod}, {it:wooldid} will switch to emplo |
| ident:F3 | medium | DOC_OMITS_CONDITION | CONFIRMED | wooldid.sthlp:1154 | 4 doc edits applied | {it:did_imputation} given that one uses the {it:att} estimand-type, uses a pooled referenc |
| perm:PERM-2 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:323 | 1 doc edits applied, 1 items held | Randomization inference assumes clusters are exchangeable in **assignment** and fails when |
| perm:PERM-3 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:444 | 2 doc edits applied, 2 items held | A caution prints when fewer than 6 clusters are ever treated, since the set of distinct as |
| pois:POIS-2 | medium | CODE_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:907 | held: code/number-moving | with {it:cluster()} at the level of {it:i} or coarser a treated cell lies inside a single  |
| se:SEB-02 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:201 | 2 doc edits applied | 3. **`robust`** is appropriate when errors are independent across observations, which is r |
| se:SEB-03 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:95 | 6 doc edits applied, 1 items held | put the unit identifier in `fe()`, and cluster at the unit level or at any level strictly  |
| se:SEB-04 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:326 | 3 doc edits applied | Intuitively, this ought to be true whenever samples are large enough that sampling variati |
| se:SEB-05 | medium | DOC_OMISSION | CONFIRMED | wooldid.sthlp:142 | 12 doc edits applied, 3 items held | The default estimates standard errors in the underlying regression and then obtains standa |
| se:SEB-06 | medium | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:1176 | 2 doc edits applied | {cmd:wooldid y i t ttre, cluster(i) att} |
| dose:CT-6 | low | DOC_INTERNAL_INCONSISTENCY | CONFIRMED | wooldid.sthlp:823 | 1 doc edits applied | draws the ame continuous treatment effects plots without confidence intervals; the att res |
| dose:CT-7 | low | DOC_ECONOMETRIC_ERROR | REFUTED | wooldid.sthlp:826 | no change | or a {it:cluster()} strictly finer than the cohort-by-period cells. |
| dose:CT-8 | low | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid.sthlp:709 | 2 doc edits applied | lowest quantile bin if their value of {it:contreat} is in the lowest {it:latticenum}-th pe |
| dose:CT-9 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:762 | 3 doc edits applied, 1 items held | If this option is specified when the always 0 condition is not met, {it:wooldid} will retu |
| es:ES-5 | low | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid.ado:11933 | 4 doc edits applied, 1 items held | See the Design conditions for every consumer of e(V_es_*) in the help, condition Two. |
| es:ES-6 | low | DISCLOSURE_GAP | REFUTED | wooldid.ado:16051 | no change | g = svsolve(M, z) |
| ident:F10 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:607 | 3 doc edits applied | large a value attenuates the estimated effects without any diagnostic, so prefer the small |
| ident:F5 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:622 | 3 doc edits applied | use of a pooled reference period is likely disadvantageous whenever the outcome variable i |
| ident:F6 | low | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid.sthlp:905 | 2 doc edits applied | taken over the post-treatment rows ({it:relyr} at or above 0), the pre-treatment rows, or  |
| ident:F7 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:69 | 1 doc edits applied, 1 items held | fixed effects in {opt fe()}, which a coarse {it:i} also needs when the panel is unbalanced |
| ident:F8 | low | DOC_ECONOMETRIC_ERROR | REFUTED | wooldid.sthlp:153 | no change | This can be useful for tests of differences in treatment efficacy across treatment arms |
| ident:F9 | low | DOC_IMPRECISION | CONFIRMED | wooldid.sthlp:1168 | 3 doc edits applied | The {it:gteffects} option returns the cohort-by-relative-period effects underlying {it:woo |
| perm:PERM-4 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid.sthlp:483 | 3 doc edits applied, 1 items held | Romano and Wolf max-/studentized/ step-down computed from the same replications, which con |
| perm:PERM-5 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:383 | 1 doc edits applied, 1 items held | It is not robust to the treated clusters being atypical in size or in error variance - rel |
| perm:PERM-7 | low | DOC_INTERNAL_INCONSISTENCY | CONFIRMED | wooldid_agent_guide.md:383 | 1 doc edits applied, 1 items held | With `contreat()` the dose stays with its unit and only onset timing is re-randomized (the |
| pois:POIS-1 | low | CODE_ERETURN_LABEL_ERROR | CONFIRMED | wooldid_agent_guide.md:889 | held: code/number-moving | *e(vce)* is *cluster*, *robust*, *unconditional*, *bootstrap*, the first word of `customvc |
| pois:POIS-4 | low | RUNTIME_MESSAGE_WORDING | CONFIRMED | wooldid.ado:1941 | 0 doc edits applied, 1 items held | Falling back to estimation with ppmlhdfe, but with any included time-trend controls, ccc_a |
| pois:POIS-5 | low | RULED_DESIGN_CHOICE | RULED | wooldid.sthlp:527 | no change | are presented in terms of differences in linear predictions (the poisson {it:predict(xb)}  |
| se:SEB-07 | low | DOC_ECONOMETRIC_ERROR | REFUTED | wooldid.sthlp:277 | no change | Each reported row uses a t distribution with degrees of freedom one less than the number o |
| se:SEB-08 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:283 | 1 doc edits applied, 2 items held | so a cluster's total mass stays proportional to its size in every replication - exactly as |
| se:SEB-09 | low | DOC_ECONOMETRIC_ERROR | RULED | wooldid.sthlp:364 | no change | Neither mode is a remedy for few treated clusters and no wild cluster bootstrap is offered |
| se:SEB-10 | low | DOC_ECONOMETRIC_ERROR | CONFIRMED | wooldid_agent_guide.md:359 | 1 doc edits applied | There the plain share stands, which is also what `bootstrap`'s own *estat bootstrap* repor |
| se:SEB-11 | low | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid_agent_guide.md:214 | 2 doc edits applied | so that a treated cohort-by-period cell contains several clusters. That third one is not c |
| se:SEB-12 | low | DOC_MISDESCRIBES_CODE | CONFIRMED | wooldid.sthlp:145 | 7 doc edits applied | P-values are then obtained from a t-distribution with the residual degrees of freedom of t |
| se:SEB-13 | low | DOC_RECIPE_OMISSION | CONFIRMED | wooldid.sthlp:1129 | 5 doc edits applied | ({it:matrix b = e(b_es_{it:type})}, {it:matrix V = e(V_es_{it:type})}, {it:ereturn post b  |
| perm:PERM-6 | none | DOC_ECONOMETRIC_ERROR | RULED | wooldid_agent_guide.md:383 | no change | These contrasts and the across-subgroup equality tests answer the sharp null of no effect  |

## 4. What changed in the documents (16f5392)

117 anchored edits, 107 applied by script with a live uniqueness check and 10 merged by hand where two
findings shared an anchor. No message text changed. Chapter 14 of the guide was not touched. Rendering
and register checks: every added help line balances its braces, no SMCL in the guide, no dash
characters added, file-wide densities unchanged. The help grew by about 1,600 words, which the author
may want trimmed in a later pass. The largest corrections, in reading order of the help:

- The default variance is documented (conventional OLS without cluster(), robust or customvce(); robust
  under ppmlhdfe and pweights; none allows for serial correlation), examples 1 to 10 are said to omit
  cluster() for brevity, and example 0 now clusters by unit with i as onset cohort (SEB-05, SEB-06).
- The cluster level is tied to the level of treatment assignment, in the cluster() entry and the remedy
  list (SEB-03); robust is said to carry the containment problem on one-row cells (SEB-02).
- The gteffects remedies no longer claim bootstrap() or permute() repair contained cells; only a finer
  cluster() does, and the analytic band is said to inherit the understatement (SEB-01, ES-3).
- The unconditional-standard-error gap is described as heterogeneity, not sample size (SEB-04).
- The pooled versus fixed reference choice is said to change the comparison group and the identifying
  assumption, with citations (F2); the did_imputation equivalence names the level of i (F3).
- The event-study horizon names are explained by relative period, the reference period having no
  column (ES-1); the anticipation over-shift sentence is corrected (F10).
- permtstat's protective condition is turned around: several clusters inside each treated cell, which
  under permute() needs permonsetcohorts; the exchangeability sentence follows the ruling (PERM-1,
  PERM-2, PERM-3, PERM-4, PERM-5, PERM-7).
- The triple-difference recipe uses att_it with i as onset cohort by group (F4).
- The lattice over-saturation paragraph says a single-valued piece borrows another bin's slope rather
  than being dropped, and separates the standard-error mechanism (CT-2); the ame's identifying
  condition is stated (CT-4); the finest-i trap's inference bullet is corrected (CT-5); the
  contreatwithin att is described as depending on the solver's normalisation where no untreated
  reference rows remain (CT-3); the lattice bin definition, the untreatedZero premise and the ccref band
  are corrected (CT-8, CT-9, CT-6).
- The e(V) recipe carries the degrees of freedom, the df note is corrected, the estat bootstrap claim is
  fixed, and the guide's mangled copy of the cluster-level paragraph is restored (SEB-13, SEB-12,
  SEB-10, SEB-11); the Bayesian bootstrap moment claim and the pretrendpower cross-reference are fixed
  (SEB-08, ES-5); the fe(unit) rationale and the gteffects comparability sentence are corrected (F7,
  F9); the single-relyr weight sentence is corrected (F6); the pooled-reference disadvantage is restated
  (F5).

## 5. Held for the author's ruling

Number-moving or message-changing items, none built. Each is put to him by multiple choice.

1. **coarsecohortcontrols() and ccc_absorb() (F1a, F1b, high).** Both build i._i#i._t#(vars) without
   centring, so every treated cell has its own slope on the uncentred control and the reported effect
   is the cell intercept at a control value of zero (the base category for an i. control). It moves
   when a constant is added to the control: probe with true ATT 2.00 gave -0.53 for c.x and -49.4 for
   c.x+100. Fix: enter the pooled controls once, and interact the cells with a copy centred at each
   treated cell's weighted mean (0 on reference rows); same column space, same fit, the intercept
   becomes the cell ATT, invariant to location, verified by hand (1.930 for x and x+100 against 1.930
   from did_imputation). Needs the poisson path and the r(3200) retry too. Options: build it (moves
   every coarsecohortcontrols()/ccc_absorb() estimate; pre-register the battery configs); document only
   (the effect is at control value zero, centre your controls); refuse uncentred controls; leave.
2. **Partial dropped dose slopes (CT-1, high).** When some treated cells lack within-cell dose variation
   and others do not, the dropped slopes enter the ame as zeros at full weight, so the ame is pulled
   toward zero by those cells' weight share (probe: ame .72 against 1.00, p near 0, only the generic
   dropped-term caution). Options: void the ame family as the all-dropped case does; renormalise over
   live cells and report the share left out; document only. The predicate must key on cells, not
   tokens (lattice runs drop dose terms legitimately).
3. **e(vce) label (POIS-1, SEB-05).** On the ppmlhdfe path and under pweights the variance is robust
   but e(vce) says oim or ols and e(vcetype) is blank. Options: post robust (label only, no number
   moves); leave.
4. **Poisson containment (POIS-2).** The one-row-cell containment screen has no ppmlhdfe arm, so a
   default poisson gteffects run on a fine i prints no caution and posts no e(gt_ncontained). Options:
   add the arm (one caution and e(gt_ncontained) appear on such runs); leave.
5. **permtstat caution trigger (PERM-1).** Widening from cluster() == i to "cluster() nests i without
   permonsetcohorts" would catch i = unit with cluster(state); it goes against the 2026-09-06 F049
   ruling and test_permute3 Q20. Options: widen; keep.
6. **Bootstrap resolution caution under bc and bca (from the 2026-09-24 docs check).** At reps() below
   40 the parse-time caution puts bootdisplay(bc) and (bca) on the 2/reps grid, which 27e277d and
   6adca1b took them off. Options: skip it under bc and bca; a bc/bca-specific sentence of his.
7. **New cautions proposed, none built:** a default caution when a lattice leaves a treated cell-by-bin
   piece without its own slope (CT-2); a caution naming the periods whose contreatwithin att rests on
   the solver's normalisation (CT-3); a caution when the analytic band is built with cluster() at or
   above i (ES-3); a premise screen when the dose is nonzero outside treated cells under untreatedZero
   (CT-9); a caution when the adversarial covariance is rank-deficient (ES-6, merging blue-sky
   hv-adv-ame-rank); the containment screen's iid case (SEB-05).
8. **Message wordings proposed on the review page** (cards m055, m056, m173, m174, m175, m026, m276,
   m034, m329; each proposal is in the card's suggested field with the reason): the gteffects
   containment and degeneracy cautions (bootstrap does not repair contained cells); the two permtstat
   cautions (several clusters per cell, not many observations); the few-treated-clusters permutation
   caution (three options); the si8ilev caution's fe() rationale; the Bayesian bootstrap moment clause;
   the rc-3200 fallback's second line (time trends never enter absorb()); the pretrendpower stale
   cross-reference; the verbose dropped-term line "If using a lattice, this is likely innocuous".
9. **Comment-only ado corrections that accompany those wordings** (PERM-1, PERM-2, PERM-3, SEB-08,
   PERM-7): built with whichever wording he picks.

## 6. Deletion candidates (for his choice; oldsyntax excluded, already gone)

| feature | burden | value | what it is | why delete | what breaks |
|---|---|---|---|---|---|
| cleanmatrices | high | low | Runs `clear mata` and `matrix drop _all` at the start and end of the call, as a recovery for sessions that reghdfe has left in a bad Mata st | The same recovery is two commands the user can type before the call, choosing what to lose, while the option also clears wooldid's own Mata kernels mid-run. That forced an inline-Mata screen (ado ~14763), cap-wrapped ker | A do-file naming cleanmatrices stops with r(198) and should type `clear mata` and `matrix drop _all` before the call instead; no estimate changes (a one-sentenc |
| bootalpha() (the size convention, and with it the whole option) | high | low | Chooses the bootstrap(bayes) Dirichlet convention: equal (the default, Rubin's cluster Bayesian bootstrap) or size (Gamma shape equal to clu | The help's own entry says size is 'for comparison against that command only; on panel data it understates the standard errors badly', so its only non-default value is wrong on the package's core data type and has produce | bootstrap(bayes) runs are unchanged, but a do-file naming bootalpha() stops with r(198), e(boot_alpha) is no longer posted, and the core config boot_bayes and r |
| suppresssgprimaryeffects | high | low | Under subgroups(), skips the per-subgroup estimates and reports only the pairwise subgroup contrasts and the equality test, as a speed optio | Since EFF3 the subgroup level margins come free from the pwcompare pass that makes the contrasts (ado 3393-3404, 3710-3720), so it now saves only the per-subgroup event-study, AME and escomp passes (measured 11.8 s to 6. | A do-file naming it stops with r(198); the same call without it returns a superset of the results, roughly 1.8x slower on subgroup runs with an event-study wind |
| ties | medium | low | Half-counts bootstrap replications exactly equal to the rep-0 comparator when forming the bias-correction constant z0, as the ties option of | The help already says it has no effect on continuous statistics, and a probe at G=5 (dta_binary, cluster = mod(i,5), reps 300) showed it is inert even in the one regime where ties should occur: the 8 identity draws land  | A do-file naming ties stops with r(198) and e(boot_ties) is no longer posted; no reported number moves (the battery config boot_bayes_sgjt drops the word from i |
| permside(left|right) | medium | low | Makes the printed permutation p-values one-sided (share at or below, or at or above, the observed statistic). | Every permute() run already posts both one-sided shares unconditionally, by coefficient name, as e(perm_p_lo) and e(perm_p_hi) (ado 9221-9227; the footer says so), and the interval stays two-sided under the option anyway | The printed p column, e(p) and so esttab can no longer carry a one-sided p (users read e(perm_p_lo)/e(perm_p_hi)), and the continuous-treatment plot cells lose  |
| permrawp | medium | low | Reports the plain absolute-value permutation share c/n without the add-one adjustment, the number ritest and Stata's permute prefix print. | It is permabsp without the +1: identical under permexact, and on sampled draws recoverable exactly as ((n+1)*p_absp - 1)/n with n from e(perm_n), while being the anti-conservative form and the only route by which a sampl | A do-file naming permrawp stops with r(198), and the ritest cross-check (newfeat/ritest_crosscheck.do), the mc_ri2 raw arm and the perm_rawp config must convert |
| permanalyticse | medium | low | Under permute(), takes the observed run out of fast mode so the analytic standard errors fill the se_analytic column and e(V), e(b_es_*) and | Only the observed run leaves fast mode (ado 6195, 6243), so every number it adds is exactly what the plain analytic run of the same specification gives, at the cost of one base fit, and permtstat already shows the analyt | Under default permute() the se_analytic column stays blank and no e(V)/e(V_es_*) is posted, so a user wanting both, or feeding honestdid, runs the analytic spec |
| altmethod_adversarial | high | medium | A modifier of adversarial that fits a line to the pre-period coefficients plus the reference and tests the post coefficients jointly against | It is a second statistic for an already niche diagnostic, threaded through every engine (a gapskip status all four call sites branch on, a separate Mata forecast branch, the permute and bootstrap statistic columns), and  | A do-file naming it stops with r(198) and its e(adversarial_*) results go (adversarial alone is unchanged), and the configs adv_altb and adv_r2001 with their ba |
| emptycellsoverride | low | low | Turns off wooldid's program-scoped `set emptycells drop`, so the underlying regression carries columns for empty interaction cells. | It changes no number (probe: e(b) and e(V) mreldif 0, 7.1 s vs 37.9 s, so about 5x slower rather than the 'roughly doubles' the help says), and the help itself says it is only of interest with verbose; a user wanting the | A do-file naming emptycellsoverride stops with r(198); nothing else changes. |
| Deprecated aliases e(perm_lb_*) and e(perm_ub_*) | low | low | Scalars carrying exactly the same two numbers as e(perm_bandlo_*) and e(perm_bandhi_*), documented in the help as deprecated aliases. | They duplicate the canonical band scalars on every run (ado 310-313), and permute() is not in the v1.3 tag, so no do-file written for a tagged release can rely on them. Deleting them reverses your 2026-09-07 keep-with-de | Code reading e(perm_lb_*)/e(perm_ub_*) gets a missing value and must switch to e(perm_bandlo_*)/e(perm_bandhi_*), and the replay footer at ado 13235 must be re- |

Considered and kept by the synthesiser: bootdisplay(bca) (Heavy (about 300 lines, a jackknife loop), but it is the only BCa route when the user has ); pretrendpower() (mdt_wald, the GLS pre-trend slope and the implied bias have no counterpart in Roth's pretr); histogramcohorteffects (For a binary treatment gteffects supersedes it, but under contreat() (where gteffects is r); latticeignoreweights (A user-built unweighted lattice would not re-bin on each bootstrap replication, which is e); cfxplotnose (It keeps the main-table SEs while skipping the expensive plot-cell margins SEs, which fast); ccc_absorb() and ccnoabsorb (Memory and estimation-stability routes bound up with the poisson paths the author ruled on); update (Five lines with no interactions, and it is the package's distribution channel; net install); e(permdist) and e(boot_percentile_mainresults) (Not pure aliases: e(permdist) holds per-estimand draws on the e(perm_stat) scale (the t-ra); customvce() (Real burden, but it is the only route to HC2/HC3 and other regress vce() forms, and its co); permtci (The author ruled 'Keep it' on 2026-09-23 (changes_review_2026-09-22.md:277).).

## 7. Where it lives

round1_results.json (49 findings, 16 verdicts, three deletion lenses and the synthesis),
round1_confirmed.txt (the 16 round-1 confirmations in full), round2_results.json (33 verdicts, 36 edit
lists with the verifiers' results and the held items), edits_applied.json (the 117 edits as applied,
with each verifier's verdict), oldsyntax_removal_notes.md. Commits: aeec885 (oldsyntax), 16f5392
(documentation corrections). Not run: the certification suite (documentation and an inert removal).

## 8. Rulings given 2026-09-25 (multiple choice, in session)

| item | ruling | state |
|---|---|---|
| coarsecohortcontrols() and ccc_absorb() report the effect at control value zero (F1a, F1b) | **Build the centred model**: each treated cell's slope on the control minus the cell's weighted mean, zero on reference rows; poisson path and the r(3200) retry included | to build |
| partial dropped dose slopes attenuate the ame (CT-1) | **Drop the cells from the estimation sample** through the estimability checks (treated cohort-period cells with no within-cell variation in the dose), except untreatedZero's all-zero cells; text approved: "Caution: estimability checks found N treated cohort-period cell(s) with no within-cell variation in `contreat', causing them to be dropped." | to build |
| e(vce) says oim or ols over a robust variance (POIS-1, SEB-05) | **Post robust** (e(vce) robust, e(vcetype) Robust) on the ppmlhdfe path and under pweights without cluster() or robust; oim stays for poisson on the unconditionalse path | to build |
| one-row-cell containment screen has no ppmlhdfe arm (POIS-2) | **Add the arm**; the vzero remedy list stops offering robust to a run that already has it | to build |
| permtstat caution trigger (PERM-1) | **Keep cluster() == i** | no change |
| bootstrap resolution caution at reps() below 40 under bootdisplay(bc) and (bca) | **Skip it under bc and bca** | to build |
| the seven ES-3 band edits beside the 2026-09-23 band ruling | **Keep the edits** | done (16f5392) |
| example 0 now clusters by unit with i the onset cohort (SEB-06) | **Keep** | done (16f5392) |
| deletion candidates (a) | **Delete**: cleanmatrices, bootalpha(), suppresssgprimaryeffects, ties, permrawp, emptycellsoverride, the e(perm_lb_*) and e(perm_ub_*) aliases. **Keep**: permside(), permanalyticse, altmethod_adversarial | to build |
| new cautions | **Draft and build** CT-2 (lattice borrowed slope, by default), CT-3 (contreatwithin normalisation), ES-3 (analytic band at cluster() >= i), ES-6 (adversarial rank deficiency, df posted), SEB-05 (iid containment: the existing caution's trigger extends to plain OLS, no new text); **not** CT-9 | texts approved, to build |
| CT-2 text | approved: "Caution: `n' treated cohort-period-by-lattice-bin piece(s) have no within-piece variation in `contreat', so their dose slopes are borrowed from another bin of the same cell and the ame includes those borrowed slopes. Coarsen the lattice or drop the affected cells." | |
| CT-3 text | approved: "Caution: under contreatwithin, `n' period(s) have no untreated reference rows (`list'), so one treated cell's level effect in each was set to zero and the att there depends on that choice. Rely on the ame results." | |
| ES-3 text | approved: "Caution: simultaneousbands is built from the analytic covariance with cluster() at the level of i or above, so the band is too narrow. bootstrap() gives a valid band." | |
| ES-6 text | approved: "Caution: the event-study covariance for the adversarial test is rank-deficient (rank `k1' of `K'), so the test uses only the directions a generalized inverse retains and its degrees of freedom fall to `k1'." | |
