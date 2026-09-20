# fineclus-cellwt-refute — skeptic pass over FCW-01..06

Rule followed: every verdict rests on a do-file I wrote from the CLAIM, on data I
generated. I read the prober's do-files for orientation only; none of them was run.

## Status
- r04_caution.do   DONE  -> FCW-04 reproduces, but severity must be cut (see below)
- r06_permgate.do  PARTIAL (Y2 arm had a coding error: cohort coded 0; rerun queued)
- r01a_straddle_real.do  RUNNING (FCW-01 on a realistic DGP)
- r02_phantom.do   RUNNING (FCW-02)
- r03_compvar.do   QUEUED (FCW-03 mechanism)
- r05_emarker.do   QUEUED (FCW-05 e() diff)
- n01_keepi_oracle.do QUEUED (NEW: is the keep-i branch itself honest?)
- n02_noclus_atti.do  QUEUED (NEW: the third deleted refusal, unclustered att_i)

## Facts established by reading the pinned build (no Stata)
- ado:6459 `if "`clus'" != "" {` ... ado:6501 `}` really does enclose the `if `reps' > 0`
  block (brace count over 6455-6505: open at rel-5, reps block rel-19..24, close rel-47).
  So FCW-04's structural claim is right.
- `ifr(`wdp_if')` at ado:394 is the RAW user `if` (c_local wdp_if `"`if'"' at ado:615), so
  bootrun's `smp` = user-if & nonmissing cluster, with no markout. FCW-02's mechanism is
  right on the source.
- wdb_nest appears at 406/407, 8969-8998, 9323, 9348 only; no ereturn. FCW-05's grep is right.
- BUT sthlp:241 DOES document the two-regime rule in the help: "Where cluster() nests i the
  copies become distinct units; where it is finer than i, or cuts across it, the copies stay
  inside their cohort". So the RULE is documented; only the per-run record is missing.

## r04 result (2026-09-17)
reps(10) throughout, one self-generated 12-county panel.
  RUN1 att att_i att_it, no cluster  -> rc 0, NO "Caution: with reps(10)"
  RUN2 same + cluster(cty)           -> rc 0, caution PRINTS
  RUN3 att only, no cluster          -> rc 0, NO caution
  RUN4 att, no cluster, bootdisplay(wald)  -> rc 0, NO caution
  RUN5 att, no cluster, boottstat    -> rc 0, NO caution
  RUN6 bootstrap(bayes), no cluster  -> rc 0, NO caution
So the gate is wider than FCW-04 said (all three caution variants + bayes).
*** REFUTATION MATERIAL: the disclosure FOOTER is NOT cluster-gated (ado:11784-11785).
RUN1/RUN3 still print "The rank-based coefficient p-values have resolution 2/10 = 0.200:
a p printed as 0 means p below that", and RUN5 prints the studentized analogue. The
substantive content of the caution therefore still reaches the user on those legs.
The ONE leg where information is genuinely lost is bootdisplay(wald): ado:11785 is
`else if "`wdbdsp'" != "wald"`, so RUN4's footer says nothing at all about resolution,
and the wald caution that would have said the joint/equality/adversarial p cannot fall
below 1/(reps+1) is the one sitting inside the cluster block.

## r06 partial result
X1 i=state (finer than the onset grouping), cluster(county), permute(shuffle) -> r459 nesting refusal
X2 same + permonsetcohorts -> r459 bijection refusal        (both as FCW-06 reported)
Y1 i=cohort (IS the onset grouping), cluster(county), no option -> r459 nesting refusal
Y2 i=cohort + permonsetcohorts -> r411 in this run because I coded never-treated cohort as 0
   and wooldid requires positive integers. RERUN with cohort in {1,2,3}.
Z  i=state, cluster(county), bootstrap(standard) -> rc 0  (bootstrap accepts what permute refuses)

## RESULTS LEDGER (2026-09-17, all from my own do-files)
- r01a (FCW-01, realistic): 6 paired ratios se(att_i) B/A = .9547 .9486 .8861 .8575 .8726 .8709,
  mean .8984; se(att) ratio mean 1.0091; b(att_i) diff exactly 0.  -> 23% is a degenerate ceiling.
- r05 (FCW-05 + FCW-01 degenerate): e() dumps byte-identical between regimes (diff empty) while
  se(att_i) ratio = .78957 (theory .7629).
- r02 (FCW-02, realistic, 3 completed pairs, rest killed for budget): se(att_i) ratios .854 .875
  .843 with se(att) and b bit-identical and e(N)=419 in both arms.
- r03 (FCW-03 mechanism): cluster(state)==i, HET se .0886/.0946/.0899/.0969 vs HOM .0044/.0041/
  .0057; predicted sqrt(.16/20)=.0894.  Composition term confirmed at full size.
- r04 (FCW-04): cluster-gating confirmed and wider than claimed; BUT the ado:11784-11785 footer is
  not cluster-gated and still states the resolution on every leg except bootdisplay(wald).
- r06 (FCW-06): X1/X2 refused as reported; Y2 (i = onset grouping, cluster finer) rc 0 with
  e(perm_irebuilt)=1.  The motivating pairing IS reachable.
- KEY DOCUMENT FACT: PREREG_FEATURE_FINECLUS_2026-09-16.md amendment 3 already retracts the
  sthlp:255 analogy that FCW-03 spends its claim correcting.  The stale copy survives in
  test/suite/mc/mc_fineclus.do:36-38.

## LATE RESULTS
- n03 (NEW): C7 arm - 14 design counties, 7 measured. header "N clusters = 7", e(df_r)=13,
  footer "t(13)"; analytic run on identical data gives e(df_r)=6. C12 arm - 2 measured
  clusters, 55 of 99 reps failed, still e(df_r)=13 and no few-cluster caution.
- n02 (NEW): MC SD att_i .21252; cluster(hosp) keep-i boot .21398/.20873 (ratios 1.007/0.982);
  unclustered obs boot .29362/.29409 (ratios 1.381/1.384).
- n04: cluster(t) cross-cutting i runs rc 0, 199/199, se ~1.5x the cluster(unit) bootstrap.
- r04b: unclustered + bootdisplay(wald) + jointtests at reps(10) prints no caution, but the
  joint tests were ABSENT in that run, so I could not exhibit a printed joint p with no floor
  disclosure. FCW-04 severity -> low.
- r03b: at cluster(state)==i, HET, se(att) = .0709192322 on the RE-ID leg (att att_i) and
  .0709192322 on the COUNT path (att alone) - BIT-IDENTICAL. The composition term is not
  injected by the re-id leg.
