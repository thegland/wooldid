# boot-reproduce -- adversarial probe notes
Agent: BR (bootstrap reproducibility / two-leg differential testing)
Method: DIFFERENTIAL TESTING between two code paths computing the same object.

## Code reading (pinned build, sha c1d8979)
- ado:378  `wdb_reid` = reid iff bootstrap(standard) AND (att_it|att_i|att_itime requested
  OR (boottstat & cluster)). So adding att_i switches the ENGINE for att too.
- ado:392  `wdb_phys` = phys iff bootstrap(standard) & boottstat & NO cluster & user weight
  & not poisson.
- ado:9291 count leg: bsample on a COPY OF THE MAP FRAME (1 row per cluster), cluster()+idcluster().
- ado:9340 reid leg: frame copy of the FULL data, keep if smp, bsample cluster()+idcluster().
  Comment claims "same-seed draws equal the count path's".
- sthlp:532 DOCUMENTS: count leg is "the same statistic as physical resampling and matches
  Stata's own bootstrap prefix replication for replication under a common seed."
- ado:6455-6478 the whole reps() caution block IS nested inside `if "`clus'" != "" {`  -> unreachable
  without cluster().  CONFIRM BY RUNNING.
- ado:6282 bootsave() HAS a comma check (RF-6, 2026-09-13).
  ado:5868 permsave() has NO comma check -> permsave(name, replace) should write "name, replace.dta".

## Probe plan
P1 leg agreement: att-only (count) vs att+att_i (reid), same rseed, clustered and unclustered.
   Compare b (must match exactly), se, and the RAW cluster multiplicity draws via bootsave w_* cols.
P2 external oracle: Stata's own `bootstrap` prefix vs count leg, same seed (documented claim).
P3 reproducibility: same cmd twice in one session; two fresh sessions; e(boot_rngstate) replay.
P4 aw=1 no-op differential: boottstat no-weight (fweight count leg) vs [aw=one] (phys leg).
P5 reps() caution reachability without cluster().
P6 bootsave/permsave: contents vs help, replace, comma handling.
P7 ties: what does it change, deterministically.
P8 reps() vs permreps() resolution.

## Status
- [ ] nothing run yet

## RESULTS LEDGER (updated live)
- probe1+cmp1: CLUSTERED count vs reid legs -> w_* cluster multiplicity draws BIT-IDENTICAL,
  att draws ~1.5e-15, att SE 1 ulp. SEEDED LEAD CLEAN under cluster().
  UNCLUSTERED: att draws differ by 1.18e-01; se 0.07137 (count) vs 0.06154 (reid) at reps(30).
- probe2: reps() caution unreachable without cluster() CONFIRMED (help promises it unconditionally).
  permsave(name, replace) writes "name, replace.dta" CONFIRMED (bootsave refuses).
  string cluster() refused r(109) everywhere. reps/permreps cross-engine all r(198).
- probe4: Stata's own bootstrap prefix matches the count leg draw for draw (2.6e-8 = float in
  Stata's saving() file).
- probe5+cmp3: aw=1 no-op leg switch (fweight count vs physical) -> draws 8e-16.
- probe3/probe3b: same-session repetition, fresh-session repetition, and the documented
  `set rngstate e(boot_rngstate)` replay are ALL bit-identical. dots/bootdisplay/bootdf no-ops.
- probe6: `ties` == Stata's _bs_sum.ado convention EXACTLY (invnorm(#{<=}/R), minus #{=}/2 under
  ties). Deterministic. Under boottstat it moves e(boot_z0) only.
- probe7+cmp3: row-order invariance on both legs (w draws bit-identical); reps(20) is exactly
  the first 20 of reps(40); keep-i regime (cluster finer than i) count vs reid clean.
- probe10: no sortseed dependence on any leg; rseed repetition identical, no-rseed differs.
- probe12: b unchanged by bootstrap(); e(V_analytic) == plain e(V) bit-identical;
  nofastcollapse under bootstrap agrees to 1e-15.
- probe9: att_i draws identical whether asked alone or with att; att draws identical when all
  four estimands are asked; `if`-restricted count vs reid clean.
  CAVEAT: cmplib's `confirm variable` accepts ABBREVIATIONS, so the "(i) ... main_att 6.01e-02"
  line in probe9.log is spurious (it matched q1main_att_i). Ignore that one line only.

## STILL RUNNING
- probe8  : unclustered count vs reid at reps(300), two seeds per leg -> is the gap MC?
- probe11c: gteffects / subgroup surfaces count vs reid (fixed comparator).
- probe13 : missing-y panel; `smp` vs e(sample) universe; prefix match; pre-existing save file.

## FINAL STATE (all probes complete)
Findings: BR-1 (permsave comma, seeded), BR-2 (reps caution unreachable unclustered, seeded),
BR-3 (unclustered legs draw different resamples at one seed, seeded question / self-measured
mechanism), BR-4 (bootstrap resamples `smp' not e(sample): phantom clusters, e(boot_df) off by
one, documented prefix match fails; SELF-FOUND).
Self-contained repro files: repro_br1.do, repro_br2.do, repro_br3.do; probe13.do + probe14.do for BR-4.
Everything else came back CLEAN -- see FINDINGS.md's negative-result section.

## BR-4 ESCALATED to HIGH (repro_br4.do, 2026-09-17)
Two arms with the IDENTICAL estimation sample (e(N)=72, e(N_clust)=9):
  phantoms in data  -> e(boot_df)=13 (G=14), NEITHER few-cluster caution printed
  phantoms dropped  -> e(boot_df)=8  (G=9),  BOTH cautions printed
t(13) vs t(8): 6.8% narrower CI on every coefficient, and the documented <10-cluster and
<6-treated-cluster cautions are suppressed by the mere presence of unusable rows.
