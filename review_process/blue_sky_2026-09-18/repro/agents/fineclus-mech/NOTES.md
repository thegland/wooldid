# fineclus-mech — running notes

TARGET: the two-regime bootstrap unit rule under finer-than-i clustering
(ado:8987-8998 detect; ado:9348-9349 apply). METHOD: code-mechanism probing.

## Mechanism as read (pinned build)
- `wooldid_bootrun` builds `smp = ((if)==1) & !missing(clusv)`  (ado:8911-8924).
  NO markout for missing y / controls / dose.  `smp` is ALSO the row set the
  physical re-id leg keeps (`qui keep if smp`, ado:9345).
- Regime detection (ado:8986-8995), only when `reid` AND `clus` are non-empty:
    gegen min(clusv) by(i) if smp  ==  gegen max(clusv) by(i) if smp
    -> wdb_nest = 1 (nesting)  : __wdbnewi = group(__wdbns, i)
    -> wdb_nest = 0 (finer/XC) : __wdbnewi = i
- `reid` fires when bootstrap(standard) AND (att_it|att_i|att_itime requested
  OR (boottstat AND cluster())).   (ado:378)
- `wdb_C` (the map-frame cluster count, -> bootdf(clusters) = G-1) is also
  counted on `smp`, i.e. on rows that need not be in e(sample).
- NOTHING records which regime ran: no display line, no e() macro.
  (compare e(perm_irebuilt) for permonsetcohorts).

## Probes
- p00_smoke.do      build assertion. PASS: `which wooldid` -> pinned build path.
- p01_regime.do     builds four variants + replicates the detection predicate.
    A  clean nesting                                    PREDNEST=1  N=768
    B  firm 1 really changes industry in 2005           PREDNEST=0  N=768
    C  firm 1 changes industry in 2005 AND y is missing
       for firm 1 from 2005 on (rows NOT in e(sample))  PREDNEST=0  N=764
    C0 y missing for firm 1 from 2005 on, no move       PREDNEST=1  N=764
    => C and C0 have IDENTICAL estimation samples but OPPOSITE regimes.
       The regime is decided on rows the estimator never sees.
- p04_regime_measure.do  measures BOOTSE for A/B/C/C0 at reps(300), same rseed,
    plus two hand-rolled oracles (draw industries w/ replacement; re-id vs keep-i).
- p02_edges.do      edge geometries (2-way cluster, cluster==i, finer, per-row,
    constant, missing cluster ids).   NOT YET RUN
- p05_finer_coverage.do  SIZE of the new finer-than-i arm (i=cohort group,
    cluster=firm, true effect 0, 120 MC draws).   NOT YET RUN

## RESULTS SO FAR (2026-09-17)
- p01: regime prediction A=1 B=0 C=0 C0=1; C and C0 have the SAME N (764).
- p04 (reps 300, rseed 4242):
    A  nest   BOOTSE(att_i)=0.14110902275043  (analytic cluster(ind) 0.04087)
    B  keep-i BOOTSE(att_i)=0.10057797358284  ratio to A = 0.7128
    C  keep-i BOOTSE(att_i)=0.10043761347134  N=764
    C0 pending
- p02 T1: two-way cluster(ind year) works analytically (SE .0298);
  bootstrap(standard) refuses it r(198) on BOTH the count path and the reid
  path -- the documented screen at ado:6427 fires. CLEAN.
- theory: under keep-i a twice-drawn cluster's cells are identical duplicates
  that the equal-cell weighting counts ONCE, so the per-cluster multiplicity
  collapses from Poisson(1) (var 1) to Bernoulli(1-1/e) (var 0.2325).  The
  composition component of Var(att_i) is therefore cut by ~4.3x, and the SE
  ratio tends to sqrt(0.2325)=0.482 as that component dominates.  Measured
  0.713 here, i.e. about half the variance was composition.

## 2026-09-17 later
- p04 C0 landed: BOOTSE(att_i)=0.14122699826312 vs C 0.10043761347134.
  IDENTICAL N (764), IDENTICAL att_i (1.01054627838245), IDENTICAL analytic
  SE (0.04030370871346).  Ratio 0.7112.  ==> FCM-2 CONFIRMED.
- p02 T2 (cluster == i): RC=0, NCLUST=96, BOOTSE(att_i)=0.09122637551236.
  The boundary works; nesting fires (min==max trivially).
- p02 T3/T3B/T4 are INVALID TESTS, not wooldid defects: in p02's panel the
  onset varies by firm INSIDE cgrp, so `wooldid y cgrp ...' is refused
  r(459) for a reason unrelated to clustering.  Redone properly in
  p11_finer.do, where onset is assigned BY GROUP.  Likewise most of
  p07_boottstat.do's cases (they all use i=cgrp on p02_base) -> superseded
  by p11.
- p02 T5 (constant cluster variable): r(459), correct and informative.
- p02 T6 (missing cluster ids): plain run N=764 NCLUST=24, consistent.

## more results
- P06 DECISIVE: att se BIT-IDENTICAL C0 vs C (0.12587695969826); att_it/
  att_i/att_itime move by 0.765/0.782/0.782.  e(sample) subset of smp
  confirmed (0 rows the other way).  NEST ON smp=0, NEST ON e(sample)=1.
- P11 A_firm: i=cgrp, cluster(firm) [finer], att+att_i, reps(100) -> RC=0,
  N=768, NCLUST=96, DFR=95, att se .0904, att_i se .0758 (analytic .0778).
  The finer-than-i arm RUNS and produces sane numbers.
- T7 (i with missing values): RC=0, N=760.  T8 (user `if' removing the
  straddling rows): RC=0, N=384.

## FINAL STATE (all probes)
DONE: p00 p01 p02(T1-T10) p04(A,B,C,C0,ORACLE_NEST,ORACLE_KEEP) p06 p08(k=0,1,5)
      p11(S0,A_firm,A_pair,A_cgrp,A_rowid,B_NOCLUS,C_BTS_*,D_SINGLETON,E_NOSUB,F)
RUNNING at report time: p08 k=20; repro_FCM2.do (validation run of the
      self-contained repro; its setup lines already reproduced).
NOT RUN: p05_finer_coverage.do (size of the finer-than-i arm, 80 MC draws x
      79 reps) -- the box carried 14-19 competing Stata processes for the whole
      session and an MC of that size could not land.  p07_boottstat.do is
      SUPERSEDED by p11 (its i=cgrp cases were invalid on p02_base).
KEY FILES: FINDINGS.md, repro_FCM2.do (self-contained), p08_dose.do
      (self-contained), p11_finer.do (self-contained), p02_edges.do.
