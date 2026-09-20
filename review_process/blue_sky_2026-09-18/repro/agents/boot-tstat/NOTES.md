# boot-tstat prober NOTES (method: recompute-from-the-dump + shape invariants)

Target: boottstat, bootdf(), bca, simultaneousbands.
Build: /private/tmp/.../scratchpad/bluesky/build (adopath ++, never `do`).

## Code map (pinned build line numbers)
- `wooldid_bootrun`  ado:8887   — smp, map frame (wdb_C = G), three draw legs
  - count path (aweight multiplicities)          ~9410
  - phys leg (unclustered + weight + boottstat)  ~9245
  - reid leg (cell-weighted estimand OR boottstat+cluster) ~9300
  - `WOOLDID_BOOTNOSUB`  ado:9367-9380 (suppresses cluster(__wdbns) substitution)
  - P7 assertion (rep clustered on resampled id) ado:9386-9420
  - BCa jackknife loop  ado:9489-9520 (count path, zero weight, G fits)
- `wooldid_bootwire` ado:9543 — rewiring; `wdbdf` from bootdf(); boottstat block ~9840
- mata `wooldid_boottp`  ado:14332
- mata `wooldid_bootbca` ado:15760
- mata `wooldid_bootband` ado:16776 (simultaneous band)
- `wooldid_bootsavemk`  ado:14871 (dump layout; se_ names truncated to 32 chars)

## Key semantics confirmed by reading
- `wdb_reid` set when bootstrap(standard) AND (cell-weighted estimand OR (boottstat & cluster)).
  So boottstat+cluster ALWAYS takes the re-identification leg. (ado:378)
- bootdf(clusters) -> df = `wdb_C` - 1 where wdb_C = # clusters on `smp`
  (`smp` = user if & !missing(cluster); NO markout on y/controls).
- boottstat: table col2 = analytic se_obs, col3 = b/se_obs, col4 = add-one rank p,
  cols 5/6 = b - qhi*se_obs, b - qlo*se_obs.

## RESULTS
### p02_recompute.do  -- CLEAN (exact)
G=20 states, i=state, cluster(state), reps=199, boottstat, esprelength(2) espostlength(2).
Recomputed n, t_obs, p, lb, ub for all 6 coefficients from d2.dta.
`RECOMP_MAXRD n=0 t=0 p=0 lb=0 ub=0 b0=0 se0=0` -- BIT-EXACT agreement with
e(boot_tstud) and e(wdidmainresults). The add-one convention 2*(min+1)/(n+1) is
confirmed on the studentized coefficient rows (not just cells).
Shape invariants held: rowtotal(w_*) == 20 on every replication; #se_ vars == #coef
vars; dump `_rep==0` == e(b) and e(boot_ase_obs) exactly.
Observation to chase: mean(se_r)/se_obs ranges 0.78-0.93 (systematically < 1).

### p03_df_nosub.do
- **A (bootdf/dead clusters): DEFECT.** 20 states, y suppressed in 3.
  `e(N_clust)=17` but `e(boot_df)=e(df_r)=19`. Same data with the 3 states
  physically dropped gives `e(boot_df)=16`. Bootstrap se 0.177585 (A) vs
  0.170171 (B) = +4.4%; 1 of 99 replications failed in A, 0 in B.
  The printed header says "N clusters = 17" and the footer says "t(19)".
- **B (WOOLDID_BOOTNOSUB): behaves as documented.** ratio_rms(se_r/se_obs)
  0.927 normally, 1.204 with the global set (x1.30; theoretical sqrt(2-1/G)=1.40).
  It is a global, so it persists across runs and across wooldid calls in a
  session; P7 is skipped while it is set. Undocumented but WOOLDID_-prefixed.
- **C: all four refusals fire with informative messages** (bc+boottstat,
  wald+boottstat, bca without cluster, bca+boottstat). CLEAN.

### p05_bca.do -- BCa CLEAN
- e(boot_accel) matches an INDEPENDENT delete-one-cluster jackknife I ran
  myself (20 separate wooldid fits, aweight 0/1): reldiff = 0 on all 6 columns.
- z0, both BCa endpoints and the dual p recomputed from e(bootdist) +
  e(boot_comparator) + e(boot_accel): exact.
- e(boot_jk_completed) == G on both G=20 and G=8.
- Estimate inside the BCa interval on every row at G=20 and G=8; no blanking
  needed, accelerations ~0.01-0.05 at G=8.
  (My first recompute printed p=1 where wooldid printed 0 on three columns --
   that was MY missing r0<=0 guard, not a wooldid defect.)

### p04_bands.do -- DEFECT (S3)
S3 = boottstat + simultaneousbands, skewed (log-normal) outcome, G=20, 6 ES horizons.
  horizon pre3: pointwise percentile-t [-9.64223368, 0.65471983]
                simultaneous band      [-6.96720703, 4.39249989]
  The 95% SIMULTANEOUS band is NARROWER than the 95% POINTWISE interval on the
  lower side (by 2.68). Cause: the pointwise interval is EQUAL-TAILED from the
  signed t-quantiles while the band is SYMMETRIC b +/- q_.95(max_h|t_rh|)*se.
  For symmetric t_r, q_.95(|t|) = q_.975(t) so the band dominates; skew breaks it.
  S1/S2 (non-studentized, same data) had no violation on that seed.

### p07_p7dead.do
- T1: boottstat + 3 dead clusters: NO spurious P7 error (rc 0, 98/99 complete).
  Same df defect: N_clust 17, boot_df 19.

### p08_legs.do -- studentization legs
- L1 unclustered boottstat (fweight leg):  mean(se_r)/se_obs 0.945-0.966  CLEAN
- L2 unclustered + aweight (phys leg):     0.953-0.971                    CLEAN
- L3 boottstat + bootstrap(bayes):         REFUSED r(198), informative
- L4 clustered boottstat (reid leg):       rms 0.915-0.972 (help claims ~0.92+-0.04) CLEAN
- L5 bootdf(0) / bootdf(-5) / bootdf(2.5): all refused r(198)
     bootdf(1e9) ACCEPTED -> e(boot_df)=e(df_r)=1e9 (absurd input, low priority)

### p09_bandskew.do -- band, 2 horizons, G=12, one high-effect cluster
- K1 (percentile display):  no violation, but razor thin: (q.975-b)/sd = 2.1702 vs c = 2.1738
- K2 (boottstat):           VIOLATION at post0. point [-4.76649500, 9.35624776]
                            band  [-4.43840483, 9.39162977]  -> band lb INSIDE point lb.
                            (2nd independent confirmation of BT-01, different design,
                             GAUSSIAN errors -- so skewed errors are not required.)
- K3 (bootdisplay(wald)):   VIOLATION on BOTH sides at BOTH horizons.
                            band crit c = 2.1738 < invt(11, .975) = 2.2010.
                            Here the pointwise interval IS symmetric, so the help's
                            documented non-nesting caveat (scoped to rank/asymmetric
                            displays) does not cover it.  The band's critical value
                            never sees bootdf().

### p06_bandrare.do -- R1/R2 no violation; espostlength(5) gives post0..post4, so
    my "rare horizon" never materialised (all reps complete). Retry needs espostlength(8).

### p07_p7dead.do T2 -- bootsave jointW_* observed row is WRONG
    e(joint_horiz) = 0.505634 ; dump jointW_att[_rep==0] = 0.0176935 ; 88/88
    replication W_r are >= that value, so recomputing the documented add-one share
    from the file gives p = 1.0 instead of 0.5056.
    Cause: wooldid_bootsavemk (ado:14915) indexes `jobs' -- one column per
    (family, PRE-HORIZON) observed COEFFICIENT -- by `jwacc's family column index.

### repro validation (2026-09-17, after the advisor check)
- repro_band_narrower.do            FAIL = 1 (1 horizon)          -- BT-01 confirmed
- repro_band_ignores_bootdf.do      e(boot_df)=7, zc=2.364624, bandcrit=2.199558,
                                    HORIZONS_WITH_BAND_INSIDE_POINTWISE=2, FAIL=1 -- BT-03 confirmed
- repro_bootsave_jointW.do          FAIL = 1; dump _rep==0 jointW_att = -0.33935852
                                    = e(b_es_att)["pre3"]        -- BT-04 confirmed
- repro_bootdf_deadcluster.do       FAIL_DF = 1 (A df 19 vs e(N_clust)-1 = 16, B df 16)
                                    SE ratio A/B at reps(499) x 3 seeds:
                                    0.97312 0.99366 1.04714 (mean 1.00464) -> the SE gap
                                    I saw at reps(99) was NOISE and is RETRACTED.
                                    A failed reps 3/2/2 of 499 across seeds.
