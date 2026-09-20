# fw-weights probe — NOTES

## PREMISE INVERSION (recorded first, 2026-09-17)
The seeded lead says "`fw` is a NEW weight type in this build". **It is not a user feature.**
`wooldid.ado:620-632`:

```
  * fw is parseable only so the bootstrap replication path can pass an fweight; a user
  * typing it directly is refused. ...
  if "`weight'" == "fweight" & "`wdrepmode'" == "" {
    dis as error "Error: wooldid does not accept fweights. Specify aweights or pweights."
    error 198
  }
  if "`weight'" == "fweight" local weight aweight
```

`fw` appears in the `syntax` line only so the internal bootstrap replication call can pass one
while `${WOOLDID_REPMODE}` is set; the user surface refuses it. The help (`sthlp:57`) advertises
only `[aw/pw=weights]`, so there is NO doc-vs-code divergence at the surface.

Consequence for the assigned METAMORPHIC method: the "fw of k == k duplicated rows" identity is
still the decisive test, but the vehicles are
 (a) `[aw=k]` vs physical duplication  — the substitute the refusal message itself prescribes;
 (b) `[pw=k]` vs physical duplication  — point estimates only;
 (c) scale invariance `[aw=k]` vs `[aw=c*k]` — ado:1654 asserts weight scale "doesn't matter";
 (d) the two internal bootstrap LEGS (fweight leg ado:9048 vs physical-expansion leg ado:379):
     with a user weight identically 1 the data is the same, only the leg changes, so every
     studentized quantity must agree. ado:9040 claims "se_fw^2/se_phys^2 = 1.0000000 on all 400
     pairs" — this is a directly testable claim.

## Audit done by reading (no probe needed)
- All seven `global WOOLDID_REPMODE 1` sites (7672, 7786, 9065, 9250, 9380, 9450, 9500) sit
  inside `nobreak { ... cap ... global WOOLDID_REPMODE }`. No leak path found by inspection.
  `wooldid_worker` consumes and clears the global at entry (ado:538-539), so even a leak would
  be consumed by the next call.
- Second `syntax [aw fw pw /]` at ado:4572 belongs to `wooldid_validate`, not a second user
  entry point.

## Probe log

### probe1_refusal.do — CLEAN
`[fw=k]` r(198) with the documented message on every parse path tried: plain, `oldsyntax`,
non-integer k. `[aw=k]` and `[pw=k]` rc 0. `aw` with 20 zero values -> N 552 (dropped, as
documented); missing -> same; negative -> r(411) (VAL-G06, known).
REPMODE leak simulation: with `global WOOLDID_REPMODE 1` set by hand, `[fw=k]` is ACCEPTED and
`level(90)` is silently ignored (e(level)=95 vs 90 on the clean run). But all seven
`global WOOLDID_REPMODE 1` sites sit inside `nobreak{...cap...}` and the worker consumes+clears
the global at entry, so **no leak path was found**. Not filed.

### probe2_core.do — POINT ESTIMATES CLEAN
`[aw=k]` (576 rows, k in 1..4) vs `expand k` (1440 rows). Max relative diff on `e(b)`:
  ATT 9.8e-16 | ALL4 9.8e-16 | ALL4 nofastcollapse 9.8e-16 | ES 5.4e-13 (maxabs 9.3e-15)
  ES nofastcollapse 5.4e-13 | GTE 9.8e-16 | SUBGROUPS 2.3e-15 | cluster(i) 9.8e-16
  cluster(fineclus) 9.8e-16 | controls(x) 0 | fe(sg) 3.1e-15 | summarystats 9.8e-16
  customweightmultiplier 9.8e-16 | **[pw=k] 9.8e-16**
All inside the project's 1e-10 indifference threshold.

### probe3_dump.do / probe4_scale.do — escomp & gteffects CLEAN
Per-column attribution of the two matrices probe2 flagged:
  e(escomp): only `nobs` differs (54 vs 126 etc. = physical rows vs weight mass).
    **sumw_att, wshare_att, wshare_att_it, ncohorts, nunits, relyr all match exactly**, and the
    scalars `escomp_refW`=276 and `escomp_trtW`=1080 are IDENTICAL across the two frames while
    refN/trtN correctly differ. This is aweight semantics, documented, not a defect.
  e(gteffects): only `se`, `p`, `lb_95ci`, `ub_95ci`, `nobs` differ. **`estimate`, `w_att`,
    `w_att_it` match**; `w_att_i`/`w_att_itime` are missing in both (estimand not requested).
WEIGHT-SCALE INVARIANCE (ado:1654's claim) HOLDS:
  [aw=k] vs [aw=1000k]:   b 3.3e-16, V 0, main table 4.0e-16
  [aw=k] vs [aw=k/7]:     b 3.3e-16, V 0, main 3.8e-16
  [pw=k] vs [pw=1000k]:   b 3.3e-16, V 9.9e-16, main 5.0e-16
  ES+jointtests x1000:    b 3.3e-13, V 1.1e-15, main 3.3e-13
  cluster(i) x1000:       b 3.3e-16, V 4.7e-15, main 2.4e-15

### probe5_legs.do
  [aw=one] vs NO WEIGHT: b/V/main all 0.000e+00 (bitwise) on att/att_i/att_it/att_itime and
    on cluster(i); 1.9e-13 with an event study. The iw-vs-aw margins plumbing is consistent.
  LEG EQUIVALENCE, unclustered boottstat, att: no weight (internal FWEIGHT count leg) vs
    [aw=one] (physical-expansion leg): b 0, main 7.2e-16, V 1.0e-15. **ado:9040's claim holds.**
  Same for att_i: main 1.4e-1. NOT the same legs: no-weight takes the REID leg (plain
    `bsample`, ado:9345) and [aw=one] takes the PHYS leg (`bsample, cluster() idcluster()`,
    ado:9248). ado:9296 itself records that "plain bsample diverges" from the cluster+
    idcluster form, so the two legs draw different samples at a common seed -> Monte Carlo,
    not arithmetic. Chased at 3 seeds x 150 reps in probe6.
  permute(shuffle) aw=k vs duplication: b 1.3e-15, main 7.5e-15 -> CLEAN.
  permute + permtstat: main 4.672e-2 - the analytic SE in the t statistic, see probe6.

### probe6_boot.do  -- the bootstrap engines under the identity
  bootstrap(standard) cluster(i) att  (count leg):  b 9.8e-16, FULL TABLE 5.3e-15  CLEAN
  bootstrap(standard) cluster(i) att_i (reid leg):  b 9.7e-16, FULL TABLE 5.2e-15  CLEAN
  bootstrap(bayes)    cluster(i) att:               b 9.8e-16, FULL TABLE 2.9e-15  CLEAN
  boottstat           cluster(i) att (phys leg):    b 9.8e-16, table 4.67e-2 -- attributable:
      se_analytic .0527255 vs .0503719 (the aweight N/df correction), t 25.76 vs 26.96,
      **p identical (.0327869)** and the studentized CI 1.21229/1.51022 vs 1.21257/1.51061
      (2e-4). Percentile-t is invariant as theory says; the only mover is the analytic SE,
      which aweights legitimately change.  NOT A DEFECT.

### NEW LEAD (found by reading, probe8/probe10 pending)
1. WEIGHT EXPRESSION. `local w `exp'' (ado:607) keeps the raw expression TEXT and pastes it
   into arithmetic at ado:1619/1620/1630/1637/1644/1650/1792/3832-3. A top-level `+`/`-`
   binds wrong: `[aw = wf+wm]` makes `w_att_it = wf + wm/W_it`. Official Stata commands
   materialise the expression first, so there is no warning. Probe: probe8_wexpr.do,
   repro: repro_wexpr.do.
2. FLOAT CELL-WEIGHT DENOMINATORS. ado:1619/1620 are `gegen `w_per_it' = total(...)` with NO
   `double`, while the escomp twin at ado:3833 IS `gegen double`. Commit 3cb0c41 ("B21:
   gen double for ringer, interactive-control, and estimand-weight copies") converted the
   four `gen `w_att*'` NUMERATORS to `gen double` and left these two DENOMINATORS alone -
   its own message scopes itself to "the four w_att* estimand weights". If gegen defaults
   to float the B21 fix is incomplete on the denominator side, and only population-scale
   weights expose it (small integer totals are exact in float). Probe: probe10_float.do.

## FINAL STATE 2026-09-17
Three findings written up in FINDINGS.md: FW-01 (weight expression x customvce/unconditionalse
-> silently wrong cell-weighted estimands, HIGH), FW-02 (pw robust covariance labelled
e(vce)="ols", MEDIUM), FW-03 (the fweight refusal's remedy is not equivalent and undocumented,
LOW). Repros: repro_wexpr.do, repro_pw_vce_label.do, repro_collapsed_admin.do - all three run
green and self-contained.
Everything else in FINDINGS.md's clean section. Probes 1-17 + 3 repros all kept with .log files.
UNTESTED: the lattice ESTIMATION arm (r(322) on both my fixtures - the weighted quantile CUT
POINTS were verified bitwise identical instead); permonsetcohorts + weights (r(459) on both
frames, fixture); bootstrap(bayes) with a user weight AND bootalpha(size); subgroups x weights
under bootstrap; weights x honestdid/pretrends handoff.
ALSO UNTESTED: bootstrap(standard) UNCLUSTERED with [aw=k] vs physical duplication. Skipped
because observation resampling draws N obs and N differs between the frames (576 vs 1440), so
a common rseed() is not comparable; the DISTRIBUTIONAL comparison that would still be
meaningful there was not run. This is the one seeded bootstrap cell left open.
