# fineclus-mech-refute — SKEPTIC notes (current)

Target: refute FCM-1..FCM-4 (two-regime bootstrap unit rule, finer-than-i clustering).
Rule followed: every piece of evidence is from MY OWN data and MY OWN do-files.
I read the prober's scripts only to know what was claimed; I never ran them as evidence.

## Own DGP (prog mkdat in r00/r01/r02/r03/r04/r06 and all three repro files)
24 counties x 4 plants = 96 plants, T=2001-2008.  onset 2004/2006/2008 by plant mod 4,
mod4==0 never treated.  ctau ~ N(1,.8^2) per county, ptau = ctau + N(0,.4^2) per plant.
cluster(county) is COARSER than i=plant.  county enters ONLY via cluster(), so moving a
plant between counties CANNOT move the point estimate.  Cleanest possible isolation.

## Files (all with .log beside them)
- r00_struct.do     structural: regime predicate on 6 variants; string cluster; missing i
- r01_boot.do       regime isolation (ARM_A/ARM_B) + dose k=0,1,8,24 + dead cluster/df
- r02_mc.do         400-draw Monte Carlo: SD_SUPER / SD_FINITE of att_i
- r03_oracle.do     hand-rolled pairs-cluster bootstrap OUTSIDE wooldid, both rules
- r04_record_df.do  e() scan + verbose (FCM-3); 15-dead-of-20 (FCM-4); weight doors
- r05_finer.do      NEW: is keep-i calibrated where a cluster is GENUINELY finer than i?
- r06_deadse.do     NEW: do dead-cluster draws move the SE?  500 reps x 3 seeds x 2 arms
- repro_FCM1.do repro_FCM2.do repro_FCM3.do repro_FCM4.do   self-contained repros

## RESULTS
### FCM-1 CONFIRMED (high)
reps 250 rseed 606060, b(att_i) bit-identical at every k:
  k=0 PREDNEST=1 se(att_i)=0.12595 se(att)=0.12820
  k=1 PREDNEST=0 se(att_i)=0.09646 se(att)=0.12844     <- 23.4% drop on ONE straddler
  k=8 PREDNEST=0 se(att_i)=0.09760 se(att)=0.13545
  k=24 PREDNEST=0 se(att_i)=0.09419 se(att)=0.13773
Ratio 0.766 == sqrt((1-p)/p), p=1-1/e = 0.763.  The claim's sqrt(0.2325)=0.482 is WRONG
(it drops the renormalisation by the realized number of included clusters).
Structural argument (decisive, does not need calibration): under a cluster COARSER than i
every cell lies wholly inside one cluster, so duplicating a cluster cannot move a cell
mean; keep-i can therefore only REMOVE variance.  Always anti-conservative.
Hand-rolled oracle (r03, 250 draws, outside wooldid): NEST 0.12717 / KEEP 0.09632 vs
wooldid 0.12570 / 0.09616 -- wooldid implements BOTH rules faithfully.
Data-implied between-county component: SD of the 24 county means of ptau = 0.5950,
/sqrt(24) = 0.1215; with estimation noise 0.049 -> ~0.131, so the nesting arm (0.126) is
about right for THIS sample.  (One dataset; stated as support, not as a coverage study.)
400-draw MC: SD(att_i - superpop ATT)=0.1797, SD(att_i - own-sample mean)=0.0491.
NOT a doc-vs-code divergence: sthlp:532 says "finer than i, OR CUTS ACROSS IT", and one
straddler does make the partitions cross-cutting.  It is a design defect.
§7's "bootstrap(standard) at cluster()==i OVER-covers" does NOT transfer: there
group(drawid,i) degenerates to the drawn id and the composition object is different.

### FCM-2 CONFIRMED (high)
DOOR 1 missing outcome: A vs B, N=764 both, b and analytic se identical to every hex
digit, se(att) identical to 3 ulps (same draws).  NEST on smp A=1 B=0; NEST on rows with
y both 1.  boot se(att_i) 0.12570 -> 0.09616.
DOOR 2 NEW (mine): missing aweight opens the same door, 0.12638 -> 0.10528.
Missing i does NOT (gegen by() gives the missing-i rows their own group): D5_PRED=1.
Overread to correct: the comment at ado:8975-8985 claims the smp placement beats the RAW
FRAME (true); it never claims smp is the estimation sample.  The sthlp:532 divergence stands.

### FCM-3 CONFIRMED (medium)
0 names containing nest/reid/regime in e(macros)+e(scalars) on both arms; verbose prints
nothing; the two tables differ only in the SE column (0.1472 vs 0.1192 at reps 40).

### FCM-4 CONFIRMED, severity UP to high
one dead cluster of 24: e(N_clust)=23, e(df_r)=e(boot_df)=23 (rule gives 22) -- tiny.
15 dead of 20 (outcome collected in 5 counties): header "N clusters = 5", e(df_r)=
e(boot_df)=19.  t(19)=2.0930 vs t(4)=2.7764 -> bootdisplay(wald) interval 25% too narrow;
e(V)+e(df_r) are what esttab/coefplot/honestdid and hand test/lincom read (the footer says
so itself).  Control with the dead rows dropped: e(df_r)=4, identical point estimates.
bootdf(clusters) is the DEFAULT under cluster() (ado:423) -- r04's NESTRUN/KEEPRUN passed
no bootdf() and still show BOOTDF=23.  Boundary: the DEFAULT bootdisplay(percentile)
table is unaffected by df.
Framing correction: sthlp:532 says "G being the number of RESAMPLING clusters", and a dead
cluster IS one here, so help and code are self-consistent; the defect is that it is a
resampling unit at all, and that e(N_clust) contradicts e(df_r) in the same output.
r06 (500 reps x 3 seeds) answers the claim's "NOT MEASURED": PAD/DROP se(att_i) ratios
0.989 / 1.119 / 1.176, mean 1.094, and the padded arm's across-seed spread is 16% against
the clean arm's 3%.  Mechanism: drawing 20 slots when only 5 can contribute makes a
replication's live-cluster count Binomial(20,.25) instead of a fixed 5.  Direction and
instability established; the +9% magnitude is INDICATIVE only (G=5, noisy).

### NEW (mine, seeded=false) — dead clusters silence BOTH few-cluster cautions
wooldid_bootctl ado:6459-6462 builds its own sample concept with no outcome markout and
every cluster-count screen hangs off it (two-cluster refusal 6464, "<10 clusters" 6479,
"<6 treated clusters" 6498).  15-dead-of-20 run: NEITHER caution prints.  Control with the
same estimation sample: BOTH print.  Severity medium; co-occurs with the wrong df above.
NOT tested and NOT claimed: whether the two-cluster refusal can be evaded.

### NEGATIVE RESULTS (clean)
- r05: where the cluster is GENUINELY finer than i (i=cgrp 12 groups, cluster=firm 96,
  effects iid across firms), the keep-i regime is CALIBRATED.  MC SD_SUPER=0.08290 vs the
  mean bootstrap se(att_i) over 12 datasets (see r05_finer.log FINER_MEAN_SEBOOT).
  §7 asked for "any finer-than-i arm that FAILS coverage": this one does not.
- String cluster variable refused informatively ("string variables not allowed in option
  cluster()", r(109)) on both the plain and the bootstrap paths.
- A straddle carried on rows with a MISSING i does not open the phantom door.

## FINAL STATE (2026-09-17) — everything below is COMPLETE, no probe left running
r00 r01 r02 r03 r04 r05 r06 r07 r08 r09 all finished; repro_FCM1/2/3/4 all executed and
all ended in their *_DONE sentinel.  Logs beside every do-file.

r06 final (500 reps x 3 seeds, se(att_i), 15-dead-of-20 vs the dead rows dropped):
  s11 0.40684/0.41155  s22 0.44773/0.39996  s33 0.47193/0.40140  mean ratio 1.094
  padded spread across seeds 16% vs clean 3%; padded loses 1-2 replications per 500.
r05 final: MC SD_SUPER 0.08290 vs mean boot se 0.08448 (n=12) vs mean analytic 0.08314
  -> the keep-i regime IS calibrated where the cluster is genuinely finer than i.
r07/r08/r09: a lead I raised and then REFUTED myself.  wooldid's analytic cluster(county)
  se(att) is bit-identical for k=1,2,8,24 straddlers while a reghdfe benchmark moves 12%,
  which looked like a VCE blind to the row-level cluster map.  r09 kills it: moving the
  same plant from 2003 not 2005, moving a different plant, moving to a brand-new 25th
  county (NCLUST=25), and half the plants straddling ALL give different numbers.  The
  k-ladder coincidence is a property of that symmetric ladder in a saturated design.
  DO NOT re-chase.  (And §1 of known_bydesign already owns analytic SEs at/above i.)
