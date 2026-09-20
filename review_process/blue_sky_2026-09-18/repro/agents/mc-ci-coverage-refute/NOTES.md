# SKEPTIC notes — refuting mc-ci-coverage findings

Build asserted: `which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado in every probe log.
All data generated here (`mkdat` inside each do-file). No repo fixtures, no prober fixtures.

## p1_quick.do / p1_quick.log   [DONE]
Own DGP: 60 units, i=cg (3 onsets 2004/2006/2008 + never), 12 periods, unit random
intercept sd .70, AR(1) rho .4 errors sd .55, WITHIN-cohort effect heterogeneity
(unit multiplier 1+N(0,.25)), cohort effects 1.00/1.15/0.85.

BLOCK A (mcci-03): one dataset, four runs. e(N)=720 and e(df_r)=59 in all four.
  none                         att=0.929543611  se=0.124672021
  esprelength(1) espostlength(1) att=0.850282631 se=0.125284747
  esprelength(2) espostlength(2) att=0.858508187 se=0.132884576
  esprelength(3) espostlength(3) att=0.894165051 se=0.151295323  (se +21%)
  NO caution anywhere in the log (grep -i "aution" -> nothing).
  e(escomp_refN) did NOT print even with ES horizons requested -> check separately.

BLOCK B (mcci-01): one balanced dataset, 60 units.
  B1 impvalues cluster(unit)  b=1.088394008 se=0.022510329 df_r=44
  B2 impvalues cluster(cg)    b=1.088394008 se=0.105582803 df_r=2   <-- textbook IM
  B3 analytic  cluster(unit)  b=1.111133012 se=0.118553096 df_r=59
  B4 analytic  robust         b=1.111133012 se=0.135143616
  B5 impvalues cluster(unit), UNBALANCED (late entry/holes, 693 rows)
                              b=0.985960369 se=0.064685882 df_r=44
  Same b in B1 and B2 (cluster-equal-weighted mean), but se ratio 4.69x and
  df 44 vs 2. NO caution.

## p2_mcimp.do  [RUNNING]  -> p2_mcimp.txt
140 draws, 3 arms (A impvalues cluster(unit)+ES post horizons, B impvalues
cluster(cg), C analytic cluster(unit)). Decisive statistic: mean(se) vs sd(b).

## p3_mcband.do [QUEUED]  -> p3_mcband.txt
250 draws, 12 units, i=cg, cluster(unit), esprelength(3) espostlength(2),
simultaneousbands level(90). Relative-period-only true effects so every
weighting gives the same horizon truth. Computes measured simultaneous
coverage, the ORACLE crit, and pointwise coverage.

## p2_mcimp RESULTS (140 draws)  [an_p2.py]
arm            n   mean(b)    sd(b)   mean(se)  se/sd   cov_att  cov_ceq  cov_MCmean  df
A  imp cl(unit) 140  0.9977   0.1376    0.0233  0.169    0.257    0.264      0.243    44
B  imp cl(cg)   140  0.9977   0.1376    0.1092  0.794    0.971    0.964      0.964     2
C  analytic     140  1.0145   0.1373    0.1252  0.912    0.957    0.950      0.936    59
ES_post0        140  0.9956   0.1278    0.0286  0.224    0.314    0.329      0.336    44
ES_post1        140  1.0023   0.1359    0.0285  0.209    0.343    0.336      0.286    44
=> mcci-01 CONFIRMED on my own DGP. The control arm B (cluster == i) is the
   discriminator the prober never ran: the IM machinery itself is FINE there.
   ES horizons under impvalues also under-cover (0.29-0.34) - prober asserted,
   I measured.

### EXACT MECHANISM (arithmetic, p1 numbers)
se(cluster(cg)) = 0.105582803 on df 2  => sd over the 3 distinct reads = 0.18287478
predicted se if those 3 values are each repeated 15 times over 45 unit-reads:
  0.18287478*sqrt(2*15/44)/sqrt(45) = 0.022510329234
observed se(cluster(unit))            = 0.022510329     ratio 0.99999999
=> the 45 per-cluster reads take EXACTLY 3 distinct values. CORRECTION to the
   prober's mechanism: the sampling component is NOT "absent entirely" - those
   3 reads carry it (arm B se/sd = 0.794). What is wrong is that each distinct
   read is counted with multiplicity, so the ttest divides by sqrt(45) and
   refers to t(44) where the honest quantities are sqrt(3) and t(2).

## p4_mech RESULTS
(A) espostlength(2) alone   : att 0.929543611 == the no-ES att, refN = 270
    esprelength(2) alone    : att 0.858508187,                  refN = 180
    esprelength+espostlength: att 0.858508187                    refN = 180
    => mcci-03's mover is esprelength() ALONE; espostlength() does not move it.
       e(escomp_refN) IS returned once ES horizons are asked for.
(B) i=unit cluster(unit) impvalues se=0.069774 df 44 vs analytic se=0.123567
    => only 0.56x on one draw: MC in p5.
(C) i=cg cluster(reg), 5 regions crossing cg: IM se/t/p/CI ALL BLANKED with
    "A collapsed IM standard error usually reflects deterministic or noiseless
     data"; this panel's idiosyncratic sd is 0.55 and unit sd 0.70.

## p3_mcband RESULTS (250 draws, G=12, level 90)  [an_p3.py an_p3b.py an_p3c.py]
SIMULTANEOUS coverage 0.788 at nominal 0.90 (3 binom SE 0.057)  => ~6 SE low
pointwise coverage per horizon 0.872-0.928, se/sd 0.96-1.16     => pointwise OK
mean reported crit 2.1536, df_r 11
DECOMPOSITION of the 0.112 deficit (400k-draw simulations on the MC's own
empirical horizon correlation):
  reported crit                       2.1536      coverage 0.788
  normal sup on the TRUE correlation  2.2112 (+2.7%)        0.804
  t(11) sup, common scale             2.4950 (+15.8%)       0.856
  t(11) sup, per-horizon scale        2.6054 (+21.0%)       0.888
  MC oracle q90                       2.7574 (+28.0%)       0.904
=> mcci-02 CONFIRMED and its OPEN QUESTION ANSWERED: the z-not-t reference is
   ~85-90% of the deficit; plugging in the estimated correlation as if known
   contributes ~2.7% of width (0.016 of coverage).

## p5_imunit RESULTS (100 draws, i=unit cluster(unit))
arm  n   mean(b)  sd(b)  mean(se) se/sd cov_att cov_ceq cov_MCmean df
I  100  1.0059  0.1249  0.0827  0.662  0.820  0.830  0.840  44
C  100  1.0233  0.1273  0.1153  0.906  0.890  0.890  0.910  59
repro_S4_imunit.do (80 draws) reproduces: imp se/sd 0.702 cov 0.850/0.863;
analytic se/sd 0.967 cov 0.900.
=> NEW FINDING skep-01. Distinct mechanism from mcci-01 (no multiplicity here;
   the 45 reads are distinct but positively correlated through the shared
   control group and period effects).

## REPRO VALIDATION (all four run clean end to end)
repro_S1_impvalues.do : ratio observed/predicted = 1.000000000000000;
   MC 60 draws -> impvalues coverage 0.283 vs analytic 0.983.
repro_S2_band.do      : 200 draws -> simultaneous coverage 0.780 at nominal .90,
   mean crit 2.1537, oracle 2.8292, ratio 0.761, pointwise 0.875-0.950.
repro_S3_esprelength.do: table reproduced, zero Caution lines in the log.
repro_S4_imunit.do    : 80 draws -> imp se/sd 0.702 cov 0.850.

## STATE: all verdicts settled, FINDINGS.md complete.
