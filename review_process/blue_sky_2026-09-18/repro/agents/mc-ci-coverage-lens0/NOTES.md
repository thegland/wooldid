# LENS0 (independent reproduction) — mcci-01-impvalues-se

## Status
- explore.do / explore.log DONE. Deterministic identity CONFIRMED at machine precision.
- mc_lens0.do running (R=120, 3 arms).

## DGP (my own, distinct from prober and skeptic)
60 units, 4 cohort groups of 15; i = cg (cohort-level FE, the canonical Wooldridge ETWFE
shape); onsets 5/7/9, cg4 never treated; T=1..12; unit random intercept sd 0.80 (NOT
absorbed by i=cg); within-cohort effect heterogeneity sd 0.30; iid noise sd 0.60;
true cohort effects 1.00/1.20/0.80.

## Realism anchor
sthlp:251 — "The intended practice is to cluster BELOW the level of i - at a level whose
units are nested inside cohorts, so that a treated cohort-by-period cell contains several
clusters." cluster(unit) with i=cg IS that geometry.

## Deterministic result (explore.log, single draw, seed 71113)
- impvalues cluster(unit): b=.90516082 se=.03535599 t=25.60 p=4.6e-28 CI=[.83390551,.97641613]
  df_r = imdf_min = imdf_max = 44.   NO caution of any kind in the log.
  The IM target (cluster-equal-weighted mean) is 1.00 — the printed 95% CI MISSES it.
- impvalues cluster(cg)   : se = .03535599*4.6904158 = .16583, df = 2.
- analytic cluster(unit)  : se = 2.9913889 x the IM cluster(unit) se.
- IDENTITY: se_unit == se_cg * sqrt(k*m*(k-1)/((G-1)*G)) with k=3,m=15,G=45
  predicted ratio .21320072, observed .21320072, rel err 3.925e-16.
  b identical across the two arms to 1.349e-15.
  => The 45 per-unit reads ARE the 3 per-cohort values each repeated 15 times.
     The defect is MULTIPLICITY (skeptic's correction is right), not an absent
     sampling component. Interval length ratio = .2132*invt(44,.975)/invt(2,.975)
     = .2132*2.0154/4.3027 = .0999 — TEN times too short.

## MC RESULT (mc_lens0.do / analyze_scope.do, R=120, my own draws)
nominal 95%; 3*binomSE = 0.0597
  ARM A impvalues cluster(unit)  n=120 mean(b)=0.99402 sd(b)=0.10438 mean(se)=0.02593
        se/sd=0.2484  meanCIlen=0.1045  cov(own MC mean)=0.342  cov(IM truth 1.00)=0.375
  ARM B impvalues cluster(cg)    n=120 mean(b)=0.99402 sd(b)=0.10438 mean(se)=0.12163
        se/sd=1.1653  meanCIlen=1.0467  cov=1.000 (conservative, df=2 -- honest)
  ARM C analytic  cluster(unit)  n=120 mean(b)=1.01600 sd(b)=0.10589 mean(se)=0.10183
        se/sd=0.9617  meanCIlen=0.4075  cov(att truth 1.02222)=0.958   <- correct
  CI_LENGTH_RATIO A/B = 0.09986  (predicted .2132*invt(44,.975)/invt(2,.975) = .0999)
b is IDENTICAL between A and B, so this is purely a variance/df defect, not an estimand one.

## SCOPE, measured (analyze_scope.log, contreat_scope.log)
Every row of the table is deflated by the SAME structural constant .21320072:
  att .21320  att(pretreat) .21320  pre2 .21320  pre1 .21320  post0 .21320  post1 .21320
Continuous treatment (contreat(dose) contreatcontrols(untreatedZero)):
  ame  se_unit .01251832 vs se_cg .0587161  -> ratio .2132003  (EXACT multiplicity constant)
  att  se_unit .02623144 vs se_cg .07403158 -> ratio .354326  (attenuated: under contreat the
       att read varies with each unit's own realized doses, so there are more than 3 distinct
       reads; still deterministic in X, still not a sampling SE)
=> the claim's contreat-ame surface, which both earlier agents had by CODE INSPECTION only,
   is now MEASURED.

## UX chain that makes this worse
On the contreat + cluster(cg) run wooldid itself prints:
  "Consider switching to use of bootstrap SEs, IM p-values, or randomization inference if
   finer clustering levels are non-viable" and "(impvalues dispersion inference, a finer-than-
   cell clustering, or conditional robust standard errors)".
A user who takes BOTH pieces of that advice (impvalues AND a finer cluster) lands exactly in
the defective cell. sthlp:251 independently calls cluster-below-i "the intended practice".

## Pre-trend consequence (one draw, analyze_scope.log)
pre1: b=0.15896 se_unit=0.013988 -> t=11.4, p~0 (spurious pre-trend rejection).
Same draw, honest arm: se_cg=0.065611, t=2.42 on df 2 -> p=0.137, not significant.

## VERDICT: CONFIRMED, critical, realistic.
