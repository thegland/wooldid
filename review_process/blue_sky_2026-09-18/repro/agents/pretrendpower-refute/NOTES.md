# SKEPTIC notes — refuting PTP-1..PTP-5 (pretrendpower)

Build asserted: `which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado, v1.4 2026-08-19. OK.

## My own design (independent of the prober's)
16 states, 4 onset cohorts of 4 states, onsets t=6,8,10 + one never-treated cohort, T=12,
state random intercepts, quadratic+cosine common year path, AR(1) rho=0.5 within state.
Call: `wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(...)`
-> e(df_r)=15, e(N_clust)=16, keff=3, pre coefficients at relyr -4,-3,-2, reference relyr = -1,
   so m = (-3,-2,-1)'. NO caution fires (grep Caution in section A: 0 hits). Matches PTP-1's setup.

## s1_orient.do / .log  (orientation, own data)
- A (gaussian, k=3): F(3,15), nissize 0.06742, mdt_wald(80)=0.627726, mdt_nis(80)=0.596476.
- B (gaussian, k=1): e(roth_nissize_att) = 0.033052124846, and
  2*normal(-invttail(15,.025)) = 0.033052124846  -> DIFF = 0.000000000000 (exact identity).
  2*ttail(15,invttail(15,.025)) = 0.050000000000.  => PTP-2 mechanism CONFIRMED exactly, on my data.
  ALSO at k=1 the two MDT columns disagree although they power the SAME test:
  mdt_wald 1.09139 vs mdt_nis 1.08237 (ratio 1.00833). mdt_nis is the optimistic one.
- C (poisson, same panel): e(df_r)=. while e(N_clust)=16; prints "chi2(3)"; mdt_wald 0.202212.
  F-branch with df=15 would be 0.230673 -> chi2 value is 0.877 of it (12.3% smaller). CONFIRMED.
  BUT `ppmlhdfe ... cluster(st)` posts e(df_r)=. and Stata's own `poisson ..., cluster(st)` posts
  e(df_r)=. -> the chi2 branch is INHERITED from the poisson family, not invented by wooldid.
  Note the NIS quadrature is INTERNALLY CONSISTENT on the no-df branch (ptpniscrit=invnormal): 0.1331.
- E (help example 7 verbatim): on MY panel `i` had only 4 values, so the pre-block was rank-deficient
  and pretrendpower was SKIPPED entirely. That is my artifact; re-running with i = the 16 states.

## Analytic result I derived (python, s2_ideal.py) — the real mechanism for PTP-1
wooldid refers W = b'S^-1 b to k*F(k,df_r) (ado:10559 uses invF(keff,dfr,.95) and npnF).
That is the KNOWN-Sigma reference. With S estimated at df_r degrees of freedom the exact
reference is Hotelling: W*(df-k+1)/(df*k) ~ F(k, df-k+1).
  k=3, df=15 : wooldid's cutoff W>9.8621 has EXACT size 0.0784 (sim 0.0794), not 0.05.
  k=3, df=79 : 0.0546.   k=3, df=279: 0.0513.   k=1: exactly 0.0500 (Hotelling collapses to F(1,df)).
And the reported MDT is a biased estimate: mdt_hat/mdt_true = sqrt(chi2(df-k+1)/df) exactly, so
  k=3, df=15 : mean 0.9132, median 0.9070, CV 0.198, 5-95% [0.627, 1.221]
  k=3, df=79 : mean 0.9841, CV 0.081.
=> the prober's "CV 23%, range .296-.925" and "median MDT 17% low" are LARGELY standard
   estimated-covariance theory, NOT an S-inversion failure specific to few clusters; and the two
   errors (liberal cutoff, deflated MDT) OFFSET at "power at reported MDT vs nominal cutoff",
   which is exactly why their own CI [.7348,.8099] includes 0.80.

## RESULTS (all runs complete)
- s3_mc_null (250 null, own DGP/seeds): wooldid joint size 0.1120; reghdfe testparm on 3 pre
  dummies, same 250 datasets, 0.0680; slope CI coverage 0.9120; se/sd 0.775;
  median mdt_hat/mdt_true 0.8919; mean c2_hat/c2_true 1.4645.  Identity check: hand-computed
  W / slope / se / mdt match e(roth_*) to 1 ulp (rel err 4.4e-16 and 0).  2nd wooldid call 0.35s.
- s4_mc_power (independent stream): null size 0.1280; median reported mdt 0.54610;
  mdt_true 0.63007 (ratio 0.8667); POWER at the median reported MDT, at wooldid's OWN nominal
  cutoff = 0.7350 (n=400, se 0.0221) -> 95% CI [0.690,0.777] EXCLUDES 0.80;
  power at mdt_true = 0.8880; size-corrected 0.5575 / 0.7880.
- s5 PART1: help example 7 verbatim with i = unit -> caution fires, diagnostics still computed
  (roth_ok=1, df_r=59, N_clust=60).  PART2 at-i null MC: size 0.1000 (n=150), coverage 0.8933,
  se/sd 0.8121 -- NOT the prober's 0.28.  Ideal-model prediction at df_r=59: 0.0562.
- s6_nis_from_null.py: realised NIS rate 0.1000 CP95 [.0658,.1441]; wooldid reports mean 0.0769;
  the SAME rectangle under mvt(15) gives 0.1094.  => the fix is the reference distribution.
- s7_hotelling.do: closed-form size table (k x df_r) + 100k-draw Mata ideal model, all in Stata.
  k=3,df=15: size 0.0770, coverage 0.9108, se/sd 0.8448, median mdt ratio 0.9083.
  k=4,df=15: size 0.1033, coverage 0.8831.   k=3,df=279: 0.0498 / 0.9488 / 0.9955.
- s8_boot.do: under bootstrap() the reference is keyed on bootdf() only -- F(3,15) at reps 39 and
  199 alike; bootdf(normal) -> chi2(3) with no correction.  roth_W/slope_se identical across
  bootdf settings.  Could NOT establish the correct reference; filed uncertain/low.

## VERDICTS: PTP-1 PARTIAL(medium) / PTP-2 CONFIRMED(medium) / PTP-3 CONFIRMED(low)
##           PTP-4 CONFIRMED(low) / PTP-5 PARTIAL(low)
## NEW: SKP-NEW-1 (known-Sigma references, medium/certain), SKP-NEW-2 (bootstrap ref, low/uncertain)

## CAVEAT recorded for the coordinator
My at-i size (0.1000, 60 units) versus the prober's (0.28, 80 units) is a large gap and I did
NOT diagnose it. I measured a DIFFERENT at-i panel; I did not vary their panel's parameters to
find what drives 0.28. Read my PTP-5 verdict as "the magnitude is design-specific and mine is
0.100", not as "0.28 is wrong".
