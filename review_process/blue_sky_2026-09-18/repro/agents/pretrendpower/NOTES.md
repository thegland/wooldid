# pretrendpower() — simulation against the claim. NOTES

WORKDIR: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/pretrendpower

## What the code actually does (read from ado:10308-10680, wooldid_ptp)
S   = e(V_es_TYPE)[pre,pre], symmetrized;  bp = e(b_es_TYPE)[pre]
m_j = relyr_j - esrelativeto        (built from the "preK" colname: -K - ref)
c2  = m' S^-1 m ;  slope = m'S^-1 bp / c2 ;  se = 1/sqrt(c2)
W   = bp' S^-1 bp ; p = Ftail(keff, df_r, W/keff)   [chi2 if no df_r]
mdt_wald(P) = sqrt( npnF(keff, df_r, invF(keff,df_r,.95), 1-P) / c2 )
bias_*(P)   = mdt * (hbar - ref)
mdt_nis(P)  = bisection on 1 - mvnormalcvqp(-c*sd, +c*sd, g*m, S) = P,
              where c = invttail(df_r,.025)  <-- t crit inside a NORMAL probability
nissize     = 1 - mvnormalcvqp(-c*sd, c*sd, 0, S)      (k<=4 only)

## VERIFIED EXACTLY (p02_algebra.do / p02_algebra.log)
Hand-recomputed slope/se/W/p/mdt_wald/bias_wald from the posted e(b_es_att),
e(V_es_att), e(df_r) reproduce e(roth_*) to 1-2 ulp.  So the displayed numbers
ARE that algebra; the open question is whether that algebra delivers the power.
-> this licenses the cheap MC harness: recompute W/p/mdt from a plain wooldid
   run (no pretrendpower) instead of paying the 7-12s NIS quadrature per rep.
   Plain run 3.16s vs 7.14s with pretrendpower(80) [p02_algebra.log timer 2/3].

## External oracle (Roth's own `pretrends` package), same b/V:
  pretrends power 0.80  -> slope 0.217610      wooldid mdt_nis(80) = 0.218494
  ratio 1.00406 == invttail(279,.025)/invnormal(.975) = 1.9684/1.95996 = 1.00431
  => the ONLY difference is that wooldid uses a t critical value where Roth uses
     normal.  At df=279 that is 0.4%.  MUST re-measure at small df (few clusters),
     where t vs normal is ~10-25% -- and where nissize is then computed as a NORMAL
     probability of a t-sized event.

## RUNS
p01_orient.do/.log  -- design, labels, hbar arithmetic verified by hand
                       (hbar=1.8333 = mean post relyr over FULL post support; bias
                        = mdt*(hbar-ref) confirmed)
p02_algebra.do/.log -- exact algebra + pretrends oracle + timing
p04_simA.do         -- RUNNING (launched 01:45). DESIGN A, iid errors, default OLS
                       vce, i=unit, k=3 pre coefficients, df_r=279.
                       arm0 A0_null.csv  delta=0   R=300
                       arm1 A1_pow80.csv delta=median(mdt80) R=300  -> expect .80
                       arm2 A2_pow50.csv delta=median(mdt50) R=300  -> expect .50
                       ~47 min total.

## PLANNED
- Design B: i=unit + cluster(unit)  (documented-oversized geometry)
- Design D: robust
- slope-CI coverage (free from the same CSVs)
- NIS size check: arm0 nisrej rate vs reported nissize (0.1226 on this design)
- esrelativeto(-2): m straddles the reference, sign handling
- anticipation(1): hbar measured on the SHIFTED onset (code says shifted; verify)
- help EXAMPLE 7 (sthlp:1546) uses cluster(i) + pretrendpower -> does the run print
  its own "not a calibrated test" caution?  (doc/UX candidate)
- bootstrap() path; non-estimable pre-period; permute() refusal

## 2026-09-17 ~02:10  vech-ordering hypothesis: CLEAN (negative result)
wooldid passes `vech(S)'` (Mata vech = COLUMN-major lower triangle) to
`mvnormalcvqp(L,U,M,V,1000)`.  Stata's help does not state the ordering, so a
row-major consumer would silently mis-assemble S for k>=3.  Decided from p01's
posted S (hex-exact) with scipy as the oracle:
   correct S, t-crit bounds -> nissize 0.12259799
   MISREAD (row-major)      -> nissize 0.09631734
   wooldid reports             0.12257446
=> ordering is CORRECT; the 2.4e-5 gap is mvnormalcvqp's 1000-point quadrature.
Script: (inline python, reproduced in p11_vech.py)

## Sharp prediction for DESIGN E (small G)
`nissize` and `mdt_nis` evaluate a NORMAL probability at a *t* critical value
(ptpniscrit = invttail(df_r,.025)).  At df_r=279 that costs 0.002 in size.
At df_r=15 (G=16 clusters) invttail = 2.1314 vs 1.96: wooldid will report
nissize ~= P(any |Z_j| > 2.1314) ~ 0.08, while the test it describes actually
rejects when |b/se| > 2.1314 with b/se ~ t(15), i.e. 5% per coefficient,
~0.12 jointly.  sthlp:1384 calls e(roth_nissize_*) "the actual size of the
any-individually-significant pre-test".  Design E measures the realised rate.

## The caution predicate (ado:1413-1444), read exactly
Fires only when every `i` sits inside ONE cluster value
(count of distinct (i,cluster) pairs == count of distinct i), and only
when neither bootstrap() nor permute() is in force.  So:
  Design A (no cluster)              -> silent
  Design B (i=unit, cluster(unit))   -> CAUTION fires
  Design C (i=coh, cluster(unit))    -> silent
  Design E (i=coh, cluster(state))   -> silent   <-- the geometry I test
The caution's own remedy text asserts "bootstrap-based clustered SEs give a
joint test that is valid under this design".  That assertion is testable and
UNTESTED (cost: each MC rep needs a full bootstrap).  Flagged as untested.

## 02:50 DESIGN E null arm (385/400) -- THE HEADLINE IS HERE
16-state policy panel, i = onset cohort (4 codes), fe(st), cluster(st), G=16,
df_r = 15, k = 3, AR(1) errors within state.  wooldid prints NO caution here
(cluster is strictly finer than i - the very remedy the at-i caution names).
  joint pre-trend F test size   : 39/385 = 0.1013  CP95 [0.0730,0.1359]  (nominal .05)
  Roth slope 95% CI coverage    : 350/385 = 0.9091 CP95 [0.8758,0.9359]  (nominal .95)
  mean reported slope se 0.1437 vs MC sd of the slope 0.1822 -> se/sd = 0.789
  NIS rejection 0.1091 vs reported nissize 0.0781
  per-coefficient |b/se|>tcrit  : .039 / .039 / .047   <-- MARGINAL inference is FINE
  reported mdt_wald(80) median 0.5463, sd 0.1228 (CV 22%)
KEY STRUCTURAL POINT: the individual pre-period t tests are correctly sized; it
is the QUADRATIC FORM b'S^-1 b and the GLS slope - the two things pretrendpower()
is entirely built on - that break, because S is a 15-df cluster-robust estimate
and both statistics invert it.

## 02:55 STATE
DONE: p01 p02 p04(A) p06 p06b p08(E) p12(poisson) p13(E extra arm) p11_vech.py
RUNNING: p09_simB.do  (Design B = i=unit + cluster(unit), the CAUTIONED at-i
         geometry; arm0 250 reps then arm1 220 reps; ~28 min from 02:54)
NOT RUN: p07_simC.do (Design C, i=cohort + cluster(unit), G=80 -- the same
         "sound" geometry as E but with many clusters; would say whether the
         PTP-1 effect is purely a few-cluster effect).  repro_power.do and
         repro_nissize.do are written but NOT yet executed (they are the
         self-contained repro artifacts; their content is a strict subset of
         p08/p06b which HAVE been run).
FINDINGS so far: PTP-1 (power shortfall on an uncautioned 16-cluster design),
PTP-2 (roth_nissize evaluates a normal probability at a t critical value),
PTP-3 (poisson takes the chi2 shortcut at G=16), PTP-4 (help example 7).

## 03:35 FINAL STATE
All designs measured.  FINDINGS.md is complete (PTP-1..PTP-5 + clean list +
untested list).  Repro artifacts EXECUTED: repro_nissize.do/.log,
repro_power.do/.log (+ repro_power_singledraw.log), repro_helpex7.do/.log,
repro_designB.do/.log (running at 03:35, 150+150 reps).
Headline: on a 16-cluster state-level DiD that wooldid does not caution about,
the joint pre-trend test has size 0.1025 and the reported mdt_wald(80) has
achieved power 0.7533 (nominal cutoff) / 0.5800 (size-corrected); the honest
80%-power trend is ~0.676 against a reported 0.546 (19% optimistic).
The same geometry at G=80 (Design C) is exact, and Design A (OLS) is exact
(0.8000 and 0.5000 on the nose), so the machinery itself is right and the
defect is the few-cluster behaviour of the quadratic form b'S^-1 b.
