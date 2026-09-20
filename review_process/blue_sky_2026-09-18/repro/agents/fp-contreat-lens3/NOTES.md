# LENS 3 (independent reproduction) — fp-contreat-B2

## What I did
Built the case a third time from the claim alone. Fixture = the claim's own shape
(12 states, 5-12 counties each, 20 years, unit-level lognormal dose with 10% time
noise, staggered onsets 8/12/15, true dose effect exactly +0.4 on treated post rows,
dose-quantile lattice(all 3)). `i` = state (12; 4 never-treated).

## Probe A (A_attribute.do / .log) — DONE. Decisive.
1. PHENOMENON REPRODUCED independently: ame = 17.976 against truth 0.400 (~45x).
   Default path: se MISSING (e(V) all missing). customvce(robust): se 0.41497.
   unconditionalse robust: se 9.8035. Same point estimate on all three.
2. ATTRIBUTION (new; neither prior agent did this): hand `regress` + `margins,
   dydx(dose) subpop(post) noestimcheck` on the identical model reproduces
   wooldid's customvce(robust) se to 6 s.f. (0.41496782 vs 0.41496882).
   => wooldid is faithfully reporting the delta-method HC1 object. NOT a wooldid
   arithmetic bug.
3. THE SAME LINEAR COMBINATION, SAME POINT ESTIMATE, other variance strings:
      vce(robust)/HC1 : se  0.4150
      classical OLS   : se 24.0194
      vce(hc3)        : se 70.8291
   The classical and HC3 SEs DO react to the ill-conditioning; HC1 is blind to it.
4. MECHANISM CONFIRMED: 228 (state x year x bin) treated groups; 83 of them have
   n<=2 (min 1, median 3). n<=2 => exactly-fitted slope => zero residuals => the
   HC1 meat contributes ~0 for those slopes while their true variance is huge.
5. WARNING AUDIT (contradicts the claim's "no path gives the reader a usable
   warning"): the DEFAULT path prints THREE cautions, incl. "estimation returned
   an all-zero variance-covariance matrix ... no standard error in this run is
   identified" and "standard errors were not identified for 2 of 2 reported
   estimates". BUT the caution MISDIAGNOSES: it says the issue is "likely related
   to ... analytic clustered SEs, with clusters occurring at a high level of
   aggregation" when there is NO cluster() here (bare robust), and its named
   remedy list ends with "conditional robust standard errors" — i.e. it steers the
   user to the one SE that is measurably the worst (0.41 on a 45x miss).
6. `margins` WITHOUT `noestimcheck` returns dydx = 0 with se 0 on all three
   variance strings: Stata's own estimability check does not pass this AME.
   wooldid passes noestimcheck. (flagged, uncertain, side observation)

## Probes B and C — RUNNING
B_fastMC.do : 25 draws, realized sampling distribution of lattice ame (fast
  reghdfe arms) + no-lattice control + does the all-zero VCV also fire under iid
  and cluster(st) + whether the prober's extra gquantiles line explains 18.091 vs
  my 17.976.
C_seMC.do : 8 draws x {customvce(robust), customvce(ols), unconditionalse robust},
  reported se + 95% coverage of the truth 0.400.

## Verdict forming
PARTIAL. Point-estimate half of the claim (systematic 45x miss) is REFUTED — the
true DGP is exactly nested in the lattice model (cell##lat##c.dose contains
cell##c.dose), the AME is a fixed linear combination of a conditionally unbiased
OLS b, so E[ame|X] = 0.400 exactly; 18.09 is one draw of a heavy-tailed variance,
as the skeptic said. Inference half CONFIRMED and sharpened: the reported HC1 ame
SE is blind to the conditioning by construction (zero-residual cell-bins), which is
why it is ~60x smaller than the classical SE on the identical estimate.

## Probes B, C, D — DONE. Results (all on the claim's own fixture, 25 draws)

### B0 — the 18.091 vs my 17.976 discrepancy is NOT a wooldid issue
The prober's loop passes `esprelength(2) espostlength(2)`; I did not. The extra
`gquantiles latu` line changes nothing (B0_WITHLATU == B0_NOLATU to 8 d.p.).
No finding.

### B2/D — the point-estimate half of the claim is REFUTED
25 independent draws, lattice(all 3), truth 0.400:
  mean -0.702, SD 6.059, median 0.206, min -26.048, max 9.769
  SE of the mean = 1.21, so the truth is 0.91 sigma away: NO detectable bias.
  median |ame - 0.400| = 1.505.
The true DGP is exactly nested in the fitted lattice model (cell##lat##c.dose
contains cell##c.dose), the AME is a fixed linear combination of a conditionally
unbiased OLS b with weights summing to one over slopes that are all truly 0.4,
so E[ame|X] = 0.400 for EVERY X. 18.09 is one draw of a heavy tail, not a miss.
No-lattice control: mean 0.369, SD 0.179 — well behaved.

### D — THE ALL-ZERO VCV IS SANDWICH-SPECIFIC  (classifies known_bydesign s7's
### "undetermined" item, whose diagnosis run "has no result")
Same data, same point estimate 2.46744112 on all three:
  lattice(all 3)             -> ame se 1.4959614   FINITE
  lattice(all 3) robust      -> ame se .           MISSING
  lattice(all 3) cluster(st) -> ame se .           MISSING
Over 25 draws: iid default returned a finite ame se in 25/25; `robust` in 0/25.
It is not occasional. It is every draw.

### D — WHICH REPORTED SE IS HONEST (25 draws, nominal 95%, truth 0.400)
  default iid (reghdfe)      : mean se 6.421, median 1.681, COVERAGE 25/25
  customvce(robust) = HC1    : mean se 0.411, median 0.417, COVERAGE 10/25 (40%)
  realized SD of the ame     : 6.059
The iid SE tracks the realized dispersion and tells the user the estimate is
worthless (rep14 ame -26.05 se 69.24; rep7 ame 9.77 se 21.79). HC1 fails at nominal 95% (10/25 = 40% coverage; the mean-ratio 14.7x and the
median-ratio ~4x are both heavy-tail-sensitive, coverage is the honest number)
and is and blind to the conditioning (range 0.371-0.643 while the
estimate ranges -26 to +10).
From C (8 draws): customvce(ols) covered 8/8, unconditionalse robust 7/8,
customvce(robust) 5/8. Only the sandwich arm fails.
ATTRIBUTION (probe A): wooldid reproduces the hand `regress`+`margins` HC1 delta
se to 6 s.f. — it is faithfully reporting that object, which is itself worthless
here because 83 of 228 (state x year x bin) treated groups have n<=2, giving
exactly-fitted slopes with zero residuals and therefore ~zero HC1 meat.

### D — NEW, NEITHER PRIOR AGENT: THE POINT ESTIMATE IS ENGINE-DEPENDENT
In 3 of 25 draws the reghdfe path and the regress path (`customvce()`) report
DIFFERENT ames; max gap 28.62. Draw 7 (seed 70007), OLS fits of the identical
model (E_zoom.log):
  reghdfe : df_r 1750, rss 3297.15031382, ame  9.76939915
  regress : df_r 1752, rss 3297.52081621, ame  0.12496007, exactly 2 coefficients
            marked omitted
So the two engines retain a column set differing by TWO near-collinear columns at
their respective collinearity tolerances: regress drops them, reghdfe keeps them
and so fits 0.011% lower rss. In 3/25 draws the margins weights load on those
disputed directions and the reported AME moves by up to 28.6 (here 78x).
(Do NOT read e(rank) as kept columns: for regress it is the rank of e(V), which
is why probe A shows rank 456 under ols but 327 under robust at the same
df_r 1744 — that drop is the zero-residual sandwich mechanism showing up a second
way, not an engine difference. Use df_r.)
`customvce()` is documented as a variance choice; here it silently changes the
reported point estimate.
Side observation, UNVERIFIED: `margins` WITHOUT `noestimcheck` returns 0 with
se 0 under all three variance strings, i.e. Stata's estimability check does not
pass this AME. margins flags non-estimable routinely under subpop() with empty
factor cells, so this is consistent with, but not evidence of, non-estimability.
The engine-dependence is the evidence.

### The caution audit (probe A) — the claim's "no warning" is REFUTED, but the
### warning that does print misdiagnoses and mis-recommends
The default path prints three cautions, including "estimation returned an
all-zero variance-covariance matrix ... no standard error in this run is
identified" and "standard errors were not identified for 2 of 2 reported
estimates". So there IS a usable warning. Two problems with it:
  (i) it blames "analytic clustered SEs, with clusters occurring at a high level
      of aggregation" when no cluster() was given (bare robust);
  (ii) its remedy list ends with "conditional robust standard errors", which
      under bare `robust` is circular advice — the user already asked for robust.
      The only in-package route to an ame SE once `robust`/`cluster()` is in play
      is `customvce(robust)`, measured here at 10/25 coverage. Following that
      advice moves the user from "no SE" (safe) to "se 0.41 on an estimate off
      by 9".
  (iii) collateral: "2 of 2" means the att SE is lost too, although the att is
      well identified and untouched by the lattice pathology (att 0.70 with
      se 0.121 under customvce; se 0.135 on the no-lattice run).

## FINAL VERDICT: PARTIAL (severity high)
