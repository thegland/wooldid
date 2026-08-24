# Theory A — locked per-case predictions for the wooldid ame-SE Monte Carlo
### Lens: cluster-robust asymptotics and incidental parameters
Author: theory agent A. Written 2026-08-24, BEFORE any MC output existed; the only empirics used are
six single-draw mechanical probes (§5; one Stata run, no Monte Carlo). Code citations are to
`wooldidinstall/wooldid.ado` (bugfixes branch). Background relied on: `review_process/inv1_v_es_ame_2026-08-24.md`,
`wooldid_old_bug_report.md` §5.2–5.4.

---

## §1 Framework

**The estimator, by FWL.** With `contreat(dose)`, wooldid fits
`reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t <contreatcontrols>) vce(...)` (L1099).
Under `esfixedbaseperiod` + `esrelativeto(-1)` on this grid (treated at t=6, t=1..12, esprelength
internally bumped 4→5 at L3998-99), every treated cohort g gets **11 saturated cells**: relative
periods −5..−2 and 0..+6 each carry a visible dummy D_gh and a visible slope S_gh = D_gh·dose. The
base period (rel −1, t=5) is recoded to (_i=0,_t=0) (L641/646): the 6 base-period obs of each treated
cohort sit OUTSIDE the cells, pooled with the controls. This is the linchpin of the whole case grid.
Because (D_gh,S_gh) saturate the {1,dose} span within every cell, OLS zeroes the within-cell moments
Σû=0, Σd·û=0, and the fitted line in every cell is the cell's own OLS line with raw slope s_gh. The
margins-visible cell slope is then exactly `γ̂_gh = s_gh − φ̂(g,h)`, where φ̂(g,h) is whatever baseline
dose slope the contreatcontrols absorb: **basic** → one pooled φ̂ (from ~1080 control obs + 90 treated
base obs); **byCohort** → cohort-own φ̂_g, which (base period!) equals **s_g,base, the raw 6-obs OLS
slope of cohort g's own base-period cell** (probe-verified to 7 digits); **byTime** → period-specific
φ̂_t from the ~90 control obs of the cell's own calendar period. The reported AME at horizon h is
b_h = Σ_g w_gh γ̂_gh (equal cells → (1/15)Σ_g), so its TRUE sampling variance always contains the
within-cell slope noise Σ w² Var(s_gh), with Var(s_gh)=σ²/Σ_cell(d−d̄)² > 0 whenever dose varies
within the cell — full-rank truth, per-horizon scale set by 6-obs cell slopes.

**What plug-in CRVE can and cannot see — the rank ladder.** The CRVE meat uses cluster scores
X̃_c'û_c. Under cluster(i)=cohort, every treated cell lies inside one cluster, and the normal
equations have already purged û of exactly the within-cell moments the cell coefficients load on: the
own-cell part of every slope score is identically zero. What survives in V(b_es_ame) is only the
variance of the SUBTRACTED baseline slope, and its survival is governed by where that slope lives
relative to cluster boundaries. **byCohort**: φ̂_g is itself cohort-g-local; α_g and φ_g's normal
equations zero the base cell's 1- and d-moments too, so treated cluster scores vanish in EVERY
direction, and since γ̂_gh = s_gh − s_g,base is an exact function of cohort-g data alone, control
clusters contribute nothing either → **V of the entire ame block = 0 exactly; rank 0**. **basic**: one
pooled ringer moment crosses clusters → V_es_ame = Var_CRVE(φ̂)·11′, **rank 1**, every ame SE =
SE(φ̂). **byTime**: φ̂_t is a per-period cross-cluster (control-group) moment, one per reported
horizon → V_es_ame ≈ diag(Var_CRVE(φ̂_t(h))) + small off-diagonals, **rank 8, full**, horizon-varying
SEs — no degeneracy signature, but the cell-slope noise component is still annihilated, so variance is
understated by the factor Var(φ̂_t)/(Var(φ̂_t)+cell noise) ≈ 0.3–0.5. One mechanism, three
signatures. No cohort-level plug-in CRVE variant can ever recover the annihilated component (the
score is identically zero; CR2/CR3 divide by I−H_gg, singular in exactly those directions). It IS
estimable by other estimator classes: cluster bootstrap of the AME functional over cohorts, or the
analytic hybrid Var(b_h) = Σw²·σ̂²_gh/Σ_cell(d−d̄)² + V̂ar_CRVE(φ̂-part) (pieces from disjoint obs,
asymptotically uncorrelated). That existence is what keeps most of these cases in bucket (iii) rather
than (iv).

**Fine VCEs and small-cell leverage.** Under cluster(unitid) each cell spans 6 unit-clusters, under
robust 6 observations; scores no longer vanish and the plug-in V is consistent with the right
structure. The finite-sample bias is downward: 2 within-cell parameters per 6 obs shrink residuals by
(1−h_ii), average 2/3, and the slope-relevant, d̃²-weighted leverage is worse for skewed
(lognormal-ish) dose because the largest |d̃| obs carries the largest h_ii — effective shrinkage
0.5–0.7 — partly offset by the (N−1)/(N−K)·G/(G−1) ≈ 1.2 multiplier (K ≈ 375 of N = 2160). Net:
per-horizon variance ~0.65–0.85 of truth, t-tests mildly oversized (~6–10% at 5%), while 8-dimensional
quadratic forms (adversarial, joint) compound the per-direction bias and stay badly oversized
(20–35%; NF-3's measured 31.5% is the anchor). On buckets: I use "correct SE" = consistent for the
true sampling variance of the posted estimator. Bucket (ii) is deliberately rejected for every
degenerate case below: the posted object is NOT a correct SE of the posted estimator (the estimator's
true variance is positive and finite), so "correct but degenerate" never applies on this grid.

---

## §2 Per-case predictions

DGP recap: 30 cohorts, 15 treated at t=6, 6 units/cohort, t=1..12, y never contains dose; truth of
every AME estimand is 0, so MC dispersion of b̂_h is the truth benchmark. "SE(φ̂)" magnitudes below
are analytic orders for THIS grid (d1: ≈.06; d2: ≈.09–.12). "Truth SE" = MC sd of the per-horizon
AME. Buckets are the user's verbatim four; refinements stated where the honest answer straddles.

| # | case | true variance structure of reported AMEs | predicted reported-SE behavior | bucket | prescribed remedy |
|---|------|------------------------------------------|-------------------------------|--------|-------------------|
| 1 | d1 · basic · cluster(i) | Full rank; per-horizon Var ≈ (1/225)Σ_g Var(s_gh) + Var(φ̂) ≈ .04–.09 (SE ≈ .20–.30); cross-horizon corr small (~.05–.2, from shared φ̂); main-ame SE ≈ .10–.13, pre-ame ≈ .12–.15. Estimator well-defined, non-degenerate. | Rank-1 V_es_ame = Var(φ̂)·11′; ALL ame SEs (8 horizons + main + pre) identical ≈ .06 (probe draw: .0590403 on all 10 rows); contrast variances exactly 0; variance understated ~8–20x; per-horizon t rejects ~35–60% at 5%; adversarial refuses (q=0). | (iii) — refined: NOT estimable by any cohort-level plug-in CRVE (score ≡ 0; CR2/CR3 singular); estimable by cluster bootstrap/analytic hybrid, and under this iid DGP also by robust/finer clustering. If "correct SE" were restricted to the CRVE class, (iv). | Code-fix: detect (every treated cell inside one cluster of cluster()) → refuse to post e(V_es_ame)/blank ame SE column with explicit message; steer to robust/unit-level cluster when iid plausible; constructive: cohort bootstrap or hybrid estimator. |
| 2 | d1 · basic · cluster(unitid) | Same truth as #1. | Full rank, horizon-varying SEs tracking truth from below: variance ~0.7–0.85 of truth (leverage, partly offset by dof multiplier); per-horizon t ~6–9% at 5%; joint/adversarial still oversized ~15–35%. | (i) — with the caveat that "normal desirable properties" holds for the SE column only; V-based 8-dim joint tests remain badly oversized (leverage bias compounds). | Optional: CR2/HC2-type or bootstrap refinement for joint tests; document. |
| 3 | d1 · basic · robust | Same truth as #1 (iid errors → robust valid). | Full rank, horizon-varying; variance ~0.65–0.8 of truth (slope-leverage on skewed dose, HC1 offset); per-horizon t ~7–10%; joint tests 20–35% (NF-3 anchor 31.5%). | (i) — same caveat as #2. | Same as #2. |
| 4 | d1 · byCohort · cluster(i) | Estimand recentred: b_h = (1/15)Σ_g(s_gh − s_g,base). Well-defined, but Var ≈ 2× basic (≈ .09–.15; SE ≈ .30–.40) and cross-horizon corr ≈ +0.5 (shared 6-obs base-cell slope); att point estimates unchanged vs basic, ame points differ. | **Total collapse, presenting as missing**: cluster-g scores ≡ 0 in every direction and γ̂ is cohort-local → ame V-block exactly 0 (rank 0); implementation converts this to b-only posting: ALL 20 displayed SEs missing "." (att rows collaterally poisoned via V entries the att Jacobian touches), caution "standard errors were not identified for 20 of 20 reported estimates, so no variance-covariance matrix e(V) was posted"; e(V_es_ame)/e(b_es_ame) NOT posted (f2ok gate L3352-59 fails: no positive-finite diagonal); expect e(*_es_att) withheld too. Reproduces every rep (structural, not draw-specific). | (iii) — refined: the missing-SE output is an ACCIDENTAL refusal; for the ame block, within the CRVE class the SE is genuinely not estimable (score ≡ 0), so a referee insisting on that class would say (iv)-and-the-refusal-is-correct; but bootstrap/hybrid estimate it fine, and the ATT-side SEs (positive true variance, CRVE-estimable) are destroyed as collateral — that part is squarely (iii). | Say-not-estimable with the RIGHT message (name the mechanism: cohort-level clustering + cohort-own baseline slope ⇒ zero cluster scores), preserve the att-side V instead of losing everything, and steer: byCohort + cluster at the cohort level is self-annihilating; use finer clustering or bootstrap. Also user-steer: byCohort silently recenters the estimand on a 6-obs base-period slope, doubling true variance. |
| 5 | d1 · byCohort · cluster(unitid) | Same truth as #4 (Var ≈ 2× basic, corr ≈ +.5). | Full rank; SEs ≈ √2× the basic-unitid SEs; reported off-diag corr ≈ +.4–.6 (probe robust draw: mean .62); variance ~0.7–0.85 of truth; joint tests oversized. | (i) — caveat as #2; plus the estimand-recentring user-steer. | Document byCohort's estimand change; joint-test caveat. |
| 6 | d1 · byCohort · robust | Same truth as #4. | As #5 (probe draw SEs ≈ .28–.31, corr ≈ .62). | (i) — caveat as #3. | As #5. |
| 7 | d1 · byTime · cluster(i) | Full rank; per-horizon Var ≈ cell noise + Var(φ̂_t) ≈ .07–.12 (SE ≈ .26–.35); corr across horizons small; each horizon carries its own calendar-period control slope. | **No degeneracy — the dangerous case**: full-rank 8×8 V (probe eig .0011–.051), horizon-VARYING SEs ≈ SE(φ̂_t) (probe .086–.199), positive contrast variances; but cell-slope noise annihilated → variance ~0.3–0.5 of truth; per-horizon t rejects ~15–30% at 5%; adversarial/joint POST p-values here and over-reject. | (iii) — the cleanest (iii) on the grid: correct SE exists (hybrid/bootstrap/finer VCE), code posts a plausible-looking but ~2–3x-understated V with no degeneracy signature. | Hardest to detect from output — remedy must be the a-priori predicate (treated cells saturated within clusters), not a rank check. Same fix menu as #1. |
| 8 | d1 · byTime · cluster(unitid) | Same truth as #7. | Full rank, tracks truth from below (~0.7–0.85); per-horizon ~6–9%. | (i) — caveat as #2. | As #2. |
| 9 | d1 · byTime · robust | Same truth as #7. | As #8, slightly lower. | (i) — caveat as #3. | As #3. |
| 10 | d2 · basic · cluster(i) | As #1 but cells are less dose-dispersed (within-cell dose var ≈ .087 vs d1's ≈ .32 avg): truth SE ≈ .32–.48; corr small. | Rank-1, all 10 ame SEs = SE(φ̂) ≈ .09–.12; variance understated ~10–20x; otherwise exactly as #1. | (iii) — as #1. | As #1. |
| 11 | d2 · basic · cluster(unitid) | As #10. | As #2 (d2's less-skewed within-cell dose → leverage bias equal or slightly milder). | (i) — caveat. | As #2. |
| 12 | d2 · basic · robust | As #10. | As #3. | (i) — caveat. | As #3. |
| 13 | d2 · byCohort · cluster(i) | Var ≈ 2× #10 (SE ≈ .45–.65), corr ≈ +.5. | Identical failure signature to #4: b-only, all SEs missing, no e(V), no F2 matrices. Within-unit dose variation changes nothing — the score identity is algebraic. | (iii) — as #4. | As #4. |
| 14 | d2 · byCohort · cluster(unitid) | As #13 truth. | As #5. | (i) — caveat. | As #5. |
| 15 | d2 · byCohort · robust | As #13 truth. | As #6. | (i) — caveat. | As #6. |
| 16 | d2 · byTime · cluster(i) | Truth ≈ cell noise (.153) + Var(φ̂_t) (≈ .11): SE ≈ .45–.55. | As #7: full rank, horizon-varying, variance ~0.35–0.5 of truth. | (iii) — as #7. | As #7. |
| 17 | d2 · byTime · cluster(unitid) | As #16. | As #8. | (i) — caveat. | As #8. |
| 18 | d2 · byTime · robust | As #16. | As #9. | (i) — caveat. | As #9. |
| 19–21 | d1 · within · {cluster(i), cluster(unitid), robust} | IDENTICAL to #1/#2/#3 in every respect: with a control group present and full support, contreatwithin touches only screens (L484/557/673; safety-off skips L1244), never the regression. Probe: mreldif(b)=mreldif(V)=0 vs basic. | Bit-identical to #1/#2/#3 (rank-1 SE(φ̂) under cluster(i), etc.). | Same buckets as #1/#2/#3. | Same as #1/#2/#3; plus doc note: contreatwithin is inert on data that has a control group — users invoking it for "within-only variation" are not getting a different estimator here. |
| 22–24 | d2 · within · same 3 | Identical to #10/#11/#12. | Bit-identical. | Same. | Same. |
| 25–27 | d3 · basic · {cluster(i), cluster(unitid), robust} | **Estimand unidentified**: dose constant within every cell ⇒ S_gh ∝ D_gh; no within-cell dose variation exists anywhere in the design's identifying span. No meaningful AME estimator exists. | VCE-INVARIANT across all three: reghdfe omits all 165 treated-cell slopes (dummies expand first, later collinear columns dropped — probe: caution "a treatment period coefficient interacted with something ... was dropped for multicollinearity ... triggered 165 times"); ame point estimates display as literal 0 with SE missing "." (10 of 20 rows), "standard errors were not identified for 10 of 20 ... e(V) contains placeholder entries"; e(b_es_ame)/e(V_es_ame) NOT posted (f2ok=0); att side intact and clean. rc = 0 — the run "succeeds". | (iv) — genuinely not estimable in principle (the point estimand itself is unidentified), and wooldid posts a results block (b=0 rows, placeholder V) anyway. Refinement: the displayed 0/(.) is an SE-shaped output of a non-estimator; the generic caution points near, but never says, the true cause. | Refusal is the only correct output: detect Var_cell(dose)=0 across treated cells pre-regression and error with "contreat() has no within-cohort-period variation; continuous-treatment AMEs are not identified — dose is collinear with the cell dummies", instead of 165 generic cautions + zero rows. |
| 28 | d1c · basic · cluster(i) | ≈ #1's truth DESPITE the cohort-period shock: the shock is absorbed by the treated-cell dummies (Var(s_gh) unchanged), and for the unit-constant, within-cohort-demeaned dose the per-(g,t) control score contribution Σ_i d̃_i·c_gt ≈ 0 (exact for control cohorts on this balanced panel, up to the small perturbation from treated base-period obs entering the ringer partialling), so even Var(φ̂) barely inflates. MC ame dispersion ≈ d1's. | Rank-1, all ame SEs = SE(φ̂) ≈ .06 (≈ unchanged vs d1) — the same collapse in the one case where cohort-level clustering is conceptually mandatory (for the att side it genuinely is). | (iii) — with the (iv)-straddle at its sharpest: here "use robust instead" is NOT a defensible generic steer (the user's clustering motive is real); within the CRVE class the ame SE is unestimable; a cohort bootstrap or the hybrid estimator is the only principled repair. I defend (iii) because those exist. | Say-not-estimable + constructive fix (bootstrap/hybrid) is the only full remedy; the steer-to-robust shortcut happens to be nearly valid for the AME functional on THIS design (the shock is annihilated for the ame by the cell dummies) but that is a design accident, not advice. |

Auxiliary locked predictions, all cluster(i) cases: displayed df = 29; under basic (rows 1, 10,
19/22-cluster(i), 28) the main-ame and pretreatment-ame SEs equal the common ES-ame SE exactly;
adversarial/joint ame p-values refuse on rows 1/10/19/22/28 (rank 1 ⇒ q=0) and on rows 4/13 (nothing
posted), but POST and over-reject on rows 7/16 and on all unitid/robust rows.

---

## §3 Sharpest discriminating predictions (what the MC should score)

1. **The rank ladder (basic=1, byCohort=0, byTime=8) under cluster(i)** — same data, same clustering,
   three contreatcontrols: basic posts one repeated positive SE (=SE(φ̂), ≈.06 d1 / ≈.10 d2) with
   corr-1 V; byCohort posts NO SEs at all (b-only, "20 of 20 not identified", no e(V_es_ame)); byTime
   posts full-rank, horizon-varying, plausible-looking SEs understated ~2–3x in variance. Only the
   incidental-parameter/score-locus account predicts this three-way split; "few-clusters bias",
   "reghdfe dof bug", or "margins delta bug" theories predict none of it.
2. **within ≡ basic bit-for-bit** (b and V, mreldif 0) on every d1/d2 draw — contreatwithin is a
   screen-relaxation flag, not a spec change, whenever a control group exists.
3. **byCohort recenters the estimand**: its ame point estimates differ from basic's on the same data
   (b_h = mean_g(s_gh − s_g,base); FWL identity probe-verified to 7 digits), its TRUE per-horizon
   variance is ≈ 2× basic's, and its true (and unitid/robust-reported) cross-horizon correlation
   jumps to ≈ +0.5 (shared base-cell slope) vs ≈ 0–.2 under basic. Its ES_att points are unchanged.
4. **d3 is VCE-invariant**: identical 0/(.) ame output and identical caution under cluster(i),
   cluster(unitid), and robust — the degeneracy lives in the design matrix, not the VCE. If d3's
   output varies with the VCE, my account of the drop mechanism is wrong.
5. **d1c ≈ d1** in both truth (MC ame dispersion) and reported behavior (rank-1, ≈ same SE(φ̂)):
   the cohort-period shock that justifies clustering is annihilated for the ame functional by the
   very saturation that breaks the CRVE. A materially larger d1c reported-or-true ame variance would
   falsify the absorption argument.
6. **Ordering of reported vs truth on the non-degenerate cases**: robust ≲ cluster(unitid) < truth
   per-horizon variance (ratios ~0.65–0.85), per-horizon t-rejections 6–10%, while 8-dim
   joint/adversarial ame tests stay at 20–35% at nominal 5% on the same runs whose att-side
   adversarial is correctly sized (NF-3 pattern).
7. **d2 truth SEs ≈ 1.5–1.8× d1's** at every horizon (smaller within-cell dose dispersion), while
   d1/d2 reported basic-cluster(i) SEs stay in the .06 / .09–.12 bands — the understatement factor
   GROWS from d1 to d2.

---

## §4 The Stata-convention ruling

`regress y i.g, vce(cluster g)` posts near-zero SEs on cluster-saturated dummies and StataCorp ships
it. **wooldid's obligation is different, and the convention is no defense here.** Four grounds.
(1) *Altitude*: in regress the pathological numbers sit on raw incidental coefficients the user can
see, next to the factor notation that produced them; caveat emptor is at least coherent. wooldid
AGGREGATES those coefficients into a headline estimate it names "ame", prints in a labeled SE/t/p/CI
column, and returns as e(V_es_ame) explicitly for downstream inference consumers (honestdid, lincom,
test, pretrendpower). A package whose entire value proposition is "we do the aggregation and
inference correctly so you don't have to" owns the statistical properties of the aggregate it
certifies; "the regression did it" is not available at that altitude. (2) *Silence and systematicity*:
the failure mode is invisible in the output (basic prints ten identical plausible SEs; byTime prints
varying ones) and arises in the package's own core recommended configuration — contreat requires
within-cohort dose variation, hence multi-unit cohorts, and cohort-level clustering is the
conventional default for cohort-level adoption. The convention Stata ships never silently rescales a
USER-FACING summary by 8–20x. (3) *Own precedent*: wooldid already accepts exactly this obligation —
it refuses the adversarial p on a degenerate V ("composite-null covariance is degenerate"), errors
(322) on dropped treatment dummies, and ships a caution when V_es rows are degenerate (L3365). Posting
rank-1 e(V_es_ame) while refusing the test that would have consumed it is the inconsistent posture.
(4) *Consumer contract*: e(V_es_ame) exists FOR honestdid/pretrendpower; feeding a rank-1 (or
silently ~2–3x understated, byTime) covariance to a partial-identification procedure produces
confidently wrong robust CIs with no diagnostic anywhere. Ruling: wooldid must implement the cheap
detection predicate (every treated cohort-by-period cell inside a single cluster of cluster()) and
then refuse-with-diagnosis: blank the ame SE column, withhold or degeneracy-flag e(V_es_ame), name
the mechanism, and steer (finer clustering / robust where defensible / bootstrap where clustering is
mandatory). What it must NOT do is either of the things the two builds have done: post the rank-1 V
silently (old and current basic path), or post nothing but a generic "SEs not identified" caution
whose true cause the user cannot recover (current byCohort path).

---

## §5 Probes run and what they verified

One Stata instance, one do-file (`/Users/hegland/.claude/jobs/640b9294/tmp/amese/thA/probeA.do`, log
`probeA.log`), six single-draw runs on the grid's own d1 DGP (seed 2010) and one d3 variant. No Monte
Carlo, no *amese*csv / mc_report* files touched.

- **RUN D, basic+cluster(i)**: all 8 ES-ame SEs = .0590402862, e(V_es_ame) an exactly constant
  matrix (max−min spread 1.16e-15), maxdiag 3.486e-3 — rank-1 confirmed on this grid's DGP; main and
  pretreatment ame SEs also .0590403 (RUN C table).
- **RUN A, byCohort+cluster(i)**: b-only posting; ALL 20 SEs missing (att rows included); caution
  "standard errors were not identified for 20 of 20 reported estimates, so no variance-covariance
  matrix e(V) was posted"; e(b_es_ame)/e(V_es_ame) not posted. **FWL identity verified**: manual
  (1/15)Σ_g(s_g,t6 − s_g,t5) from 30 per-cell regressions = 0.15328824 vs displayed post0 ame
  0.1532882. The att-rows-missing observation is what pins the implementation layer: exact-zero slope
  variances are converted (reghdfe/margins) into missing V entries that poison every functional whose
  Jacobian touches a slope column — by contrast d3 (slopes OMITTED rather than degenerate) leaves the
  att side intact, confirming the poisoning runs through present-but-degenerate columns.
- **RUN B, byTime+cluster(i)**: full-rank 8×8 (eigenvalues .00109–.0506), horizon-varying SEs
  .086–.199, max |offdiag corr| .565 — no collapse, as predicted by the cross-cluster locus of φ̂_t.
- **RUN C, within+cluster(i)**: mreldif vs basic = 0.0e+00 for both e(b_es_ame) and e(V_es_ame) —
  within ≡ basic exactly.
- **RUN E, byCohort+robust**: SEs .280–.310 (≈ √2× basic-robust scale), mean off-diag corr .620 —
  the shared-base-cell covariance structure appears in a correctly-sighted VCE as predicted.
- **RUN F, d3+basic+cluster(i)**: multicollinearity caution triggered exactly 165 times (= 15 treated
  cohorts × 11 cells — every treated-cell slope), ame rows displayed as estimate 0 / SE "." (10 of 20
  flagged), att rows clean, e(b_es_ame) not posted, rc 0. Pre-registered alternative (dummies dropped
  → error 322) is REJECTED by the probe; slopes-omitted is locked as the d3 prediction, and its
  VCE-invariance (P4) follows since the drop precedes the VCE.

Corrections the probes forced on the pre-probe theory, disclosed for scoring honesty: (a) I predicted
byCohort+cluster(i) would post literal-zero SEs; the machine-zero V is instead surfaced as MISSING
SEs/b-only posting (the structural rank-0 claim stands; the presentation layer differs). (b) I flagged
a possible per-cohort slope drop under byCohort; the base period sitting outside the cells (L641)
removes the collinearity, and no drop occurs — probe confirms no drop caution and exact FWL identity.
