# Theory B — Design-based / finite-sample predictions for the wooldid ame-SE case grid

Author: theorist B (design-based / finite-sample sampling-theory lens). LOCKED before any MC
results. Sources: wooldidinstall/wooldid.ado (read-only; L906-919, L1017-1035, L1099, L1785ff,
L2186ff, L3302-3371, L5571-5650), review_process/inv1_v_es_ame_2026-08-24.md,
wooldid_old_bug_report.md §5.2-5.4. Verification probes (single datasets, deterministic
structure checks only, NO Monte Carlo): §5.

## §1 Framework: what is actually random, and what each variance estimator is FOR

**The fitted model** (all cases; verbose-verified in INV-1): every treated (cohort g, relative
period h ≠ base) cell carries a visible dummy theta_gh AND a visible cell dose slope gamma_gh;
the baseline dose control is absorbed and invisible to margins ("ringer"):
basic/within → one pooled slope phi; byCohort → per-cohort phi_g (all 30 cohorts); byTime →
per-calendar-period phi_t (12). With esprelength(4) espostlength(4) esfixedbaseperiod on this
DGP (onset t=6, t=1..12), EVERY treated relative period except the base (-1) is saturated
(-5..-2, 0..+6); the 8 reported ame horizons are relyr {-5,-4,-3,-2,0,1,2,3} (probe-confirmed;
validate bumps esprelength by 1 under esf).

**Exact finite-sample identities** (fixed X, any error draw; OLS normal equations):

1. Within every saturated cell: sum(uhat)=0 and sum(d·uhat)=0. Hence the cell's fitted line is
   its own OLS line: gammahat_gh + (baseline slope evaluated for that cell) = s_gh, the raw
   within-cell slope of y on dose. n_cell = 6, Var(s_gh) = sigma^2 / (6·Var_cell(d)) > 0
   whenever dose varies within the cell (d1, d2; NOT d3).
2. basic/within: gammahat_gh = s_gh − phihat; phihat is the pooled FE-purged dose slope over
   the UNsaturated observations (15 control cohorts × 12 periods × 6 obs + 15 treated base
   cells) — a many-cluster, many-obs quantity with comparatively small variance.
3. byCohort: for a treated cohort g the ONLY unsaturated cell is its base cell, and phi_g's +
   alpha_g's normal equations saturate that cell too ⇒ phihat_g = s_{g,base} EXACTLY, so
   gammahat_gh = s_gh − s_{g,base}. AME_h = mean_g s_gh − mean_g s_{g,base}: a functional of
   cohort-g data ONLY, with within-cell-demeaned influence weights in both cells.
4. byTime: gammahat_gh = s_gh − phihat_{t(h)}; with a common onset, horizon ↔ calendar period
   is a bijection, and phihat_t (t≠5) is identified purely off the 15 CONTROL cohorts' cells
   at period t (~90 obs). The reported horizons never touch phi_5 (the base period).
5. within-arm ≡ basic: contreatwithin only relaxes sample screens (wooldid.ado L484/557/673)
   and forces safety off (L4015-4027); on these DGPs (controls present in every period) the
   estimation sample and specification are unchanged. Probe: mreldif(b)=mreldif(V)=0 exactly.

**True (design-based) covariance of the 8 AMEs** (iid errors; d1c adds cell shocks — see row):
Var(AME_h) = (1/15^2)·Σ_g Var(s_gh) + Var(baseline term) − 2Cov(·); the s-part and the
baseline part are essentially independent (disjoint observations, exactly so under byCohort/
byTime for reported horizons). Full rank always (the s_gh carry independent within-cell noise).
d1 specifics: dose is unit-constant, so each cohort applies the SAME within-cell weight vector
at every horizon ⇒ (fixed X) per-horizon true variances are EQUAL across horizons, s_gh ⊥
s_gh'; byCohort's cross-horizon truth correlation = 1/2 exactly (shared base-cell noise);
basic's ≈ Var(phihat)/Var(AME_h) ≈ .05–.2. d2: same structure, variances vary somewhat by
horizon (dose redrawn each period), byCohort corr ≈ 1/2 approximately.

**Estimand vs plug-in.** cluster(c) CRVE's estimand: the sampling variance of the estimator
under independence across c-clusters, arbitrary dependence within — well-defined and positive
in every d1/d2/d1c case, for every c considered. HC robust's estimand: variance under
observation independence — correct for d1/d2/d3's iid errors, WRONG for d1c. The plug-in CRVE
computes scores from OLS residuals; any coefficient whose influence vector is (a) supported
inside single clusters and (b) orthogonal to the residuals within those clusters by the normal
equations has an IDENTICALLY ZERO plug-in score: its own-noise variance component is
annihilated — the incidental-parameter/cluster-saturation pathology. No cohort-level CRVE
variant repairs it (CR2/CR3 divide by (I−H_gg), singular in exactly those directions); a
cohort-level BOOTSTRAP of the margins functional, or the analytic hybrid
Var(AME_h) = (1/G1^2)Σ_g s2_gh/(6·Var_cell(d)) + Var_CRVE(baseline term), does estimate the
estimand. That distinction — "estimand fine, plug-in collapses" (all cluster(i) rows, d1/d2/
d1c) vs "point estimand itself undefined" (d3 only) — organizes the whole grid.

**Consequences by contreatcontrols under cluster(i)** (each cell ⊂ one cluster):
- basic/within: the gamma block of e(V) collapses EXACTLY to Var_cl(phihat)·11' (the only
  surviving channel is the shared absorbed ringer, which spans control clusters). Rank 1,
  every ame SE (main, pretreatment, all 8 horizons — 10 numbers) = SE_cl(phihat), all
  cross-horizon contrasts get 0 reported variance. Probe: diag .0015519795 constant, corr
  ≡ 1, eig {.0124, ~1e-18}; adversarial ame REFUSES (q=0).
- byCohort: gammahat_gh = s_gh − s_{g,base} has NO cross-cluster channel at all ⇒ its plug-in
  CRVE variance is EXACTLY ZERO (not just deficient). Observed: reghdfe/margins SEs come back
  unidentified — wooldid displays ALL 20 SEs (att AND ame) missing, posts b-only e() (no
  e(V)), withholds e(b_es_ame)/e(V_es_ame) (f2ok gate, L3352-3358), and BOTH adversarial
  tests refuse. Total inference blackout.
- byTime: gammahat_gh = s_gh − phihat_{t(h)}; the phi_t's live on 15 control clusters ⇒
  nonzero scores, horizon-specific. Reported V_es_ame[h,h'] = Cov_CRVE(phihat_t(h),
  phihat_t(h')): FULL RANK, unequal diagonals, generic correlations — but it contains ZERO of
  the s-part, i.e., none of the sampling noise of the treatment-cell slopes being averaged.
  Probe: diag .0050–.0236, corr −.66..+.44, eig .038..0022, adversarial ame POSTS
  (p = .00011 on a true-null draw — the oversizing shows through immediately). This is the
  dangerous case: numerically healthy-looking, structurally blind, undetectable by rank or
  equal-SE checks — only the structural predicate (every treated cell inside one cluster of
  cluster()) catches it.

**cluster(unitid) / robust** (d1, d2): cells span 6 unit-clusters / 6 observations, all scores
alive; the estimand is valid (iid errors) and the plug-in is consistent. Residual defect is
generic small-cell leverage: per-cell slopes fit 2 parameters on 6 obs, so CR1/HC1 variance is
biased down by roughly the within-cell leverage (~1/3 before the (N−1)/(N−K)·G/(G−1)
multiplier partially offsets), a fixed proportional bias that does NOT vanish as G grows and
that 6–8-dimensional quadratic forms amplify (NF-3's measured 31.5% adversarial-ame / 22%
joint at nominal 5% under robust; my two independent null draws posted p=.024, .0027).

**d3 (no within-cell dose variation):** the within-cell AME estimand is 0/0 — NOT identified,
by ANY estimator, at this design: every D_gh·dose column is exactly collinear with D_gh
(and with cohort-constant dose the ringer itself is collinear with cohort FE). What a correct
program should DO: refuse the ame quantities loudly, naming the cause (dose has no
within-cell variation; wooldid's cell-saturated contreat design identifies dose effects only
off within-cell contrasts) and steer the user to a different design/estimator class
(cross-cohort dose-level comparisons under much stronger identifying assumptions). What
wooldid DOES (probe, both cohort-constant and cell-constant-time-varying dose, all VCEs):
reghdfe drops the interaction columns (## expansion order), wooldid prints the generic
wordcount>3 multicollinearity caution, then displays every ame row as estimate 0 with SE
missing ("standard errors were not identified for 10 of 20"), withholds e(V_es_ame), refuses
adversarial-ame. The missing SE is the right instinct; the fabricated "0" point estimate in a
labeled ame row, with only a generic caution, is the residual defect.

## §2 Per-case table (28 cases)

Buckets (verbatim): (i) correct SEs with normal desirable properties; (ii) correct SEs but
totally undesirable properties; (iii) a correct SE is possible but the code accidentally posts
something incorrect; (iv) a correct SE is genuinely not estimable in principle but we post
some sort of SE anyway.

| case | what's random / true finite-sample variance | estimator's estimand + plug-in behavior | bucket | remedy |
|---|---|---|---|---|
| d1 basic cluster(i) | iid noise through s_gh (equal true per-horizon var across h, near-diagonal corr) + small phihat part; full rank; point estimator well-defined | estimand (cohort-independence variance) fine and positive; plug-in CRVE annihilates every s_gh component: V_es_ame = Var_cl(phihat)·11' exactly; all 10 ame SEs = SE_cl(phihat); contrasts 0; reported/true per-horizon var ≈ .15–.35 (SE ratio ~.4–.6); adversarial ame refuses (never posts) | (iii) | code-fix: structural detect (cells ⊂ clusters) → blank/flag ame SE column + suppress or caution e(V_es_ame); provide cohort bootstrap of the margins functional or the analytic hybrid (within-cell Var(s_gh) + CRVE Var(phihat)); NOT estimable by any cohort-level plug-in CRVE variant (CR2/CR3 singular) |
| d1 basic cluster(unitid) | same truth | estimand valid (iid); plug-in full rank, consistent; ~5–15% downward SE bias from CR1 leverage in 6-obs cells; per-horizon coverage ~91–94%; adversarial-ame oversized (~15–35% at 5%) while att-side ~5–8% | (i) | (refinement) optional CR2/HC3-style fix + sthlp caveat that many-dim ame joint tests oversize; per-horizon SEs are as correct as standard practice |
| d1 basic robust | same | same as unitid up to tiny differences (dose unit-constant makes them nearly identical here); estimand correct for iid errors | (i) | as above |
| d1 byCohort cluster(i) | AME_h = mean_g(s_gh − s_g,base); true per-horizon var ≈ 2× basic's s-part; cross-horizon truth corr = 1/2 exactly (d1) | estimand fine; plug-in CRVE variance of every gamma is EXACTLY ZERO (influence supported in single clusters, orthogonal to residuals) → observed total blackout: all 20 SEs missing, e(V) not posted, e(b_es_ame)/e(V_es_ame) withheld, both adversarial refuse | (iii) | code-fix: same structural detect; same bootstrap/hybrid (hybrid: Var(s_gh)+Var(s_g,base) within-cell pieces); the blackout at least fails loudly — keep the honesty, add the named mechanism + remedy pointer |
| d1 byCohort cluster(unitid) | same truth | valid estimand; full rank; SEs ≈ sqrt2× basic's per-horizon; reported corr ≈ 1/2; leverage caveat doubled (2 six-obs cells per cohort in play) | (i) | as basic-unitid |
| d1 byCohort robust | same | same; probe: corr .38–.54, full rank, posts adversarial p | (i) | as basic-robust |
| d1 byTime cluster(i) | AME_h = mean_g s_gh − phihat_t(h); true var = s-part + control-cell phi_t part (comparable magnitudes → roughly 2× basic per horizon) | estimand fine; plug-in reports ONLY Cov_CRVE(phi_t(h), phi_t(h')): full rank, unequal diags, but ~40–65% of true per-horizon var (SE ratio ~.65–.8), contrast variances similarly blind to the s-part; adversarial ame POSTS and over-rejects (predict 25–60% at 5%) | (iii) | code-fix: MUST be the structural predicate (numerical rank/equal-SE detectors pass this case); then bootstrap/hybrid as above |
| d1 byTime cluster(unitid) | same truth | valid; full rank; mild leverage bias | (i) | as basic-unitid |
| d1 byTime robust | same | same | (i) | as basic-robust |
| d2 basic cluster(i) | as d1-basic but per-horizon true variances vary by horizon (dose redrawn per period) | identical collapse: rank-1 = Var_cl(phihat)·11', all 10 SEs equal, contrasts 0, adversarial refuses (INV-1's repro (a) IS this case) | (iii) | as d1 counterpart |
| d2 basic cluster(unitid) | same | valid, full rank, leverage caveat | (i) | as d1 counterpart |
| d2 basic robust | same | same | (i) | as d1 counterpart |
| d2 byCohort cluster(i) | as d1-byCohort, corr ≈ 1/2 approx | identical EXACT-ZERO blackout (argument uses no d1 structure) | (iii) | as d1 counterpart |
| d2 byCohort cluster(unitid) | same | valid | (i) | — |
| d2 byCohort robust | same | valid | (i) | — |
| d2 byTime cluster(i) | as d1-byTime | full-rank phi-only V; posts oversized adversarial p | (iii) | as d1 counterpart |
| d2 byTime cluster(unitid) | same | valid | (i) | — |
| d2 byTime robust | same | valid | (i) | — |
| d1 within cluster(i) | ≡ d1 basic cluster(i) — bit-identical b and V (probe: mreldif 0) | ≡ | (iii) | ≡ basic row; note contreatwithin also disables the safety screens that would catch degenerate variants |
| d1 within cluster(unitid) | ≡ d1 basic unitid | ≡ | (i) | — |
| d1 within robust | ≡ d1 basic robust | ≡ | (i) | — |
| d2 within cluster(i) | ≡ d2 basic cluster(i) | ≡ | (iii) | — |
| d2 within cluster(unitid) | ≡ d2 basic unitid | ≡ | (i) | — |
| d2 within robust | ≡ d2 basic robust | ≡ | (i) | — |
| d3 basic cluster(i) | NO within-cell dose variation: the within-cell AME point estimand is unidentified (0/0); the posted statistic is the constant 0 | no variance estimator has a well-defined target for the LABELED estimand — the estimand itself is undefined at this design (the one grid case where this, not plug-in collapse, is the verdict); observed: interactions dropped, ame rows display 0 with missing SE, e(V_es_ame) withheld, generic caution only, VCE-invariant | (iv) | say-not-estimable: hard refusal of the ame block naming the cause (no within-cell dose variation) + misspecification-steer to a cross-cohort dose design outside this machinery; refinement: the missing SE is arguably "correct" for the degenerate constant actually computed ((ii)-flavored), but relative to the labeled ame estimand this is (iv); the fabricated 0 estimates must go |
| d3 basic cluster(unitid) | same | same, VCE-invariant | (iv) | same |
| d3 basic robust | same | same, VCE-invariant | (iv) | same |
| d1c basic cluster(i) | cell shock nu_gt is ABSORBED by the treated cell dummies ⇒ the s-part of the truth is UNCHANGED from d1; only Var(phihat) inflates (via control cells, where no dummies absorb nu) | cluster(i)'s ESTIMAND is the right one (within-cohort dependence real; robust's estimand is wrong here — it misses the within-cell covariance in the phihat channel); plug-in collapse identical to d1: rank-1 Var_cl(phihat)·11' with a LARGER common SE than d1's; per-horizon reported/true ratio improves as sigma_nu grows (both →Var(phihat)) but contrast variances are still exactly 0 and the joint tests still refuse — the irreducible signature | (iii) | code-fix, emphatically NOT a steer to robust: refuse the plug-in numbers with the structural caution and provide the cohort bootstrap (consistent for BOTH parts under cohort independence) or the hybrid (within-cell Var(s_gh) + CRVE Var(phihat), which is exactly right here since nu never enters s_gh); when the conceptually-right VCE is computationally degenerate the program must fail loudly with the design-based alternative, never print the collapsed numbers and never silently downgrade the dependence assumption |

## §3 Sharpest discriminating predictions (confirm vs falsify)

P1 — Three-way cluster(i) signature on the SAME design (d1 and d2), the core fingerprint of
the score-support account:
  basic/within → all 10 ame SEs one identical positive number (= SE_cl(phihat)); e(V_es_ame)
  posted, corr ≡ 1, rank 1; e(adversarial_pval_ame) missing in ~100% of reps.
  byCohort → all-missing SE blackout (20 of 20), NO e(V), NO e(b_es_ame)/e(V_es_ame), BOTH
  adversarial refused, in ~100% of reps.
  byTime → full-rank 8×8, horizon-VARYING SEs, adversarial-ame POSTS in ~100% of reps and
  rejects far above 5% (predict 25–60%).
  Falsifiers: byCohort showing equal-positive SEs like basic; byTime showing rank-1 or equal
  SEs; basic posting an ame p.
P2 — within ≡ basic bit-identity: per rep, identical b and V (mreldif 0) on d1 and d2 for all
  three VCEs. Any systematic difference falsifies my reading of contreatwithin as
  screen-relaxation only.
P3 — cluster(unitid) anywhere (and robust on d1/d2): full-rank V; per-horizon reported-SE /
  MC-truth-SD ratio in [0.80, 1.00] (leverage bias), per-horizon coverage 90–95%; adversarial-
  ame size 15–35% at 5% while the SAME reps' att-side adversarial sits ~5–8%; d1 vs d2
  qualitatively identical; byCohort-variant per-horizon true and reported variances ≈ 2× the
  basic variant's with cross-horizon corr ≈ .5 (exactly .5 in d1's truth).
P4 — d1c (cluster(i) DEFENDED despite the collapse): same rank-1/equal-SE/refusal signature
  as d1-basic-cluster(i), with a strictly LARGER common SE than the matched d1 case; MC-truth
  per-horizon variance minus d1's differs only through the phihat channel (cell shocks
  absorbed by dummies — the s-part contribution is unchanged); reported contrast variance
  still ~0 vs positive truth. The right program response is refuse+bootstrap/hybrid, not
  robust: robust on this DGP would understate the phi-channel and its estimand embeds a false
  independence assumption even where its number is accidentally closer.
P5 — d3: every rep, every VCE: ame point estimates exactly 0 with missing SEs, e(V_es_ame)
  absent, adversarial-ame refused, identical across the three VCEs (VCE-invariance is the
  tell that nothing about dose is being estimated at all).
P6 — MC-truth shape checks of the framework itself: in d1 the true per-horizon variances are
  ~equal across all 8 horizons (unit-constant dose ⇒ same weights every period) while in d2
  they vary; basic's truth cross-horizon correlations are small (.05–.2) vs byCohort's ~.5.

## §4 Stata-convention ruling

Official Stata also posts degenerate CRVE SEs on cluster-saturated dummies — but wooldid's
obligation is HIGHER, for three reasons. (1) Stata's output is a transparent primitive: the
degenerate SE sits next to the coefficient it belongs to, and Stata claims nothing about
aggregates. wooldid engineers the saturation deliberately and hides the baseline slope from
margins by construction (the ringer, L898-906) — the user cannot see that the "ame" SE is the
SE of an absorbed nuisance coefficient. (2) wooldid relabels the aggregate as an estimate row
("ame", "ES_ame pre2") with an SE column: a representation that the number measures the
sampling variability of the labeled quantity. It provably does not: the reported variance
estimator's estimand EXCLUDES a variance component (the within-cell slope noise) that
wooldid's own point estimator provably POSSESSES. A package may not post such a pair without
flagging it. (3) The F2 handoff: e(V_es_ame) exists precisely so honestdid/lincom/test can
consume (b, V) as a coherent joint object; posting a rank-1 V there converts a display blemish
into silently wrong downstream confidence sets — a contract plain Stata never enters.
wooldid already concedes the principle elsewhere: the e(V)-preserve caution, the f2ok
degenerate-row gate, and the adversarial degeneracy refusal all exist to avoid shipping
degenerate inference silently. The rank-1 case is the one gap: equal-positive collapsed SEs
currently print with NO caution (probe RUN1) and e(V_es_ame) posts. The convention defense
("regress does it too") fails exactly where wooldid adds its value — curation and aggregation
imply ownership of the coherence of what is posted. Minimum obligation: the cheap structural
predicate (every treated (g,relyr) cell inside one cluster of cluster()) triggering a named
caution + suppressed/flagged e(V_es_ame); the constructive obligation: bootstrap/hybrid.

## §5 Probes run (single-dataset deterministic structure checks; no Monte Carlo)

Scratch: /Users/hegland/.claude/jobs/640b9294/tmp/amese/thB/ (probeB1.do/.log, probeB2,
probeB3, probeB4). d1-style draw: 30 cohorts × 6 units × 12 periods, even cohorts treated
at t=6, unit-level dose ~ N(1,1), y = cohortFE + .05t + N(0,1), no dose effect.

- probeB1 (correct build, loaded via `do wooldidinstall/wooldid.ado`): RUN1 basic+cluster(i):
  V_es_ame ≡ .0015519795 (all 64 entries), corr ≡ 1, eig {.0124, ~1e-18}, all 10 ame SEs
  .0393952, adversarial-ame refused. RUN2 byCohort+cluster(i): 20/20 SEs missing, b-only e(),
  F2 withheld, BOTH adversarial refused. RUN3 byTime+cluster(i): full-rank (eig .038→.0022),
  diag .0050–.0236, corr −.66..+.44, adversarial-ame POSTED p=.00011 (true null).
  RUN4 within+cluster(i): mreldif vs basic = 0 for b AND V. RUN5 byCohort+robust: full rank,
  diag .0223–.0338, corr .38–.54, p=.013.
- probeB4: d3a (cohort-constant dose) and d3b (cell-constant, time-varying) × cluster(i)/
  robust: rc 0; interaction columns dropped; ame rows 0 with missing SEs ("10 of 20");
  e(V_es_ame) absent; adversarial-ame refused; VCE-invariant. d1 basic robust: full rank,
  diag .0131–.0298, p=.0243. d1 basic cluster(unitid): full rank, diag .0136–.0297, p=.0027.
- HARNESS HAZARD (found via probeB1/B2/B3 first attempts): a stale wooldid.ado (2024-07-31)
  in the PLUS adopath SHADOWS `adopath +` (which appends at lowest priority). The stale build
  lacks F2 entirely and POSTS g-inverse adversarial p-values on the degenerate V (observed
  p=.00002 on basic+cluster(i)). Any MC driver must load the repo build via `do` (or
  `adopath ++`) and can assert build identity by checking e(b_es_ame) existence.
