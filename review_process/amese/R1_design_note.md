# R1+R4+R3 implementation design note (author-ruled 2026-08-24)

Author ruling: implement R1 and what naturally follows (R4, R3); SKIP R2 (no
bootstrap SEs in wooldid, so none here); R1 must not interfere inappropriately with
randomization inference or IM p-values; a useful explanatory warning must print when
it fires; consider generalizable guards against posting SEs that do not represent
the labeled quantity. This note fixes the scope BEFORE implementation.

## Scope boundaries (load-bearing — verifiers must not extend them)

1. The refusal is AME-SIDE ONLY. The predicate is TRUE on wooldid's bread-and-butter
   binary configuration (i = unit id, cluster(i), singleton cells), so an att-side
   refusal or warning would fire on essentially every standard run. The att-side
   understatement (G1a, NF-1) remains a documented ticket (R5/R8 + standing NF-1),
   NOT part of this change.
2. Predicate (referee-adjudicated ANY semantics): the run uses a cluster-type VCE
   (cluster() or customvce cluster) AND contreat() is active AND ANY treated
   (_i,_t) cell in the realized estimation sample lies wholly inside a single
   cluster of that cluster variable. Computed once, post-markout, never on the
   fastnose path (permute reps are fastnose).
3. On fire: (a) every ame-family SE/t/p/CI (main ame, pretreatment ame, all ES_ame
   horizons, per-subgroup ame rows) is routed through the EXISTING per-estimate
   not-identified machinery (47defa0 sentinel placeholders + the existing "standard
   errors were not identified for X of Y" caution plumbing) — no parallel mechanism;
   (b) e(V_es_ame*) withheld while e(b_es_ame*) STILL POSTS (splitting the F2
   single gate — points are valid, the garbage V is not); (c) R3 companion: the
   ame-side V_es-consuming tests are suppressed with the same explanation —
   adversarial ame p, joint_Choriz*, pretrendpower on ame blocks. ATT-side tests
   unchanged this chunk. (d) The explanatory warning (below) prints once.
4. Non-interference requirements (probe-verified before landing):
   - permute(): rep runs are fastnose (no V) — the predicate code must be
     unreachable there; observed-run warning prints once; perm_* outputs and
     permtstat (att-scale) untouched.
   - impvalues: IM SEs are cross-cohort dispersion — the estimator class that
     SURVIVES the pathology. Never suppressed by R1. If the probe confirms
     impvalues x contreat produces valid ame IM results, the warning steers to it;
     otherwise the steer is omitted.
   - customvce(cluster v): predicate uses v. unconditionalse: probe determines
     whether the pathology persists (expected yes) and the predicate applies.
   - poisson/ppmlhdfe leg gets the same treatment (probe).

## R4 (blackout repair)

Detection is UNCONDITIONAL: after the internal estimation, if the stored e(V) is
identically zero (the Stata `ereturn post` screen signature — asymmetric/highly
singular V zeroed WHOLE; wooldid's qui swallows the only warning), print an
explanatory caution naming the mechanism. Restoration of the att block via the
explicit-dummies regress fallback is GATED on a feasibility check (absorbed-levels
count — the fallback materializes FEs and would hit r(907) at SI-8 scale); when
infeasible, the message covers it ("att variances were destroyed by Stata's variance
screen; re-run with robust or a finer cluster() to recover them"). Restored att SEs
under a firing R1 predicate carry the G1a caveat, never certified-valid.

## Warning text (drafts; final wording in the commit)

R1 fire: "Caution: with cluster(<v>), every observation of <k> treated cohort-by-
period cell(s) lies within a single cluster, so the cluster-robust variance cannot
see within-cell variation in <dose>. The continuous-treatment (ame) standard errors
this would produce do not measure the sampling variability of the ame estimates
(they collapse to the variance of the absorbed baseline dose control), so wooldid
has set the ame standard errors, tests, and e(V_es_ame*) to missing. Point
estimates are unaffected and remain reported. To obtain ame inference, use
vce(robust) [valid when errors are independent across observations; not under
within-unit serial dependence] or cluster() on a variable finer than the
cohort-by-period cells[; impvalues provides cluster-level inference via
cross-cohort dispersion — steer included only if probe-verified]."
V==0 fire: names the Stata screen, states what was lost, and the recovery options.

## Generalizable guards (capped at two, both CAUTIONS not refusals)

G-a: V==0 detector (= R4's detection; also unswallows the reghdfe warning).
G-b: identical-SE collapse signature: if >=3 reported SEs within one estimand
family are pairwise identical at reldif <= 1e-12 UNDER A CLUSTER-TYPE VCE, print a
generic degenerate-VCE caution. Gated on clustered VCE because balanced unclustered
designs produce legitimately identical SEs (fixture-verified). This is the
belt-and-braces layer for unknown paths; the structural predicate is primary
(rank/eigen guards are known-spoofable and are NOT used).

## Sequencing collision: T-M2 / cont_sg

R3 suppresses joint_Choriz* wherever the predicate fires. If dta_cont's cells are
cluster-contained (fire-list probe confirms), the T-M2 defect sites stop executing
on cont_sg — the SOLE baseline pinning the defective "." — and cont_sg.csv gets
re-blessed in this wave regardless. Therefore: BANK the T-M2 mandated evidence
(>=3-dataset console + mid-run r(p)/r(b) matrix dumps, per
review_process/mini_assessment_2026-08-23.md §2b) BEFORE any ado change, archive
under review_process/, and surface to the author the option of bundling the
[1,2]->[1,1] fix into this same re-bless (not fixed without their ruling).

## Verification signature (preregistered before the battery runs)

- battery: EXACTLY the fire-list configs differ (re-blessed in-commit, justified
  row by row); every other config bit-identical at rtol 0.
- oldnew: new LARGE ledger rows scoped to the fire-list; 0 stale / 0 unscoped after.
- cert: a new test_*.do pinning — fires-with-message on the demo design; the ANY-
  semantics acid test on the committed referee mixed-cluster data (refA.dta);
  byTime-now-refuses; robust/cluster(unit) unchanged; permute end-to-end clean with
  one warning; fastnose no-fire; impvalues behavior per probe; V==0 message; b_es
  posts without V_es on fire.
- then a 3-verifier wave on the IMPLEMENTATION (the diagnosis is already reviewed).

## Bookkeeping on completion

demo_vesame_pathology.do gains a post-fix header note (Exhibit A output changes);
R1/R3/R4 move to wooldid_solved_issues.md leaving R5-R8 standing; the sthlp gains
the warning text + design-conditions paragraph (the R8 subset that naturally
follows); standing-issues + memory updated.
