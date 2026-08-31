# wooldid standing issues — the bug-fixing work list

Derived 2026-08-24 from `wooldid_old_bug_report.md` (the 2026-08-23 report of record
plus its §5 new-feature-validation extension). This file lists ONLY what remains to
be addressed; nothing here is fixed yet. (Status annotations added 2026-08-25:
items marked IMPLEMENTED/DONE/ADDED are complete and await only author review;
items marked DIAGNOSED/REPORTED BACK carry proposals awaiting an author ruling.) Excluded by construction: everything in §1
(fixed old bugs), §3 (intended changes), §5.1 (validated-clean surfaces), the pure
records of §2 (items 4, 5, 8: environment pins, the corrected-T1 reading, parity
pins), and the author-ruled ignores from the mini assessment (T-M1; T-M3 waits for
the efficiency pass). Evidence pointers cite the report section, the INV-1 mechanism
report (`review_process/inv1_v_es_ame_2026-08-24.md`), the mini-assessment protocol
(`review_process/mini_assessment_2026-08-23.md` §2b), and the runnable demo
(`review_process/demo_vesame_pathology.do`).

---

## Working principle on effort allocation (author-stated 2026-08-27)

**Do not spend substantial effort, or accept material code complexity, to improve behavior
that only occurs on datasets so pathological that a good result is impossible anyway.**

The author's reasoning, recorded because it should govern future triage: users build
specifications up incrementally and debug by subtraction. Someone who hits a failure on
`att subgroups(g)` will almost certainly have run the same command without `subgroups()`
first, and when something breaks the standard move is to drop options one at a time until
it works. So on a genuinely degenerate dataset the user diagnoses the situation quickly on
their own, and the marginal value of a more polished message in that corner is small —
often smaller than the cost of the code that produces it.

SI-7 is the worked example. The crash was real and its message was cryptic, but it fired
only where the variance-covariance matrix had already collapsed, so nothing usable was
being lost. The fix was kept because it turned out to be cheap and inert (fires only where
the old code aborted; battery bit-identical; one fixture's diff), not because the polish
justified the investigation that preceded it. Had it required invasive changes or ongoing
maintenance, the correct call would have been to leave it alone.

Applying this: before opening an investigation, ask whether any dataset reaching the defect
could have produced a trustworthy answer. If not, cap the effort - document the behavior
and move on. Keep such a fix only when it is safe, small and self-contained; otherwise
toss it.

## A. The V_es / AME inference pathology family — FULL DIAGNOSIS (2026-08-24 grid)

Diagnosed by author directive via a preregistered 26-case x 256-rep MC-truth grid
(designs x contreatcontrols x VCE; PREREG_AMESE_2026-08-24.md, bands locked before
results), a deterministic structure pass, an att-side supplement (2 x 258, paired
seeds), TWO locked-prediction theory reviews (incidental-parameters; design-based),
a blind third theorist, a code-path audit, and an adversarial referee whose required
edits are incorporated below. Artifacts: review_process/amese/ (theoryA_locked,
theoryB_locked, blind_theorist, code_audit, referee, mc_results) + the INV-1 report.
Bucket definitions applied under ONE lens throughout (the referee's mandate): a
"correct SE" is one that tracks the sampling variability of the LABELED quantity
under the design; where a within-method reading would change the bucket, it is noted.

### A.1 Final classification (the author's four buckets)

| grp | design x options | measured (SE/truth ratios etc.) | bucket | verdict detail |
|---|---|---|---|---|
| G1 | cluster(i) + basic (unit-dose or time-dose data); also within-arm on unit-dose | e(V_es_ame) exact rank 1; every ame SE = SE_cl(pooled phi-hat); ratios: ES horizons ~0.20, pretreat ~0.40, main ~0.52 (ratio^2 tracks #cells averaged); SILENT | (iii) | The AME estimand is fine and a correct SE is estimable (R2) — but not by plug-in CRVE here (cell scores identically zero). The posted number is a CORRECT SE of the WRONG object (the absorbed baseline dose slope) mis-attributed to every ame quantity — the named nearest-(ii) miss; under a strict within-requested-method reading this is (iv). |
| G1a | ATT side under cluster(i), same designs (att supplement, paired seeds) | att ratios 0.66-0.75 across all 10 att quantities; robust att 0.93-1.04 | (iii) | The pathology is NOT ame-only: cohort-clustered att SEs are ~25-30% understated (milder than ame's ~80% because att cells lack the slope channel). This is the microfoundation of NF-1's att-side 20%+ joint-test rejections. |
| G2 | cluster(i) + byTime | FULL-RANK V, horizon-varying SEs — but each horizon SE = SE_cl(phi_t-hat) to all displayed digits (confirmed on two shapes/seeds); ratios 0.38-0.69 | (iii) | The full rank is SPOOFED by the period-varying baseline block; within-cell noise still 100% annihilated. Kills any rank/eigenvalue guard design. |
| G3 | cluster(i) + byCohort | total blackout: reghdfe e(V) arrives bitwise ZERO (Stata's `ereturn post` screen zeroes a nonsymmetric/highly-singular V WHOLE); "not identified 20 of 20" caution; no e(V_es_ame), no e(b_es_ame) | (iii) | ame side = G1. The att collateral: att variances are computable via the regress path but are themselves G1a-understated (~0.7) — so the R4 repair restores them only WITH the R1 caveat, never as certified-valid. The honest caution is the right instinct; losing att inference entirely is not. |
| G4 | all-treated staggered + cluster(i) (within-style) | stochastic blackout 95/256 reps (37%; 3 of 4 frozen realizations); posting reps are G1 rank-1 (all SEs = SE_cl(phi-hat) exactly); main/pre ratio ~1.0 is a COMPOSITIONAL COINCIDENCE, ES ~0.61-0.76 | (iii) | All-or-nothing numerical behavior on an over-parameterized design; R1's structural predicate makes handling deterministic. |
| G5 | cluster(unitid) or robust, d1/d2, all contreatcontrols | ratios 0.68-0.97 (typical 0.72-0.91); full rank, horizon-varying; one case at the WRONG boundary | (iii) [relabeled from "(i) with caveat" on referee ruling: the preregistered band itself excludes "approximately correct"] | Textbook CR1/HC1 small-cell leverage bias (6-obs cells, 2 within-cell params): ~8-13pt under-coverage at the margins. Displayed-SE remedy is R5 (caveat + opt-in correction); the joint tests built on these V (NF-3: 20-35% at 5%, measured in the NF-1/NF-3 MCs, not this grid) require R3. att under robust is essentially correct (0.93-1.04). |
| G6 | all-treated + cluster(unit)/robust | ratios 0.86-1.11 — inside the preregistered band | (i) pointwise | Compositional improvement (truth weight shifts onto the well-estimated baseline block), not a repair of the leverage bias; joint tests remain uncertified. Shape-qualified. |
| G7 | dose constant within cohort-period cells (d3), any VCE | deterministic: 165 collinearity drops; ame rows shown as structural-ZERO point estimates with missing SEs behind generic cautions | (iv) | The within-cell AME estimand is genuinely unidentified, and no option rescues it (lattice only refines cells — audited). Correct output: a NAMED refusal of the ame block + steer (R6). wooldid rightly posts no SE; wrongly presents fabricated 0 estimates unflagged. |
| G8 | contreatwithin on data WITH controls | bitwise identical to basic (estimation-inert) | n/a | A user requesting within-identification is silently not getting it — screen or document (R7). |
| G9 | cluster(i) + cell-correlated errors (d1c) | = G1 (rank-1, ratios 0.18-0.46) with knife-edge stochastic blackouts (33/256 MC reps ~13%; 1 of 16 frozen realizations); robust survives at 0.75-1.05 because treated-cell dummies absorb the cell-level shocks | (iii) | Sharpest normative case: cluster(i) is the CONCEPTUALLY right level yet numerically worst. Default remedy = a valid cluster-level estimator (R2); a steer to robust is CONDITIONAL on the dependence model (valid for pure cell-level shocks, invalid under unit-serial dependence) and must be worded so. |

Bucket (ii) is EMPTY on this grid; the two nearest misses are named in G1 and G3.
Old wooldid shares every pathology above (same spec/VCE paths) and additionally
posted the invalid g-inverse adversarial-ame p that current wooldid refuses.

### A.2 Theory scorecard (locked predictions vs measurement)

Confirmed by all three theory reads: the cluster(i) three-way signature (rank-1
silent / byCohort blackout / byTime full-rank-posted), within==basic inertness, d3
unidentifiability, robust/unit-cluster leverage-bias direction. Quantitative scoring:
theory A's cu/rb band 0.65-0.85 vs measured 0.68-0.97 (best); theory B's 0.8-1.0
slightly high-biased; both under-predicted the cluster(i) severity (measured ES ~0.2).
Resolved disagreements (referee-mandated, now explicit): the G5 bucket goes to the
blind theorist's (iii) over the draft's "(i)-with-caveat"; the d1c robust-steer goes
to theory B's position with the referee's conditionality refinement. Measured but
unpredicted by anyone: the ereturn-post V-zeroing blackout channel; the compositional
main/pre ~1.0 coincidence; the partial-containment false negative (below).

### A.3 Remedy set R1-R8 (renamed to avoid feature-ID collisions), with response type and blast radius

| id | remedy | response type | blast radius |
|---|---|---|---|
| R1 (P1) | STRUCTURAL detect + refuse: predicate = ANY treated (i, relyr) cell wholly inside one cluster of cluster() ["any", not "every" — the referee DEMONSTRATED a partial-containment false negative: mixed cluster variable, full-rank V, no caution, spurious p=1.2e-08 on a true null]. On firing: suppress plug-in ame SEs, e(V_es_ame), and the V_es-based ame tests; print a design-specific caution naming mechanism + remedies. Belt-and-braces: also screen margins r(table) for the blanking behaviors. NEVER guard on rank/eigenvalues (G2). | say-not-estimable | wooldid.ado validate/worker; battery cont_* pins re-bless (cluster(i) configs change output); oldnew ledger rows; sthlp. LARGE — the flagship change. |
| R2 (P3) | Replacement estimator, opt-in first: cluster bootstrap of the margins functional, or analytic cross-cohort dispersion of per-cohort AMEs (IM-style). Audit confirms ingredients are in hand at the margins site (cell slopes in e(b), weights in r(Jacobian), Var_cl(phi-hat)). Context: the F1/boottest wild-bootstrap route was previously ruled infeasible for the stored model — this is a different, self-contained bootstrap of the margins functional. | code-fix (new option) | new option + newfeat validation; no re-bless until made default. |
| R3 (P2) | Joint tests on V_es (adversarial, jointtests, roth prep, NIS size, AND the pretrendpower mdt columns — input validity): refuse under the R1 predicate; elsewhere (small-cell cu/rb) bootstrap/CR3 the V or refuse above a dimension threshold. 20-35% rejection at nominal 5% is not caveat material. | code-fix | advstat/jointtests/ptp paths; battery adv_*/esf_adv/ptp_two pins re-bless; NF-1/NF-2/NF-3 tickets close together. |
| R4 (P1, gated) | Blackout repair: detect e(V)==0 after the reghdfe call (the ereturn-post zeroing signature; unswallow the warning wooldid's qui hides), restore the ATT block only, via the existing regress fallback path — under the R1 predicate the restored att SEs carry the G1a caveat (~0.7 understatement), and the ame side stays refused (referee: an ungated repair would resurrect the rank-1 ame SEs the zeroing usefully suppressed). | code-fix | worker post-estimation block; G3/G4/G9 blackout behavior changes; cautions re-worded. |
| R5 (P4) | Small-cell leverage bias (G5): sthlp caveat + printed note when cells are small; consider opt-in CR2/HC3-style corrections; default unchanged. Magnitudes are shape-qualified (single 30x6x12 shape measured; ratios should approach 1 with larger cells). | doc + opt-in code | sthlp; optional new sub-option. |
| R6 (P4) | d3-class: NAMED refusal of the ame block ("dose has no within-cell variation; ame not identified") + steer, instead of structural-zero rows behind generic cautions. The existing catastrophe counter already detects the condition. | say-not-estimable + steer | one screen + message; ext/patho configs re-pin where they hit it. |
| R7 (P4) | contreatwithin-with-controls: screen or document the estimation-inertness. Separate sub-ticket: the byCohort+ame POINT estimand is normalization-dependent across solvers (audited: reghdfe .723 vs regress .554, att invariant) — decide what that estimand IS before offering SEs for it. | misspecification-steer / design decision | screens + sthlp; the byCohort sub-ticket may alter estimates (author decision). |



| R8 (P4) | sthlp design-conditions rewrite for ALL V_es consumers (honestdid/pretrendpower handoffs included), folding the singleton-cell O16 caveat; plus a user note that a stale 2024-07-31 wooldid.ado in the PLUS adopath can shadow interactive sessions. | doc | sthlp only. |

Priority: P1 = R1 + gated R4 (stops silent-wrong and repairs the blackout collateral
in one coherent behavior change); P2 = R3 (test validity); P3 = R2 (restores
inference the refusals remove); P4 = R5-R8. R1/R3/R4 are one design decision
executed together — they share the predicate, the messages, and one re-bless wave.

### A.4 Method limitations (recorded)

Single panel shape (30 cohorts x 6 units x 12 periods): G5/G6 magnitudes, blackout
rates, and the d1c-robust performance are shape-sensitive; bucket assignments and
mechanisms are not. T=256 gives ~4.4% SD noise (no verdict is within noise of
flipping). Joint-test rejection rates cited from the NF-1/NF-3 MCs, not this grid.
Harness: expand/bysort tie order after a wooldid call makes individual reps
non-reproducible standalone for ame quantities (distributionally harmless — unit
labels exchangeable; att invariance audited); the structure pass's "1x1 rank 0"
lines mean e(V_es_ame) ABSENT (matrix-coercion artifact, interpretation unaffected).

### A.5 The author's coarse() matched-pair example (2026-08-24 addendum)

The author's bonus example (contreat(dose) contreatcontrols(basic) vs
coarse(c.dose), same data, robust / cluster(unit) / cluster(i)) was run verbatim
and probed (review_process/amese/probes/ bonus.do/bonus2.do). It introduces no
revision to the classification — it is a clean single-dataset demonstration of the
diagnosis's central split, and it sharpens three answers:

1. WHAT IS GOING ON: the two cluster(i) legs are the two faces of ONE mechanism.
   coarse(c.dose) puts a visible dose slope in EVERY (_i,_t) cell (traced: 
   `i._i#i._t#(c.dose)`, incl. the pooled-control cells) with no absorbed pooled
   ringer; its cluster(i) CRVE is so thoroughly score-dead that Stata's
   ereturn-post screen zeroes it WHOLE (probe: reghdfe e(V) bitwise zero on all
   577 diagonals) -> wooldid's honest "not identified 4 of 4" blackout = regime
   G3/G4. contreat's cluster(i) CRVE survives the same screen because the absorbed
   pooled ringer routes control-cohort dose variation into one clean variance
   block -> the silent regime G1: ameSE = 59.04 on this draw = SE_cl(phi-hat)
   identically for every ame quantity (rank 1 of 2), and the "normal-looking" att
   SE (143 vs 170 robust) is the G1a ~0.7-understated att, not a healthy one.
   Which face a design lands on is numerical happenstance (G4 tripped 37% of
   reps), not a substantive difference.
2. NOT MISSPECIFICATION, AND THE byTime APPEARANCE IS AN ILLUSION: the DGP is
   null and the AME estimand is well-defined — nothing here is model
   misspecification; it is the variance estimator's score support. The
   observation that byTime-style contreat controls "don't show issues" is exactly
   the G2 trap: byTime+cluster(i) posts full-rank, horizon-varying, plausible-
   looking SEs that equal SE_cl(phi_t-hat) to all displayed digits and measure
   0.38-0.69 of truth. It LOOKS the healthiest and is the most dangerous member
   of the family.
3. DOES coarse() IMPLY NO SEs SHOULD BE PRODUCED FOR contreat UNDER cluster(i)?
   For the plug-in CRVE: YES — coarse's blackout is the honest answer for both
   legs, which is precisely remedy R1 (make contreat's cluster(i) behavior match
   the refusal, deterministically and with a design-specific message). In
   principle: NO — a correct cohort-level SE exists (R2: cluster bootstrap of the
   margins functional, or cross-cohort dispersion of per-cohort AMEs), so the
   refusal should name the alternative. On the "what would the design-based
   econometricians say" question: the Abadie-Athey-Imbens-Wooldridge clustering
   framework points the same way — if cohort-level dependence/assignment is the
   reason to cluster, valid variance must come from BETWEEN-cohort variation of
   cohort-level statistics, which the plug-in CRVE on a cell-saturated regression
   cannot see; and with ~15 treated cohorts any cohort-level inference is
   few-cluster inference, where bootstrap-t/IM-style dispersion methods are the
   standard recommendation regardless.
   Note on the example's middle pair: `cluster(unit)` there is unit = 1..6 within
   cohort, i.e. SIX cross-cutting clusters — cells span clusters so no
   annihilation (hence "sensible"-looking), but 6-cluster CRVE inference has its
   own small-G caveat; the grid's cluster(unitid) (180 unique units) is the
   better-behaved analogue and sits in G5 (~0.7-0.97 with leverage bias).

### A.6 IMPLEMENTATION RECORD (2026-08-25): R1 + R3 + R4 LANDED; R2 skipped by ruling

The implementation pool is named the **Group A Patches** — the canonical commit
manifest with tiers and the revert procedure is
`review_process/amese/GROUP_A_PATCHES.md` (author is testing these in extra cases;
revert stays a clean single step).

Author ruling implemented across commits 7054878 (core), ea69544 (R4 narrowing),
c184ea0 (cert convention), the V-wave fix commits, and the AMESE ledger commits:
- R1: structural cell-containment predicate (ANY treated cohort-by-period cell wholly
  inside one cluster of cluster()/customvce-cluster incl. Stata abbreviations), fired
  only with contreat and never under fastnose/impvalues/unconditionalse. On fire: an
  explanatory caution names the mechanism and the three valid alternatives
  (impvalues — probe-verified valid dispersion inference; conditional robust; finer
  clustering); every ame-family se/t/p/ci cell is scrubbed to missing across main/
  pretreatment/ES/subgroup/pairwise/cfx surfaces (plots drop CI bands via noci);
  e(V_es_ame*) withheld while e(b_es_ame*) posts; ame-side V_es tests skipped
  (adversarial ame, joint_Choriz*, SG-equality ame chain; ptp ame keys self-gate).
- R4: unconditional detection of Stata's ereturn-post V-zeroing (with an explanatory
  caution; the swallowed reghdfe warning surfaced); an att-rescue refit through the
  explicit-dummies regress path, gated NARROW: firing predicate + single cluster() +
  no customvce/poisson/fe()/ccc_absorb + maxvar feasibility. The former byCohort
  20-of-20 blackouts now return att inference.
- Generic guard: the identical-SE collapse signature caution (>=3 identical SEs with
  DIFFERING point estimates under a cluster-type VCE).
- R2 (bootstrap/dispersion estimator): SKIPPED per author ruling.
Verification: 10-scenario cert suite (fires incl. the referee's partial-containment
acid data; silent under robust/cluster(unit)/impvalues/fastnose/binary); battery
34/34 rtol-0 with the 8 fire-list configs re-blessed; full oldnew 450/450 with 26
scoped AMESE ledger rows, 0 unexplained/stale/unscoped; a 3-verifier wave
(review_process/amese/verifier1-3) whose findings were all fixed or ruled:
cfx-ame leak fixed; refit narrowed off fe()/ccc_absorb (bare-factor linearization);
glevelsof capture (rc-920 hard-abort at scale); sescan differing-estimates
signature (two demonstrated false-positive classes cured); cluster-abbreviation
parse; ledger reason/bound amendments. T-M2's mandated evidence was banked PRE-R3
(review_process/tm2_evidence/, commit a603f38) since R3 suppresses the defect sites
on the pinning config.

STILL OPEN in family A after this implementation:
- SI-A1 (was NF-1b/att side): att-side understatement under cohort-level clustering
  (G1a: 0.66-0.75) and the att joint tests' oversizing — deliberately NOT touched
  (standard-usage surface; needs its own ruling). V1's DQ6 (pwsgs att rows with
  se ~8.7e-09 / p ~1e-188 uncautioned) corroborates.
  - **JOINT-TEST SIDE ADDRESSED 2026-08-25 (author-ruled three-part fix).**
    (1)+(2) sthlp: the jointtests entry now documents both regimes — continuous
    case correctly sized only with many observations per i-t cell (with the
    lattice()/contreatpoly rider), and incompatibility with analytic clustered
    SEs at the cohort-i level or above (steer: bootstrap clustering or
    randomization inference); the safety() paragraph lists the new check.
    (3) code: a pre-estimation screen in the worker — jointtests + any cluster
    variable (or cluster-type customvce) in which `i` is nested → error 198
    under default safety with the mechanism + permute()/bootstrap steer +
    override mention; under safety(off) the same text prints as a caution and
    estimation proceeds with the p-values computed. Gated on `safetyinput`
    (contreatwithin's internal safety overwrite does not demote the error);
    nestedness is a per-cluster-variable gegen tag-count on the original `i`;
    fires on ANY member of a multi-way cluster(); validate-stage contract
    errors (fastnose, impvalues×jointtests) still fire first by construction.
    Cert: test/test_jointclus.do (8 scenarios incl. two-way and customvce).
    Suite: trig_jtclus pins the refusal 0→198; eight cluster(i)+jointtests
    configs (basic, sg_attit, sg_attitime, cont_sg, sg_missif, sg_degpre,
    cov_kitchen, adm_base) and six battery configs moved to safety(off) so
    their estimation surfaces and scoped ledger rows stay alive under the
    caution path. STILL OPEN in SI-A1: the pointwise att SE understatement
    (0.66–0.75) and the adversarial test (author: set aside for now).
    VALIDATED 2026-08-25: full suite ALL COMPONENTS PASS first run (454/454;
    0 unexplained / 0 stale; cert 24/24 incl. test_jointclus; battery
    IDENTICAL after the cmdline-only re-bless; trig_jtclus judged
    rc-by-design 0→198; the two TM2 rows also fired on this run).
  - **ADVERSARIAL SIDE MEASURED 2026-08-26 (AB campaign; report of record
    `test/suite/mc/REPORT_ADVBOUND_2026-08-26.md`, prereg
    `test/suite/PREREG_ADVBOUND_2026-08-25.md`).** 24 configs / 17,056 reps,
    preregistered with locked theory predictions (25 of 26 inside band).
    RULING: binary `adversarial` = KEEP, GATED — it is correctly sized at
    4.80% (6-unit cells, unit clustering), 4.60% (robust), 5.00% (20-unit
    cells) and 5.60% (AR(1) with matching clustering); power 98.8-100%. The
    governing parameter is CLUSTERS PER CELL, not cluster count or cell size
    (the campaign's decisive result: 2 clusters/cell = 14.50% FAIL even with
    60 clusters and 6-unit cells). Fails: containment 56.35%, cohort-period
    shocks 12.30% at only 3.8% of error variance, 20 clusters 17.46%,
    staggering mildly (7.20%, chi2-flagged). Robust under AR(1) is INERT not
    safe (0.00%, chi2 1072). `altmethod_adversarial` is uniformly better at
    EVERY boundary and certifies once where standard fails (2-unit cells,
    6.94% vs 9.10% paired) — but at the exact band edge (35/504), so that one
    claim is provisional pending M>=2000. ame-side adversarial = WITHDRAW
    (9.92% at best geometry); and the clustering advice INVERTS between sides
    (fine clustering rescues att, harms ame: 48.41% vs robust 31.55%).
    Ship recommendations in the report; no code written pending author ruling.
  - **CLOSED 2026-08-27 by author ruling — DOCUMENTATION ONLY, no runtime gate.**
    The author ruled the adversarial line tests properly addressed by an sthlp
    note that they are particularly sensitive to small cluster counts and tend
    to be severely oversized when each i-t cell does not contain MANY clusters;
    explicitly NO runtime tests for cluster counts. Implemented in the
    adversarial() entry, grounded in the measured campaign (>=6 clusters per
    cell near nominal; ~3x oversized at 2 clusters per cell; order-of-magnitude
    at 1; cell size immaterial once the cluster condition is met and no repair
    when it is not; small overall cluster counts degrade independently). Two
    further conditions were included in the same paragraph and are flagged to
    the author as ADDITIONS beyond the literal ruling, easily struck: the
    continuous/ame version additionally needs cells of order 50 observations,
    and a cohort-by-period common shock is invisible to every analytic VCE
    (worsening with cell size), with a steer to permute()/bootstrap in both
    cases. Campaign artifacts retained: REPORT_ADVBOUND_2026-08-26.md,
    PREREG_ADVBOUND_2026-08-25.md, adversarial_calibration.html, advb_grid.json,
    and the mc_advb/run_advb/analyze_advb harness.
  - **CLOSED 2026-08-25 by author ruling (pointwise att side).** No caution for
    the standard binary run: the AME case is covered by the containment screen,
    the joint-test case by the SI-A1 screen, and the author's judgment is that
    the residual pointwise exposure collapses into cases where clustered SEs
    are known-problematic anyway (e.g. small cluster counts). The targeted
    binary pointwise MC was offered and declined. SI-A1 is now closed except
    for adversarial, which the author is taking up separately.
- SI-A2 (= R7): byCohort ame POINT estimands are solver-normalization-dependent
  (audited reghdfe .723 vs regress .554; V1 demonstrated O(1) shifts under the R4
  refit — the caution now discloses it). PARTIALLY GUARDED 2026-08-25 by the Group B
  screen (7e61af6, review_process/amese/GROUP_B_PATCHES.md): zero uncovered-pool
  dose variation now REFUSES r459 naming the cohort; thin variation cautions with a
  regtol() steer. Remaining open: the estimand-definition question itself (the
  well-identified byCohort ame is fine; what the near-threshold decomposition MEANS
  is still an author call) and contreatwithin-on-controls inertness.
  - **CLOSED 2026-08-25 by author ruling:** the imposition of identification
    requirements (the Group B screen) implicitly resolves this; closure covers
    both the estimand-definition question and the G8 contreatwithin-inertness
    observation that lived inside it.
- SI-A3 (= R5): small-cell leverage bias under robust/cluster(unit) (0.68-0.97) —
  sthlp caveat + opt-in CR2/HC3 consideration.
  - **CLOSED 2026-08-25 by author ruling:** one-sentence sthlp warning added at
    the contreat() entry (i-t cells need many observations relative to the
    within-cell parameters estimated; the jointtests entry already carries the
    joint-test version). No opt-in correction: intuitive to most researchers.
- SI-A4 (= R6): d3-class named refusal (unidentified ame presented as 0-rows).
  - **IMPLEMENTED 2026-08-25 per author ruling.** Detection rides the existing
    multicollinearity drop screen (the audit's own fix-site): a compact scan of
    e(b) colnames counts live vs omitted treated-cell dose terms — omitted
    continuous interactions carry the `co.` marker (a naive `c.`+`o.` match
    never fires; caught on the frozen dataset) — and a data-shape predicate is
    provably WRONG (probe of record: every suite continuous dataset has zero
    within-cell dose variation, because ordinary usage is singleton cells where
    the ame is identified through the slope parameterization; only the DROPS
    mark the d3 class). On all-dropped + none-live: a named caution (mechanism
    + coarser-i / drop-contreat steer + safety(off) escape note), then the ame
    block is voided — the R1 scrub machinery reused with its column floor
    lowered to include the POINT column (fabricated structural zeros must not
    print), and e(b_es_ame*) withheld along with V under void (not trusting
    f2ok=0 across the class). att untouched by construction. untreatedZero
    exception implemented structurally: the scan scopes to targetflag==1
    terms, so expected zero-dose reference periods never count (cert C2 proves
    no false fire under DEFAULT safety). Verification of record on
    review_process/amese/probes/d3_r1.dta: old = silent fabricated ame zeros;
    new = named caution + ame rows missing; att rows bit-identical. Cert:
    test_amescreen C1/C2 (ALL_PASS); suite: trig_d3void + SIA4-d3 ledger row.
    untreatedZero SE answer (author asked): under cohort clustering, cont_uz
    sits in the AMESE R4 family — OLD suffered a total ereturn-post V blackout
    (b-only); NEW's containment+R4 machinery detects, refits, and revives att
    inference with the ame side cautioned (rows AMESE-R4-l/r/p-cont_uz).
    Under robust/fine clustering untreatedZero inherits only the standard
    small-cell condition now documented by SI-A3's sentence; nothing
    untreatedZero-specific beyond that in any campaign measurement.
- SI-A5: non-containment V-zero blackouts (e.g. lattice+cluster(unitid) 4-of-4) are
  pre-existing shared behavior, now with an explanatory caution but NO refit (out of
  the narrow gate by design); revisit if the author wants rescue there too.
  - **CLOSED 2026-08-25 by author ruling:** no rescue needed; the explanatory
    caution stands as the final behavior.
- TESTING RULE (author-clarified 2026-08-25, supersedes the earlier deletion
  recommendation): the wooldid copies in the Stata PLUS adopath and at
  /Users/hegland/Documents/Research/wooldidinstall/ (outside the repo) are the LIVE
  external version BY DESIGN — they mirror what outside users access and what Stata
  loads by default, and update only when the repo's edits go live. They must NOT be
  modified. Consequence for all testing in this repo: never trust default adopath
  resolution — always load the development build explicitly via
  `do <repo>/wooldidinstall/wooldid.ado` (this bit twice during the Group A wave:
  agent probes and a cert-suite path slip silently ran the live build).

### Legacy ticket mapping

SI-1 (adversarial size) and SI-2/SI-3 (jointtests/ptp scalars) are consequences of
G1/G1a/G2/G5 and close via R3 (+R1). SI-4 (demo issue) = G1's silent rank-1 regime —
the demo file remains the canonical reproduction; its fix is R1 (+R2 for restored
inference). SI-5 (ame oversized under every VCE) = G1+G5 jointly; closes via R1+R3+R5.

## B. Functional defects on both builds (invisible to the old-vs-new comparison)

- **SI-6 — `ccc_absorb()` broken on both builds.** The regression completes, then a
  bare "invalid syntax" r198 in the margins phase (`cov_ccc_absorb`). Standing
  defect since REPORT_2026-08-18 ticket 2; the T-M1 ccnoabsorb corner note is
  adjacent context. §2.1.
  - **REPORTED BACK 2026-08-25 (author challenged; investigation confirms the
    author on substance).** The author's cohort-level example runs clean on both
    builds (ATT 1.42 vs true 1.5; `ccc_absorb(dose)` and `coarse(c.dose)`
    point-identical). The failing config's minimal ingredient is `ccc_absorb()` +
    **unit-level `i`**: `absorb(i._i#i._t#c.(x))` puts a control slope in every
    treated cohort-period cell, and with unit-level `i` each such cell is a single
    observation — reghdfe drops all treated rows as singletons, so the spec is
    unestimable **by construction**, not an estimator bug. The defect is only the
    delivery: with no treated dummies in e(b) the horizon scan leaves `relyrlist`
    empty and `matrix ... = ()'` throws the bare r198 (the dropped-dummy screen
    misses it because the levels are absent, not `o.`-marked). PROPOSED (awaiting
    ruling, not implemented): an informative treated-sample-survival error
    steering to coarser `i`/fe()/coarsecohortcontrols, an empty-`relyrlist` guard,
    and a passing cohort-level `ccc_absorb` suite config (the option has zero
    passing coverage today). Full report:
    `review_process/si_reports/si6_ccc_absorb_report.md`. No T-M1 interaction.
  - **CLOSED 2026-08-25 by author ruling.** The sthlp ccc_absorb() entry now
    carries the author's exact caution ("be extra careful to ensure that
    specifications featuring ccc_absorb() are identified; not all cases of
    identification failure in this case are bound to be automatically
    detected"). No screens or other tests to be implemented: the name
    coarsecohortcontrols itself flags that the option is for coarse cohorts
    (`i` containing multiple observations).
- **SI-7 — `o.`-omission-marker naming fragility.** `ext_contract10_rcs_factor_fe`
  is green-by-pin (0 -> 198, "5o. invalid name") — the pin does NOT certify the
  behavior as intended; the colleagues' Stata 17 old build hit the same error, our
  Stata 19 old build escapes by omission-pattern luck. Related: the dupcol entry's
  `o.`-marker truerow-strip trigger (report §1.11) is its own future corner. §2.2.
  - **DIAGNOSED 2026-08-25 (author asked "concretely, what fails?"; report
    `review_process/si_reports/si7_omission_marker_report.md`).**
    (a) Exact failure: the subgroup event-study block rebuilds relative periods
    by string-stripping the colnames of `margins, over(_t) subpop(&sg)` r(table);
    when a subgroup×horizon margin is NOT estimable its column comes back
    omission-marked (`5bno.__000003`), the sanitizer strips `bn` but not `o`,
    and `local relyrinj = 5o. -6 -1` (worker L2546) throws the bare r198. Not a
    dropped regression coefficient — e(b) is full, zero `o.` flags. (b) The old
    build's rc 0 is NOT good behavior: its returned tables are degenerate
    (all-zero clustered VCV, every SE/t/p/CI missing, pre-trend "estimates" =
    machine zeros ~9e-16, sg1's pretreatment row silently absent) — survival by
    accident because its margins name happens to parse as a number. (c) Minimal
    failing set: `att` + `subgroups()` + any ES horizon whose subgroup-level
    margin is inestimable — an ordinary spec (e.g. `subgroups(g)` with
    `fe(i.g##i.t)`, which absorbs subgroup pre-trends by construction, or a
    subgroup simply lacking a horizon). safety-on cannot catch it today: the
    deliberate omission screen reads e(b), not margins output. (d) Regression
    vs luck: the fragile chain is byte-identical in both builds; the flip is
    numeric — B21's float→double w_att (3cb0c41) exposed a knife-edge that a
    one-line toggle reproduces in BOTH directions (new+float → rc 0, old+double
    → rc 198); untouched old also crashes on the unweighted variant and on
    Stata 17. B21 exposed the crash, it did not create it. (e) Boundary:
    wooldid's DELIBERATE omission-pattern exploitation lives in reads of the
    regression's e(b) (o.-excluding regexes, the safety screen); the ACCIDENT
    is that margins-over() colnames can also carry omission markers, and the
    subgroup ES path is the one place such a name flows into an arithmetic
    evaluator instead of a filter. PROPOSED (awaiting ruling, none
    implemented): strip/guard `o` markers in the sg-horizon parser and report
    omitted horizons with a caution (mirroring the pooled path); or derive sg
    relyrs from e(b) names like the pooled path; optionally warn on fe()
    interacting the subgroups variable with time; re-label the suite pin
    (post-fix expectation = rc 0 with the sg1 pre-period reported as omitted,
    NOT restoration of the old degenerate table).
- **SI-8 — unit-level `i` at realistic scale dies late with cryptic r(907).**
  Cell-dummy grid exceeds maxvar (~2,000 counties x periods); correct usage is
  cohort-level `i`, but the failure is late and unexplained. Ticket: a pre-flight
  size check with an informative message steering to coarser `i`. §2.3 (ticket T-A).
  - **IMPLEMENTED 2026-08-25 per author ruling (083ea89).** Under default safety:
    >100 distinct `i` values → caution steering to fe() + coarser `i` (cohorts /
    onset groups); estimated variable demand > c(maxvar) → early informative
    r(907) with the same steer. Both skipped under safety(off). Cert:
    `test/test_si89.do` T3–T5 (SI89_ALL_PASS); sthlp safety() paragraph synced.
- **SI-9 — `poisexpresults` without `poisson` silently accepted on both builds.**
  The last member of the REPORT_2026-08-18 ticket-4 validation-gap family (the
  other two closed by b42b71f). Candidate r198 screen. §2.7.
  - **IMPLEMENTED 2026-08-25 per author ruling (083ea89).** Option-contract
    r(198): "poisexpresults requires poisson to have any effect". Suite pin
    `ext_contract10_poisexp_without_poisson` rc_new 0→198 (8f7e830); positive
    poisexpresults+poisson path re-checked in cert T2. sthlp synced.
- **SI-10 — T-M2: continuous-SG joint pre-trend p stored as literal missing.**
  The sg-continuous overjoint `r(p)` is 1x1 with the real p at [1,1] (instrumented:
  sg1 6.209e-06, sg2 4.109e-07); the `[1,2]` read silently discards it and
  `e(joint_Choriz_s*)` stores "." where significant p-values exist. AUTHOR-RULED
  DEFERRAL WITH MANDATED PROTOCOL (mini-assessment §2b): before fixing, dump the
  full wooldid console output plus the mid-run r(p)/r(b) matrices at the two sites
  from >= 3 different datasets incl. rich non-degenerate data; then [1,2] -> [1,1]
  at the two sites (~L2685/L2696 vintage), re-bless `cont_sg.csv` (sole baseline
  pinning the "."), and document the old-vs-new suite diff explicitly (old shares
  the bug; the comparison alone can never flag it). §2.6.
  - **DEMO DELIVERED 2026-08-25 (author-requested; d0f2625):**
    `review_process/demo_tm2_sgjoint.do` — planted sg-differentiated dose
    pre-trends; att sg joint p's real (sg3 = 5.7e-05), every ame twin "."
    (Exhibit A); the 1x1-continuous vs 1x2-binary r(p) geometry reproduced on
    plain `reg` + `margins` (Exhibit B).
  - **FIXED 2026-08-25 (author-ordered).** Protocol honored in full: (1) the
    mandated evidence was banked pre-fix at `review_process/tm2_evidence/`
    (3+ datasets, instrumented + stock, pre-R3 vintage — R3 makes the sites
    unreachable on the cluster(i) pinning config); (2) [1,2] → [1,1] at the two
    sg-continuous sites (`joint_Choriz_s#` / `joint_Choriz_i_s#`), matching the
    pooled-continuous reads that were already [1,1]; the binary twins'
    1×2-geometry [1,2] reads untouched; (3) the old-vs-new diff is now
    explicitly pinned per the author's instruction: new suite config
    `trig_tm2_sgjoint` (deterministic planted sg pre-trends) + ledger rows
    TM2-sg-p / TM2-ext-p, with the EXPECTED direction recorded — old posts
    p=., new posts a real p; a new-side regression to p=. leaves no diff, the
    scoped row goes STALE, and the run FAILS (the mechanical deviation flag).
    Probe evidence on the exact pinned config: old ame twins all "."; new ame
    p's real and recovering the planted structure (sg1 .81 no-trend, sg2
    1.4e-13, sg3 1.4e-36) with att twins bit-identical across builds (%21x).
    Re-bless is MOOT: `test/baseline/cont_sg.csv` contains zero joint_Choriz
    entries post-Group-A (the containment gate skips the block under
    cluster(i)), and the battery's three contreat+subgroups configs all use
    cluster(i) — battery stays bit-identical. Suite NOT rerun per author
    instruction; next full run expects the two TM2 rows to fire (on
    trig_tm2_sgjoint and the ext_c2_* grid, whose 24 configs captured the
    "." on both builds pre-fix).
- **SI-11 — `allsgssame`-ame joint tests missing-valued on both builds.** The
  continuous SG-equality tests ("all SGs have the same ame"/"...pre-treatment ame")
  are missing in both probe cases AND pinned missing in the blessed cont_sg
  baseline, where the clean case expects rejection. Separate site family from SI-10
  (E1 chain machinery), undiagnosed — needs its own diagnosis first. §2.6 adjacent.
  - **DIAGNOSED 2026-08-25 (author-directed extensive diagnosis; report
    `review_process/si_reports/si11_allsgssame_diagnosis.md`). NOT a code defect —
    it is the family-A (G1) cluster-in-cell covariance degeneracy, already
    handled by Group A.** Direct answer to the author's clustering question:
    **turn clustering off and the missingness disappears** — under robust, no
    vce, and cluster(unit) the SG-equality ame tests return real, sane values
    on BOTH builds (a truly-false null rejected at 6e-41; a true null p=.667).
    Mechanism: under cohort-level cluster(i)+contreat the posted margins
    covariance of the per-SG AMEs is exactly rank 1 (measured: every e(V)
    entry identical), so the equality contrast has exactly zero estimated
    variance; `test` drops the zero-variance constraints (r(drop)=1, r(df)=0)
    and the overjoint fallback returns r(p)[1,1]="." — or, when rank detection
    only half-fires, silent garbage (5.2e-80 / 1.03e-109, bit-matched to the
    banked pins). The chain machinery itself (wooldid_sgchain, rank guard, the
    continuous [1,1] read) verified correct — no T-M2-class indexing issue
    here (r(p) confirmed 1x1). The att twin "works" only because att cells
    lack the annihilated slope channel (near-rank-1; ties to open SI-A1). At
    HEAD, the amescrub/R3 gate correctly skips the chain under containment
    with the explanatory caution; the ungated `unconditionalse` carve-out
    probed clean on both builds. PROPOSED (awaiting ruling, none implemented):
    F0 close as diagnosed + one sthlp sentence; F1 optional hardening — on
    r(drop)==1 skip the overjoint fallback (it posts residual-rank garbage)
    with a caution; F2 optional tidy — suppress "."-valued locals at the
    ereturn site and cap the uncapped r(p)[1,1] read to match its pre twin.
    No valid plug-in number exists under containment; bootstrap/IM dispersion
    (the approved future R2 work) remains the only principled restoration.
  - **CLOSED 2026-08-25 by author ruling: solved by the Group A work.** F0/F1/F2
    all declined — no sthlp change, no hardening. No further action.

## C. Documentation-grade actions (no estimator change)

- **SI-12 — sthlp caveat bundle.** (i) Singleton-cell designs: with unit-level `i`,
  treated residuals vanish and V_es carries control-cluster noise only (SE ratios
  0.84-1.07 vs the pooled event study; b/att unaffected at 1e-15) — caveat on the
  honestdid/pretrendpower handoff entries. (ii) The family-A design-conditions
  language (= remedy (a); lands with whatever family-A choice is made). (iii) A
  note near the honestdid example that Rambachan-Roth C-LF robust sets can be
  legitimately EMPTY on data with strong pre-trends (observed; not a handoff
  defect). §5.5, oracle O16.
  - **Author ruling 2026-08-25:** (i) demo + commentary delivered —
    `review_process/demo_singleton_ves.do` (d0f2625): same common-onset data as
    unit-level `i` vs pooled cohort; att reldif 4.2e-16, ES SE ratios 0.84–1.07.
    When it arises: any design where `i` is a unit-level identifier, i.e.
    wooldid's most ordinary binary usage — each treated cell is fit exactly by
    its own dummy, so the clustered V_es is built from control-side residual
    variation only; it matters mainly for the honestdid/pretrendpower handoffs,
    which receive a structurally different covariance than the pooled-cohort
    formulation of the same event study. (ii) family-A design-conditions sthlp
    language: DONE 2026-08-25 per author ruling — the continuous-treatments
    caution section now opens with the author's framing (valid clustered SEs
    cannot be obtained analytically for the ame-family quantities when
    clustering at the cohort level or above, necessitating bootstrapping or
    other methods). (iii) skipped
    per author.
- **SI-13 — upstream courtesy report.** honestdid v1.1.2 crashes (mata assert 3498
  in `_flciFindWorstCaseBiasGivenH`, parallel AND sequential paths) on
  binary_mixed_imbalanced delta(sd) while delta(rm) runs clean on the same b/V.
  Reproducible; report to the package maintainer. External action only. §5.5.
  - **DONE 2026-08-25:** author is handling the upstream report directly; no
    further action in this repo.

## C2. New defect found in passing (bootstrap R2 review, 2026-08-30)

- **SI-18 — CLOSED 2026-08-31, no action needed: a DEV-HARNESS ARTIFACT, not a user
  bug.** No .mlib is shipped and none is needed; the mata stays inline in wooldid.ado.
  Detail below.** Measured: with wooldid loaded via the adopath, as any installed user has it,
  every case returns rc 0 - gteffects+cleanmatrices in one run, bootstrap and permute
  after a cleanmatrices run. The r(3499) failures reproduce ONLY when the package is
  loaded with `do "<path>/wooldid.ado"`, which is how every test suite, harness and agent
  brief in this project loads it, and how all the original evidence was gathered.
  Mechanism, isolated in a ten-line reproduction: a missing mata function makes Stata
  re-read the defining ado from the adopath, re-running its mata block; a `do`-loaded
  program has no such association. CONSEQUENCE FOR SI-17: its refusal blocks a
  combination that works for real users (pre-SI-17 build, adopath-loaded,
  permute+cleanmatrices = rc 0) - see the SI-17 entry. Original report retained below.** Found by the commit-E adversarial
  verifier (finding F2, measured). SI-17 refuses `permute` + `cleanmatrices` on the same
  command line, but `cleanmatrices` runs `clear mata`, which destroys the ado's own
  compiled mata functions for the rest of the session — and because the ado is already
  loaded and cached, its mata block never re-runs. Measured in one clean session:
  `bootstrap(bayes)` rc 0, then a plain `cleanmatrices` run rc 0, then the SAME bootstrap
  command → r(3499) `wooldid_bootdw() not found`; a subsequent `permute()` → r(3499) too.
  Both failures land after a fully successful estimation.

  **SCOPE CORRECTION (measured 2026-08-31, after the first filing).** This is not only a
  cross-run problem, and the existing refusals do not contain it. `cleanmatrices` clears
  mata early in the worker (~L287), while `gteffects` calls `wooldid_gtx()` much later in
  the same worker (~L1564) — so **`gteffects` + `cleanmatrices` fails r(3499)
  `wooldid_gtx() not found` on one command line**, with no refusal anywhere. Fresh-session
  controls: `gteffects` alone rc 0, `cleanmatrices` alone rc 0. Refusing combinations one
  at a time therefore does not scale — it would need a third refusal now and a new one for
  every future mata consumer, whereas a structural fix covers all of them at once.

  **MECHANISM EVIDENCE for the candidate fixes (all measured):**
  - Functions living in a **mata library survive `clear mata` and reload on demand**
    (probe: built an mlib, cleared mata, called the function — rc 0, correct value).
  - **Re-running the ado from inside one of its own programs FAILS r(110)**, with or
    without the `cap prog drop` pattern this file uses — a program cannot be dropped or
    redefined while it is executing. So "re-`run` the ado to recompile" does not work.
  - **Re-sourcing a mata-only file from inside a running program DOES work** (rc 0), and
    the failure is detectable first (`cap mata:` a probe call returns 3499).
  - Mata object names **are** enumerable by logging `mata: mata describe` and parsing,
    so a selective drop is mechanically possible though it means parsing a log listing. The do-no-harm verifier
  independently hit the same class on the OLD build (a permute config with no
  `cleanmatrices` at all died r(3499) because an earlier config in the session had used
  it), so this predates the bootstrap work. User-side recovery today is `discard` or
  re-running the ado.
  CANDIDATE FIXES, reordered by the evidence above:
  1. **Ship the mata in a library** (`lwooldid.mlib` built by `mata mlib create/add/index`,
     added to the `.pkg`). `clear mata` then cannot reach wooldid's functions at all —
     they reload on demand — so `cleanmatrices` keeps its full documented brute-force
     meaning and every consumer (gteffects, permute, bootstrap, anything later) is immune
     by construction. Structurally the cleanest. Costs: a build step in the release
     process, a `.pkg` entry, and a decision about the dev workflow, since `do
     wooldid.ado` currently defines the functions inline and inline definitions would
     shadow (or be shadowed by) the library copies.
  2. **Split the mata into its own file and re-source it after clearing.** `cleanmatrices`
     does `clear mata` then immediately re-sources `wooldid_mata.do` (via `findfile`), so
     third-party junk goes and ours comes back — exactly the option's documented intent
     ("clear out ... content in mata" left by reghdfe). Proven to work from inside a
     running program. Smallest behavioral change; costs a second shipped file.
  3. **Selective drop**: enumerate mata objects and drop only those outside the
     `wooldid_*` namespace. Same end state as 2 without a second file, but it means
     capturing and parsing `mata describe` output — the most fragile of the three.
  4. **Document only**: state in the help that `cleanmatrices` ends the session's ability
     to use gteffects/permute/bootstrap until wooldid is reloaded, and leave the behavior.
     Cheapest, but leaves a silent-looking r(3499) after a successful estimation.
  Not implemented pending your ruling: 1 is a packaging change and 2 changes the shipped
  file set, so neither is mine to choose unilaterally.

  **RESEARCHED 2026-08-31 — see `review_process/mlib_decision_2026-08-31.md`.** The
  library fix is PROVEN to cure all four failing cases (each r3499 -> rc 0, with
  wooldid.ado unchanged), but three complications were measured that make the prebuilt
  form a decision rather than a formality: (i) an mlib records the Stata version that
  compiled each function, and an older Stata IGNORES the member with a silent r(3499) —
  no `version()` option exists and `version 16` does not help — so a library built here
  reaches Stata 19 users only, while the package declares 16.0; (ii) after `clear mata`
  resolution falls to whatever library the adopath holds, rc 0 and no message, so a stale
  or older-release library silently shadows the running source — the exact dev-vs-shipped
  divergence the chunk set out to prevent; (iii) reghdfe and ftools, the closest
  comparators, deliberately ship no prebuilt library and recompile when c(stata_version)
  changes. A fourth design emerged and is now recommended: have `cleanmatrices` build the
  library from the in-memory functions immediately BEFORE clearing, which dissolves all
  three at once (built by the running Stata, from the running source, no binary shipped).
  Build script at `build/build_mlib.py` and cert suite at `test/test_mlib.do` are ready
  and work; neither is wired into anything pending the choice.

- **SI-17 — REVERTED 2026-08-31 by author ruling; CLOSED.** The refusal is gone and
  `permute()`/`bootstrap()` with `cleanmatrices` return rc 0 again (measured, adopath-
  loaded). Rationale below.

- **SI-17 (history) — the refusal that was shipped and then reverted.** The refusal was ruled and shipped on evidence gathered under `do`-loading.
  Measured 2026-08-31: on the pre-SI-17 build loaded via the adopath - i.e. as an
  installed user has it - `permute(shuffle) permreps(20) cleanmatrices` returns **rc 0**.
  So the shipped refusal removes a working user configuration, and the same applies to the
  `bootstrap() + cleanmatrices` refusal added later. RECOMMEND revisiting both: either
  revert them, or narrow them to the `do`-loaded case if that can be detected. Not changed
  pending the author's ruling, since they were author-ruled in. Original entry follows.

- **SI-17 (as ruled and shipped) — author-ruled 2026-08-31; fix (a), refuse the combination.**
  `permute()` now refuses `cleanmatrices` at parse with a named message, mirroring the
  bootstrap refusal; `test_permute2.do` carries S17a-S17d (rc, message text, and both
  options still working alone), and the pre-fix build was confirmed to reach the
  r(3499) failure on the same command line. Original report follows.

- **SI-17 (original) — permute() + cleanmatrices kills the permutation engine after a
  successful estimation (r3499).** MEASURED at HEAD 4b08962 by the bootstrap plan's blind R2
  reviewer (probe P2, artifacts in the session scratchpad r2/): `permute(shuffle)
  permreps(20) cleanmatrices` dies r(3499) "wooldid_permx() not found" - the worker-end
  `clear mata` that cleanmatrices requests destroys the ado's own compiled mata
  functions, so the shell's post-worker permutation engine cannot run (control without
  cleanmatrices: p=.238, rc 0). Pre-existing since F6; not introduced by the bootstrap
  work. The bootstrap feature refuses bootstrap+cleanmatrices v1 for the same reason
  (plan amendment A16). CANDIDATE FIXES for author graduation: (a) refuse
  permute+cleanmatrices (cheap, mirrors bootstrap); (b) scope the cleanmatrices
  `clear mata` away from the ado's own namespace (e.g. mata drop of user objects only);
  (c) defer the mata clear until after the shell epilogue. Not implemented pending your
  ruling.

## E. Feature tickets (author-ruled, not yet implemented)

- **SI-19 — per-horizon permutation p-values (author-ruled 2026-08-31: "we absolutely
  should compute per horizon p values for the permutation tests", as a separate ticket).**
  Today `permute()` reports statistics for the POOLED estimands only: `e(wdid_permest)` is
  e.g. "att att_it" and `e(permdist)` carries exactly one column per estimand, with no
  `ES_` entries anywhere (verified 2026-08-31). Yet every replication already computes the
  full event-study horizon vector and discards it - measured at 16% of permute wall time,
  and 43% under `permtstat` where each horizon's delta-method SE is computed and thrown
  away too. So this feature costs almost nothing to add: the work is already being done.
  WHY IT MATTERS: randomization inference exists precisely for the few-cluster setting
  where analytic pre-trend tests are unreliable, so per-horizon RI p-values are wanted
  exactly where they are hardest to get any other way.
  IMPLEMENTATION NOTE: each replication's `e(b)` already carries the complete named
  `ES_att_pre#/post#` vector, and `wooldid_bootrun` already performs exactly this
  name-aligned accumulation (stacking rep vectors by column name with per-column
  availability counts) - so the pattern to copy is in the file. Open questions for the
  implementation: the e() naming for per-horizon results, the display surface, and whether
  the 5th-95th null band and the add-one p-value carry the same interpretation per horizon
  as they do pooled.
  CONSEQUENCE FOR THE EFFICIENCY REVIEW: the "permute computes a discarded event study"
  finding (16-43%) is hereby RETIRED as an efficiency item - the computation is not waste,
  it is the input to this feature, and must NOT be gated off.

## D. Suite-side residuals (test infrastructure, not wooldid defects)

- **SI-14 — NFADV unpost-path coverage.** The pooled-reference/degenerate-q
  adversarial unpost paths are ledgered only for `ext_*` configs; no core/adm
  config exercises them. Add one. §4 gap 1.
  - **ADDED 2026-08-25 per author ruling (8f7e830):** `cov_advpool` (adversarial
    without esfixedbaseperiod → the pooled-reference unpost path) + 3 scoped
    ledger rows (NFADV-core-p/-l/-ln) mirroring the ext_* family.
- **SI-15 — `cov_icbad` under-triggers B21(b).** Engineered-collinearity levels
  ~1e6 are float-representable; a ~1e8 variant would pin the class old-side. §4
  gap 2.
  - **ADDED 2026-08-25 per author ruling (8f7e830):** `cov_icbad2` (icbad8 prep,
    levels ×100 → ~1.02e8 > 2^24) + 5 wildcard rows and the README
    accepted-blind-spot listing (G10B precedent). Both new configs validate on
    the full suite run now in flight; stale-row amendments, if any, follow the
    AMESE-R4-b precedent.

---

Sequencing note for the chunk: family A is one design decision executed five ways —
settle the remedy choice first (it determines whether SI-1..SI-5 are one patch + one
re-bless or documentation-only); SI-10 has its protocol pre-scripted and touches two
lines + one baseline; SI-11 needs diagnosis before any fix; B-items are independent
of each other and of A. Every wooldid.ado change re-runs the affected suite subsets
per the established protocol before the campaign state is re-claimed.
