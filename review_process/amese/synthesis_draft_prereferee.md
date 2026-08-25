# DRAFT synthesis — AME-SE diagnosis classification (for adversarial referee attack)

Evidence base: structure pass (28+7 cases), MC grid (26 x 256, preregistered bands),
locked theory A (incidental-parameters) + B (design-based), blind theorist
(probe-anchored), code-path auditor (pending — d1c/byCohort channels + fix sites).
User's four buckets applied per case group.

## Classification table

| group (cases) | MC evidence | bucket | refinement |
|---|---|---|---|
| G1. cohort-level cluster(i) + basic/within-unit-dose (d1_basic_ci, d2_basic_ci, d5_within_ci, d4_within_ci posting reps, d1c_basic_ci posting reps) | rank-1 V; every SE = SE_cl(pooled phi-hat); ratios: ES horizons ~0.2, pre ~0.4, main ~0.5 (ratio^2 tracks #cells averaged); posted SILENTLY | **(iii)** | The estimand is fine; a correct SE is estimable (cross-cohort dispersion of per-cohort AMEs; cluster bootstrap of the margins functional) — but NOT by plug-in CRVE on this regression (scores identically zero; within the requested method a referee could argue (iv)). The posted number is a real SE of the WRONG object (the absorbed baseline slope), attributed to every ame quantity. |
| G2. cluster(i) + byTime (d1_bytime_ci, d2_bytime_ci) | full-rank V but each horizon SE = SE_cl(phi_t-hat) to 7 digits (blind probe); ratios 0.38-0.69 | **(iii)** | Same annihilation with a period-varying baseline block replacing the rank signature — DEFEATS any rank/eigenvalue guard. Fix predicates must be structural (design-based), not numeric. |
| G3. cluster(i) + byCohort (d1_bycoh_ci, d2_bycoh_ci) | total blackout: "not identified 20 of 20", no e(V), no e(V_es_ame) | **(iii) — and the blackout itself contains a repairable collateral defect** | AUDITED: reghdfe's stored e(V) arrives bitwise ZERO because Stata's built-in `ereturn post` screen zeroes the WHOLE matrix when it is nonsymmetric/highly singular beyond tolerance — destroying the VALID att variances in transit (plain regress on identical data posts att SEs .15-.21; the ame slope block is the genuinely-annihilated part). So: ame side = G1's (iii); att side = pure code-repairable collateral (detect e(V)==0 post-reghdfe, fall back to wooldid's own regress path L1102-1108 or repair/split the V). The honest caution is right; losing att inference is not. |
| G4. all-treated + cluster(i) (d4_within_ci) | 95/256 reps blackout (audited: 3 of 4 frozen realizations trip the same ereturn-post zeroing; margins additionally blanks att rows on a near-zero horizon), posting reps rank-1 with all SEs = SE_cl(phi) exactly; main/pre ratio ~1.0 is compositional coincidence, ES ~0.7 | **(iii)** | Stochastic, numerically-triggered all-or-nothing behavior on an over-parameterized design; the structural predicate (F1) makes the handling deterministic. |
| G5. cluster(unitid)/robust on d1/d2, all contreatcontrols (12 cases incl. d2_bycoh_cu at band edge; d1c_basic_rb) | ratios 0.68-0.91 — consistent ~15-25% understatement; full rank, horizon-varying | **(i) with a material caveat** (displayed SEs) / **(iii) for the joint tests** | Textbook CR1/HC1 small-cell leverage bias (6-obs cells, 2 within-cell params; slope-weighted leverage ~0.45). Displayed-SE harm = 8-12pt under-coverage (caveat + opt-in CR2/HC3); but the 6-8-dim joint tests built on these V's over-reject 20-35% — that crosses into required-change territory. |
| G6. all-treated + cluster(unit)/robust (d4/d5 cu, rb) | ratios 0.86-1.11 — inside band | **(i)** | Pointwise only; does not certify joint tests. The improvement is compositional (truth weight moves onto the well-estimated baseline block), not a repair of the leverage bias. |
| G7. dose constant within cohorts (d3, all VCEs — structure-only) | deterministic: 165 collinearity drops, ame rows 0 with missing SEs, generic cautions | **(iv), currently semi-handled** | The within-cell AME estimand is genuinely unidentified. Wooldid does not fabricate an SE (good) but fabricates 0-point-estimate rows behind generic cautions; correct output is a NAMED refusal ("dose has no within-cell variation; ame not identified; consider cohort-level designs/steer"). |
| G8. within-option on data WITH controls (d1/d2_within_*, structure-only) | bitwise = basic | (not an SE bucket) | contreatwithin is estimation-inert when controls exist — doc gap or missing screen; user asking for within-variation identification is not getting it. |
| G9. d1c cluster(i) (cell-correlated shocks) | AUDIT CORRECTION to the structure pass: NOT a systematic blackout regime — 15/16 realizations are ordinary G1 rank-1; the blackout is a 1-in-16 y-driven numerical knife-edge through the same ereturn-post zeroing. MC: WRONG 0.18-0.46 with ~13% blacked-out quantities; robust FLAG 0.75-1.05 because the treated-cell dummies absorb the cell shocks | **(iii)** (= G1 with stochastic G3-style wipeouts) | Normative (theory B, upheld): cluster(i) is the CONCEPTUALLY right level; the remedy must be a valid cluster-level estimator (bootstrap/dispersion), NEVER a steer to robust (robust survives here only because the shocks are cell-level; unit-serial dependence would break it). |

Cross-cutting: bucket (ii) is EMPTY — no posted number is a correct SE of the posted
estimator anywhere on the grid. Old wooldid shares every pathology (same spec/VCE),
plus posted the invalid g-inverse adversarial-ame p that new refuses.

## Fix recommendations (for the standing-issues file; NO code changes now)

F1 (P1). STRUCTURAL DETECTION + REFUSAL: predicate = every treated (i,relyr) cell
lies inside a single cluster of the specified cluster() (cheap; catches basic,
byCohort, byTime, within, all-treated — INCLUDING the spoofed-full-rank G2). On
detection: suppress all plug-in ame SEs and e(V_es_ame), print a design-specific
caution naming the mechanism and the remedies. Never guard on rank/eigenvalues.
F2 (P1b). NAMED REPLACEMENT ESTIMATOR (opt-in or default under the predicate):
cluster bootstrap of the margins functional, or the analytic cross-cohort dispersion
of per-cohort AMEs (IM-style). Blast radius: new option + docs; battery/oldnew
untouched until made default.
F3 (P2). JOINT TESTS on V_es (adversarial, jointtests, roth prep, NIS size): under
the predicate — refuse; elsewhere (cu/rb small-cell) — bootstrap/CR3 the V or refuse
above a dimension threshold; 20-35% at 5% is not caveat material. Touches advstat/
jointtests/ptp paths; battery adv_* pins re-bless.
F4 (P3). Make the G4 stochastic blackout deterministic via F1's predicate.
F5 (P4). Small-cell leverage bias (G5): sthlp caveat + a printed note when cells
are small; opt-in CR2/HC3-style correction consideration; default unchanged.
F6. G3/G4 blackouts: detect e(V)==0 after reghdfe (the ereturn-post zeroing
signature), repair via the regress fallback or a split post that preserves the valid
att block; keep the caution, make it design-specific. G7: named unidentified-ame
refusal instead of structural-zero rows sitting unflagged in e(b) (the L1244-1332
catastrophe counter already detects the condition). G8: screen or doc for
estimation-inert contreatwithin-with-controls.
F6b. byCohort+ame estimand normalization-dependence: separate sub-ticket — decide
what the byCohort ame estimand IS before offering SEs for it.
F7. Remove/replace the stale-adopath hazard note + sthlp design-conditions rewrite
(remedy (a) from §5.2), folding the F2 sthlp caveat for e(V_es_ame) consumers.

## Auditor findings folded in (all audited, review_process/amese/code_audit.md)

- The blackout channel is Stata's `ereturn post` V-zeroing screen (bracketed
  empirically: asymmetry in (1e-12,1e-6] or off-diagonal-carried indefiniteness
  zeroes e(V) wholesale; PSD rank-deficiency passes) — partially-valid VCEs die
  whole. wooldid's `qui` swallows the only warning. One-line e(V)==0 check after
  the reghdfe call (L1099) + the existing regress fallback path = the repair site.
- Every wooldid SE detector tests only MISSINGNESS (L3307-3309, L3352-3358): a
  present-but-garbage V (G1/G2) passes silently; an absent V (G3/G4) cautions.
  That is the silent-vs-caution asymmetry's exact location.
- F2 posts b_es and V_es inside ONE gate: blackout consumers get neither (safe);
  G1/G2 consumers get b WITH garbage V — the real handoff gap.
- byCohort ame POINT estimates are normalization-dependent across solvers
  (reghdfe .723 vs regress .554 on identical data; att invariant) — a distinct
  estimand-definition question for byCohort+ame, flagged as its own sub-ticket.
- Correct-SE ingredients for the replacement estimator are already in hand at the
  margins site (cell slopes in e(b), weights in r(Jacobian), Var_cl(phihat)).
- Harness limitations (documented, not validity-threatening): expand/bysort tie
  order after a wooldid run makes individual MC reps non-reproducible standalone
  for ame quantities (the DGP distribution is unchanged — unit labels exchangeable);
  the structure harness's `matrix V = e(...)` coercion made "absent" print as
  "1x1 rank 0" (interpretation already treated them as absent).
