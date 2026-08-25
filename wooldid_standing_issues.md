# wooldid standing issues — the bug-fixing work list

Derived 2026-08-24 from `wooldid_old_bug_report.md` (the 2026-08-23 report of record
plus its §5 new-feature-validation extension). This file lists ONLY what remains to
be addressed; nothing here is fixed yet. Excluded by construction: everything in §1
(fixed old bugs), §3 (intended changes), §5.1 (validated-clean surfaces), the pure
records of §2 (items 4, 5, 8: environment pins, the corrected-T1 reading, parity
pins), and the author-ruled ignores from the mini assessment (T-M1; T-M3 waits for
the efficiency pass). Evidence pointers cite the report section, the INV-1 mechanism
report (`review_process/inv1_v_es_ame_2026-08-24.md`), the mini-assessment protocol
(`review_process/mini_assessment_2026-08-23.md` §2b), and the runnable demo
(`review_process/demo_vesame_pathology.do`).

---

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
- **SI-7 — `o.`-omission-marker naming fragility.** `ext_contract10_rcs_factor_fe`
  is green-by-pin (0 -> 198, "5o. invalid name") — the pin does NOT certify the
  behavior as intended; the colleagues' Stata 17 old build hit the same error, our
  Stata 19 old build escapes by omission-pattern luck. Related: the dupcol entry's
  `o.`-marker truerow-strip trigger (report §1.11) is its own future corner. §2.2.
- **SI-8 — unit-level `i` at realistic scale dies late with cryptic r(907).**
  Cell-dummy grid exceeds maxvar (~2,000 counties x periods); correct usage is
  cohort-level `i`, but the failure is late and unexplained. Ticket: a pre-flight
  size check with an informative message steering to coarser `i`. §2.3 (ticket T-A).
- **SI-9 — `poisexpresults` without `poisson` silently accepted on both builds.**
  The last member of the REPORT_2026-08-18 ticket-4 validation-gap family (the
  other two closed by b42b71f). Candidate r198 screen. §2.7.
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
- **SI-11 — `allsgssame`-ame joint tests missing-valued on both builds.** The
  continuous SG-equality tests ("all SGs have the same ame"/"...pre-treatment ame")
  are missing in both probe cases AND pinned missing in the blessed cont_sg
  baseline, where the clean case expects rejection. Separate site family from SI-10
  (E1 chain machinery), undiagnosed — needs its own diagnosis first. §2.6 adjacent.

## C. Documentation-grade actions (no estimator change)

- **SI-12 — sthlp caveat bundle.** (i) Singleton-cell designs: with unit-level `i`,
  treated residuals vanish and V_es carries control-cluster noise only (SE ratios
  0.84-1.07 vs the pooled event study; b/att unaffected at 1e-15) — caveat on the
  honestdid/pretrendpower handoff entries. (ii) The family-A design-conditions
  language (= remedy (a); lands with whatever family-A choice is made). (iii) A
  note near the honestdid example that Rambachan-Roth C-LF robust sets can be
  legitimately EMPTY on data with strong pre-trends (observed; not a handoff
  defect). §5.5, oracle O16.
- **SI-13 — upstream courtesy report.** honestdid v1.1.2 crashes (mata assert 3498
  in `_flciFindWorstCaseBiasGivenH`, parallel AND sequential paths) on
  binary_mixed_imbalanced delta(sd) while delta(rm) runs clean on the same b/V.
  Reproducible; report to the package maintainer. External action only. §5.5.

## D. Suite-side residuals (test infrastructure, not wooldid defects)

- **SI-14 — NFADV unpost-path coverage.** The pooled-reference/degenerate-q
  adversarial unpost paths are ledgered only for `ext_*` configs; no core/adm
  config exercises them. Add one. §4 gap 1.
- **SI-15 — `cov_icbad` under-triggers B21(b).** Engineered-collinearity levels
  ~1e6 are float-representable; a ~1e8 variant would pin the class old-side. §4
  gap 2.

---

Sequencing note for the chunk: family A is one design decision executed five ways —
settle the remedy choice first (it determines whether SI-1..SI-5 are one patch + one
re-bless or documentation-only); SI-10 has its protocol pre-scripted and touches two
lines + one baseline; SI-11 needs diagnosis before any fix; B-items are independent
of each other and of A. Every wooldid.ado change re-runs the affected suite subsets
per the established protocol before the campaign state is re-claimed.
