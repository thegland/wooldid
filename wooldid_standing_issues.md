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

## A. The V_es inference pathology family (one root cause, one remedy menu)

Root cause (mechanism-proven, §5.2-§5.4 + INV-1): saturated cohort-by-relative-period
cells inside single clusters zero the CRVE scores, understating the event-study
covariance DIRECTIONALLY — cross-horizon contrast directions carry 7-20x less
reported than true variance while diagonals look only 1-2x off. Every
V_es-consuming multi-dimensional test inherits over-rejection. The test statistic
code itself is sound (robust-VCE att-side sizes at 5.75%).

Shared remedy menu (§5.2, rising ambition): (a) sthlp design-conditions rewrite
(when are these p-values trustworthy: clustering finer than the cells, or
robust/unclustered; narrow windows do NOT rescue multi-unit designs); (b) detect
cell-saturated clustering (cheap predicate: every treated (i, relyr) cell inside a
single cluster of the specified cluster()) and print a specific caution on every
V_es-consuming p-value and on the honestdid handoff; (c) cluster-bootstrap or
CR3/jackknife V for the ES block; (d) eigenvalue-floor regularization inside
wooldid_advstat (defensible but ad hoc). Any fix here moves real results: battery
adv_*/esf_adv pins and the oldnew adversarial surfaces re-bless with documentation.

- **SI-1 (NF-1a) — adversarial line test oversized under cohort-level cluster().**
  Flat-null rejection at nominal 5%: 20.5-23.2% (G=40 singleton, incl. AR(1) and
  planted-line nulls), 48.4% (G=40, 8-unit cohorts), 44.8% (G=100 multi-unit),
  24.8% (2/2 window multi-unit), 11.7% (G=200 singleton); altmethod 11.2-12.6%;
  conservative (2.4%) at G=10. Evidence: §5.2 rate table, M=1000/config MC.
- **SI-2 (NF-1b) — `jointtests` pre-trend p (`joint_horiz`) oversized.** Same V,
  old feature: 11.3-36.0% at nominal 5% across the G=40 arms; **65.7% at G=10**.
  The sthlp's existing few-clusters caution radically understates this. §5.2.
- **SI-3 (NF-1c) — `pretrendpower` reported p/size scalars inherit the pathology.**
  `e(roth_prep)` rejects 14.0% at nominal 5% under the null on the saturated
  design; reported `e(roth_nissize)` .1325 vs empirical .2132. The power claims
  (mdt columns) themselves validated and need no change. §5.2.
- **SI-4 (NF-2, demo issue) — all ame SEs identical + rank-1 `e(V_es_ame)` under
  cohort-level cluster().** Every displayed ame SE (main, pretreatment, every ES
  horizon) equals the SE of the absorbed baseline dose slope; e(V_es_ame) is
  EXACTLY rank 1, so every cross-horizon contrast carries zero reported variance
  while true variances are ~5x larger and full rank; honestdid/lincom/test ingest
  the matrix silently. Old wooldid shares the V and additionally POSTED an invalid
  g-inverse adversarial-ame p (.6147, verified in probe8_old.log) where new
  correctly refuses — so the refusal half is already fixed; the standing defect is
  computing/posting that V (and those SEs) without warning. The F2 sthlp entry
  needs the caveat regardless of remedy choice. Runnable reproduction:
  `review_process/demo_vesame_pathology.do` (identical SEs 0.059040286 x8,
  correlation 1.0, eigenvalues .0279 vs 6e-17, contrast variance -2.3e-17 against
  an estimate difference of 0.491; robust rerun: SEs 0.166-0.211, full rank).
  Evidence: §5.3 + INV-1 report.
- **SI-5 (NF-3) — ame-side tests oversized under EVERY tested VCE, incl. robust.**
  On iid multi-unit data with `robust`: adversarial-ame 31.5%, `joint_Choriz`
  22.0% at nominal 5%, while the att side is simultaneously correct (5.75%) —
  small-cell dose-slope variance bias amplified by the 6-8 dimensional quadratic
  forms. Continuous-treatment ES inference currently has NO configuration with
  trustworthy joint p-values. §5.4.

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
