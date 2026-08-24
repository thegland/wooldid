# Preregistration — new-feature statistical validation (2026-08-24)

Written and committed BEFORE any Monte-Carlo run fires (f239262 precedent). Scope: the
work-file "New Features" chunk — external/statistical validation of the feature surface
with no old-build analogue: randomization inference (`permute()` family), the
adversarial line tests (author-directed inclusion), `pretrendpower()`, `fastnose`,
the F2 event-study b/V returns (honestdid handoff), and the reghdfe common-onset
oracle. The MC suite is a NEW component (`mc/`) deliberately excluded from
`run_full.sh all` (too costly for the standing campaign run).

Phase-0 measured costs (probe0-probe4, this machine): RI B=99 ≈ 7-10 s/rep;
permexact(20) ≈ 1.4 s; adversarial binary G=40 ≈ 1.3 s; continuous ≈ 6 s;
pretrendpower+jointtests ≈ 3.2 s. Pool cap: 8 concurrent Stata instances.

## Global protocol

- Every MC rep sets `set seed <base+rep>` immediately before data generation, and the
  seed is recorded in the rep's output row, so any rep reproduces standalone
  (permute consumes the session RNG).
- Drivers load wooldid.ado once per shard process; per-rep wooldid calls are `cap noi`
  with rc recorded. Output = one CSV row per rep.
- Analyzer: exact binomial acceptance regions (no normal approximation at the margins),
  chi-square GOF on the achievable p-grid as the secondary uniformity test
  (10 equal-probability bins; flag at chi2 p < 0.01).
- Any post-hoc change to a band, config, or rule below is a documented DEVIATION in the
  final report, with the reason.

## M1 — RI size calibration (`permute()`)

Design: per rep, draw the TRUE assignment from the config's mechanism, generate the
outcome under the sharp null (assignment never enters y), run wooldid with the SAME
mechanism as the permutation scheme, collect `e(perm_p_att)`, `e(perm_failed)`,
`e(perm_lb_att)`, `e(perm_ub_att)`, `e(perm_comparator)[1,1]`, rc, seed.

Exactness being tested (verified against the code, wooldid.ado ~L4769-4794):
tail count is inclusive (|s| >= |obs|); sampled p = (1+nex)/(nuse+1); exact
p = nex/nuse. With (obs, draws) i.i.d. from the randomization distribution
conditional on y, P(p <= k/(B_eff+1)) = k/(B_eff+1) exactly (sampled; ties only
conservative — assignment spaces chosen >= ~3000 so tie prob <= ~3%); exact mode:
P(p <= k/M_enum) = k/M_enum.

Configs (M = 400 MC reps each; B = permreps(99) unless noted; panels 10 periods,
outcome = cluster FE N(0,1) + noise as stated; NO treatment effect anywhere):

| id | design | mechanism | errors | notes |
|---|---|---|---|---|
| R1 | 12 clusters, onsets {5,7,never}×{3,3,6} | shuffle | iid N(0,1) | space 18480 |
| R2 | as R1 | shuffle | AR(1) rho=.7 within cluster | cluster-asymptotics-hostile; RI must stay exact |
| R3 | as R1 | shuffle + permtstat + permtci | sd multiplied 1..4 across clusters | t-scale path; CI coverage secondary |
| R4 | 16 clusters, 2 strata of 8, onset-3-of-8 per stratum | shuffle + permstrata | iid | space 56^2=3136 |
| R5 | 10 clusters, per-cluster probs {5:.3, 7:.3, never:.4} | probabilities | iid | space 3^10; all-treated draws fail (~0.6%) — conditioning argument documented; expect nfail>0 occasionally |
| R6 | 6 clusters {4,never}×{3,3}, 8 periods | shuffle + permexact | iid | 20-assignment enumeration; cleanest exactness test |
| R7 | as R1 + analytic weights [aw] (cluster-level weights 1..4) | shuffle | iid | weighted estimator, same exchangeability |

PRIMARY acceptance per config: empirical rejection COUNT at alpha = .05 and .10 inside
the central 95% binomial acceptance region for Binomial(M=400, alpha), normal
approximation rounded outward — these integer counts are authoritative:
alpha=.05 -> accept counts [11, 29]; alpha=.10 -> accept counts [28, 52].
SECONDARY: chi-square grid-uniformity (flag <0.01); null-band check: count of reps with
comparator inside [perm_lb, perm_ub] within [348, 372] (Binomial(400, .90));
R3 permtci: count of CIs covering 0 within [371, 389] (Binomial(400, .95)).
Secondary flags are findings to investigate, not automatic FAILs. `perm_failed` must be 0 on all R1-R4/R6/R7 reps (any nonzero =
investigate); R5 failures documented and expected rare.

## M2 — adversarial line test size calibration

Per rep: null DGP (below), run wooldid `adversarial` (or `altmethod_adversarial`),
collect `e(adversarial_pval_att)` (+`_ame` where posted), `e(joint_horiz)` secondary,
rc, seed. Panels: G clusters × 12 periods, unit FE + time FE + N(0,1) noise; treated
half at onset 6; `att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod`.

Null-fidelity PRETEST (required before each config's MC fires): one deterministic
large-N tiny-noise run planting the config's null effect; the reported ES coefficients
must lie on the intended line to reldif <= 1e-4 (guards against endcap-pooling
artifacts silently breaking the null). If a pretest fails, the config's panel/window is
redesigned so no pooling occurs, and the change is a documented amendment.

Configs (M = 1000 binary / 500 others):

| id | G | null | method |
|---|---|---|---|
| A1 | 40 | flat zero effect | adversarial |
| A2 | 40 | planted line: effect = 0.15*(relyr+1) all periods (zero at ref -1) | adversarial |
| A3a | 40 | flat | altmethod_adversarial |
| A3b | 40 | planted line (A2) | altmethod_adversarial |
| A5 | 10 | flat | adversarial — documented degradation, NOT pass/fail |
| A6 | 40 | flat; errors AR(1) rho=.7 within cluster | adversarial |
| A4 | continuous-treatment ame p — DESIGN PENDING the V_es_ame rank-1 investigation (INV-1, running). If no sound full-rank configuration exists, that is itself a reportable finding and A4 is replaced by a continuous-data att-side config. |

PRIMARY (A1, A2, A3a, A3b, A6): the test is asymptotic (F with cluster-dof), so:
PASS if the rejection count at alpha=.05 is inside [36, 64] for M=1000 (central 95%
binomial region, outward-rounded; for M=500 configs use [15, 35]); MARGINAL/documented
if the rate is in [.025, .085]; FAIL beyond (real mis-sizing). At alpha=.10:
PASS counts [81, 119] (M=1000; M=500: [36, 64]); MARGINAL rates [.070, .135].
A5: report only.
SECONDARY: chi-square uniformity; e(joint_horiz) sizing in A1/A6 (documented).

## M3 — pretrendpower validation

M3a (deterministic, with the honestdid grid): feed the SAME `e(b_es_att)`,
`e(V_es_att)`, `e(numpre_att)` to Roth's `pretrends` package (v0.5.1, installed
2026-08-24). Checks: (i) package slope estimate == `e(roth_slope_att)` (they fit the
same GLS line; reldif <= 1e-6 expected — investigate beyond); (ii) package's
slope-for-power-50/80 vs `mdt_nis` rows: |reldif| <= 5e-2 PASS (wooldid uses exact
quadrature, the package a normal-critical-value approximation of the same convention);
(iii) package power evaluated AT wooldid's mdt_nis in [.78, .82].

M3b (MC, M = 400 + 100-rep pilot): fixed binary design G=40 (A1 shape), window
verified pool-free by the null-fidelity pretest. Pilot at slope 0 records the
median mdt_wald / mdt_nis and `e(roth_nissize_att)`. Then:
- arm W: plant pre-period slope = pilot-median mdt_wald; empirical rejection rate of
  the pre-trend Wald test (`e(roth_prep_att)` <= .05) in [.72, .87] (binomial ±2se
  plus mdt sampling-variation slack).
- arm N: plant slope = pilot-median mdt_nis; empirical rate of "any pre-period
  |t| > 1.96" (recomputed from e(b_es_att)/e(V_es_att) diagonals) in [.72, .87].
- arm 0 (slope 0, reuse pilot + 300 more reps): empirical NIS rate within the exact
  binomial 95% region around the mean reported `e(roth_nissize_att)` (validates the
  size-exceeds-nominal claim quantitatively).
Planting protocol note: the planted "slope" is a linear PRE-period trend in the
treated outcome matched to the mdt definition (per-relative-period trend through the
reference); the post-period carries the same line (irrelevant to pre-tests).

## M4 — fastnose parity grid (deterministic)

For each grid config: run wooldid full and with `fastnose`, compare every shared e(b)
entry (and e(gteffects) cell-by-cell where present): PASS iff |diff| <= 1e-12 OR
reldif <= 1e-10 (documented fast-path divergence is ~1e-15 abs / ~1e-13 rel worst;
this gives >=100x headroom while catching real divergence). Worst diff reported per
config. Grid: >= 20 configs spanning binary/continuous, aw/pw, subgroups, ES windows,
esfixedbaseperiod, anticipation, poisson, gteffects, att_it/att_i estimands, lattice.
Refusal contract: fastnose × {unconditionalse, adversarial, jointtests,
pretrendpower(), impvalues, permtstat, permtci} each refused loudly — rc pinned (198
expected) per pair.

## M5 — honestdid grid (deterministic)

>= 10 cases spanning: dta_binary, the admin panel, rich suite data, subgroup ES
(b_es_att_s#), anticipation(1), asymmetric windows, aw weights; each under delta(rm)
mvec(.5(.5)2) AND delta(sd) mvec(.5(.5)2). Rules per case: honestdid rc 0; Original
row reproduces the wooldid post0 CI at reldif <= 1e-6; every robust row finite with
lb < ub; delta(rm) CI length monotone nondecreasing in M (jitter tolerance 1e-8).
Breakdown M recorded for the report. The ame handoff joins the grid ONLY as INV-1
allows. User-calibrated qualitative standard: honestdid-side p/CI agreement within
~5e-3 absolute would already count as "pretty similar"; our rules are stricter.

## M6 — reghdfe common-onset oracle

Prototype protocol delegated (running): common-onset TWFE event study via reghdfe must
reproduce e(b_es_att) (reldif <= 1e-8) and e(V_es_att) (reldif <= 1e-6) and the att
aggregation (reldif <= 1e-8). If the prototype fails to match, the oracle is dropped
with the mismatch documented — the honestdid grid + identity checks still cover the
surface. If it matches: variants aw + asymmetric window (+ poisson only if cheap).

## AMENDMENT A4 (2026-08-24, post-INV-1, before any A4 rep fired — documented deviation)

INV-1 resolved as a DEFECT-class finding: under cohort-level cluster(i), the CRVE
scores of every cell-specific dose-slope coefficient are identically zero (each
saturated treated cell sits inside one cluster), so e(V_es_ame) collapses to exactly
rank 1 = Var(phi_baseline)*11' and NO cohort-clustered contreat config can ever
exercise a well-posed adversarial-ame p (the refusal is correct given that V; the MC
probe verified the true covariance is full rank and ~5x larger). Sound config classes
(probe-verified): unclustered/robust, or clustering strictly finer than the
cohort-by-period cells. A4 is therefore fixed as:

- A4: 30 cohorts x 6 units x 12 periods, dose varying within cohort-period cells,
  flat null, iid errors, `robust` VCE, adversarial + jointtests; M = 504; collect
  e(adversarial_pval_att), e(adversarial_pval_ame), e(joint_horiz), e(joint_Choriz).
  PRIMARY bands: the M=500 counts of §M2 applied to BOTH the att and the ame p
  (alpha=.05 -> [15,35]; alpha=.10 -> [36,64]; marginal corridors as in §M2).

Driver: mc_adv2.do (separate from mc_adv.do so the already-running A1-A6 family is
untouched). The cohort-clustered ame pathology itself is documented in the campaign
bug report, not fixed this chunk (author directive: bug-fixing follows later).

## Outcome classification

Each check lands in: PASS | MARGINAL (documented, secondary flag or corridor) |
FAIL->FIX (real defect; fix via the established protocol: targeted cert + battery +
affected oldnew subsets rerun) | FAIL->DOCUMENT (real issue, no in-scope fix; recorded
in wooldid_old_bug_report.md per the chunk directive). The MC record (seeds, CSVs,
analyzer report) is preserved under test/suite/mc/ with the summary committed.
