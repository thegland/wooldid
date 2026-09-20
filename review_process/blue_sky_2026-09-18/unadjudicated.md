# Open questions this review could not settle

**Status line, dated 2026-09-18.** These are stated as questions, not findings. Each names what
would settle it. Nothing here should be acted on as if it were established.

## 1. Is the bootstrap `e(V)` rank deficiency a defect or an intrinsic property? — REFRAMED

Recon measured `reps(8)` with 20 coefficients giving `rank(e(V)) = 7`, and a `test` over three
event-study coefficients returning p = 1.5e-14 on df 3 with no constraint dropped. A probe cell then
escalated this to "rank-deficient at EVERY `reps()`". The refuter for that cell was lost to the
usage limit, so the coordinator checked it directly:

| reps | K | rank(e(V)) | min eigenvalue | `test` p (3 ES coefs) |
|---|---|---|---|---|
| 8 | 20 | 7 | -1.35e-16 | 1.05e-06 |
| 40 | 20 | **10** | -4.96e-17 | 7.06e-07 |
| 200 | 20 | (run exceeded its deadline) | — | — |

At `reps(40)` the rank is **10**, not `min(R-1, K) = 20`. So the deficiency is **not** a low-reps
artifact. The likely explanation is structural: `att`, `att_i`, `att_it` and `att_itime` are four
different weighted averages of one underlying set of cell effects, so the 20 reported quantities
live in a lower-dimensional space and a genuine covariance among them **must** be singular.

If that is right, the framing flips: the bootstrap path is reporting an *honest* rank-deficient
covariance, while the analytic path hides the same structure behind a covariance that is
`diag(se^2)` with event-study blocks spliced in (PS-62, ruled 2026-09-07) — full rank only because
its off-diagonals are zero by construction. **Neither supports an arbitrary joint `test` across
estimand families**, which is the real user-facing statement, and which the help already gestures at.

**What would settle it:** the `reps(200)` arm, plus a size study — run `test` under a true null at
several `reps()` and see whether the rejection rate is nominal. Until then this is a reframing, not
a defect.

## 2. Do `bootstrap()` and `permute()` keep the coverage gain at a finer-than-`i` cluster?

The finer-than-`i` MC covered the **analytic** path only, and the headline recipe is validated only
there. The `bootstrap(standard)` keep-`i` branch (`wdb_nest = 0`, `ado:9349`) and the
`permonsetcohorts` rebuild are one day old and were never coverage-tested. This matters because the
`SMP` mechanism gives concrete reasons to expect trouble: the regime flag is global over a per-unit
property, and the resampling universe is not `e(sample)`.

**What would settle it:** a `bootstrap(standard)` coverage MC in the keep-`i` regime at `reps(99)`,
roughly an hour per 100 draws on this hardware. Note the open, unruled observation pointing the
other way — `bootstrap(standard)` at `cluster() == i` *over*-covers `att_i` at 0.9933.

## 3. Is `esfixedbaseperiod` a lever on the analytic SE understatement?

One probe cell claims the understatement is estimator-specific and that `esfixedbaseperiod` fixes
it. Its refuter was lost to the usage limit. If true it would change the author's remedy entirely —
a reference-period choice rather than a clustering choice. **Single-source; do not act on it.**
**What would settle it:** rerun the 2026-09-14 coverage MC with the reference convention as the arm.

## 4. The 45x lattice AME point miss

Honestly marked `uncertain` by its finder, and the finder gave the right reason: a large point miss
paired with a small *conditional* robust SE could be a legitimate conditional-variance statement
about an ill-conditioned design rather than a standard-error defect. What is not in doubt is that the
default path delivers neither a standard error nor a caution naming within-bin dose variation.
**What would settle it:** decide whether the conditional robust variance is the right object for an
AME aggregated over near-singular within-bin slopes.

## 5. Coverage at other cluster counts

Every coverage number in this review is at G = 40. Small-G is exactly where cluster-robust intervals
break, and G = 20 is an ordinary US-state design after dropping non-adopters. The 2026-09-14 MC was
scoped to G = 40 by ruling and a G = 200 follow-up was ruled SKIP.

## 6. What the usage limit cost

Two refuter cells (`ereturn-contract`, `mc-boot-perm-size`) and most third-lens escalations did not
run. Concretely: 15 findings carry NO_VERDICT and are single-source; the `e()` contract cell and the
bootstrap/permute **size** study are the two areas with the weakest verification in this report.
The size study matters — it is the one measurement that would say whether the replication engines
deliver their nominal level under a sharp null.
