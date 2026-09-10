# Held for the author, 2026-09-09

Items reached during the batch implementation of the 2026-09-07 rulings where carrying the
ruling out as written would contradict another standing ruling. Nothing below is built.

## PS-2 — "the stored covariance and degrees of freedom follow `bootdisplay()`"

**What the code actually does today** (checked 2026-09-09, ado 8837-8856, 5971-5975, 6058-6062,
8278-8284):

| stored entry | value | set by |
|---|---|---|
| `e(V)` | the bootstrap replication covariance, every display | `wooldid_bootwire`, `repost V` |
| `e(V_analytic)` | the analytic covariance, hidden twin | same |
| `e(df_r)` | `.` unless `bootdf()` was given; `nclus-1` under `bootdf(clusters)`, else the number given | `local wdbdf = .` then the three `bootdf` branches |
| `e(df_r_analytic)` | the analytic residual df, hidden twin, only when `e(df_r)` is blanked | same |

**The first half of the ruling has no lever.** `e(df_r)` is already missing by default under
every display, precisely so that a hand-run `test` or `lincom` goes to chi-squared/z and
reproduces the printed numbers. It becomes finite only when the user asks for it with
`bootdf()`. So "the degrees of freedom follow `bootdisplay()`" describes what already happens,
except in the case where the user has explicitly overridden it - and overriding the reference
distribution is what `bootdf()` is for.

**The second half runs into two standing rulings.** Two readings, both blocked:

1. *Make `e(V)` match the printed p and interval under a rank display.* The percentile, BC and
   BCa families produce an interval and a rank p, not a covariance. The only way to store them
   as a variance is to back a standard error out of the interval width - which the author
   overruled himself on (2026-09-03/04 display rulings) in favour of reporting the bootstrap
   standard deviation, on the ground that an interval-implied standard error is not a standard
   error. Note also that the printed **standard-error column** under the rank displays already
   IS the bootstrap standard deviation, so `e(V)` and the printed SE column already agree; what
   diverges is only the p and the interval.

2. *Under `boottstat`, store the analytic covariance, since the printed SE column is
   `se_analytic`.* This is implementable and is the only reading with real content. It
   contradicts the footer at ado 9800, which states that `e(V)` and `e(df_r)` carry the
   standard-deviation family whichever family is displayed, and it inverts the
   `V_analytic` hidden-twin design (the bootstrap covariance would become the twin).

**Recommendation.** Ruling 2 is the only one with a lever, and taking it means rewriting the
9800 footer and swapping the twin. That is a contract change, not a bug fix, so it is held.
If the intent was the weaker thing the practice sweep proposed - name `esttab`, `coefplot` and
`honestdid` in the footer beside `test` and `lincom` - that is a one-line documentation change
and can ship immediately.

**Ruling wanted:** (a) swap `e(V)` to analytic under `boottstat` and rewrite the footer; (b)
name the exporters in the footer and leave the stored family alone; (c) both.

## BS-15 — the degeneracy tolerance is a judgment call, not a ruled number

The ruling says "when the subgroup pairwise covariance is numerically rank-deficient,
blank the pairwise SEs and both equality tests with a named caution". It fixes no cut.
Built and shipped at 1 - 1e-2 on the implied correlation between the two subgroups'
estimates. The measurement, on `test/data/dta_binary.dta`, the fixture the bug sweep used:

| run | main att row | pretreatment att row |
|---|---|---|
| `cluster(i)` | rho .99940 | rho .99715 |
| `robust` | rho .51141 | rho .49779 |

On the clustered run the contrast standard errors are 1/44 and 1/20 of their robust
counterparts - the same pathology on both rows. A cut at 1 - 1e-3 was tried first and
caught only the main row, leaving the pre-treatment row reporting a standard error twenty
times too small, which is the defect the ruling names. 1 - 1e-2 catches both and still
sits half the range away from a healthy design.

Across the 34-configuration battery it fires on seven, all of them subgroup runs, at rho
between .9990 and 1.0000. The equality p-values it withheld were 1.1e-12, 3.3e-88,
9.2e-46 and 4.2e-25 - the residue-driven false certainty the ruling is about. Nothing
outside a subgroup run moved.

The cut lives in one local, `wdsgdtol` in `wooldid_sgdeg`. **Ruling wanted only if the
author wants it elsewhere**; the work does not wait on it.
