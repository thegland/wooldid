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
