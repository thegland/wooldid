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

## PS-62 — "zero the sentinel row and column": the diagonal is kept, not zeroed

**The ruling.** *"Zero the sentinel row and column. Where a row's variance is unidentified,
zero its row and column rather than carrying the -1e100 sentinel into a matrix with
off-diagonals, which the posting screen can zero wholesale."*

**Built:** the off-diagonals of a sentinel row and column are left at zero, and the -1e100
diagonal cell is left exactly as the partial-identification rescue posted it.

**Why the literal reading was not taken.** The reason the ruling gives for zeroing is a real
hazard — Stata's `ereturn post` screen does zero a nonsymmetric or highly singular VCE
wholesale, which is what `wooldid_amescreen , mode(vzero)` exists to detect. It does not arise
on this path. The splice reaches `e(V)` through `ereturn repost V=`, not `ereturn post`, and
`repost` was measured (probe at `~/.claude/jobs/be9644d0/tmp/w1/probe_ps62.do`, 2026-09-10) to
accept verbatim both an indefinite matrix and one carrying a -1e100 diagonal beside live
off-diagonals: the posted matrix came back unchanged and `ereturn display` still rendered that
row's standard error as missing. With the off-diagonals of such a row held at zero, the matrix
the splice posts has exactly the block structure the run already ships today.

**What the literal reading would cost.** Zeroing the diagonal turns a variance of -1e100 into
0, and a standard error that prints as `.` today would print as `0` in `ereturn display`,
`esttab` and `lincom` — a fabricated certainty in place of an honest absence, and a reported
number moving on a fix that was not supposed to move any. That is the one thing the campaign's
standing instruction rules out, so the diagonal was kept and the question is put here.

**To reverse it** is a two-line change in `wooldid_vspl` (the Mata kernel): also write 0 into
`V[a, a]` for every row whose diagonal is not positive-finite. The comment block above the
function records this decision at the site.

## PS-45 — the bootstrap simultaneous band studentizes by the coefficient-draw SD, even under `boottstat`

**What was built.** With `bootstrap()`, `simultaneousbands` forms the band from the replication
accumulator the joint test already uses: for each listwise-complete replication `r` it computes
`M_r = max_h |A[r, h] - O[h]| / sd_h`, where `sd_h` is the standard deviation of the
coefficient draws at horizon `h`, and the critical value is the `level()` quantile of `M_r`.
The band is then `O[h] +/- c * sd_h`.

**The question.** Under `boottstat` the run already carries each replication's own analytic
standard error in `e(bootdist_se)`, and every other reported quantity on that run is
studentized by it rather than by the draw SD. A percentile-t band would use `se_r[h]` in the
denominator of `M_r` and `se_obs[h]` in the endpoints, matching the rest of the `boottstat`
surface. What is shipped instead is one construction for both, so a `boottstat` user gets a
band whose scale is the coefficient-draw SD while the pointwise intervals beside it are
percentile-t.

**Why it was built this way.** The ruling asked for the band to come from "one pass over the
accumulator that already feeds the joint test", which is the coefficient accumulator; the
`permstepdown` procedure shipped in the same batch studentizes the same way; and a single
construction means the band does not change meaning between bootstrap displays. It is also the
only form available without `boottstat`, so the alternative is two constructions rather than a
replacement.

**Cost of the alternative.** Small and contained: `wooldid_bootband` would take the replication
standard-error matrix as a fourth argument, skip replications with a missing or non-positive
`se_r`, and use the observed analytic standard errors for the endpoints — roughly the shape the
studentized coefficient tail already has. The `esband_engine` marker would need a third value
so the store says which construction produced the band.

**No number moves unless the author rules for the change**: `simultaneousbands` is off by
default, so this affects only runs that ask for it.

## PS-46 — the measurement is in, and it points away from the instrument the ruling named

**The ruling.** *"Conditional. First establish whether the permutation fast-mode default
actually breaks esttab-class table tools. If it does, post an ARTIFICIAL covariance that lets
those tools run and that is transparently false to anyone who inspects it (an all-zero matrix
being the example the author gave), rather than leaving nothing posted. If it does not break
them, leave as is."*

**Answer to the condition: yes, it breaks them — silently.** Full report and the raw probe
record at `review_process/research_2026-09-10/permute_fastmode_esttab.md`. `esttab` never
refuses. It prints the coefficient, an **empty** standard-error row, and then a footer
asserting *"Standard errors in parentheses"* with a full star legend — and writes that same
assertion into `.tex` and `.csv` files that outlive the session. `estout`'s default
`cells(b)` is clean; Stata's own `estimates table` and `estat vce` refuse honestly at r(321)
/r(322) because `e(properties)` is `b`.

**But the named instrument fails the ruling's own test.** The all-zero covariance is
"transparently false" on exactly two surfaces — `matrix list e(V)` and `estat vce`, both
commands a suspicious analyst runs deliberately. Everywhere else it is invisible or worse.
Measured:

- `esttab ., ci` goes from a blank bracket row to `[0.824,0.824]`, `[-0.0457,-0.0457]` — not
  a transparently false number but a plausible, infinitely precise one, written to file.
- `estimates table S, b se` goes from an honest r(321) refusal to a printed `0` in a
  standard-error row.
- `di _se[main_att]` goes from `_se not allowed when e(V) is not present` to `0` at rc 0.
- `r(table)` goes from an honest 1x6 b-only matrix to a full 9x6 carrying `df 79` and
  `crit 1.9904502` for statistics that do not exist — which downstream scripts read.
- On the `esttab` surfaces the remedy was meant to protect, it converts a blank cell to `(.)`
  and leaves the false footer and the star legend exactly as they were. No grade improves.

**A route that was not available when the ruling was made.** `estout.ado` (3156-3172) checks
for posted `e(se)`, `e(p)`, `e(ci_l)`, `e(ci_u)` matrices and uses them **ahead of** anything
derived from `e(V)`. Posting those — and no covariance at all — makes `esttab` print the
**real permutation numbers**, while `test`, `lincom`, `_se[]`, `estat vce` and
`estimates table` go on refusing honestly, because `e(V)` is still absent. It fabricates
nothing and it fixes exactly the tools that were lying. Measured working:
`esttab ., se` → `(0.130)`, `(0.0767)`; `esttab ., p` → `(0.095)`, `(0.857)`.

**Why this is put to you rather than built.** It moves a printed number, and one input has no
designated source. `e(p)` is obvious — `e(perm_p)` is already posted and is what the tables
print. `e(se)` is not: `wooldid_permwire`'s own header says the standard-error column is left
blank under the fast-mode default **by design**, and the probe used `e(perm_se_shift)` as an
arbitrary stand-in. So:

| option | what `esttab` would then print |
|---|---|
| **A. the ruling as written** (all-zero `e(V)`) | `(.)` under the same false footer; fabricated zero-width CIs; honest refusals elsewhere become `0` |
| **B1.** post `e(p)` only | `esttab ., p` correct; `esttab ., se` still blank under a false footer |
| **B2.** post `e(p)` and `e(se)` from `e(perm_se_shift)` | both correct, but commits to calling the shift standard error *the* standard error of a fast-mode run |
| **B3.** B2 plus `e(ci_l)`/`e(ci_u)` from `e(perm_ci)` | also fixes the CI path — the one place A actively fabricates |
| **C. leave as is** | the silent misreport stands; only the help changes |

**One thing is being fixed now regardless of the ruling, because it is simply false.**
`wooldid.sthlp` says that under fast mode *"`lincom`, `esttab` standard errors, and anything
else needing a variance-covariance matrix will refuse."* Measured, `esttab` does not refuse;
it prints an empty standard-error row beneath a footer claiming standard errors are there.
That sentence is corrected in the P8 documentation sweep whichever way PS-46 goes.

**Untested, and it is the one thing worth testing before B ships:** every `e(perm_p)` in the
fixture was above .095, so no star was ever due and the star path was never exercised.
`coefplot` was also not probed, and under A it would draw zero-width whiskers.
