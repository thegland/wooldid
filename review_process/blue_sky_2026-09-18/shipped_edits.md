# What this review actually changed in the package — 2026-09-19

Everything else in this directory is identification and investigation, per the original instruction.
These are the only edits to `wooldidinstall/`, all uncommitted in the working tree, all made on the
author's explicit instruction of 2026-09-18.

## 1. `wooldid.sthlp` — two prose edits (+1,756 bytes, SMCL-balanced, verified)

**(a) The `impvalues` entry** now says not to pair `impvalues` with a `cluster()` finer than `i` or
one that cuts across it, with the measured coverage figures and the reason: the per-cluster estimates
are read off a single fitted model, so clusters sharing a level of `i` are re-weightings of the same
coefficients. It states plainly that this is the opposite of what the default analytic path wants and
that the two cannot be satisfied at once.

**(b) "What level to cluster at"** now gives the recommended shape directly: `i` at the treatment
onset cohort, `cluster()` at the unit, and `fe()` naming the unit — and says why the `fe()` term is
not optional on an unbalanced panel, where the unabsorbed permanent heterogeneity costs precision
always and moves the point estimate whenever the holes are related to treatment.

## 2. `wooldid.ado` — the `impvalues` geometry caution (+49 lines at ado:1392)

A caution, not a refusal. The author's condition was that a guard here be *very computationally
cheap even on large datasets*; this is two `gtools` tag passes over the treated rows only, and it
runs only under `impvalues`. Measured on the shipping form at 500 `i`-levels x 50 clusters:
**0.100 s at 200k rows, 0.131 s at 1M, 0.785 s at 5M.**

### The predicate was WRONG as first written, and the review caught it

The first version counted distinct `i` and distinct `(i, cluster)` over all of `ifstatement`. Both
`IMP` reviewers falsified that. All six of the Ibragimov-Müller `over(cluster)` reads are `margins`
subpops restricted to the treated event window (`ado:2042, 2100, 2219, 2320, 2375, 2715`), so a
never-treated row contributes to **no** cluster's read. Counting those rows makes the screen fire on a
sound geometry: the common idiom that pools every never-treated unit into one `i` code —
`gen cohid = cond(ttre == ., 1, ttre - 2000)`, which is this project's own fixture idiom — gives that
single `i` level a row in every cluster.

**Reproduced directly, three builds, separate Stata processes** (`scratchpad/impvfix/`). Design: three
treated cohorts each confined to one region, never-treated units spread across all regions, so the
treated levels nest in `cluster(region)` while the pooled never-treated level does not.

| case | original predicate | corrected predicate |
|---|---|---|
| **A** nested among treated, never-treated pooled — *sound, must be silent* | **FIRES — false positive**, and miscounts "4 treated levels" (it counted the never-treated pool as treated) | **silent** ✓ |
| **B** `cluster()` strictly finer than `i` — *must caution* | fires ("4 treated levels … 36 groups") | fires, correctly: "the **3** treated levels of `i` are split across 18 `cluster()` groups" ✓ |
| **C** `cluster() == i` — *must be silent* | silent ✓ | silent ✓ |

The counts are now taken over `ifstatement & targetflag == 1`. `targetflag` is the whole event window,
pre and post — the union of all six subpops — so it must **not** be narrowed to `t >= ttre`, which
would drop the pre-window reads. `targetflag` is created at `ado:1009` and the guard sits at `ado:1392`
in the same program, so it is in scope.

### Two further corrections the review forced into the message text

- The original text asserted flatly that refining `cluster()` below `i` "shrinks it further with no
  change to the data". **False on a measured class**: going from 8 cluster levels to 48 on the
  reviewer's `r1` design shrinks the standard error by exactly zero — it is bit-identical. The
  message now says the standard error *can* shrink and the interval *can* cover below its level.
- Nesting is **sufficient** for the per-cluster reads to be linearly independent, **not necessary**.
  A panel where a couple of units change cluster mid-panel (migration, a crosswalk revision, a firm
  changing owner) straddles without meaningfully damaging the inference — measured se ratio **0.968**
  against the geometry the screen permits — and still trips the screen. The message now says so, and
  reports the straddle magnitude so a user can tell 2-of-24 from 45-of-3. This false-positive class is
  ordinary data, and it is the reason this is a caution and not a refusal.

### Why a caution and not a refusal, against the precedent

This is the `IMP` proposal's Variant B. It does not contradict the **T4** ruling of 2026-09-15, which
declined a guard about the *estimand* under cross-cutting clusters; this is about the *standard
error*. The closest precedent, **PS-64**, ruled **for** a caution in exactly this class. The proposal
argues for escalating to a refusal (Variant A) and the case is real, but Variant A also refuses the
mover class above, and both reviewers found `test/test_fastcollapse.do` F15 and `adm_impv` would
abort under it. Variant B leaves both green. **The escalation is the author's to rule on.**

## Verification state

- **Bit-identity of the shipped build against pinned `c1d8979`** — measured on the *rewritten* guard,
  not just the first version. Three arms in separate Stata processes: one where the caution **fires**
  (`cluster()` finer than `i`), one where it is silent (`cluster() == i`), and one with no
  `impvalues` at all. `e(b)`, `e(V)`, `e(wdidmainresults)` every cell at `%21x`, plus `e(N)`,
  `e(df_r)`, `e(r2)`, `e(level)`, compared as strings: **51 lines per build, zero differences, all
  three arms rc 0.** The guard is display-only, including on the arm where it prints.
- `python3 test/check_validate_sync.py` — passes (worker/validate syntax lines byte-identical;
  c_local closure holds, 46 assigned / 102 returned; every syntax option covered).
- `wooldid.sthlp` SMCL braces balanced: 4,341 open, 4,341 close.
- The guard is silent on every non-`impvalues` path by construction.

## Not verified by me: the caution will fire on two of the project's own fixtures

Both `IMP` reviewers established that `adm_impv` and `test/test_fastcollapse.do` **F15** cross among
**treated** cohorts (cohid 9/12/15 are national across 48 states), so the `targetflag` predicate fix
does **not** spare them — the caution will print on both. Reviewer analysis says both stay green
under a caution rather than a refusal: `degeneracy_expected.csv` ships no `ncaution` column and
`adm_impv` has no row in it, a caution is not a stem, and F15's `wooldid` calls are bare `qui`. **I
did not run the suite to confirm this** — the review was instructed to stay outside it — so treat
"stays green" as reviewer analysis, not a measurement. Worth knowing before the next battery run
prints something new.

## One process note, recorded because it cost a run

The first attempt to test this put both builds in **one** Stata session and switched with
`adopath ++`. Stata had already loaded `wooldid_worker` from the first build and did not reload it,
so both arms silently ran the same code and the guard appeared not to fire at all. This is the hazard
`review_process/agent_brief.md` documents under "Comparing two builds" — one process per build, dump
to disk, compare offline. The results above use separate processes.
