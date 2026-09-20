# SKEPTIC's verdicts on hv-bootrank-01 .. -04  (agent: refuter)

Build: `adopath ++ .../scratchpad/bluesky/build`, asserted in every log by `which wooldid`.
Stata 19, c(MP)=1. Data: MY OWN generator `gen_skep.do` (36 states x 4 counties x 2006-2019,
5 single-state early adopters + three multi-state cohorts + 13 never-treated, 3 regions,
AR(1) county noise, NO pre-trend). The project's test suite was not used. Every number below
comes from a log in this directory; none of the original prober's do-files were run.

---

## hv-bootrank-01  -> CONFIRMED (low)
Own repro: `r1b_post.do` (+ `r1c_repost.do`), `r2_pilot.do`, `r5_b60.do`.

* `ereturn post` has no rank screen. r1b case A: a symmetric `a*a'` matrix, exact rank 1 of 3,
  posts rc 0 and comes back verbatim (V11 .0025, V12 .0015, V33 .0049, rank 1). Case F: a
  bootstrap-shaped `dev'dev/(B-1)` from 3 draws over 5 columns, rank 2 of 5, posts rc 0 and
  `test c1 c2 c3 c4 c5` returns df 2 after printing "Constraint N dropped".
  => the seeded claim's "the repost bypass is what lets a rank-deficient e(V) through" is false,
  exactly as the finding says.
* rank(e(V)) = e(boot_N_V) - 1 whenever that binds: r5 seed 5101 N_V=13 rank 12 of 67;
  seed 5102 N_V=12 rank 11; r4 N_V=4 rank 3. The finding's arithmetic reproduces on my data.
* `test` announces dropped constraints, so "no constraint dropped" at q <= rank is correct.

No correction needed. Severity stays low: nothing is silently wrong here.

## hv-bootrank-02  -> CONFIRMED (low), with one refinement the prober did not have
Own repro: `r1b_post.do` (post cases, top level) and `r1c_repost.do` (repost cases, inside an
eclass program -- `ereturn repost` is rc 152 outside one, which is a harness fact, not a defect).

| case                              | `ereturn post`        | `ereturn repost`      |
|-----------------------------------|-----------------------|-----------------------|
| symmetric PSD, exact rank 1 of 3  | rc 0, verbatim        | rc 0, verbatim        |
| V[1,1] = -1e100 sentinel          | **rc 0, VERBATIM**    | rc 0, verbatim        |
| nonsymmetric (off-diag + .01)     | rc 0, e(V) ALL ZERO   | rc 0, e(V) ALL ZERO   |
| **missing value on the diagonal** | **rc 504**            | **rc 504**            |
| finite negative diagonal          | rc 0, lincom se = .   | (same family)         |

The code's stated justification -- "`ereturn repost V=' does NOT re-run `ereturn post''s
positive-definiteness screen - MEASURED, a 2x2 with V[1,1] = -1e100 reposts rc 0 and comes back
verbatim - and the -1e100 sentinel rows ... would otherwise trip it" (ado:10046-10049) -- does
not hold on this Stata: post accepts the identical matrix. My r504 row is the refinement: the
one screen `post` does apply (missing values) `repost` applies too, so there is no screening
difference between them at all to justify the comment. The real reason `repost` is right is
that `ereturn post` wipes the rest of e() and moves its source matrices -- a sound choice, a
wrong justification. "Zeroes a nonsymmetric OR highly singular VCE" is half right: nonsymmetric
yes, singular no.

## hv-bootrank-04  -> CONFIRMED (low)
Own repro: `r2_pilot.do`. Verbatim from r2_pilot.log line 482, in a run where e(boot_N_V) = 6:
`  e(V) holds the bootstrap covariance (complete-replication listwise, e(boot_N_V) replications;`
The token, not the value. Two lines above it, `e(boot_completed)` interpolates correctly, so the
house style is not in doubt. This is the disclosure that would have mitigated finding -03, and
it is a one-token fix.

---

## hv-bootrank-03  -> PARTIAL (severity stays high, but for a different reason)
Own repros: `r2_pilot.do` (reps 30), `r5_b60.do` (four seeds at reps 60), `r3_big.do`
(reps 300 + a 400-block study), `r4_collapse.do` (six subgroups), `r6_df.do` (Monte Carlo).

### WHAT REPRODUCES (all of it, on my own panel and my own arithmetic)
* The two moment passes are what the finding says. BOTH sides are locked against wooldid's
  own output: my per-column recompute of `main_att` in the pilot is 0.078807 and the printed
  `e(wdidmainresults)` se cell is 0.078807; my listwise recompute equals `e(V)`'s diagonal
  with MAXRELDIFF_DIAG = 0. I recomputed BOTH from `e(bootdist)`:
  per-column SD over that column's available draws, listwise complete-case covariance over
  the non-excluded columns. `MAXRELDIFF_DIAG = 0` against `e(V)` in every run, and my
  per-column availability counts equal `e(boot_navail)` exactly.
* e(boot_N_V) really does collapse: **6 of 30** (pilot, 3 subgroups), **13/12/21/14 of 60**
  (four seeds), **80 of 300**, and **4 of 60** with six subgroups. rank(e(V)) = N_V - 1
  whenever that binds.
* The silence reproduces. At reps(60) no caution fires: the reps() caution needs 2/reps > .05
  i.e. reps < 40, and the low-availability caution needs a column below 50% (my thinnest was
  188/300 = 63%). The only caution my six-subgroup run printed was about identical analytic
  SEs. And the one line that could name the count prints the token (finding -04).
* The exported numbers differ from the printed ones, materially, and `esttab`/`lincom`/`test`
  do carry the exported ones: r4, printed main_att se .07589904 vs esttab (0.073936).

### WHAT DOES NOT REPRODUCE: the magnitude and the "systematic, downward" characterization
The finding's headline is `main_att` exported/printed = **0.504**, "median .634", "biased
DOWNWARD", "systematic rather than noisy". Measured:

1. **At reps(300), where N_lw = 80 and the sampling noise in an SD is only ~8%, the median
   ratio is 0.936**, not 0.5-0.63. Full-availability columns 0.944; the columns that actually
   lose draws 0.845; `main_att` 0.845. So the systematic selection component is about
   **-6% in SE for the median column and -15% for the affected ones** - real, one-directional,
   and an order of magnitude smaller than claimed. Triangulating the three instruments - the
   reps(300) run (median .936 at N_lw 80), the 400-block study (median of block medians .898)
   and the four independent reps(60) seeds (.832/.876/.910/.885) - against an MCAR baseline of
   ~.98 puts the systematic component at **-6% to -11% in SE for the median column** and about
   **-15% for main_att and for the columns that lose draws**.
2. **The selection is real** - I confirm that half. Against 1000 random 80-row subsets of each
   column's own available rows, 22 of 67 columns sit at p_lower <= .05 (expected 3.4);
   `main_att` p_lower = .010. The mechanism the finding describes (conditioning on the
   single-cluster cohorts being drawn deletes the extreme resamples) is the right mechanism.
   It is simply not worth a factor of two.
3. **The rest is noise, and it is not one-directional.** A 400-block study on the same 300
   draws (each block a random 60 draws, both passes redone inside the block, median N_lw 16 -
   i.e. the exact regime the finding was measured in) gives for `main_att`:
   min .3164, p05 .5833, med .8298, p75 .9296, p95 1.0428, max 1.1679; 12% of blocks EXCEED 1.
   **The finding's 0.504 sits between the 2nd and 5th percentile of that distribution.**
   (The 400 blocks overlap heavily - 60 drawn from 300 - so that percentile is indicative; the
   clean corroboration is the four independent seeds, whose main_att ratios are .944, .896 and
   the pilot's 1.3205.)
   Independently: my pilot's `main_att` ratio was **1.3205** (exported 32% WIDER than printed)
   and my six-subgroup run's `ES_att_post0` was **1.32**. A run-to-run draw, not a factor.
4. **The count evidence is not what it looks like.** "65 of 84 ratios below 1" and "17 of 19
   pooled rows" are not 84 or 19 observations: every column's exported SD rests on the SAME
   listwise rows, so they move together. The unit of observation is the RUN. My four seeds
   give run-level median ratios .832 / .876 / .910 / .885 - a consistent but modest tilt.

### THE CORRECTION: the systematic harm is the DEGREES OF FREEDOM, not the level
e(V) carries N_V - 1 degrees of freedom while `test`, `lincom`, `esttab` and `coefplot` refer
it to `e(df_r)` = G - 1 = 35 (verified: r4's `lincom` CI is exactly b +/- t(.975,35)*se, and
its `test` prints F(2,35)). That mismatch IS systematic and IS one-directional.
`r6_df.do`, pure Monte Carlo, nominal 5% / nominal 95%:

|  N_V | df_r | q=1  | q=2  | q=3  | q=6  | CI coverage |
|------|------|------|------|------|------|-------------|
|    4 |   35 | .137 | .312 |      |      | **.8647**   |
|    6 |   35 | .099 | .187 | .332 |      | .9032       |
|   14 |   29 | .062 | .082 | .114 | .294 | .9374       |
|   30 |   35 | .052 | .056 | .069 | .108 | .9473       |
|   60 |   35 | .049 | .047 | .050 | .061 | .9536       |

The last row is what the user would get if e(V) rested on the replications the printed SEs
rest on. The finding put this in its last paragraph as "not a size demonstration here"; it is
the whole of the systematic damage, and at N_V = 4 it is a nominal-95% interval that covers
86.5% and a nominal-5% two-coefficient test that rejects 31% of the time.
A live instance: r4's `test ES_att_pre2 ES_att_pre3` returns **p = .00136** on a panel whose
pre-trend is exactly zero by construction.

### ONE POINT OF FAIRNESS
"Total silence" is slightly overstated. The footer does say "complete-replication listwise",
and it does name esttab, coefplot and honestdid as tools that export the stored family. What
it cannot say is the COUNT - and only because of finding -04.

### WHAT IS STILL TRUE AND STILL SERIOUS
A user who writes reps(60) gets a covariance built on 4-21 of them, exported without the count,
to every tool that reads e(V). Severity **high** stands; the claim needs restating as
"noisy and anti-conservative", not "a factor of two narrower".

---

# NEW FINDINGS (mine, found while refuting)

## skep-sgexport-01 [high/certain] `subgroups()` silently rescales every exported standard
## error in e(V) - identical estimates, identical printed table, different export
Repro: `r7_threeses.do`. Two wooldid calls on the same panel, the same `rseed(6201)`, the same
reps(60); the second differs ONLY in dropping `subgroups(reg)`. Every estimate and every
printed standard error agrees to the eight digits the tables print (both Main Results blocks
in r7_threeses.log are identical line for line). `e(boot_N_V)` is **52 without
subgroups and 16 with**, and `sqrt(diag(e(V)))` - the number `esttab`, `coefplot`, `lincom`
and `test` use - moves accordingly:

| coefficient  | printed (both runs) | e(V) no subgroups | e(V) with subgroups(reg) |
|--------------|---------------------|-------------------|--------------------------|
| ES_att_pre6  | 0.111737            | 0.117609          | **0.082039  (-30%)**     |
| ES_att_post8 | 0.179295            | 0.152829          | **0.107046  (-30%)**     |
| ES_att_post9 | 0.150121            | 0.150121          | **0.096468  (-36%)**     |
| ES_att_post0 | 0.082181            | 0.078058          | 0.064258   (-18%)        |
| main_att     | 0.076951            | 0.076059          | 0.077015   (+1%)         |

Requesting subgroup effects is not a statement about the main event study, and it changes
nothing the user can see; it narrows the exported covariance of the main event study by up to
a third, because the extra columns shrink the listwise complete-case intersection. This is the
same root as hv-bootrank-03, but it is the version a user can actually walk into: nothing in
the run differs except an option that is documented as adding rows to the output.

## skep-threese-02 [medium/certain] One coefficient, three standard errors, and the footer
## points `honestdid` at the wrong one
Repro: `r7_threeses.do`. In the subgroups run, `ES_att_post9` carries:
printed 0.150121 | `sqrt(diag(e(V)))` 0.096468 | `sqrt(diag(e(V_es_att)))` 0.150121 - a 56%
spread between the two stored matrices for the same coefficient. `e(V_es_att)` is a listwise
covariance over that family's columns ONLY, so it keeps more replications; measured, run A's
`e(V_es_att)` equals run B's whole `e(V)` to the six digits printed for every ES coefficient,
which pins the mechanism. `e(V_boot)` is bit-identical to `e(V)` (the finding's claim, confirmed).
The bootstrap footer says "the tools that read e(V) without being asked: esttab, coefplot and
honestdid export the standard-deviation family" - but the package's own honestdid recipe
(`sthlp:1560-1567`, verbatim: `matrix b = e(b_es_att)` / `matrix V = e(V_es_att)` /
`honestdid, numpre(...) b(b) vcov(V) ...`; and `ado:10068` calls e(V_es_*) "the honestdid
entry point") feeds the OTHER matrix. Whichever a user follows, the footer's premise that those three tools see
one covariance is false on a subgroups run. Secondary note for a maintainer: `sthlp:1318`
("since 2026-09-07 the same blocks are spliced into e(V)") and the PS-62 description of e(V)
as carrying the event-study blocks are written for the analytic state; under bootstrap
e(wdid_vcvstate) reads `bootstrap` and e(V)'s ES block is a DIFFERENT listwise covariance from
e(V_es_att). The pointer to e(wdid_vcvstate) makes that defensible, but the two matrices
disagreeing by 56% on the same coefficient is worth a sentence somewhere.

## skep-dfmismatch-03 [high/certain] e(V) carries N_V-1 df and is referred to e(df_r)=G-1;
## on a plausible spec N_V is 4 of 60 and a nominal-95% interval covers 86.5%
Repro: `r4_collapse.do` (the run) and `r6_df.do` (the size/coverage Monte Carlo).
Same 36-state panel with SIX subgroups (think census divisions) instead of three, reps(60),
cluster(state): completed 60, failed 0, **e(boot_N_V) = 4**, K = 107, rank(e(V)) = 3,
e(df_r) = 35, and the whole run prints no caution about any of it. `lincom main_att` returns
exactly b +/- t(.975,35)*se, and `test ES_att_pre2 ES_att_pre3` returns F(2,35), p = **.00136**
on a panel built with exactly zero pre-trend. The Monte Carlo (table in the -03 verdict above)
gives, for (N_V=4, df_r=35): nominal-95% coverage .8647 and nominal-5% q=2 size .312.
wooldid's OWN joint and equality tests are unharmed (rank tests on a fixed metric: p = .73,
.95 in the same run) - the damage is confined to the hand-run and exported path, which is
exactly the path the footer tells the user to expect numbers from.
A cheap disclosure fixes the diagnosis: print the realized e(boot_N_V) (finding -04) and, when
it is far below e(boot_completed), say that e(df_r) does not describe the covariance.

---

# NEGATIVE RESULTS (attacked, found clean)
* **`ereturn post` is not the barrier the code says it is, but nothing user-visible is wrong
  because of it.** I tried five distinct matrix pathologies through both post and repost
  (r1b, r1c). They behave identically in all five. No wooldid number depends on the choice.
* **wooldid's own bootstrap inference is not affected by the listwise covariance.** The joint
  pre-trend and across-subgroup equality p-values are add-one rank shares on a FIXED metric,
  so a noisy or rank-deficient metric does not bias them: on the N_V = 4 run they returned
  .73 and .95 where a Wald on the same e(V) returned .0014 on a true null. The printed
  coefficient table (SD, percentile, BC families) is likewise built from the per-column
  draws, not from e(V).
* **The per-column pass and e(boot_navail) are exactly right.** My independent recomputation
  from e(bootdist) matched e(V)'s diagonal to the last bit (MAXRELDIFF_DIAG = 0) and matched
  e(boot_navail) exactly, in every run (reps 30, 60 x4, 300). There is no arithmetic error in
  wooldid_bootfin; the defect is entirely in which rows each pass uses and in what is
  disclosed about it.
* **e(V_boot) is bit-identical to e(V)** - checked coefficient by coefficient in r7. The
  finding's claim that the documented fallback carries the same matrix is correct.
* **No sentinel leakage.** Across every run above, no `-1e100` reached e(V) (all diagonals
  positive, including the N_V = 4 run), so the "observable failure" path never fired and the
  narrow/noisy numbers arrived looking perfectly ordinary.

---

# ONE ROOT, ONE FIX
hv-bootrank-03, skep-sgexport-01 and skep-dfmismatch-03 are three consequences of the single
listwise rule at `ado:14536-14576`; they should not be counted as three independent defects.
A fix that makes e(V) rest on the same replications as the printed standard errors (or that
discloses e(boot_N_V) and refuses to call e(df_r) the covariance's df when it is not) removes
all three. skep-threese-02 is separate: it is about two stored matrices disagreeing.
