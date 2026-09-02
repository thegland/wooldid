# Cell collapse and subgroups(): can the gate come off? (2026-09-01)

Assessment requested by the author alongside the default-on switch (`111ecbb`): "assess
if there is a way to make [the cell collapse] compatible with subgroups." Assessment only;
nothing in the ado changes here. Companion to `cell_collapse_2026-08-31.md` and
`cell_collapse_followup_2026-09-01.md`.

**Status of the switch itself.** The full suite (`test/suite/run_full.sh`) on `111ecbb`
finished after that commit was written: oldnew 473/473 compared, `unexplained=0`,
`prov_fail=0`; cert 34/34; newfeat 7/7; battery IDENTICAL; coverage OK - "ALL COMPONENTS
PASS". The default-on polarity is verified end to end.

**Headline.** Yes, there is a way, and the evidence says the risk is narrower than the
gate assumes. Across 522 collapsed-vs-uncollapsed subgroup runs (472 suite configs plus a
50-config adversarial round built to break the claim), every difference in *reported
structure* sits on a run whose uncollapsed output is itself a function of row order: the
same cells flip between "missing" and "se = 2^-32" when the input rows are merely
permuted, with no collapse anywhere. Two routes keep the cardinal rule (never change a
reported number) intact: a **partial collapse** (exact, ~8 lines, 1.45x on the timing
case) and a **collapse-then-verify** (per-read fallback to the full data when a
degeneracy signal fires; ~40-60 lines, the full 5.6x on non-degenerate runs). A latent
rowname bug in the subgroup event-study assembly (mislabel or r(503)) is reachable on the
full-data path too and should be fixed regardless of what is decided here.

---

## 1. Why subgroups() was gated: the three reads

The collapse changes the order in which margins sums cell contributions, so any quantity
that is *mathematically zero* comes out as bitwise `0` in one arm and as residue
(`~1e-20` for a variance, `~1e-16` for a point estimate) in the other. Outside the
subgroup section nothing reads a number for exact zero, so the residue is last-place
noise - the property the default-on decision accepted and the help file quantifies.
The subgroup section has three reads that branch on exact zero:

1. **The omission marker** (`wooldid.ado` ~2789). Under `noestimcheck`, margins tags a
   subgroup-by-horizon margin `o.` iff its `b` AND `V` are bitwise zero (probe result,
   this session). SI-7 reads that marker and drops the horizon from that subgroup's event
   study with a Caution (2800-2855). A residue on either side of zero flips the drop.
2. **`r(table_vs)` cells** (2986, 3030). A pairwise-subgroup contrast whose variance is
   bitwise zero reports se/t/p as missing; a residue variance reports `se = 2^-31..2^-32`,
   `t ~ 1e8`, `p ~ 1e-85`.
3. **The chain-test rank guard** (3011-3024; pre twin 3078-3089). `test` decides the rank
   of the chain-contrast VCE with a tolerance; on a near-singular VCE the decision flips
   between "full rank, report `r(p)`" and "fall back to `contrast(overjoint)`", and the
   two engines disagree wildly on a singular direction (`1e-83` vs `.`). The chain
   builder itself (`wooldid_sgchain`, 5096) also skips any posted subgroup level that
   carries the `o.` marker, so a level flip changes which subgroups the chain compares.

Everything else in the section (`glevelsof` of `_t` within a subgroup, the EFF3 level
grab, the rowname pool) is invariant to the collapse because the collapse keeps one row
per `(_i, _t, targetflag, subgroup)` cell and so preserves every set of levels.

## 2. Where the two arms actually differ

Harness: the HEAD ado with line 4621 (`subgroups()` -> `fcbad`) neutralized, run as arm C
(collapse engages) and arm U (`fullprecisionslow`) on the same configs, every e() scalar,
macro, and matrix cell dumped in `%21x`, then a structural differ (rc, varlist, e(sample)
N, matrix inventory/dims/names, scalar names, Cautions, present-vs-missing cells) plus a
tier-2 pass for any cell moving more than 1e-9 relative.

**Suite sweep (D).** 472 suite configs, 255 with the collapse engaged. 3 structural
differences, all `survey_rcs.dta` with `subgroups(group)` + `fe(group_time)` /
`fe(i.group##i.t)` (the fixed effect *is* the cell, so every cell contrast has zero
variance) or the shock-subtracted exact-fit outcome ("standard errors were not
identified"; no e(V)).

**Adversarial round.** Two lenses, each told the claim and told to refute it on
*non-degenerate* data, each free to write its own configs. `wfc` (weights, fixed effects,
controls, `if`, poisson, `oostreatedcontrols`, `anticipation`, `timetrends`,
`customweightmultiplier`, person-level fe, mixed grain, admin subsamples, `ttyear`,
`covars`, singleton, imbalanced, deep pre): 25 configs, 22 engaged, 3 errored in both
arms identically (2 x r(198) poisson refusals, 1 x r(2001)). `est` (`subgroups(ttre)`,
`subgroups(cohort)`, `att_i`/`att_itime` alone, `suppresssgprimaryeffects`,
`gteffects`+`histogramcohorteffects`+`makeplots`, `cleanmatrices`, `esfixedbaseperiod
esrelativeto(-3)`, a bayes bootstrap `reps(10)`, never-treated-only subgroup via `if`,
tiny, `bigcodes`, `deeppre`, `oostreatedcontrols(1)` on a no-post shape, a straddling
subgroup): 25 configs, 25 engaged, **0 structural differences**; no matrix cell moved
more than 1e-9 relative, and the largest move the lens found by decoding the dumps itself
was 1.5e-11 relative on a confidence bound formed by cancellation (9e-16 absolute).
Between the two lenses, 36 engaged runs were adjudicated non-degenerate (well-identified
SEs, no degeneracy Caution) and every one of them matched to the last place.

`wfc` produced 4 structural differences. All four are read (2) or read (3) flips on
contrasts that are degenerate *in both arms*:

| config | dataset / spec | the differing cells |
|---|---|---|
| `wfc_bb_pw_fecat` | binary_balanced, `subgroups(sg_three) fe(categorical_control) cluster(cohort)` | `pwsgs[10,*]` `pre_att_it:2vs1`: C missing, U `se=3.3e-10, t=-4.4e7, p=1e-79` |
| `wfc_bb_antic` | same + `anticipation(1)`, aw | `pwsgs[12,*]`: C `se=4.7e-10, t=-3.2e8, p=4e-89`, U missing; `joint_pre_att_it_allsgssame` C `2e-85`, U `.` |
| `wfc_mg_clcoarse` | binary_mixed_grain, same spec | `pwsgs[3,*]` and `[4,*]` flip in opposite directions (`se=1.6e-10` / `3.3e-10`); `joint_pre_att_allsgssame` C `.`, U `2e-83` |
| `wfc_rcs_esf_if` | survey_rcs, `subgroups(group) cluster(state)`, `if t>1`, fixed base | `pwsgs[4,*]` `pre_att_itime:2vs1`: `b=-8.9e-16`, C `se=2.3e-10`, U missing; `att:2vs1` is `b=1, se=.` in both arms |

Every cell that differs has `se <= 2^-31` in whichever arm reports it - a variance of
one to a few ULPs (`2^-62`, `2^-63`, `3 x 2^-64`) on a point estimate of `0.1` - and the
*same* matrix carries identical `se = 2^-32`, `p = 1e-90` rows in both arms. Why these
specs: `cluster(cohort)` with the unit *being* the cohort (the `binary_*` datasets) makes
each cluster a union of saturated cells, so the cluster score in every treated-cell
direction is a sum of residuals that vanishes identically and the clustered VCE is zero
in those directions - visibly so: the by-subgroup `att_it` SEs on these runs are
identical to the pooled `att_it` SEs, i.e. the between-subgroup differences carry
exactly zero clustered variance. wooldid already says so: both arms print "jointtests was combined with analytic clustered SEs clustered at
the level of the cohort variable i or above ... joint test p-values are severely
oversized" and "4 reported standard errors appear to be identical under a cluster-type
VCE ... failure of analytic clustered standard error estimation", and the configuration
is refused without `safety(off)`.

**The row-permutation control - the decisive piece.** `wfc_bb_pw_fecat` re-run on the
committed HEAD ado with `fullprecisionslow` (no collapse anywhere, `e(wdid_fcol)=0`) four
times: the original row order and three seeded random permutations of the rows.
`e(pwsgs)` rows 2, 5, 7, 8, 9, 10, 11, 12 each flip between `se = 2^-31 / 2^-32` and
missing across the four orders (row 8 `att_it:3vs1`: `2^-31, ., ., .`; row 5
`pre_att:3vs1`: `2^-32, 2^-32, 2.6e-10, .`). `e(joint_att_it_allsgssame)` is `1.8e-83`
in the original order and `.` in all three permutations; `e(joint_pre_att_allsgssame)`
is `7.1e-90, 7.1e-90, 3.0e-87, 3.9e-11`. So the differences the collapse produces on
these runs are the same differences a `sort` produces, on output wooldid has already
labelled a failed VCE. This is the subgroup-section form of the follow-up report's
scale-free statement: **the collapse perturbs the answer no more than permuting the
input rows does** - here including the branch decisions, because the uncollapsed branch
decisions are themselves order-dependent on exactly the runs where they differ.

**No non-degenerate divergence was found**: not in 472 + 50 configs, not by either lens
whose brief was to find one. The tier-2 report on the clean runs shows only last-place
digit moves (`1e-15`) in p-value macros, the accepted residue.

## 3. The latent rowname bug (fix regardless)

`relyrlist_s<sg>` (2871-2879) is reset per `horizitem`, so after the event-study-by-
subgroup loop it holds the *last* estimand's horizon list (att_it when both att and
att_it run). The rowname filter at 3563 then screens **both** the `ES_att:` and the
`ES_att_it:` rows of the stacked `wdid_sg<sg>` (att rows appended 3522, att_it rows 3526)
with that one list. If att and att_it keep different horizon sets for a subgroup - which
happens exactly when the SI-7 marker drop fires for one estimand and not the other - the
name list is too short (Stata pads by repeating the last name: silent mislabel, the same
hazard the comment at 3544-3550 documents for the non-ES rows) or too long (r(503)).
This is the r(503) recorded in `4b5679c`; it is not a collapse defect - it is reachable
on the full data whenever the two estimands' marker drops differ - and the collapse
merely made it easier to hit. Fix: keep `relyrlist_s<sg>_<horizitem>` per estimand and
filter `ES_att:` rows with the att list and `ES_att_it:` rows with the att_it list. A few
lines; deserves a cert case that drops a horizon from one estimand only.

## 4. The prize

Timing case: `dta_s_admin` (29,667 rows), `wooldid y cohid year ttre, att att_it
cluster(state) subgroups(region) esprelength(2) espostlength(2) jointtests safety(off)`,
one Stata, warm cache:

| build | `e(wdid_fcol)` | seconds |
|---|---|---|
| HEAD (subgroups gated: nothing collapses) | 0 | 19.5 |
| full collapse (gate removed) | 1 | 3.5 (5.6x) |
| partial collapse (main block collapsed, subgroup section on full data) | 1 | 13.5 (1.45x) |

The subgroup section carries most of a subgroup run's margins calls (per estimand x per
subgroup x per horizon window, plus the pairwise and chain passes), which is why the
partial route keeps so little of the prize.

## 5. Routes

**A. Partial collapse.** `restore` before the subgroup section when `subgroups()` is
given, skip the final restore. Prototype exists and timed above; ~8 lines. The subgroup
section runs on the full data exactly as today, so every subgroup number and every
branch is unchanged by construction; the main block moves at the accepted last-place
residue, as every non-subgroup run already does. Cardinal rule: kept. Prize: 1.45x.

**B. Pre-detect degeneracy, then collapse.** Family (i) (all-zero / unposted VCE) is
already detected before the section runs (`vzero`). Family (ii) - a singular *direction*
of a clustered VCE - is a property of each contrast, not of `e(V)` as a whole (a
full-rank test on `e(V)` would flag every clustered run with more cells than clusters
and forfeit the prize on exactly the large runs that need it). Not sufficient alone.

**P. Collapse-then-verify (recommended if the prize is wanted).** Run the section
collapsed; at each of the three reads, compute a degeneracy signal *before* reading, and
if it fires, `restore, preserve`, re-issue that margins call on the full data, read from
that, and stay on the full data for the rest of the section. Signals:

- read (1): any `o.` marker, or any margin with `V_kk <= 1e-9 * max diag(V)` in the
  subgroup-by-horizon call;
- read (2): any `o.` marker or `V_kk <= 1e-9 * max diag` among the posted subgroup
  levels, any missing cell in `r(table_vs)`, or any pairwise contrast with variance
  `<= 1e-9 * max diag` (one signal for the level and pairwise reads: they are the same
  margins call);
- read (3): `test` reporting `r(df) < sgchainn - 1`, or the chain-contrast VCE with a
  minimum/maximum eigenvalue ratio `<= 1e-9` (`matrix symeigen` on a matrix of order
  #subgroups - 1; negligible cost).

Why this keeps the cardinal rule: every branch flip needs a quantity that is bitwise
zero on the full data. The collapse residue is `~1e-16` relative to the terms summed, so
the collapsed value of such a quantity is `<= 1e-16` of its scale - three to seven orders
below the `1e-9` band - and the signal fires. When it fires, the numbers reported are the
full-data computation, i.e. today's output for that row order (the same guarantee
`fullprecisionslow` gives). When it does not fire, no read was within `1e-9` of a branch
boundary in either arm, and the only difference is the last-place residue the default-on
decision already accepted. The residual risk is a signal I have not thought of; the
adversarial round found no divergence outside the three reads. Cost: a flagged run pays
the collapsed prefix plus the full-data remainder (bounded by ~1.2x today); an unflagged
run gets the full 5.6x. The per-read form avoids the two things a whole-section retry
would need - snapshotting the accumulators the section appends to (`f2keys`,
`sgmatricestodrop`, `relyrlist_s*`, `joint_*`) and buffering Cautions printed by a
discarded attempt - because the read decision is made once, on the call that survives.
Estimate ~40-60 lines across five sites (2741/2789, 2986, 3011, 3030, 3078). Cert: the
three degenerate suite shapes must reproduce HEAD byte for byte with the collapse on,
and a non-degenerate subgroup run must show `e(wdid_fcol)=1` with no fallback taken.

**I. Tolerance the reads themselves.** Treat a variance below `1e-9 * max diag` as "not
estimable" in *both* arms (drop the horizon with the Caution; report the contrast as
missing; take the fallback). This makes the subgroup output order-stable and honest on
degenerate runs (today's `se = 2^-32, p = 1e-89` would become `.`), and with it the
collapse could engage unconditionally. But it changes reported numbers on the full-data
path for the degenerate class - the author's call, and a separate ticket from this one.

**Z. Drop the gate and document.** After the section-3 fix the worst case is a Caution
line and a `. / 2^-32` cell flipping on runs wooldid already labels a failed VCE. The
row-permutation control shows those outputs are not stable today either. Still a
reported-number change on an observable path; not recommended without I.

## 6. Recommendation

1. Fix the section-3 rowname bug now (independent of the collapse; removes the crash
   and the silent mislabel from the degenerate class on every path).
2. If the subgroup prize matters - it is the difference between 19.5 s and 3.5 s on
   30k rows, and subgroup runs are among wooldid's slowest - implement P. If the
   effort-cap doctrine (no machinery for defects that fire only on hopeless data)
   outweighs it, implement A: exact, eight lines, a third of the prize.
3. Consider I on its own merits: it is the only route that makes the subgroup output on
   a degenerate VCE reproducible under `sort`, with or without the collapse.

## 7. Artifacts

Session scratch (not committed): the probe build (`adv/wooldid_probe.ado`), the paired
runner and differ (`adv/runpair.sh`, `adv/differ.py`), the two lens config files
(`adv/cfg_wfc.csv`, `adv/cfg_est.csv`) with their `RESULT.txt`, `%21x` dumps and Stata
logs, the suite sweep (`estscope/D/RESULT.txt`), the timing builds and logs
(`partial/`), and the row-permutation control (`permctl/perm.do`, `perm.log`). The
harness is the suite's `mkrun` with `--extraopt fullprecisionslow` for arm U, so any of
it can be recreated from the repo's suite tooling.

## 8. Addendum (2026-09-01, after the ruling): route P built, probed, pinned

**Ruling.** The author chose route P with two amendments: the collapse runs with
`subgroups()` by default behind the per-read screen, a flagged read falls back to the full
data under the default safety setting, and under `safety(off)` there is no fallback - the
run proceeds on the collapsed rows and prints a caution. Then an adversarial probe for
cases the screen might miss, before any full assessment run.

**What was built.** The `subgroups()` term left the eligibility gate (`fcbad`, ~4687; the
seven substantive exclusions stay, and an explicit `fastcollapse` with one of them is
still `r(198)`). A per-read screen `wooldid_fcdeg` (~5193) runs before each subgroup read
while `fcsg` is armed: read 1 (`levels`), the two overjoint pre-trend calls (`joint`),
the read-2 pwcompare call, and the pre-window twin. It hands `r(V)` to Mata
`wooldid_fcdegx(1e-9, levels)`, which looks at the subgroup half of the diagonal and,
for `n >= 2` non-level reads, at the eigenvalues of the adjacent-difference contrast
VCE, and returns a reason when any of: an `o.` marker in the `r(b)` names; the largest
level variance missing or `<= 0`; a level variance `<= 1e-9 x` the largest; the largest
contrast eigenvalue `<= 1e-9 x dmax`; the smallest contrast eigenvalue `<= 1e-9 x
max(dmax, largest eigenvalue)`; or a missing/`<= 1e-9 x dmax` diagonal in `r(V_vs)`.
On a reason: default safety does `restore, preserve`, clears `fcnoes`, disarms the
screen, and re-issues that margins call on the full data (`fcrerun`), so the rest of
the section runs uncollapsed; `safety(off)` prints the caution once and continues. The
outcome is exposed as `e(wdid_fcguard)` in {empty, `fallback`, `caution`}. A fallback
in the observed run also sends `bootstrap()` replications to the full data (the reps run
under `fastnose`, where the screen has no `r(V)` to read). The section-3 rowname bug was
fixed as recommended.

**Verification.** `test/test_fastcollapse.do` F8-F10 (`FASTCOLLAPSE_ALL_PASS`, 1:54):
F8 the `subgroups()+fe()+jointtests` command that hard-errored `r(503)` collapsed
(section 3) now completes collapsed under `safety(off)` with the caution; F9 the
degenerate `survey_rcs` command (cell shock removed, `fe(i.group##i.t)`, singular in the
subgroup-contrast direction) falls back under default safety with every subgroup matrix
read after the fallback bit-identical to `fullprecisionslow` (`mreldif` exactly 0), stays
collapsed under `safety(off)` with `e(wdid_fcguard)=caution` and `rc 0`, and two healthy
subgroup runs collapse with the guard empty; F10 pins the three probe-driven fixes
below. `test/check_validate_sync.py` OK. Suite pin `fcol_sg` moved from `r(198)` to
`rc 0`, verified on its `adm_sub` subsample (guard empty, `fcol 1`). Timing on the
shipped build: the cohort-level admin subgroup section reads 3-5x faster collapsed
(2.1 s vs 9.7 s for the reads), guard empty - the prize of section 4 survives the screen.

**The probe.** Seven angles, one subagent and one Stata each, 420 distinct subgroup
configurations run paired (arm C: default, collapse engaged; arm U: `fullprecisionslow`)
through the suite's `mkrun` with `%21x` dumps of every `e()` matrix, scalar and macro,
and a verdict per config: PASS (guard empty, C agrees with U to residue), MISS (guard
empty, some cell moves `> 1e-9` relative or the structure differs), FLAG_SEFAIL (flagged
on a run wooldid itself labels SE-degenerate), FLAG_HEALTHYVCE (flagged with no wooldid
caution - false alarm or a catch wooldid's caution misses), FALLBACK_DIV (fallback taken
but C still differs), CAUTION_RC (`safety(off)` changes the return code). Angles: A1
scale heterogeneity and heteroskedastic-but-healthy VCEs; A2 tolerance-boundary hunting
(contrast variances placed at `1e-7`..`1e-11` of `dmax` by weights and rescaling); A3 the
`o.` / non-estimable path; A4 bootstrap, permutation and `if` interactions; A5
`cleanmatrices`, `survey_rcs` cluster variants; A6 subgroup structures (subgroup = i,
finer than i, straddling units, single-unit levels, one-level subgroup, missing subgroup
values); A7 the eligibility gate and the reference-only pre window. Raw tallies over
the real build: PASS 235, FLAG_SEFAIL 85, FLAG_HEALTHYVCE 32, MISS 5, FALLBACK_DIV 4
(all four fixed below), CAUTION_RC 0; the remainder are BOTH_ERR (refusals in both arms),
NOTENGAGED (the eligibility gate, as designed, plus the `cleanmatrices` harness artefact
below) and NODUMP (harness). The `dg_*` rows in `a3_5`/`a3_6` and the `a6_d2` rows are
runs of a diagnostic build with the screen disabled, made to read the raw variance
ratios; their MISS verdicts are the unguarded baseline, and the same configurations
flag on the real build (`a3_9`, `a6_5`, `a6_6`).

**Fixed from the probe (each pinned in F10):**

1. *Singular direction anchored to the level scale* (A2 `b_q1e3`, `b_q5e4`:
   `binary_balanced`, `[pw]` with two units at `1 + 1e-3` and `1 + 2e-3`,
   `subgroups(sg_three) cluster(cohort)`). The contrast VCE was rank one with its
   largest eigenvalue itself `~1e-8` of `dmax`, so the contrast-relative eigenvalue test
   passed while the collapsed joint p came out `1.5e-82` against `3.3e-83` on the full
   data. The smallest-eigenvalue floor is now `1e-9 x max(dmax, largest eigenvalue)`.
   (F10b: fallback, `e(joint_att_allsgssame)` and the pre twin equal to full data.)
2. *Bootstrap replications follow the fallback* (A4 `a4b_rcs_std5`, `a4b_rcs_if_std5`:
   `survey_rcs`, `cluster(state) bootstrap(standard) reps(5)`). The observed run fell
   back, the reps did not (nothing for the screen to read under `fastnose`), and the
   bootstrapped subgroup-equality test came back "unavailable" against `p = 0` on the
   full data. The shell now appends `fullprecisionslow` to the reps when
   `e(wdid_fcguard) == "fallback"`. (F10c: fallback, bootstrap `e(V)` and `pwsgs` equal.)
3. *Reference-only pre window* (A7: an early cohort used as a subgroup so that its
   `esprelength(2)` window contains only the reference period, `jointtests`). The
   overjoint call has a single `over(_t)` level with zero variance by construction,
   nothing is tested (`r(p)` absent, `joint_horiz_s2` empty in both arms), and the screen
   fell back for nothing - a realistic data situation. The overjoint sites now exit when
   `r(p)` is absent: nothing reported, nothing to protect. (F10a: guard empty,
   `joint_horiz_s1` equal, `joint_horiz_s2` empty in both arms.)
4. *Dead `o.` predicate* (A3). The marker lives in the `r(b)`/`r(table)` column names,
   never in `r(V)`'s (0 of 148 screened calls carried it there); the read now uses
   `st_matrixcolstripe("r(b)")`. In practice the all-zero-VCV branch fires first on
   these runs (`a3_9`: fallback / caution as expected), so this closes a gap rather
   than a measured miss.

**Document-only.** (a) *Scale heterogeneity* (A1): the read-2 level test is relative to
the largest level variance across subgroups, so a subgroup whose outcome scale is below
about `3.75e-5` of another's (variance ratio `1e-9`) is flagged on a healthy VCE; every
such run was bit-identical between arms, so the cost is the lost speedup under default
safety and an overstated caution under `safety(off)`. Not worth a per-subgroup
normalisation; the help file names the route and the caution text was reworded to say
what wooldid treats the finding as rather than to assert a failed VCE. (b) A6
`a6_sgfine_rob_off`: a subgroup strictly finer than `i` (only reachable under
`safety(off)`; default safety refuses `r(459)`) has structurally duplicated margins
within each unit, so the contrast VCE is singular by construction and the caution
prints although nothing wrong is reported; bit-identical. Low priority.

**Residue MISSes (5, immaterial).** A1 `sc_rcs_g1_1em4_rob` (`e(joint_horiz_s1)` moves
`2.3e-9` relative on `p ~ 1.2e-6`), `sb_bb_shared_1e6_conv` (`pwsgs[2,1]` moves `5.1e-9`
on a t of `1e-6`); A2 `r_nm5`, `r_n5e6`, `r_n2e5` (between-subgroup contrast t-statistics
moving `6.7e-8`, `9.55e-8`, `8.95e-9` relative, on contrasts whose variance sits within
a factor of ten of the `1e-9` band). No point estimate, standard error, reported set or
return code moved in any of them; these are the amplified last-place residue the
default-on decision already accepted, at the edge of the band by construction.

**Justified catches.** FLAG_HEALTHYVCE verdicts that were not the scale artefact were all
runs where wooldid's own identical-SE caution stays silent but the contrast variance is
exactly zero on the full data too: robust VCE on saturated designs (`dta_s_deeppre`
adversarial, `dta_tiny jointtests`), `cluster(t)` without `fe` on `survey_rcs`, few
coarse clusters, weight-driven contrast degeneracy (`d_pm4`..`d_pm6`), and rescaled
flaggers (`k_rcs_mul8`) that defeat the caution's absolute tolerance. The dominant
FLAG_SEFAIL route is a cluster variable that nests the cohort variable, which makes
between-subgroup contrasts of identical (cohort, period) composition exactly
zero-variance.

**Harness artefact, no change.** A5/A6 `cleanmatrices` after a subgroup run returned
`r(3499)` (`wooldid_fcdegx() not found`) only when the ado is loaded with `do`; from the
adopath Stata re-reads the ado when a Mata function is missing and all seven variants
run `rc 0` (`cm_adopath.log`). Same mechanism as the SI-17/SI-18 note at ado ~5559.

**Off-angle observations for the author (both arms, not caused by the collapse):**
wooldid's identical-SE scan is cluster-only and uses an absolute tolerance (three SEs
identical to `1e-12`), so it misses robust singular directions and is defeated by
rescaling; under `fastnose` the SI-7 horizon drop can never fire (no `r(V)` to read),
so bootstrap and permutation replications cannot see a non-estimable subgroup horizon;
`pwsgs` treats a zero-variance contrast inconsistently across estimands (`att` row se
missing, `att_it` row se `1.6e-10`); a one-level subgroup covering all treated units
returned `att 1.385` against the main `1.406` (unverified beyond the observation);
reghdfe `vce(robust)` on a saturated `fe(sgt)` design returns an all-zero VCE when one
subgroup is rescaled by `<= 1e-4`; an explicit `fastcollapse` that falls back is silent
without `verbose`.

**Revised certification.** The `r(503)` regression completes (F8), the measured
degenerate shape reproduces the full-data subgroup output byte for byte after the
fallback and completes with the caution under `safety(off)` (F9), healthy subgroup runs
engage with no fallback (F9, `fcol_sg`), the three probe-driven refinements hold (F10), and
across 420 probed subgroup configurations no result read while the collapse was engaged
differed from its full-data value by more than amplified rounding. The full assessment
run should add the battery's subgroup rows under default safety (expect
`e(wdid_fcguard)` empty on every row wooldid does not already caution) and re-time the
slow subgroup rows against the pre-ruling numbers.

**Artifacts (session scratch, not committed).** `adv2/a*_*/RESULT.txt` and dumps for
every probe run, `adv2/runpair.sh` (paired runner and differ), `probe/a{1,2,4,6,7}/REPORT.md`
(A3 and A5 reported through their RESULT files), `diag/wooldid_diag.ado` (the
screen-disabled diagnostic build, `FCDIAG|` lines), `fix2/verify_fixes.do`, `fix2/v5b.do`
and `fix2/pin_fcol_sg.do` with logs, `cm/cm_adopath.do` and log.
