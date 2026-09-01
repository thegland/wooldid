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
