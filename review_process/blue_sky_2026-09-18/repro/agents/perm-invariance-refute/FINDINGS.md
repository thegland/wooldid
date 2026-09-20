# SKEPTIC REPORT — PI-1 (phantom cluster in the permute() design)

Agent: skeptic on "permute() metamorphic invariance at fixed seed".
All evidence in this directory. Every number below is from a do-file written from
PI-1's CLAIM, on MY OWN data; PI-1's own script was read but never run.

---

## VERDICT ON PI-1 — **PARTIAL**, severity revised **medium -> low**

### What reproduces, and harder than PI-1 showed it
`e1_exact.do` / `e1_exact.log`. 8 states x 3 counties x 8 years, state 1 = 2004,
state 2 = 2006, rest never, state 8 = the phantom. Options
`att cluster(state) permute(shuffle) permexact`. **permexact, so there is NO RNG
at all** — PI-1's contrast was 30 sampled replications on datasets with 10 vs 9
clusters, whose shuffle streams necessarily differ, so a p gap at 30 draws was
not by itself evidence that the reference DISTRIBUTION moved. Under exact
enumeration it is:

| arm                                  | e(N) | N_clust | npossible | exact p   | 95% shift CI      |
|--------------------------------------|------|---------|-----------|-----------|-------------------|
| A phantom, y missing on all its rows | 168  | 7       |    56     | .07142857 | [-.5312, .8277]   |
| B phantom rows dropped               | 168  | 7       |    42     | .04761905 | [-.5130, .5876]   |
| C phantom cluster value -> missing   | 168  | 7       |    42     | .04761905 | [-.5130, .5876]   |
| D `if state != 8`                    | 168  | 7       |    42     | .04761905 | [-.5130, .5876]   |
| E phantom killed by a missing CONTROL| 168  | 7       |    56     | .07142857 | [-.4588, .9437]   |
| E' same, phantom rows dropped        | 168  | 7       |    42     | .04761905 | [-.4571, .5168]   |
| J phantom is the SOLE 2008 adopter   | 168  | 7       |   336     | .07142857 | [-.5312, .8277]   |
| J' same, phantom rows dropped        | 168  | 7       |    42     | .04761905 | [-.5130, .5876]   |

`att` is bit-identical across every arm (`+1.a7c7591ced275X-001`, and
`+1.e331bfc86f07cX-001` in the controls arm). B == C == D bit-for-bit.
The p **crosses .05** (.0714 vs .0476). PI-1's numbers are right, and arm E is a
new instance PI-1 did not have: the killing variable does not have to be the
outcome, a missing CONTROL does it too.

### Why it is nevertheless not the defect PI-1 describes
**wooldid.sthlp lines 430-431, inside the `permute()` entry, specifies exactly
this behaviour, and names missing outcomes as its example:**

> "Requires cluster() with a single variable - the randomization unit (**the
> randomization design consists of the clusters having at least one row where the
> *if* expression evaluates to exactly 1 with a nonmissing cluster value; rows
> *wooldid* later removes by its deterministic markout, e.g. missing outcomes or
> zero weights, still count toward the design and its possible-assignment count**)"

That single sentence supplies, verbatim: the design predicate PI-1 reports
(`(if)==1 & !missing(cluster)`), the asymmetry between the three exclusion routes
(only those two remove a cluster), the rule that markout-dead rows stay in, the
outcome as the named example, and the fact that the **possible-assignment count**
is what moves. PI-1's EXPECTED rests on "sthlp:61 states only that the `if` rule
applies identically ... nothing says a missing outcome is treated differently
from a missing cluster id". That is false: PI-1 cited sthlp:61 and missed
sthlp:430. The footer also prints the design size at run time ("56 of 56
replications completed; 56 possible assignments"), which a reader can compare
against e(N_clust).

So there is no doc-vs-code divergence (brief value #5), no silently wrong number
(#1), no failed inference claim relative to what is documented (#3). What is left
is a methodological disagreement with a deliberate, documented design choice:
PI-1 thinks the randomization population should be the estimable design, wooldid
says it is the assignable one. That is a legitimate thing to put to the author —
a unit that could have been assigned treatment is arguably in the design whether
or not its outcome was collected — but it is a design question, not a defect, and
"medium" overstates it. **Severity low.** Honest concession: arm J is the sharpest
form and PI-1 is right to single it out — the observed statistic is computed on a
sample with no 2008 cohort while 300 of the 336 reference statistics come from
samples that have one. It is still exactly what sthlp:430-431 specifies.
The mechanism PI-1 states is correct; only its novelty and severity are wrong.

---

## NEW FINDINGS (mine, not seeded)

### SK-1 (low, UX) — permute()'s nesting refusal is diagnosed on rows absent from e(sample)
`e4b_nesting.do` / `.log` (and `e4_nesting.do`, whose V1a line is spoiled by an
unrelated `dots()` parse error; its V1d carries the comparison correctly).

8 states x 3 counties. State 3 adopts in 2006 but its county 9 is a merge miss
kept in the panel skeleton: `ttre == .` AND `y == .` on every one of its 8 rows.
Onset is constant within every `i`, so wooldid's general within-`i` screen
(ado:772) does not fire.
 - V1a without `permute()`: runs clean, e(N)=184, att `+1.8da894f951c6eX-001`
 - V1c drop county 9's rows: e(N)=184, `permute(shuffle) permexact` runs, npossible=168
 - V1d keep vs drop, same run: `identical=1`, N 184 both
 - V1b keep them and add `permute(shuffle) permexact`: **r(459)** — "permute()
   requires clusters to be subsets of treatment-onset groups ..., but 24
   in-sample observation(s) sit in clusters that mix onset values or mix treated
   and never-treated rows."  V2b: same refusal when the dead county carries a
   stale onset (2010) instead of a missing one.

**Honest reading, after checking the help against my own result:** the refusal
itself is the conjunction of two documented sentences — sthlp:430 puts
markout-dead rows in the design, sthlp:431 requires every design cluster to be
onset-homogeneous — and a never-treated (or 2010) county inside a 2006-adopting
state genuinely IS an onset coding inconsistent with cluster-level assignment.
That the estimation path never sees those rows does not make the coding right.
So this is **not** an undocumented consequence and **not** a wrong number.

What survives is the diagnostic, and it is a UX defect (brief value #6), **low**:
the message counts 24 observations when 16 of them are perfectly consistent live
rows and only 8 are responsible; the 8 cannot be found through `e(sample)`,
`e(N)` or `e(N_clust)`, all of which are identical to the version that runs; and
the message names no remedy for them (its remedies are about the onset coding,
not about rows the estimator already discarded). A user who runs the same
specification without `permute()`, sees it work, adds `permute()` and is told 24
observations mix onsets has nothing in the output to point at county 9.

### SK-2 — bootstrap()'s few-cluster and few-TREATED-cluster cautions are suppressed by a cluster with no estimable row
`e6_bootscreens.do` / `.log`, `e6b_controls.do` / `.log`.

`wooldid_validate` counts clusters over `smp2 = ((if)==1) & !missing(cluster)`
(ado:6461) — the same phantom-inclusive population — and NOT over the estimation
sample. Measured, with the estimation sample held bit-identical across each pair:

| arm | data | e(N) | e(N_clust) | caution printed? |
|-----|------|------|-----------|------------------|
| A1  | 9 real clusters + 1 phantom never-treated | 216 | 9 | **no** "only 9 clusters are available for the bootstrap" |
| A2  | the phantom's rows dropped                | 216 | 9 | yes, "only 9 clusters are available for the bootstrap" |
| B1  | 5 real treated + 1 phantom treated (y missing)   | 264 | 11 | **no** few-treated caution |
| B2  | the phantom's rows dropped                       | 264 | 11 | yes, "only 5 treated cluster(s) are available for resampling" |
| D1  | 5 real treated + 1 phantom treated (CONTROL missing) | 264 | 11 | **no** few-treated caution |
| D2  | the phantom's rows dropped                           | 264 | 11 | yes, "only 5 treated cluster(s) are available for resampling" |

B/D is the damaging half: the few-treated-cluster caution is the one warning that
tells a user the pairs cluster bootstrap is badly sized in their design and points
them at `permute()` or `boottest` (sthlp:532, MacKinnon and Webb 2017/2018). A
user with exactly five estimable treated clusters and one dead one is told
nothing, and `e(N_clust)` = 11 gives no hint either. It bites only when the count
sits on the threshold (10 total / 6 treated), which is why this is medium and not
high — but a threshold caution that a phantom can push over is exactly the case
the threshold exists for. sthlp:430-431 documents the phantom-inclusive
population for **permute()**'s randomization design and its possible-assignment
count; nothing extends it to bootstrap()'s screens.

Measured and NOT a finding (reported for honesty): the same predicate governs
`bootstrap() requires at least 2 clusters in the estimation sample`, whose own
message names the estimation sample while the predicate does not. With 1 real +
1 phantom cluster that screen IS bypassed and both cautions print the inflated
count ("only 2 clusters are available"), but a downstream guard using the right
population catches it: `Error: cluster(state) leaves 1 cluster(s) in the
estimation sample` r(459) (arm C1). Message wording only; no live consequence.

### SK-3 — NEGATIVE RESULT: bootstrap() does NOT inherit PI-1's defect
`e5_boot.do` / `.log`. PI-1 extends its finding to `bootstrap()` from code
reading alone ("wooldid_bootrun uses the identical predicate (ado:8912-8923), so
bootstrap() shares the design population"). The predicate is indeed shared, but
the consequence is not. 9 clusters (4 treated) + 1 phantom vs the phantom's rows
dropped, `att cluster(state) reps(400)`, three seeds each:

| mode                 | mean bootstrap SE, phantom kept | phantom dropped | ratio |
|----------------------|---------------------------------|-----------------|-------|
| `bootstrap(bayes)`   | 0.072886                        | 0.072812        | 1.0010|
| `bootstrap(standard)`| 0.089842                        | 0.090492        | 0.9928|

Both inside Monte-Carlo noise (about 3.5% per run at reps(400), about 2% on a
three-run mean). The reason is structural, not luck: `bootstrap(bayes)` draws
Dirichlet cluster weights and the normalized sub-vector over the real clusters is
again Dirichlet over their own shapes, while the `att` leg of
`bootstrap(standard)` is cluster frequency reweighting — and aweights are
scale-free, so in both modes only the DIRECTION of the weight vector over the
clusters that actually carry rows enters the estimate, and that direction has the
same distribution with or without the phantom. Only the RNG stream moves. So
**Physical-resample legs, measured too** (`e7_physleg.do` / `.log`), because the
scale-invariance argument does not cover them: `bootstrap(standard)` with the
cell-weighted `att_i`, which the help says physically resamples and
re-identifies, gives mean SE 0.082352 with the phantom and 0.082811 without
(ratio 0.9945, three seeds each); `boottstat`, which takes the physical leg for a
different reason, gives p = 0.050 / 0.040 with the phantom and 0.040 / 0.050
without at reps(200), with A-arm intervals [.245, 1.182] and [.310, 1.217]
against B-arm [.378, 1.123] and [.278, 1.211] — the seed-to-seed spread WITHIN an
arm is as large as the spread between arms, so no structural movement is
detectable there either (reps(200) is coarse; this arm is weaker evidence than
the att_i one). So PI-1's bootstrap extension should be struck: `bootstrap()` shares the predicate
and shares SK-2's screen problem, but it does NOT report a different standard
error, p-value or interval because of a phantom cluster.

---

## WHAT I ATTACKED AND FOUND CLEAN
- The three exclusion routes really are 1-vs-2 and B/C/D agree bit-for-bit, with
  no hidden fourth behaviour (E1 arms B, C, D).
- The point estimate and `e(N)`/`e(N_clust)` never move with a phantom present,
  in eight arms across three panels — so there is no silent estimation-side
  contamination, only the design-side one the help describes.
- `e(perm_failed)` and `e(perm_blank)` stay 0 in every arm, including arm J where
  the reference set mixes two cohort structures: the replications are not
  breaking, they are answering a different question, exactly as PI-1 said.
- `e(perm_npossible)` matches the exact combinatorial count for the documented
  population in every arm (56 = 8!/(1!1!6!), 42 = 7!/(1!1!5!), 336 = 8!/5!,
  168 = 8!/(1!2!5!)), so the count is not itself wrong.
- bootstrap SEs and intervals are invariant to a phantom cluster on the
  observation-weighted `att` leg of both modes (SK-3), and on the
  physical-resample legs measured in `e7_physleg.do` (see below).
