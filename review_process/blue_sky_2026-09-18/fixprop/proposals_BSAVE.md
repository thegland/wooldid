# Fix proposals — mechanism `BSAVE` (bootsave()/permsave())

Three findings, **two fix sites**, three proposals (one of them a rider the author can rule
on separately).

| finding | proposal | moves a reported number? |
|---|---|---|
| `BT-04`, `hvr-bootsave-jointwobs-01` (one defect, two agents) | **BSAVE-1** | no |
| `BR-4` — resampling universe | **BSAVE-2** | **yes**, on three named regimes (bit-identical everywhere else) |
| `BR-4` — the two silenced cautions | **BSAVE-3** (rider on BSAVE-2) | no number; changes *where* two cautions print |

Everything below was **built and measured**. My patched copy of the ado is at
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-BSAVE/patched/wooldid.ado`
(13 edits; unified diff against the pinned build at
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-BSAVE/patch.diff`).
The repo was not touched.

---

## PROPOSAL BSAVE-1 — keep the observed joint/adversarial statistics and let `bootsave()` write them

**Covers:** `BT-04` and `hvr-bootsave-jointwobs-01` (the same defect; the second also
establishes the `advW` twin).

### Root cause, restated precisely

`wooldid_bootsavemk` (ado:14918-14933) writes the `_rep = 0` row of each `jointW_<family>`
column as `JO[1,k]`, where `JO = st_matrix(jobs)` and `k` indexes the columns of
`JW = e(boot_jointW)`. `JW` has **one column per test FAMILY** (`att`, `att_it`, `ame`,
`ame_it` and their `_s#` twins). `jobs` is `wdbjobs`, built by `wooldid_bootjobs`
(ado:14666-14695) with **one column per (family, PRE-HORIZON) pair**, holding
`e(b_es_<fam>)[1,j]`. The two column sets have nothing to do with one another, so row 0
gets an observed event-study **coefficient**. Same for `advW` with `aobs = wdbaobs`.

The fix is not in the Mata writer. The value the file needs already exists: `wooldid_bootjnt`
computes exactly `W_obs` per family at ado:6903 and ado:6929 (`wooldid_bootjwobs`) and at
ado:6980 (`wooldid_bootadvobs`) — and **throws it away into a local**. The `fixedmetric`
leg of the same program already keeps its twin and posts `e(boot_jointWobs)` (ado:6872);
the `studentized` leg does not. So `boottstat` is the only engine where wooldid knows
`W_obs` and stores it nowhere.

### The exact edit (5 hunks, `wooldid_bootjnt` + the shell)

**(a)** ado:6892, immediately after `if wdb_jtok local wdb_jfk : colnames jwacc`, insert:

```stata
    tempname jwobs
    if `wdb_jtok' {
      matrix `jwobs' = J(1, colsof(`jwacc'), .)
      matrix colnames `jwobs' = `wdb_jfk'
    }
```

**(b)** ado:6903-6904 and ado:6929-6930 — at **both** `wooldid_bootjwobs` call sites, as the
first statement inside the `if wobs < .` block:

```stata
            matrix `jwobs'[1, `jk'] = `wobs'
```

(the local `jk` is already `: list posof "`fk'" in wdb_jfk`, i.e. the column index of `jwacc`,
so the alignment with `e(boot_jointW)` is by construction, not by luck.)

**(c)** ado:6946, beside the existing `cap ereturn matrix boot_jointW = jwacc`:

```stata
      cap ereturn matrix boot_jointWobs = `jwobs'
```

**(d)** the adversarial twin, three lines in the same shape: after ado:6962
(`if wdb_advok local wdb_acols : colnames awacc`) add a `J(1, colsof(awacc), .)` scratch
matrix named `awobs` with colnames `wdb_acols`; at ado:6980-6981 add
`matrix awobs[1, ak] = wobs` as the first statement inside `if wobs < .`; at ado:7007 add
`cap ereturn matrix boot_advWobs = awobs` beside the `boot_advW` posting.
(`e(perm_advWobs)` at ado:8639 is the naming precedent.)

**(e)** shell ado:437, the `wooldid_bootsave` call — two argument swaps:

```
- jobs(`=cond("`boottstat'" != "", "`wdbjobs'",  "")')
+ jobs(`=cond("`boottstat'" != "", "e(boot_jointWobs)", "")')
- aobs(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "`wdbaobs'",  "")')
+ aobs(`=cond("`boottstat'" != "" & "`adversarial'`altmethod_adversarial'" != "", "e(boot_advWobs)", "")')
```

Line 394's `jobsname()`/`aobsname()` stay exactly as they are: `wooldid_bootrun` and
`wooldid_bootjnt` still need the observed **coefficient** vectors, which is what those two
arguments are for.

`wooldid_bootsavemk` is **not** touched. Its `cols(JO) >= k ? JO[1,k] : .` guard now does
the right thing: a family whose `W_obs` could not be computed leaves a missing in row 0
rather than a coefficient.

### Help edit that must ship with it

`sthlp:1340` currently reads "... with the replication statistics in `e(boot_jointW)` /
`e(boot_sgW)` and the observed ones in `e(boot_jointWobs)` / `e(boot_sgWobs)` **(the latter
two without boottstat)**". After the fix `e(boot_jointWobs)` exists under `boottstat` too,
so the parenthetical becomes "(`e(boot_sgWobs)` without `boottstat`)", and `e(boot_advWobs)`
joins the `boottstat` list beside `e(boot_advW)`. `sthlp:557` (the `bootsave()` entry) and
the writer's comment at ado:14905-14906 need **no** change — they already describe the fixed
behaviour; today they are simply false.

### Runtime cost

**Zero measurable.** Two `1 x K` matrices (`K` = number of test families, typically 1-8, at
most 4 x (1 + #subgroups)) and one `matrix[1,jk] =` store per family. No pass over the data,
nothing that scales with N or with `reps()`. `colsof()` on a matrix already in memory.

### Does it move a reported number?

**No.** The p-values are computed from the local `wobs` and are untouched; I only keep a
copy. Measured: on the boottstat fixture below, pinned and patched give
`e(joint_horiz)` = `+1.2762762762762X-001` — **bit-identical** — and every replication
statistic in `e(boot_jointW)` is unchanged (`repW min/max` identical in both logs).
What changes is only the content of two columns' `_rep = 0` cell in the saved `.dta`, plus
two new `e()` matrices.

### Blast radius

- `bootsave()` files from `boottstat` runs: the `jointW_*` and `advW` row-0 cells change from
  a coefficient to the statistic. Anyone who built a workflow on the old (wrong) value would
  see it change — but that value was never usable, and it was negative, which a quadratic
  form cannot be.
- Two new `e()` matrices under `boottstat`. `e(boot_jointWobs)` already exists on the other
  leg, so no consumer keys on its **absence** as a `boottstat` marker (I grepped: the only
  references are ado:6872, ado:8444, ado:8639 — three postings, no reads; `e(boot_jt_family)`
  is what the display keys on, ado:11797).
- Nothing else reads `wdbjobs`/`wdbaobs` at the bootsave call site.

### Alternatives rejected

1. **Fix it inside `wooldid_bootsavemk` by name-matching `JO`'s colstripe to `JW`'s.**
   Impossible: no correct value exists anywhere in `wdbjobs`. It is a coefficient vector by
   design, and `wooldid_bootjwobs` is what turns it into `W_obs`.
2. **Write `.` into row 0 and amend the help to say the observed statistic is not in the
   file.** Cheapest possible fix, removes the silently-wrong number, and I would accept it —
   but it retires a documented capability ("so the draws can be recentred from the file
   alone") when the value is already computed four lines away. Keep as the fallback if the
   author does not want two new `e()` matrices.
3. **Recompute `W_obs` inside `wooldid_bootsavemk`.** Would duplicate the whole
   `wooldid_bootjwobs` guard chain (rank, positive-definiteness, name mapping) in a second
   place, which is exactly how the two column sets drifted apart in the first place.

### Test that proves it (no oracle needed)

`probes/prop-BSAVE/t_patched.do`, ARM E — an **internal-consistency invariant** that the
program's own output must satisfy:

```
p_recomputed = (1 + #{r>0 : jointW_att[r] >= jointW_att[_rep==0]}) / (#{r>0 : jointW_att[r] < .} + 1)
p_recomputed == e(joint_horiz)      exactly, in %21x
jointW_*[_rep==0] >= 0              (a Wald statistic cannot be negative)
```

Measured, 60 units x 11 years, `cluster(unit) jointtests esprelength(2) espostlength(2)
bootstrap(standard) reps(25) boottstat bootsave(bj)`:

| | pinned | patched |
|---|---|---|
| `jointW_att[_rep==0]` | `-1.040992fa3653dX-004` (= `e(b_es_att)[pre2]`, **negative**) | `+1.c65d1e611887aX+000` |
| p recomputed from the file | `1.000000000` | `+1.2762762762762X-001` |
| `e(joint_horiz)` reported | `+1.2762762762762X-001` | `+1.2762762762762X-001` |
| MATCH | **0** | **1** |

and the `advW` twin (`t2_*.do` leg L4, `adversarial esfixedbaseperiod`, 60 units, 15% missing y):

| | pinned | patched |
|---|---|---|
| `advW[_rep==0]` | `-1.38a70d07e21b1X-003` (**negative**; = an observed pre-horizon coefficient) | a positive Wald statistic |
| adversarial p recomputed from the file | `1.000000000` | equals the reported p |
| `e(adversarial_pval_att)` reported | `+1.3b13b13b13b14X-005` | `+1.3b13b13b13b14X-005` (**unchanged**) |
| MATCH | **0** | **1** |

The same MATCH=1 result is reproduced independently on the `boottstat` + `cluster()` arm of
`t3_patched.do` (`jointW_att[_rep==0]` = `+1.70fc6a4b208b7X+004`, p recomputed = p reported
= 0.1, on **both** arms of that fixture), so it is not a property of one dataset.

**The per-subgroup columns were tested too** (`t9_sg.do` / `o9_sg.txt`, `subgroups(sg)` +
`boottstat` + `bootsave`), because that is where the positional alignment could most plausibly
slip. `e(boot_jointWobs)` comes back `1 x 3` with colnames `att att_s0 att_s1` and values
3.2220796 / 6.3993803 / 2.1790015 — all positive — and every one of the three recomputes
exactly:

| column | p recomputed from the file | p reported | MATCH |
|---|---|---|---|
| `jointW_att` | `+1.4000000000000X-002` | `e(joint_horiz)` `+1.4000000000000X-002` | **1** |
| `jointW_att_s0` | `+1.8000000000000X-002` | `e(joint_horiz_s0)` `+1.8000000000000X-002` | **1** |
| `jointW_att_s1` | `+1.e000000000000X-001` | `e(joint_horiz_s1)` `+1.e000000000000X-001` | **1** |

Logs: `probes/prop-BSAVE/t_pinned.log`, `t_patched.log`, `o2_pinned_run1.txt`, `o2_patched_run1.txt`, `o3_pinned.txt`, `o3_patched.txt`.

---

## PROPOSAL BSAVE-2 — the resampling universe is `e(sample)`, not the user-`if` screen

**Covers:** `BR-4`, the draws-and-degrees-of-freedom half.

### Root cause

`wooldid_bootrun` builds its universe at ado:8909-8923:

```stata
  * sample concept (A7): user-if sample + nonmissing cluster (permrun byte-parallel)
  tempvar smp
  if `"`ifr'"' != "" { qui gen byte `smp' = ((`ifr') == 1) ... }
  else qui gen byte `smp' = 1
  ...
  qui replace `smp' = 0 if missing(`clusv')
```

Everything downstream keys on the `smp` marker: the map frame (`frame put ... if smp`, ado:8939)
and hence `wdb_C` -> `wdb_r_nclus` -> `bootdf(clusters)` -> `e(boot_df)` and `e(df_r)`; the
`w_<cluster>` columns of `bootsave()`; the `wdb_nest` predicate; and the two
`qui keep if smp` statements on the physical-resample legs (ado:9247, ado:9344). Rows and whole
clusters that the **estimation** dropped — a unit whose outcome is never observed, a unit
missing a control, a singleton absorbed group, a period outside the ES window — stay in it.

### The exact edit (one line + its comment)

ado:8923, immediately after the existing `qui replace smp = 0 if missing(clusv)`:

**And rewrite the now-false comment at ado:8909**, which currently reads
`* sample concept (A7): user-if sample + nonmissing cluster (permrun byte-parallel)` — after
this edit the two engines are no longer byte-parallel and the line contradicts the code four
lines below it:

```stata
  * sample concept: the REALIZED ESTIMATION SAMPLE e(sample), intersected with nonmissing
  * cluster. This used to be the user-if screen and was byte-parallel with wooldid_permrun's
  * (ado:7479-7486); since BR-4 it is not - permrun still builds the user-if screen
```

```stata
  * BR-4: the resampling universe is the REALIZED ESTIMATION SAMPLE, not the user-if
  * screen. e() here is still the observed run's post (the shell `est sto's it and does
  * not clear it; the rep-0 comparator call below is the first thing that overwrites e()),
  * so e(sample) is exactly the row marker the worker posted. Rows and whole clusters the
  * estimation dropped - missing outcome, missing control, a singleton absorbed group,
  * outside the ES window - must not be resampling units: counted, they inflate G and so
  * bootdf(clusters) = G-1 and e(boot_df); drawn, they make the realized multiset of
  * CONTRIBUTING clusters random, which breaks the replication-for-replication identity
  * with Stata's own bootstrap prefix that the help promises
  qui replace `smp' = 0 if e(sample) == 0
```

**Precondition verified before writing a line of this** (`probes/prop-BSAVE/p1_esample.log`):
`e(sample)` is a valid row marker **in the caller's frame** after a wooldid estimation,
despite the worker's inline `preserve` and `wooldid_fcol`'s `keep if`. On a 112-row file with
40 phantom rows: `e(N)` = 72, `count if e(sample)` = 72, `count if e(sample) & phantom` = 0,
`count if e(sample)==0 & !phantom` = 0. Holds on the default collapse path, under
`nofastcollapse`, with `controls()` dropping rows for a missing covariate (`e(N)` = 66,
`count if e(sample)` = 66, zero in-sample rows with a missing `y` or `x1`), and after a
bootstrap run. **Not measured under `poisson`**: `bootstrap(standard)` + `poisson` is allowed
and the worker reads its `e(sample)` from ppmlhdfe through the same `esamp` construction at
ado:1510-1511, so it is structurally covered, but I did not run that arm.

**If the author wants the dependency made explicit rather than read out of live `e()`** (a
defensible preference in this codebase), the variant is: at ado:366, beside
the existing `qui est sto wdbobs`, add a `tempvar wdbes` and `qui gen byte wdbes = e(sample)`,
pass a new `esvar(wdbes)` argument on the ado:394 call, and write
`qui replace smp = 0 if esvar != 1` in bootrun.
Same cost, larger diff, no reliance on `e()` still being live inside the subprogram. I tested
the one-line form; the pass-down form is mechanically equivalent.

### Runtime cost

**One `replace` pass over N rows**, once per run (not per replication) — the same order as
the `replace` on the line above it, which already runs unconditionally.

Measured, `probes/prop-BSAVE/t5_cost.do` / `o5_cost.txt`, **1,000,000 rows, 40,000 clusters**:

| operation | seconds |
|---|---|
| **the one added line** (`replace smp = 0 if e(sample) == 0`) | **0.041** |
| reference: one `gegen ... , by()` the program already runs on the same rows | 0.059 |

So the added guard is **cheaper than a single `gegen` wooldid already performs**, and it runs
once against a bootstrap that refits the whole model `reps()` times. It comfortably satisfies
the author's "no guards that are expensive on large datasets" constraint.

It is strictly work-*saving* downstream: the physical-resample legs' `keep if smp` now
drops rows earlier, so each replication's nested `wooldid` call sees fewer rows (it was going
to screen them out itself anyway).

### Does it move a reported number? **Yes — and exactly where it should**

Three regimes, all measured:

1. **Complete data (`smp` == `e(sample)`): bit-identical.** Measured pinned vs patched on the
   same seed: `att` + `att_i` under `bootstrap(standard)` (ARM C), `bootstrap(bayes)` (ARM D),
   and `boottstat jointtests` (ARM E) — every one of 20 `e(bootdist)` rows in `%21x`,
   `e(boot_df)`, `e(N_clust)` and `e(joint_horiz)` **identical**. Any blessed baseline built
   on data with no dropped cluster does not move.
2. **A whole cluster outside `e(sample)`: everything bootstrap moves.** `G` changes, so
   `bsample` draws a different multiset and every replication coefficient changes;
   `e(boot_df)` = `G`-1 and `e(df_r)` change; the Wald p/CI, the percentile p/CI (the
   *draws* moved, not just the reference distribution), the joint tests and the `w_*`
   column set of `bootsave()` all change. **Name `e(boot_accel)` separately in the ruling**:
   `bootdisplay(bca)`'s acceleration comes from a delete-one-**cluster** jackknife of G extra
   fits, so on the pinned build each cluster outside `e(sample)` contributes a degenerate
   delete-nothing fit to that jackknife. It moves by a different mechanism than the
   resampling families do, and `e(boot_ci_bca)` / `e(boot_p_bca)` move with it. (I ran the
   `bootdisplay(bca)` cell in the option matrix — rc = 0, identical on complete data — but
   did **not** measure its jackknife on a phantom-cluster fixture.) **This needs a re-bless on any baseline whose
   fixture has such a cluster, and it is the author's ruling.** Measured, 14 units of which
   5 have no usable outcome at all: `e(boot_df)` 13 -> 8 with `e(N)` = 72 and
   `e(N_clust)` = 9 unchanged; the Wald critical value moves from t(13,.975) = 2.160 to
   t(8,.975) = 2.306, **6.8% wider**.
3. **Partial row loss, no whole cluster lost** — measured on a 60-unit panel with 15% of `y`
   missing at random (`e(N)` = 404 of 480 rows; `e(boot_df)` = 59 on both builds):

   | leg | all 20 `e(bootdist)` rows, pinned vs patched |
   |---|---|
   | `bootstrap(standard)` + `cluster(unit)` | **bit-identical** |
   | `bootstrap(bayes) bootalpha(equal)` + `cluster(unit)` | **bit-identical** |
   | `bootstrap(bayes) bootalpha(size)` + `cluster(unit)` | **all 20 differ** |
   | `bootstrap(standard)`, **no `cluster()`** | **all 20 differ** |

   Two separate movers, each needing its own line in the ruling:
   - **`bootalpha(size)`**, because the Gamma shape and divisor `n_c` is the cluster's
     in-sample **row count** (ado:8941-8949) and that count changes. Its `bootsave()` `_rep = 0`
     row moves with it: `1/480` -> `1/404`. That direction is a **fix to a documented-vs-actual
     divergence** — `sthlp:557` says that cell holds "1/*N*, *N* the number of observations in
     **the estimation sample**", `e(N)` is 404, and the pinned build wrote `1/480`. (The ado's
     own comment at ado:8942 calls it the "in-sample ROW COUNT"; today the two readings differ,
     after the fix they agree.)
   - **The unclustered leg**, where the resampling unit is the observation itself, so every row
     with a missing outcome was a resampling unit. `e(boot_df)` does not move there (the
     unclustered default is `bootdf(normal)`, which posts none), but every draw does.

### Blast radius

Checked on an 8-cell option matrix (`t6_pinned.do` / `t6_patched.do`; `o6_pinned.txt` vs
`o6_patched.txt`): no `cluster()`; `subgroups()` + `jointtests` + `boottstat` + `bootsave()`;
`bootstrap(bayes)` + `subgroups()` + `att_i`; `adversarial` without `boottstat`; the
`fixedmetric` joint leg + `bootsave()`; a plain analytic run; and `bootdisplay(bca)`.
Every cell returns `rc` = 0 on the patched build and zero Stata errors in the log. Pinned vs
patched, seven of the eight are identical on `e(N)`, `e(boot_df)` and `e(joint_horiz)`; the
eighth is the `esfixedbaseperiod` arm above, which is the fix working. The only other
difference anywhere in the matrix is `e(boot_jointWobs)` newly present on the
`subgroups()` + `boottstat` cell — that is BSAVE-1.
`permute()` is **bit-identical** between the two builds (all 20 `e(permdist)` rows,
`t7_pinned.do` / `t7_patched.do`), as the diff says it must be: no permutation program is
touched.

`wooldid_bootrun` only. Every consumer reads the restricted universe through `smp`,
`wdb_C`, `wdb_nrows` or the map frame, so nothing needs a second edit. It does **not** touch
`permute()`; see the open question below. It cannot fire outside `bootstrap()`.

**What a re-bless has to cover**, in the order the measurements above put it:
any fixture with a whole cluster outside `e(sample)` (everything bootstrap moves); any
`bootalpha(size)` fixture with *any* dropped row; any **unclustered** `bootstrap()` fixture
with any dropped row. A `bootstrap(standard)` or `bootstrap(bayes) bootalpha(equal)` fixture
with `cluster()` and no whole-cluster loss is **bit-identical** and does not need re-blessing.

Bonus, free: the `wdb_C < 2` refusal at ado:8956-8959 currently passes on *1 real + 1 phantom*
cluster and the run proceeds to resample a universe with one contributing cluster. After the
fix it refuses, with the message it already has.

### Alternatives rejected

1. **Fix only the count** — use `e(N_clust)` for `bootdf(clusters)` and for the cautions, and
   leave the draws alone. Rejected: the refuter measured the *draws* moving (real-cluster mass
   varying 22-36 instead of fixed at 30), `sthlp:532` promises replication-for-replication
   identity with Stata's own prefix, and I measured that promise broken
   (MAXRELDIF = **1.0678e-01** on every one of 40 replications, `o2_pinned_run1.txt` leg L1).
   The draw fix is the same one line. Also `e(N_clust)` is a *minimum* under multiway
   clustering while `wdb_C` is the intersection — they are not interchangeable in general
   (multiway is refused under `bootstrap()`, but the identity would be a trap for the next
   editor).
2. **Refuse the run when `smp` has clusters outside `e(sample)`.** A unit with no usable
   outcome in the file is an ordinary administrative shape, not a user mistake. Refusing it
   would be a new gate on realistic data; the project's own T2/T3 note (ado:8960-8978) is
   explicit that it prefers a rule to an exclusion.
3. **Caution instead of fixing.** Leaves `e(boot_df)` wrong and the help's Stata-prefix claim
   false. A caution is what you write when the number is defensible; this one is not.

### Test that proves it — a **metamorphic identity that needs no oracle**

`probes/prop-BSAVE/t_pinned.do` / `t_patched.do`. Build one file; make ARM A the file with
`k` units whose outcome is missing on every row, ARM B the same file with those units
physically deleted. The two have the **same estimation sample by construction** (same
`e(N)`, same `e(N_clust)`, same point estimate). Under a common `rseed()` the bootstrap
must therefore be the same object.

14 units, 5 of which (3 treated, 2 never) have `y` missing on every row; `e(N)` = 72,
`e(N_clust)` = 9 in **both** arms.

| comparison | pinned | patched |
|---|---|---|
| ARM A vs ARM B, `e(boot_df)` | 13 vs 8 | **8 vs 8** |
| ARM A vs ARM B, all 20 `e(bootdist)` rows in `%21x` | **all 20 differ** | **bit-identical** |
| pinned ARM B vs patched ARM A | — | **bit-identical** |
| few-cluster cautions printed | ARM B only (1 each) | both arms (2 each) |
| complete-data arms C/D/E, pinned vs patched | — | **bit-identical** |

Repeated on the two code paths the first table does not exercise
(`probes/prop-BSAVE/t3_pinned.do` / `t3_patched.do`, same fixture):

| leg | comparison | pinned | patched |
|---|---|---|---|
| `att_i` — the physical re-identification leg (`wdb_reid`, `keep if smp` at ado:9247) | ARM A vs ARM B, `e(boot_df)` | 13 vs 8 | **8 vs 8** |
| | ARM A vs ARM B, 20 `e(bootdist)` rows | **differ** | **bit-identical** |
| | pinned ARM B vs patched ARM A | — | **bit-identical** |
| `boottstat` + `cluster()` + `jointtests` | **`e(joint_horiz)` reported** | **0.2 (A) vs 0.1 (B)** | **0.1 (A) = 0.1 (B)** |
| | ARM A vs ARM B, 20 `e(bootdist)` rows | **differ** | **bit-identical** |
| | `bootsave()` `w_*` columns written | 14 (A) vs 9 (B) | 9 (A) = 9 (B) |
| | `jointW_att[_rep==0]`, p recomputed == reported | MATCH **0** in both arms | MATCH **1** in both arms |

That `e(joint_horiz)` row is the concrete user harm of `BR-4`: **the reported joint
pre-trend p-value is twice as large** purely because five units with no usable outcome were
left in the data file. Nothing about the estimation sample differs between the arms.

**And it fires without any engineered missingness at all.** In the option-matrix smoke run
(`t6_pinned.do` / `t6_patched.do`, 40 units x 11 years, 8% of `y` missing at random — no unit
all-missing), the arm
`att cluster(unit) jointtests adversarial esfixedbaseperiod esprelength(2) espostlength(2)
bootstrap(standard) reps(15)` gives, on **identical** `e(N)` = 392:

| | pinned | patched |
|---|---|---|
| `e(boot_df)` | 39 (G = 40 — every unit in the file) | **37** (G = 38) |
| `e(joint_horiz)` | 0.5625 | **0.75** |

Diagnosed directly (`t8_s4why.do` / `t8_s4why.log`, pinned build, same spec without
`bootstrap()`): 40 units in the file, **40** units with at least one non-missing `y`, and
only **38** units with any row in `e(sample)`. The two that vanish are units 14 and 21, which
have **10 and 9 rows of perfectly good outcome data each** and are dropped wholesale by the
estimation's own screen — not by missing data, which is exactly why the `!missing(y)`
approximation rejected below would not have caught them. The pinned build resamples 40
clusters; the estimation used 38.

This is the ordinary case, not a contrived one: any option that narrows the event-study
window can strand a cohort, and the resampling universe never hears about it.

Second, independent test — **cross-engine**, which is what `sthlp:532` actually claims.
`bootstrap attb = _b[main_att], reps(40) seed(S) cluster(unit) idcluster(newc): wooldid ...
fastnose` against `wooldid ..., bootstrap(standard) reps(40) rseed(S) bootsave(w1)`, merged
on `_rep`, on a 60-unit panel with 15% missing `y` plus one all-missing unit:

| | pinned | patched |
|---|---|---|
| max relative difference over 40 replications | **1.0678e-01** | **2.6110e-08** |
| `e(boot_df)` | 59 | 58 |

2.6110e-08 **is** the floor — it is the `float` precision of Stata's own `saving()` file, i.e.
an exact match. The pinned build differs from the prefix on every one of the 40 replications;
the patched build matches it on all 40. That is `sthlp:532`'s promise, restored.

---

## PROPOSAL BSAVE-3 (rider) — move the two few-cluster cautions to where the count is real

**Covers:** `BR-4`, the silenced-cautions half. Separable from BSAVE-2; rule on it separately.

### Why a separate proposal

The two cautions live in `wooldid_bootctl` (ado:6479, ado:6498), which runs inside
`wooldid_validate` at **parse time** — before any estimation, so `e(sample)` does not exist
there. They are computed on `smp2` (ado:6460-6462), the same user-`if`-and-nonmissing-cluster
screen. That is why, in the measured arm, a run with 9 real clusters and 4 treated clusters
printed **neither** caution.

If BSAVE-2 ships alone, the situation gets *incoherent* rather than better: the caution
would be computed on 14 clusters while `e(boot_df)` says 8.

### The exact edit

**Delete** from `wooldid_bootctl`: ado:6479 (the `<10` caution line) and the
`if "`ttrevar'" != "" { ... }` block at ado:6491-6500 (the `<6 treated` caution). Keep the
`r(max) < 2` refusal at ado:6465-6468 — the help says explicitly that it "fires when the
options are screened, before the estimation runs", and it should keep doing so.

**Add** to `wooldid_bootrun`, immediately after the `wdb_C < 2` refusal at ado:8956-8959:

```stata
  if "`clus'" != "" & `wdb_C' < 10 dis as text "Caution: only `wdb_C' clusters are available for the bootstrap; bootstrap inference is unreliable with very few clusters."
  if "`clus'" != "" & "`ttrevar'" != "" {
    qui count if `smp' & `wdbtrtc' == 1 & `wdbctag' == 1
    if r(N) < 6 dis as text "Caution: only " r(N) " treated cluster(s) are available for resampling. What is implemented is the " cond("`mode'" == "bayes", "Bayesian bootstrap, which shares the pairs cluster bootstrap's first two moments", "pairs cluster bootstrap, which keeps only the replications that draw at least one treated cluster, so the retained distribution is conditional on that event") "; with few treated clusters the literature finds this bootstrap badly sized (MacKinnon and Webb 2017, 2018), and boottstat is not a cure for it. Consider permute() - randomization inference over the treatment assignment - or the wild cluster bootstrap-t (boottest, on the stage-1 regression outside wooldid), which wooldid does not offer."
  }
```

with the two flags built on the same `smp` marker, just before the `frame put` at ado:8939:

```stata
  tempvar wdbtrtr wdbtrtc wdbctag
  if "`clus'" != "" & "`ttrevar'" != "" {
    qui gen byte `wdbtrtr' = (`ttrevar' < .) if `smp'
    qui gegen byte `wdbtrtc' = max(`wdbtrtr') if `smp', by(`clusv')
    qui gegen byte `wdbctag' = tag(`clusv') if `smp'
  }
```

The message text is copied byte-for-byte. Two substitutions only: the local `bootstrap`
becomes `mode`, which is what `bootrun` calls the same option, and `r(max)` becomes `wdb_C`.

### Runtime cost — **net zero; the work is moved, not added**

This is a *relocation*. The four operations that compute the treated-cluster count
(`gen byte` + `gegen max` + `gegen tag` + `count`) are deleted from `wooldid_bootctl` and
recreated in `wooldid_bootrun` over the same rows — in fact over **fewer** rows, since after
BSAVE-2 `smp` is a subset of bootctl's `smp2`. The `<10` caution stops costing anything at
all: it reads `wdb_C`, which is already in hand, where bootctl needed a `gstats sum` to get
`r(max)`. Both blocks stay gated on `cluster()` being present. `bootctl`'s
`gegen group()` + `gstats sum` are **kept** — the `<2` refusal still needs them.

Measured on a 1,000,000-row / 40,000-cluster frame (`t5_cost.do` / `o5_cost.txt`):

| block | seconds |
|---|---|
| the four moved operations, as run in `wooldid_bootrun` | 0.183 |
| the whole old `bootctl` block including the `gegen group()` + `gstats sum` that stay | 0.237 |
| reference: one `gegen ... , by()` the program already runs | 0.059 |

So the honest number is **0.183 s moved from one program to another on a million rows**, plus
a small saving from dropping a `gstats sum` off the `<10` path. Nothing is added.

### Does it move a reported number? No

No estimate, no p-value, no `e()` scalar. What changes is **where in the output** two
cautions print: today at parse time, after the fix after the option screens and immediately
before the replication dots. That is a display-order change and the author's call — the RF-2
note at ado:6469-6472 shows he has positioned these deliberately before.

### Alternatives rejected

1. **Leave them at parse time and approximate `e(sample)` with `!missing(y)`.** Guesses at
   the screen; would still miss a cluster dropped for a missing control or a singleton, and
   would then disagree with `e(boot_df)` in the other direction.
2. **Keep the parse-time pair and add a second post-estimation pair.** Double-prints on every
   ordinary run.
3. **Do nothing and document it.** The help sells both cautions unconditionally
   (`sthlp:532`: "A caution fires when fewer than 6 treated clusters are available for
   resampling (beside the existing caution at fewer than 10 clusters in total)"). Documenting
   "unless some units have no usable outcome" is not a sentence anyone should have to write.

### Test

The ARM A / ARM B table above, checked on the numerals rather than on occurrence counts:

| grep | pinned (`t_pinned.log`) | patched (`t_patched.log`) |
|---|---|---|
| `only [0-9]* clusters are available` | `only 9 clusters are available` (**once** — ARM B only) | `only 9 clusters are available` **twice** (both arms) |
| `only [0-9]* treated cluster` | `only 4 treated cluster` (**once**) | `only 4 treated cluster` **twice** |

9 and 4 are the true counts in both arms, so the patched build prints the same numerals on the
arm the pinned build says nothing about. Grep for those **short** fragments, not whole
sentences — the batch log hard-wraps at 80 columns.

---

## Not proposed — and why

- **Nothing here argues for a help-only fix.** All three findings are code defects with
  cheap, testable code fixes. The two help edits named above are riders on code changes,
  not substitutes for them.

## Open questions for the author / next agent

1. **`wooldid_permrun` has the identical sample concept** at ado:7479-7486 — `tempvar smp`
   (7479), `qui gen byte smp = ((ifr) == 1)` (7482), `qui replace smp = 0 if missing(clus)`
   (7486) — the construction `bootrun`'s own comment calls "permrun byte-parallel". Whether the permutation **assignment**
   universe should also be `e(sample)` is a different question from the bootstrap's — an
   onset drawn onto a cluster that contributes nothing still consumes a draw and changes the
   distribution over the clusters that do contribute — but it is the same line of code and
   it belongs to mechanism `PERM`, not `BSAVE`. I did not probe it.
2. **The `reps()` p-value-resolution caution at ado:6473-6481 sits INSIDE the
   `if "clus" != "" {` block** that opens at ado:6459 and closes at ado:6501, so it cannot
   print on an unclustered `bootstrap()` run. Noticed while reading the caution block for
   BSAVE-3, and then measured on the pinned build: the smoke matrix's no-`cluster()` arm
   (`bootstrap(standard) reps(10)`, grid 2/10 = 0.20, well above the 0.05 threshold) printed
   it **0 times**, against 4 occurrences across the clustered arms of the same log
   (`t6_pinned.log`). Separate from every BSAVE finding, not fixed here — a one-line
   de-indent — and not otherwise probed; I have not checked whether it was deliberate.
3. `e(boot_sgWobs)` (the across-subgroup equality tests' observed statistics, ado:7153) has
   the same studentized-leg gap as `e(boot_jointWobs)`. `bootsave()` does not write `sgW`
   columns, so it is not part of this defect, but the same three-line treatment would close it.
