# Fix proposals — mechanism `PERM` (permute())

Author: prop-PERM agent, 2026-09-19.
Workdir (every probe, log and the patched ado):
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-PERM/`
Patched copy of the ado (**mine, the repo was never written to**): `probes/prop-PERM/mybuild/wooldid.ado`,
generated from the pinned build by `probes/prop-PERM/mkpatch.py`. Line numbers are on the
pinned build (`bluesky/build/wooldid.ado`, sha c1d8979, 17021 lines).

| # | findings | what it is | moves a reported number? | risk | tested? |
|---|---|---|---|---|---|
| **FIX-PERM-1** | `fineclus-perm-FP1` | refuse (safety(off): caution) a covariate that is provably inert in the observed fit when `permonsetcohorts` rebuilds `i` | only by converting a currently-printed, invalid `e(perm_p_*)` into an error | medium (new refusal) | **yes, patched and run** |
| **FIX-PERM-2a** | `SK-FC-C` | sthlp:743 + ado:1443 stop promising that "cluster strictly finer" fixes the joint test | no | low | text only |
| **FIX-PERM-2b** | `SK-FC-C` | widen the jointtests caution predicate from binary nesting to *clusters per treated cell* | no (caution only) | medium (fire rate) | **yes, patched and run**; threshold **not** settled |
| **FIX-PERM-3** | `SK-1` (the permute-nesting one) | make the r(459) nesting refusal locatable | no | none | **yes, patched and run** |

Two things I am **not** proposing, with reasons, are in "Do not fix" at the end.

---

# FIX-PERM-1 — the inert-covariate screen under the `permonsetcohorts` rebuild

**Finding:** `fineclus-perm-FP1` (high; prober PARTIAL, refuter re-grounded it on
`coarsecohortcontrols()`).

## What is actually wrong

`wooldid_permrun` substitutes the rebuilt `i` into the **varlist only** (`local repcall
\`y' \`repi' \`t'`, ado:7641). The replication option string built in the shell
(ado:125-163) still names the user's **original** variables in `fe()`, `controls()`,
`coarsecohortcontrols()`, `ccc_absorb()` and `interactivecontrols()`.

A covariate that is constant within the observed `i` is therefore:

* **exactly collinear** with the `i` fixed effects in the **observed** fit — silently
  dropped, and provably unable to move the point estimate; and
* **not** collinear inside a replication, because the rebuilt `i` *regroups* clusters drawn
  from different observed cohorts. The term goes live there and absorbs part of the
  between-observed-cohort variation the randomization test is generating.

The observed statistic and the reference statistics then stop being the same function of the
data, which is the one thing a randomization test may not give up. Nothing is printed.

### Measured on my own data, on the pinned build

Panel: 6 states (2 never-treated, 2 adopting at t=6, 2 at t=8), T=10, `i = coh` (the onset
grouping), `cluster(state)` strictly finer, `permute(shuffle) permexact permonsetcohorts`,
M = 90 enumerated assignments. `cohbase = 10*coh + .5*coh^2` is a cohort-level control;
`a_s` is a state-level one. Log: `p4_cost.log` (part C), `p2_perm.log`.

With a real effect (tau = 1.5), `b` printed in `%21x`:

| arm | `b` (hex) | permutation `p` | null SD |
|---|---|---|---|
| baseline | `+1.4a413b7465a7aX+000` | **0.02222222** | 0.51258333 |
| `+ coarsecohortcontrols(c.cohbase)` | `+1.4a413b7465a7aX+000` **bitwise identical** | **0.26666667** | 1.52686900 |
| `+ controls(cohbase)` | bitwise identical | 0.02222222 | 0.51258333 |

**A significant result (p = .022, the exact floor 2/90) becomes insignificant (p = .267, a
12x move) because the user added a control that cannot change the point estimate by a single
bit.** The null SD inflates 2.98x. On a sharp-null draw (`p2_perm.log`) the same option
leaves `b` bitwise identical and inflates the null SD 3.79x (0.27410173 -> 1.03768030),
which matches the refuter's independent median of 3.85 over 24 sharp-null datasets.

**Negative control, same build, same option:** at `cluster(coh) == i` there is no rebuild
(`e(perm_irebuilt) = 0`) and the option leaves `p` and the null SD **bit-identical**
(1.00000000 / 0.08711505 in both arms). So this is the rebuild, not the option.

### Which option slots are actually in the inert class (measured, `p1_inert.log`)

| slot with a cohort-constant covariate | observed `e(b)` | verdict |
|---|---|---|
| `controls(cohbase)` | **bitwise identical**, `e(r2)` identical | inert -> screen |
| `coarsecohortcontrols(c.cohbase)` | **bitwise identical**, `e(r2)` identical | inert -> screen |
| `fe(wave)`, `wave = (coh<=2)` | identical to 3 ulps, `e(r2)` identical | inert -> screen |
| `ccc_absorb(cohbase)` | r(322) — already refused by an existing safety screen | screen anyway (cheap, and the refusal is not guaranteed on every design) |
| `interactivecontrols(cohbase)` | **moves**: `-1.2f6d…X-005` vs `+1.c5b…X-006`, `e(r2)` .5798 vs .5677 | **LIVE — must not be screened** |

`interactivecontrols()` enters as `i.t#c.x` (the `interinject` at ado:1171), which is outside
the `i` fixed-effect span, so a cohort-constant `x` is estimable in the observed fit too and
the inert-here/live-there asymmetry does not arise. Screening it would be a false positive on
every run. This is the one scoping decision in the proposal and it is measured, not assumed.

### Which slot actually distorts the reference distribution — be precise about this

The screened slots are not equally guilty, and my own data says so:

| slot, cohort-constant covariate | observed `b` | permutation `p` | null SD |
|---|---|---|---|
| `coarsecohortcontrols(c.cohbase)` | bitwise identical | **.0222 -> .2667** | **x2.98** (and x3.79 on the sharp-null draw) |
| `controls(cohbase)` | bitwise identical | .0222 -> .0222 (**bit-identical**) | **bit-identical** |

On my design `controls()` left the whole permutation distribution bitwise unchanged; the
refuter, on his design, saw it change 208 of 210 reference statistics while moving `p` by at
most one grid step, in **both** directions. So:

* **`coarsecohortcontrols()` (and by construction `ccc_absorb()`) is the slot with the
  measured, large, one-directional distortion.** The reason is structural and is in the help:
  the option adds the control **and its interaction with every `Tz_it` dummy** (sthlp:1061-1070),
  so with a cohort-constant control the *whole block* is collinear in the observed fit and the
  *whole block* goes live at once inside a replication. It is also the option a user reaches
  for when they have a cohort-level control — precisely the inert class.
* **For `controls()` and `fe()` the case for refusing rests on observed-fit inertness alone**,
  not on a measured distortion: the term provably cannot move the point estimate, so deleting
  it costs the user nothing, and the refuter's measurement shows it *can* perturb the reference
  set on some designs. I am not claiming a large effect there and the proposal text should not
  be read as claiming one.

**Variant the author may prefer:** refuse on `coarsecohortcontrols()`/`ccc_absorb()` and
merely caution on `controls()`/`fe()`. That is a two-line change (carry the slot name
alongside each variable in `wdpinert` and branch on it). I did not make it the default
because the argument for refusing is the same in both — the term is inert in the observed
fit, so refusing it costs nothing — and a split rule is harder to document than a single one.

### The false-positive route is closed (`p7_window.do` / `p7_window.log`)

"Constant within `i`" is not literally the same as "collinear with the absorbed `_i` fixed
effect", because `_i` is set to `i` only on in-window relative periods and to 0 on the
reference pool (ado:928-929). If the reference pool held out-of-window rows from two cohorts
with different values of the covariate, the term would be estimable there and the screen
would refuse a live term. I built exactly that: T=20, onsets 14 and 17, then
`esprelength(3)` and `esprelength(5)`, so several pre-treatment years of both treated
cohorts fall outside the window and into the reference pool.

| window | `controls(cohbase)` | `coarsecohortcontrols(c.cohbase)` | `fe(wave)` |
|---|---|---|---|
| `esprelength(3) espostlength(2)` | bitwise identical to baseline | bitwise identical | bitwise identical |
| `esprelength(5) espostlength(2)` | bitwise identical | bitwise identical | bitwise identical |
| no window option (T=20) | bitwise identical | bitwise identical | identical to 1 ulp |

`e(r2)` is identical to eight decimals in every cell. So on the default path "constant within
`i` on the design sample" does imply inert in the observed fit even with a truncated
pre-window, and the predicate is safe to refuse on. I have not found a design where it is
not; see Open questions for the residual.

## The exact edit — four hunks

### Hunk 1 — `wooldid_permctl`'s syntax line, ado:5807

```
-  syntax , [permute(string) ... rawp(string) onc(string) replace]
+  syntax , [permute(string) ... rawp(string) onc(string) yv(varname) pfe(string asis) pctl(string asis) pccc(string asis) pcca(string asis) replace]
```

(`yv()` is used by FIX-PERM-3 below, not by this one.)

### Hunk 2 — the call site in `wooldid_validate`, ado:4927

```
-    ... si(`safety') psave(`permsave') level(`level') onc(`permonsetcohorts') `replace'
+    ... si(`safety') psave(`permsave') level(`level') onc(`permonsetcohorts') yv(`y') pfe(`fe') pctl(`controls') pccc(`"`coarsecohortcontrols'"') pcca(`ccc_absorb') `replace'
```

`y`, `fe`, `controls`, `coarsecohortcontrols` and `ccc_absorb` are all live locals at that
point (ado:4874 and the S5 syntax line at ado:4572).

### Hunk 3 — a new block in `wooldid_permctl`, inserted **immediately after** the
existing cluster-must-nest-`i` block, i.e. after ado:6081 (`  }` closing
`if "\`clus'" != "\`ivar'" & "\`onc'" == "" { ... }`)

This is deliberately written as the mirror of the block above it: same guard shape, same
gegen idiom, opposite `onc` condition.

```stata
  * ---- FP1 (blue-sky review 2026-09-18): the mirror of the block above, for the arm
  * the T1 relaxation opened. Under permonsetcohorts with a cluster STRICTLY FINER than
  * i, wooldid_permrun substitutes the rebuilt i into the VARLIST only; the rep option
  * string still names the user's original variables in fe(), controls(),
  * coarsecohortcontrols() and ccc_absorb(). A covariate that is constant within i on the
  * design sample is therefore exactly collinear with the i fixed effects in the OBSERVED
  * fit - it is dropped, and provably cannot move the point estimate - while inside a
  * replication the rebuilt i REGROUPS clusters drawn from different observed cohorts, the
  * term stops being collinear, goes live, and absorbs part of the between-observed-cohort
  * variation the randomization test is generating. The observed statistic and the
  * reference statistics are then not the same function of the data, which is the one
  * thing a randomization test may not give up.
  * MEASURED on the pinned build (probes/prop-PERM/p4_cost.log): coarsecohortcontrols()
  * with a cohort-level control leaves e(b) BITWISE identical and moves the permutation
  * p from .0222 to .2667 while inflating the null standard deviation 3.0x; an independent
  * refuter measured a median null-SD ratio of 3.85 over 24 sharp-null datasets, never
  * below 1. The bootstrap's A10 convention is NOT the same case and stays as it is: its
  * re-identification REFINES the partition, so a covariate constant within i stays inside
  * the new i's span and stays dropped.
  * NOT SCREENED: interactivecontrols(). Its terms enter as i.t#c.x (the interinject at
  * the head of the worker's control block), which is outside the i fixed-effect span, so
  * a cohort-constant x is live in the observed fit too and the asymmetry does not arise
  * (measured: e(b) and e(r2) both move when it is added).
  * COST: nothing at all unless permonsetcohorts is combined with a cluster that is not
  * the i variable. Then one group + two by-i extrema passes to learn whether the rebuild
  * fires at all (the same three gegens the block above runs in the no-onc arm), and one
  * tag pass per screened covariate. Measured at 200,000 rows with two covariates: 0.156
  * seconds, once, on a run that is about to fit the model 90-1000+ times.
  if "`clus'" != "`ivar'" & "`onc'" != "" {
    local wdpcvl
    foreach wdptok in `pfe' `pctl' `pccc' `pcca' {
      local wdptok : subinstr local wdptok "#" " ", all
      local wdptok : subinstr local wdptok "(" " ", all
      local wdptok : subinstr local wdptok ")" " ", all
      foreach wdppc of local wdptok {
        local wdpv = substr("`wdppc'", strrpos("`wdppc'", ".") + 1, .)
        cap confirm numeric variable `wdpv'
        if _rc == 0 & strpos(" `wdpcvl' ", " `wdpv' ") == 0 local wdpcvl `wdpcvl' `wdpv'
      }
    }
    if "`wdpcvl'" != "" {
      * does the rebuild actually fire? byte-parallel with wooldid_permrun's own test.
      * at cluster() == i as a partition the rebuild is a no-op, the replication keeps the
      * observed i, and an i-keyed covariate stays collinear in every replication - there
      * is nothing to say there, and the negative control must find nothing
      tempvar wdprg wdprmn wdprmx
      qui gegen double `wdprg'  = group(`clus') if `insmp'
      qui gegen double `wdprmn' = min(`wdprg') if `insmp', by(`ivar')
      qui gegen double `wdprmx' = max(`wdprg') if `insmp', by(`ivar')
      qui count if `insmp' & `wdprmn' != `wdprmx'
      local wdpreb = (r(N) > 0)
      qui drop `wdprg' `wdprmn' `wdprmx'
      if `wdpreb' {
        tempvar wdpti wdptv
        qui gegen byte `wdpti' = tag(`ivar') if `insmp'
        qui count if `wdpti' == 1
        local wdpnlev = r(N)
        qui drop `wdpti'
        local wdpinert
        foreach wdpv of local wdpcvl {
          qui gegen byte `wdptv' = tag(`ivar' `wdpv') if `insmp'
          qui count if `wdptv' == 1
          if r(N) == `wdpnlev' local wdpinert `wdpinert' `wdpv'
          qui drop `wdptv'
        }
        if "`wdpinert'" != "" {
          local wdpwhat = cond(wordcount("`wdpinert'") > 1, "are", "is")
          if "`si'" != "off" {
            dis as error "Error: under permonsetcohorts with a cluster() finer than `ivar', every replication rebuilds `ivar' from the onsets it drew, but the covariate specification names `wdpinert', which `wdpwhat' constant within `ivar' on the permutation design sample."
            dis as error "Such a term is exactly collinear with the `ivar' fixed effects in the observed fit - it is dropped there and cannot change the reported estimate - and is NOT collinear inside a replication, where the rebuilt `ivar' pools clusters from different observed cohorts. It goes live there and absorbs part of the between-cohort variation the randomization test generates, so the observed statistic and the reference statistics stop being the same function of the data and the permutation p-value is not a test of the statistic the table reports."
            dis as error "Drop the term from fe()/controls()/coarsecohortcontrols()/ccc_absorb() - it contributes nothing to the point estimate - or drop permonsetcohorts and cluster at or above `ivar'. A covariate that varies WITHIN `ivar' is unaffected and needs no change. To compute the permutation distribution anyway, specify safety(off)."
            error 459
          }
          dis as text "Caution: `wdpinert' `wdpwhat' constant within `ivar', so under permonsetcohorts' per-replication rebuild of `ivar' the term is dropped from the observed fit (it cannot change the reported estimate) and estimable inside the replications, where it absorbs part of the between-cohort variation the randomization test generates. The reported permutation p-value and null standard deviation are therefore not those of the statistic the table reports; measured inflation of the null standard deviation on the coarsecohortcontrols() fixtures was a factor of roughly 3 to 4. safety(off) is what allows this to run."
        }
      }
    }
  }
```

### Hunk 4 — help

`sthlp:1240` (the engine-compatibility row) and `sthlp:431`/`sthlp:68` each describe the
rebuild; one sentence belongs with them:

> Because the rebuild **regroups** `i` rather than refining it, a covariate that is constant
> within `i` — a cohort-level control, an adoption-wave fixed effect, `i` itself — behaves
> differently in a replication from the way it behaves in the observed fit: it is collinear
> with the `i` fixed effects and dropped in the observed fit, and estimable once the
> replication has regrouped the cohorts. `wooldid` refuses such a term under
> `permonsetcohorts` with a finer `cluster()` and names it; `safety(off)` computes it anyway.
> This is not the same situation as `bootstrap()`'s `idcluster` convention (sthlp:532), where
> re-identification **refines** the partition and the term stays dropped in every replication.

## Measured before/after on the patched build (`p3_patched.log`)

| arm | pinned build | patched build |
|---|---|---|
| baseline | p .88888889, SD 0.27410173, `b` `+1.c5b78511a2600X-006` | **bit-identical** |
| `controls(cohbase)` | runs silently | **r(459)**, message names `cohbase` |
| `coarsecohortcontrols(c.cohbase)` | runs silently, SD 1.03768030 | **r(459)** |
| same `+ safety(off)` | — | runs, Caution printed, p .86666667 / SD 1.03768030 = **exactly the pinned number** |
| `coarsecohortcontrols(c.a_s)` (state-level, live) | p .22222222 SD 0.86554478 | **identical, nothing fires** |
| `cluster(coh)` (== `i`, no rebuild) `+ ccc(c.cohbase)` | p 1.0 SD 0.08711505 | **identical, nothing fires** |
| `interactivecontrols(cohbase)` | p .64444444 SD 0.46641343 | **identical, nothing fires** |
| no `permute()` at all `+ controls(cohbase)` | `b` `+1.c5b…X-006` | **identical** |
| `permute()` without `permonsetcohorts` | p 1.0 SD 0.08711505 | **identical** |

Nine arms; the four that must be untouched are untouched bit for bit, the one that must be
refused is refused, and `safety(off)` reproduces today's number exactly.

## Cost

Measured on 200,000 rows / 2,000 clusters / `i` with 20 levels (`p4_cost.log` part A):

| step | seconds |
|---|---|
| does the rebuild fire (1 `gegen group` + 2 by-`i` extrema + `count`) | **0.111** |
| `gegen tag(i)` + `count` | 0.012 |
| `gegen tag(i x)` + `count`, per screened covariate | 0.015 – 0.018 |
| **total, two covariates** | **0.156 s, once** |
| for scale: one `fastnose` wooldid fit on a 16,000-row panel | 1.43 s |

It is gated twice — `permonsetcohorts` present **and** `cluster()` not the `i` variable — so
**every run that does not use `permonsetcohorts` pays exactly zero**, and a
`permonsetcohorts` run at `cluster() == i` pays zero too. A run that does reach it is about
to fit the model between 90 and several thousand times. This is the cheap end of the guards
already in `wooldid_permctl`.

## Does it move a reported number?

**No number changes on any run that keeps running.** Every arm above that is not refused is
bit-identical. What changes is that a set of runs that today print an invalid
`e(perm_p_*)`/`e(perm_sd_*)` now print an error instead. That is a behaviour change and needs
the author's ruling on the tier:

* **as shipped in my patch:** `error 459` under default safety, caution under `safety(off)`.
  That is the `bootstrap(bayes)` + cell-weighted-estimand pattern (sthlp:532) and the
  `contreatwithin` pattern.
* **one-line alternative:** delete the `if "\`si'" != "off" { ... error 459 }` wrapper and
  always print the caution. Then **nothing is refused and no number moves at all**, and the
  user is merely told. My recommendation is the refusal, because the option provably cannot
  change the estimate — refusing costs the user nothing but a deleted token — but I would
  not argue hard against the caution-only form, and it is the zero-risk choice.

## Blast radius

`wooldid_permctl` is reached only from `wooldid_validate` when `permute()` is given, and the
new block is additionally gated on `permonsetcohorts` + a non-`i` `cluster()`. That
combination is **three days old** (T1, 2026-09-16) and its assessment suite was deferred, so
there is no installed base to break. `wooldid_validate` runs again inside the permute shell
(ado:7577) and inside every replication's re-parse; the block is idempotent and, inside a
replication, `permonsetcohorts` is not in the rep option string (ado:125-163 never emits it),
so it does not fire R times.

Interactions checked in the patch run: `subgroups()` is already refused with
`permonsetcohorts` (ado:6028-6031), so there is no compound path; `permstrata()`, `permprobs()`
and `permtstat` do not touch the covariate slots.

## Alternatives rejected

1. **Auto-drop the inert terms from the rep option string `wdp_ro` (ado:125-163).** This is
   the *statistically* correct repair — it restores "observed and reference are the same
   function", and because the term is dropped in the observed fit anyway the observed estimate
   does not move. Rejected as the default because (a) it silently edits the model the user
   typed, (b) `e(perm_p_*)` moves on every affected run with nothing in the table to say why,
   so it needs a re-bless and a ruling of its own, and (c) it has a sample trap: a covariate
   that is constant within `i` but **missing** on some rows changes `e(sample)` through
   markout, so removing it from the replications would make the reps' sample differ from the
   observed one — the exact class of defect the `SMP` mechanism is full of. Worth offering to
   the author as a **second** option behind the refusal; if he wants it, the predicate must be
   tightened to "constant within `i` **and** never missing on the design sample".
2. **Re-key the covariate through the same onset->code map `wooldid_permirb` uses.** Inside a
   replication the re-keyed covariate is constant within the *rebuilt* `i`, hence collinear
   again, hence dropped — so this is alternative 1 with a Mata call added. Rejected.
3. **Detect it after the fact from the replications** (accumulate the stage-1 rank per
   replication and compare with the observed). Rejected: it reports only after every
   replication has run, it needs a new `e()` posting in the worker plus a harvest in
   `wooldid_permrun`, and it is not specific — rank legitimately moves when a drawn assignment
   empties a cell.
4. **Widen the bootstrap's A10 tripwire (ado:409-415) to the same predicate.** Rejected, and I
   think deliberately: under `bootstrap()` re-identification **refines** the partition
   (`group(__wdbns, i)`, ado:9323-9349), so a covariate constant within the old `i` is constant
   within the new one, stays collinear and stays dropped. There is no asymmetry to warn about.
   The refuter's observation that A10's bare-token test misses `controls(cohbase)` is correct
   and is already documented at sthlp:532 ("The check compares bare tokens"), and it is
   harmless there. **Do not change A10.**

## The test that would prove it correct (no oracle needed)

Three invariants, all of which my patch run already satisfies:

1. **Negative control by geometry.** The same option at `cluster() == i` (rebuild a no-op)
   must leave `e(perm_p_*)` and `e(perm_sd_*)` **bit-identical** to the pre-patch build.
   (Q6 vs P5: 1.00000000 / 0.08711505 both.)
2. **`safety(off)` identity.** The `safety(off)` arm must reproduce the pre-patch numbers
   bitwise, because the patch adds a screen and changes no computation. (Q4 vs P2:
   .86666667 / 1.03768030 both.)
3. **Inertness is the premise, so assert the premise.** For every variable the screen names,
   a plain run with and without the option must give a **bitwise identical** `e(b)`. This is
   the invariant that makes the refusal safe: if it ever failed, the screen would be refusing
   a term that does real work. `p1_inert.do` is that test; it is cheap enough to live in a
   cert leg.

A fourth, stronger test I did not run and would recommend for the suite: under `permexact`,
dump `e(permdist_b)` with and without the inert term after proving the enumerations coincide
(the `permsave` onset blocks match row for row). The difference at the identity assignment
must be exactly zero, and non-zero at every other assignment — that is the mechanism, stated
as an invariant. The refuter did exactly this (208 of 210 reference statistics changed,
difference at the identity assignment bitwise zero).

---

# FIX-PERM-2 — `SK-FC-C`: the help names a remedy that does not deliver a sized test

**Finding:** `SK-FC-C` (high, REFUTER_FOUND; two reviewer passes downgraded the magnitude to
"16-18%" and "roughly 7x" against the original 61.7%).

## The fix site and what is and is not in dispute

Not in dispute, and the part that decides the fix: **sthlp:743 and the runtime caution text at
ado:1443 both end by naming "cluster strictly finer than the cohort-by-period cells" as a
remedy that gives "a joint test that is valid under this design", and in that geometry the
caution predicate at ado:1421-1432 is silent by construction.** The predicate compares
`count(tag(i))` with `count(tag(i clus))`: it is the strictly binary question "does every `i`
lie wholly inside one cluster", so it cannot see 2, 3, 5 or 10 clusters per cell. In dispute
is only how oversized the test is there (61.7% vs 16-18% at nominal 5%).

**The direction is not in dispute, and the direction is what the text fix needs.** The help
already contains the right idea in the wrong scope: the *second* reliability condition in the
same paragraph says designs "with only a handful of observations per cohort-period cell yield
oversized joint tests **even under robust or finely clustered SEs**". That is this failure,
with clusters-per-cell scarce instead of observations-per-cell.

**This does not contradict `headline.md`.** Pointwise `att` coverage at finer-than-`i` is
fine (SE/SD 0.970, coverage 0.935 over ~1,500 draws). It is the 4-degree-of-freedom joint
quadratic form that breaks, and it breaks for the standard reason: a cluster-robust covariance
estimated from few effective clusters per tested direction is downward biased, and the bias
compounds through `V^{-1}` across horizons in a way it does not in a single `se`.

## FIX-PERM-2a (recommended, no number moves): qualify the remedy

**sthlp:743**, current text of the clause:

> ... use bootstrap-based clustered standard errors or randomization inference (see
> {it:permute()}), or cluster strictly finer than the cohort-by-period cells.

becomes

> ... use bootstrap-based clustered standard errors or randomization inference (see
> {it:permute()}). Clustering strictly finer than the cohort-by-period cell removes the
> degeneracy, but it is not by itself a cure for the joint test's size: the joint statistic
> inverts a cluster-robust covariance over several horizons at once, and that covariance is
> downward biased when a treated cohort-by-period cell contains only a handful of clusters,
> so the test stays oversized there. Cluster finer **and** check that each treated cell holds
> many clusters; the pointwise {it:att} standard errors do not carry this restriction.

and the same replacement in the runtime caution string at **ado:1443** (the sentence ends
"..., as does clustering strictly finer than the cohort-by-period cells." -> "Clustering
strictly finer than the cohort-by-period cells removes the degeneracy but does not by itself
size the joint test, which stays oversized when a treated cell holds only a few clusters.").

**Cost:** zero. **Moves a number:** no. **Blast radius:** two strings. This is the part of
the proposal I would ship regardless of what happens to 2b, because as written the help sends
a user from a documented, cautioned failure into an undocumented, uncautioned one.

## FIX-PERM-2b (needs a threshold ruling): widen the caution predicate

Inside the existing `if ("\`jointtests'" != "" | "\`pretrendpower'" != "") & "\`bootstrap'\`permute'" == ""`
block (ado:1414-1447), after the `jtcoarse` branch, add a second branch that fires when the
binary one did **not**:

```stata
    * SK-FC-C: the binary nesting test above is 'every i wholly inside one cluster', so it
    * is silent at 2, 3, 5 or 10 clusters per treated cell - which is the geometry its own
    * remedy sentence sends users to, and where the joint test is still oversized because the
    * cluster-robust covariance the quadratic form inverts is estimated from that many
    * clusters per tested direction. Only evaluated when the binary test did NOT fire, so a
    * cautioned run pays nothing extra.
    if "`jtcoarse'" == "" & "`jtclusvars'" != "" & wordcount("`jtclusvars'") == 1 {
      tempvar jtct jtnc
      qui gegen byte `jtct' = tag(`_i' `_t' `jtclusvars') if `ifstatement' & `_i' != 0
      qui gegen double `jtnc' = total(`jtct') if `ifstatement' & `_i' != 0, by(`_i' `_t')
      qui summarize `jtnc' if `ifstatement' & `_i' != 0, meanonly
      local jtmincl = r(min)
      if `jtmincl' < `jtminclthresh' {   // jtminclthresh set to 6 in the tested patch
        dis as text "Caution: `=cond("`jointtests'" != "", "jointtests", "pretrendpower()")' with analytic clustered SEs, and the smallest treated `i'-by-period cell contains only `jtmincl' distinct `jtclusvars' cluster(s). Clustering finer than `i' removes the degeneracy the caution above describes, but the joint statistic inverts a cluster-robust covariance over several horizons at once and that covariance is downward biased when a cell holds few clusters, so the joint test is still oversized here. The pointwise standard errors are not affected. Use bootstrap-based clustered SEs or permute() for a joint test, or a design with more clusters per treated cell."
      }
    }
```

`_i`, `_t` and `targetflag` are all in scope at ado:1414 (`_i` is declared at ado:793 and filled at ado:928-929;
the si8 block at ado:1361 already uses exactly this `gegen tag(\`_i' \`_t')` idiom on the same
sample, which is the precedent for the cost).

**Cost:** two `gegen` passes plus one `summarize` over the estimation sample, on runs that ask
for `jointtests` or `pretrendpower()` and are **not** already cautioned. Extrapolating the
measured `gegen tag` cost at 200,000 rows (0.012-0.018 s) plus a `total(), by()` pass, this is
about **0.05-0.08 s at 200,000 rows**, once. It is the same order as the two `gegen tag`
passes per cluster variable the block already runs unconditionally.

**Moves a number:** no — it is a caution only, printed text.

**What is not settled, and why I am not choosing the threshold.** The evidence is two rungs
from the original probe (clusters per cell 3 -> 0.617 rejection at nominal 5%; 10 -> 0.108),
both disputed by two reviewers who measured 16-18% at the same geometry, plus my own rung
below. A threshold picked from a disputed pair will either fire on most ordinary state-level
designs (if set at 10+) or miss the regime it exists for (if set at 3). **This is the
author's parameter.** What I can say is that the evidence rules out *no* threshold: every
measurement of the finer-than-`i` arm, including the most conservative reviewer's, is at
least 2x nominal.

### My own rung (`p5_joint.do` / `p5_joint.log`)

I ran the size study on my own DGP rather than repeat a disputed one. 40 states throughout
(G = 40 fixed), T = 16, `y = state random effect + 0.15*yr + iid N(0, .6)` — **no effect and
no pre-trend, so the joint pre-trend null is exactly true**. `i = coh` with 40/K cohorts of K
states each (cohort 1 never-treated), `cluster(state)`, `jointtests esprelength(4)
espostlength(3)`, p read from `e(joint_horiz)` (the scalar behind the printed "p-value for
joint test that all pre-treatment att horizon estimates are 0"). Clusters per treated
cohort-by-period cell verified directly: exactly K. 100 draws per rung.

| arm | clusters per treated cell | rejection at nominal 5% | at nominal 10% |
|---|---|---|---|
| `i`=coh, `cluster(state)` — **the help's remedy geometry** | 4 | **0.100** | 0.150 |
| `i`=coh, `cluster(state)` | 10 | 0.060 | 0.080 |
| `i`=state, `cluster(state)` (at `i`, the **cautioned** case) | 1 | *arm did not finish inside its deadline; killed. `p5_joint.log`* | — |

The at-`i` positive control (arm C in `p5_joint.do`) was killed at its deadline — `i` = state gives 40 cohorts and the run is an order of magnitude slower — so I have no harness validation against the documented 11-36%; the two finer-than-`i` rungs stand on their own. Binomial SE at 100 draws is about 0.03 at these rates, so the 4-cluster rung is about 2x
nominal and the 10-cluster rung is indistinguishable from nominal. **This is monotone in
clusters-per-cell, which is what the mechanism predicts, and it is much milder than the
probe's 0.617 — closer to the two reviewers' 16-18%.** My DGP is benign (a state random
effect plus homoskedastic iid noise, no serial correlation), which is the obvious reason my
number is at the low end; an adversarial error structure is where the larger numbers come
from. I record it as a third independent point, not as a refutation of either side.

**What it means for the threshold.** Three independent measurements now exist for the
finer-than-`i` geometry and they disagree by a factor of 6 on magnitude while agreeing on
sign and on the direction of the gradient. A threshold of **6 distinct clusters in the
smallest treated cell** would have fired on my oversized rung and stayed silent on my
nominal one; it is also the number the package already uses for "few treated clusters" under
`bootstrap()` (sthlp:532). I offer 6 as a defensible starting point, not as a measured
optimum, and the honest statement in the caution is "oversized here too", not a rate.

**The predicate as written was patched in and run** (`mybuild2`, `p6_2b.do`, `p6_2b.log`),
with `jtminclthresh` set to 6:

| arm | result |
|---|---|
| `i`=coh, 4 clusters per cell | **new caution fires**, text names "only 4 distinct state cluster(s)" |
| `i`=coh, 10 clusters per cell | silent |
| `i`=state (at `i`) | **only the existing caution fires**; the new branch is not reached |
| no `jointtests` | nothing |

and its measured cost on a 200,000-row panel (`gegen tag` + `gegen total, by()` +
`summarize`) is **0.042 s**, once, on runs that ask for `jointtests`/`pretrendpower()` and
are not already cautioned.


## Rejected for SK-FC-C

* **Changing the joint statistic** (a small-sample F reference, a `G-q` correction, a
  Bell-McCaffrey / CR2 covariance). That is a new estimator, it moves every reported joint
  p-value on every clustered run, and the project has ruled twice (PS-64/BS-9) that the
  statistic is computed as-is with a caution. Out of scope for a review fix.
* **Refusing `jointtests` in the finer-than-`i` geometry.** PS-64's ruling already went the
  other way on the at-or-above-`i` case for a good reason ("a hard error sent users into a
  documented dead end"); a refusal here would be worse, since this is the geometry the help
  recommends for everything else.
* **Making the existing binary predicate non-binary** by replacing it rather than adding to
  it. Rejected: the at-or-above-`i` caution says something different and correct (degenerate
  covariance, rank-deficient in the probe's own measurement: rank 3/4, cond 6.2e15) and must
  keep firing with its own text.

---

# FIX-PERM-3 — `SK-1`: make the `permute()` nesting refusal locatable

**Finding:** the `permute()`-nesting `SK-1` (severity **low** in `findings_raw.json`; see
Open questions on the id collision). The refusal itself is correct and the finding says so —
sthlp:430 puts markout-dead rows in the design and sthlp:431 requires onset-homogeneous
clusters. What is wrong is only the diagnostic.

## The exact edit — ado:5977-5981

```
-  qui count if `insmp' & ((`ntr' > 0 & `nnv' > 0) | (`ntr' > 0 & `tmin' != `tmax'))
-  if r(N) > 0 {
-    dis as error "Error: permute() requires clusters to be subsets of treatment-onset groups (levels of the timing variable, with never-treated its own group), but `r(N)' in-sample observation(s) sit in clusters that mix onset values or mix treated and never-treated rows."
-    error 459
-  }
+  tempvar wdpmix
+  qui gen byte `wdpmix' = `insmp' & ((`ntr' > 0 & `nnv' > 0) | (`ntr' > 0 & `tmin' != `tmax'))
+  qui count if `wdpmix'
+  if r(N) > 0 {
+    * SK-1 (blue-sky review 2026-09-18): the refusal is right - sthlp:430 puts markout-dead
+    * rows in the permutation DESIGN and sthlp:431 requires every design cluster to be
+    * onset-homogeneous - but the diagnostic was not locatable. The count is the whole
+    * offending cluster's rows, most of which are ordinary live rows, and the rows actually
+    * responsible are typically OUTSIDE e(sample), so e(N), e(N_clust) and the point
+    * estimate of the run that works are all identical to the run that is refused. Name the
+    * smallest offending cluster (the idiom the within-i screen at the head of the worker
+    * already uses), count the offending CLUSTERS, and say that the screen runs on the wider
+    * design sample. Everything here is on the error path.
+    local wdpnmixr = r(N)
+    tempvar wdpmixc wdpmixt
+    qui gegen byte `wdpmixt' = tag(`clus') if `wdpmix'
+    qui count if `wdpmixt' == 1
+    local wdpnmixc = r(N)
+    qui gen double `wdpmixc' = `clus' if `wdpmix'
+    qui gstats sum `wdpmixc', meanonly
+    local wdpmixid = strtrim(string(r(min), "%18.0g"))
+    local wdpdead 0
+    if "`yv'" != "" {
+      qui count if `wdpmix' & missing(`yv')
+      local wdpdead = r(N)
+    }
+    dis as error "Error: permute() requires clusters to be subsets of treatment-onset groups (levels of the timing variable, with never-treated its own group), but `wdpnmixr' in-sample observation(s) in `wdpnmixc' cluster(s) of `clus' mix onset values or mix treated and never-treated rows (smallest offending `clus': `wdpmixid')."
+    dis as error "This screen runs on the permutation DESIGN sample - every row your if-restriction keeps that has a `clus' value - which is wider than the estimation sample: a row the estimator itself discards for a missing outcome, weight or control still carries an onset here and still counts toward the design and its assignment count. `wdpdead' of the offending observation(s) have a missing `yv'. That is why a specification that runs without permute() can be refused the moment permute() is added."
+    dis as error "Common causes: a unit kept in the panel skeleton by a merge miss, with a missing (or stale) onset inside an otherwise homogeneous cluster; an onset recorded only from the onset period onward; or a cluster that genuinely pools units adopting at different times, in which case permute() must cluster at or above the level at which onset is assigned. Restrict the sample with if, fill in the onset on every row of the offending cluster, or cluster at the assignment level."
+    error 459
+  }
+  qui drop `wdpmix'
```

plus `yv(varname)` on the syntax line (Hunk 1 above) and `yv(\`y')` at the call site (Hunk 2).

## Measured (`p3_patched.log`, blocks S1/S2)

The finding's own design: 8 states x 3 counties x 8 years, state 3 adopts in 2006, its county
9 is a merge miss with `ttre == .` and `y == .` on all 8 of its rows. `i` = county,
`cluster(st)`.

* S1, without `permute()`: runs clean, `e(N)` = 184, `e(N_clust)` = 8.
* S2, with `permute(shuffle) permexact`: r(459), and the message now reads

  > ... but **24** in-sample observation(s) in **1** cluster(s) of `st` mix onset values or mix
  > treated and never-treated rows (**smallest offending `st`: 3**).
  > This screen runs on the permutation DESIGN sample ... **8 of the offending observation(s)
  > have a missing `yy`.** That is why a specification that runs without `permute()` can be
  > refused the moment `permute()` is added.

That is exactly the three things the finding said the user could not get from anywhere: *which*
cluster, that the design sample is wider than `e(sample)`, and *how many* of the counted rows
are the dead ones.

## Cost

The added `gegen tag` / `gstats sum` / `count` are all inside `if r(N) > 0`, i.e. on the error
path. The only success-path change is `count if \`wdpmix'` replacing `count if <expression>`,
plus one `gen byte`. Measured on 200,000 rows: **0.008 s** (`p4_cost.log` timer 5).

## Moves a number / blast radius / alternatives

Moves nothing; `error 459` and the condition are byte-identical. The message text is the only
change a user sees. Blast radius: one refusal message reached only under `permute()`.

Rejected alternative: **evaluate the screen on `e(sample)` instead of the design sample.** That
would make the specification run, but it changes the permutation *design* and hence
`e(perm_npossible)`, and it contradicts sthlp:430, which deliberately puts markout-dead rows in
the design. The finding agrees. This is a message defect, not a screen defect.

Rejected alternative: **listing every offending cluster.** `wooldid`'s house idiom (ado:772,
the within-`i` screen) names the *smallest* offending value, and a message that could print
hundreds of ids is worse. The cluster count plus the smallest id gives the user the entry point.

---

# Do not fix

1. **The bootstrap's A10 predicate (ado:409-415).** Its bare-token narrowness is real, is
   documented (sthlp:532), and is **harmless**: bootstrap re-identification refines the
   partition, so a covariate constant within `i` stays collinear and stays dropped in every
   replication. Widening it would add a caution that is never about anything. See
   FIX-PERM-1's rejected alternative 4.
2. **`interactivecontrols()` under `permonsetcohorts`.** Measured live in the observed fit
   (`e(b)` and `e(r2)` both move), so the inert-here/live-there asymmetry does not arise and
   screening it would be a false positive on every run. There *is* a residual conceptual point
   — a covariate keyed to the observed onset is a function of the assignment, so the
   randomization null is arguably not sharp for it either — but that is a modelling statement
   about the user's specification, not an implementation defect, and nothing in the
   implementation is inconsistent there.
3. **`ccc_absorb()` with a cohort-constant control** is, on the design I tried, already caught
   by the existing multicollinearity safety screen with r(322). It is still in FIX-PERM-1's
   screened list, because the r(322) is a downstream accident rather than a guarantee, and
   including it costs one more `gegen tag`.

---

# Open questions / what I could not test

1. **The tier for FIX-PERM-1: refusal or caution?** My patch refuses under default safety.
   The caution-only form is a one-line deletion and moves no number at all. **Author's ruling.**
2. **The threshold in FIX-PERM-2b.** Two rungs from the original probe, two reviewer
   re-measurements that disagree with them by 4x, and one rung of my own. Not enough to pick a
   number. What would settle it: the size study at clusters-per-cell in {2, 3, 4, 5, 8, 10, 20}
   holding G fixed, ~200 draws per rung, which is a few hours on this hardware.
3. **The residual false-positive risk for FIX-PERM-1's predicate is now small but not zero.**
   `p7_window.do` closed the route I was most worried about (a truncated pre-window putting two
   cohorts' rows into the `_i = 0` reference pool: the term is *still* dropped bitwise at
   `esprelength(3)` and `(5)`). Not tested: the non-reghdfe engines (`unconditionalse`,
   `customvce()`, `poisson`), where the fixed effects are explicit dummies rather than
   absorbed, and `esfixedbaseperiod`, which changes what the reference pool is. If the author
   wants belt and braces, the exact predicate is "constant within `_i`", which is computable —
   but `_i` is a worker tempvar and `wooldid_permctl` runs at parse time, so that version of the
   screen has to move into the worker, after `wooldid_reltime`.
4. **The id collision on `SK-1`.** Three different agents used the slug `SK-1`. In
   `findings_raw.json` the permute-nesting one is `severity_final: low`; four of the six
   reviewer verdicts printed under `SK-1` in the mechanism brief belong to the `level(99.9)`
   and `bootsave()` findings, which are not `PERM` surfaces. FIX-PERM-3 addresses the
   permute-nesting one only. The `level(99.9)` finding (CONFIRMED, high) is a genuinely
   separate defect and needs its own proposal under whatever mechanism owns `LEVEL`.
5. **`bootstrap()` and `permute()` coverage at a finer-than-`i` cluster** remains untested
   (unadjudicated.md §2). FIX-PERM-1 removes one way the permutation reference distribution
   can be wrong; it says nothing about whether the rebuild delivers nominal size when no
   covariate is involved. That measurement is the one that would validate T1 as a whole.
