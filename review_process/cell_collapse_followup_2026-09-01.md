# Cell collapse, follow-up: option cases, mechanics, and a break-even rule (2026-09-01)

Second investigation of the collapse-to-cells strategy, answering the author's specific
questions. Eight agents: two mechanics probes, four option cases, a break-even pass, and
an adversarial verifier. Companion to `cell_collapse_2026-08-31.md`, whose safe gate and
excluded vectors (`interactivecontrols()`, non-`dydx` operators, `contreat()`,
`unconditionalse`) are assumed and were not re-litigated.

**Headline: every case tested came back clean, the mechanics questions have decisive
answers, and a single `att` already pays for the collapse. The adversarial pass did not
break the method - it broke the measuring instrument and the break-even formula.**

---

## 1. The four option cases

### (a) Fast mode - yes, collapse still pays, but roughly half as much

`fastnose` removes **68-73% of the per-row cost** (2.393 -> 0.653 microseconds/row at the
main site) but not all of it: a surviving O(N) pass still forms the point estimate, and
that is exactly what the collapse deletes. `tCOLL` stays flat in N in both modes.

| N | ratio without fastnose | ratio with fastnose |
|---|---|---|
| 29k | 4.07x | 2.82x |
| 116k | 10.57x | 5.69x |
| 465k | 40.93x | 20.50x |

End-to-end, collapse on top of fastnose removes **16.1%** at 29k and **18.3%** at 465k -
near the ceiling, since margins is only 26-30% of a fast run. Break-even is under one
margins call at every size >= 29k.

**The bonus finding is the interesting one: collapse-only beats fastnose-only while
keeping full delta-method standard errors** (0.544s vs 0.585s at 29k; 5.605s vs 6.471s at
465k, projected). If that holds up, the collapse offers fast-mode speed *without* asking
the user to give up inference - which is a better trade than `fastnose` itself.

An expectation of mine was refuted: I predicted high-K compression would vanish under
`nose`. It got *worse* (2.48x -> 1.57x), because at p=683 `tCOLL` is itself 0.18-0.78s.
The fast-mode gain is strongest at large N and small p.

### (b) Sample weights at the unit-within-cohort level - clean

Random lognormal weights, one draw per unit, propagated within unit; verified in-band as
genuinely unit-level and unit-within-cohort. **Stressor confirmed real**: within-cell
weight variation in **60 of 60 cells**, weight range 0.00107 to 814.34 (763,822:1).

Verdict **ULPCLASS** for both `aweight` and `pweight`, at both the ATT and event-study
sites: worst b relative 7.6e-15, V 6.0e-16, se 1.6e-16. That is roughly **20x cleaner than
the precedent**, because the random weights moved the worst cell off the
cancellation-dominated pre-trends and onto the large ATT coefficient.

Two things make this a real test rather than a vacuous pass:
- **Negative controls with teeth.** Replacing the per-cell summed weight with `firstnm(w)`
  moves b by 5.3e-01 and with `count(w)` by 7.9e-02 - twelve to thirteen orders above the
  7.6e-15 signal. Per-cell weight *summing* is doing load-bearing work.
- **A mechanism control that finally returned non-null.** Re-running with only the row
  order permuted makes FULL vs COLL **exactly bit-exact**, while row order alone moves
  pristine wooldid by 3.9e-14. The collapse's effect is smaller than reordering the input,
  and vanishes entirely in one ordering.

Side answers: weights do **not** change the cell structure or collapse ratio (the weight is
payload, not key). And weights make the collapse pay off **more** - they inflate the full
margins call 3-7x while the collapsed call stays flat.

### (c) Subgroups spanning multiple cohorts - clean, and the literal case is refused

The premise needed splitting, because wooldid distinguishes two structures:

- **Supported spanning** (subgroup constant within each `i`, one level covering several
  cohorts): legal, and **ULPCLASS** in every configuration - weighted, unweighted, ragged,
  at both the level sites and the pwcompare sites.
- **The subgroup varying *within* an `i`-by-`t` cell**: **refused at parse time, rc 459**,
  by the F13 containment guard at lines 836-848. The author's literal configuration is a
  `safety(off)` regime. It was tested there too, and is also correct.

**TERMINOLOGY CORRECTION (author-flagged).** An earlier draft of this section described the
requirement as "the subgroup must be nested in `i`". That is backwards. The guard requires
the subgroup to be CONSTANT WITHIN each `i`, so each `i` maps to exactly one subgroup level:
**`i` nests inside the subgroup, and each `i` is a subset of a subgroup.** The subgroup is a
*coarsening* of `i`. The guard's own remedy says so - it tells the user to define `i` at the
FINER level (`egen newi = group(i subgroups)`). The analysis was never affected, because the
load-bearing statement used throughout is the correct one: the subgroup is a *function of*
`_i` within any subgroup subpop.

**The key-split worry does not materialise.** Inside any subgroup subpop the rows are
in-window treated rows, so `_i` equals `i`; F13 forces the subgroup constant within `i`;
therefore the subgroup is a *function of* `_i` there and cannot split a cell. Measured
key-split ratio **1.000 at every site** across all supported configurations. Only the
`safety(off)` regime splits (ratio 3.0), and it is exactly correct there.

Weights showed no detectable interaction: differences moved in both directions with no
pattern, the quantum unchanged. 928 subpop-ratio measurements, minimum 69.25, **zero
degenerate sites**.

**SPEED, measured afterwards to close this gap** (the case-C agent measured correctness
only). Spanning subgroup, `cluster(cohid)`, `jointtests`, N=29,667, collapse-once/restore-once,
3 interleaved reps: **2.744s -> 1.092s, a 60.2% saving (2.51x)**, with the collapse costing
0.024s and preserve/restore ~0. `NONVACUITY_FAIL=0`; a `noesample`-only control arm was
bit-identical to pristine, so the movement is the collapse and not the option.

**Do the subgroup and subgroup-test estimates come back identical?** Not bit-identical -
ULP-class, as everywhere else. At rtol 0 there are 16 difference groups, spread over
`b_es_att_s2`, `pwsgs`, `wdidmainresults` and both `wdidsgresults_*` matrices, all in the
last digits. At rtol 1e-13 only **two** survive, and both are the per-subgroup joint-test
p-values, which wooldid stores as MACROS (text at display precision) rather than as matrix
cells: `joint_horiz_s1` .0107726416627541 -> .0107726416627534 (~6.5e-14 relative) and
`joint_horiz_s2` 1.89593725727e-07 -> 1.89593725733e-07 (~3.2e-11 relative, tail
amplification on a p of 1.9e-07). Because they are stored as strings they are flagged at any
tolerance; nothing about them is structurally different.

Two operational findings worth keeping:
- `jointtests` is *required* to reach the per-subgroup joint-test sites at all, because
  `suppressjointtests` defaults on - they are otherwise dead code.
- The standalone subgroup level sites are dead code whenever more than one subgroup is
  realised, because the EFF3 change (commit `73ece5c`) defers them to the pwcompare pass.

### (d) Coarse cohort controls - the source prediction confirmed

Predicted safe from the source: line 1035 builds `i._i#i._t#(ccc)` with **no
`targetflag`**, unlike `interactivecontrols()` at line 1058. So the ccc contribution
appears identically in `xb(tf=1)` and `xb(tf=0)` and subtracts out.

**Confirmed, ULPCLASS on every configuration** - continuous ccc, categorical ccc,
ccc combined with plain `controls()`, and an amplified variant where ccc is a dominant
driver of the outcome. Worst b relative 3.9e-15.

The verification is unusually tight:
- **The option was proven live**, not silently ignored: column-name audit gives
  `k_ccc` 62-155 nonzero coefficients and `k_ccc_x_tf = 0` in every config, with
  `k_tf = 31` as a positive control that the scanner can see targetflag at all. Cross-check:
  ATT differs with ccc (0.6501) versus plain controls (0.6498), so the fit really changed.
- **Within-cell variation was 100%**: the ccc variables vary in 24/24, 12/12 and 21/21
  subpop cells, so the collapsed row's arbitrary ccc value could have mattered.
- **A falsification probe** corrupts ccc in the retained cell rows and re-runs: shifting the
  ccc contribution by ~60x the ATT moved the answer by at most **7.1e-15**. Had ccc entered
  the derivative, the AME would have moved by ~31.

Exactly 11 nonzero b differences across four configs and three sites, **every one equal to
2^-49 to the last bit**.

**SPEED (per call, directional, from the case-D runs):** the collapse is **4.2x to 15.5x**
faster at every site of every ccc configuration - continuous 12.9/4.2/6.7x, categorical
9.1/5.9/6.8x, mixed-with-controls 9.9/5.7/15.5x, amplified 7.9/4.3/5.6x.
**End-to-end, measured afterwards** (`coarsecohortcontrols(c.x1 c.x2)`, robust, jointtests,
N=29,667, 3 interleaved reps): **2.225s -> 1.454s, a 34.7% saving (1.53x)**. Correctness on
that same run: **IDENTICAL at rtol 1e-13** (4 last-digit groups at rtol 0), with the
`noesample`-only control bit-identical to pristine.

---

## 2. How the collapse should be done (the author's mechanics questions)

### Is a literal collapse necessary? No - tag one row per cell, and take the *first*

**A literal `gcollapse` is not required. Tagging one row per cell is bit-for-bit identical
to it - provided you tag the first row - and it is cheaper.**

But the author's follow-up ("or does it suffice to select one *random* observation?") has a
sharper answer: **which row you keep does move reported numbers**, by exactly one ULP of
the linear prediction. Keeping the first row is bit-exact against a full collapse; keeping
the last or a random row is ULP-class different.

The mechanism was measured, not asserted. margins forms `dydx(targetflag)` as
`xb(tf=1) - xb(tf=0)`. Both predictions carry the *same* control partial sum `S` from the
retained row, so the controls cancel **algebraically** but the two roundings do not cancel
**bitwise**: `(S+c1) - (S+c0)` leaves a residue of order `ULP(|xb|)`. The retained row sets
`S`. A tripwire confirms the channel: with `controls()` removed from the model, all six
build variants are mutually bit-exact even with the controls still varying within cells.

**Recommended build (V4)**, bit-exact against a full `gcollapse` and 1.3-2.5x cheaper:

```
gcollapse (sum) w = wvar (min) pick = _n, by(<key>) merge
keep if _n == pick
```

A full `gcollapse`'s cost tracks the **column count** (0.040s -> 0.109s when 40 columns are
added at fixed N), because it aggregates every carried column for nothing.

**One caveat worth keeping.** A literal collapse can express *per-column aggregators* -
that is how the `contreat()` breakage was repaired, carrying the dose as a weighted cell
mean. No tag-one-row build can express that. Inside the present gate the carried values are
provably irrelevant and the aggregator is dead weight; if the gate is ever widened, the
literal shape is needed again.

### Restoring the data: preserve/restore wins decisively, and frames are disqualified

Measured on 2,000,000 rows / 320 MB:

| approach | cost | verdict |
|---|---|---|
| **preserve + restore** | **0.008s + <0.0005s** | **keep** |
| correct frames pattern | 0.107s (13x) | disqualified |
| tempfile save/use | 0.124s (15x) | disqualified |

`preserve` moves 336 MB in 0.009s (~37 GB/s) - an in-memory copy, not a disk write (no file
appears in the temp directory). `restore` is a pointer swap. It is **not** copy-on-write:
an in-place `replace` on 320 MB costs the same with and without an outstanding preserve, so
no hidden cost is deferred into the collapse. The prior session's suspiciously cheap
"0.000-0.001s" was **verified, not an artifact**.

On correctness, the earlier `r(322)` hazard is **inherent, and not frames-specific**.
Running `margins` in another frame re-points the estimation results at that frame and
leaves the original marker dead - and *every* restore leg except Stata's own `restore`
kills it, including a plain tempfile round trip. There is a repair the earlier pass missed
(`_estimates hold NAME, copy` / `_estimates unhold`), which genuinely works, but it is an
**undocumented** internal - shipping on it is a maintenance bet.

The decider is the abort path. A program that collapses and then errors leaves:
- preserve: caller sees the **full data** (auto-restored, scoped to the program)
- tempfile: caller sees **62 rows** - leaked
- frames: caller sees **62 rows** plus an orphan frame - leaked

### The architecture that matters most: collapse once, restore once

Because the by-key is identical at every margins site, one collapse can serve the whole
run. This was built and measured, not just proposed - and **nothing between the margins
sites needs the full data**, verified empirically rather than by inspection (the in-region
frame/svmat/histogram/ttest blocks and the six `est restore` calls all ran clean, and
wooldid's final `ereturn post ... esample()` matched the pristine sample exactly).

| config | N | pristine | collapse-once | saving |
|---|---|---:|---:|---:|
| admin, att+ES+jt | 29,667 | 0.989s | 0.592s | **40.1%** |
| admin x4 | 118,668 | 3.120s | 1.501s | **51.9%** |
| admin x16 | 474,672 | 11.659s | 5.083s | **56.4%** |
| admin + subgroups | 29,667 | 2.142s | 1.013s | **52.7%** |
| high-K (K=765, p=683) | 150,000 | 107.6s | 99.0s | **8.0%** |

The restore leg is free at every N; the collapse is the only real cost, and this pays it
once. `preserve`/`restore` must sit inline in the worker's own body - a helper that issues
`preserve` unwinds at its own `return`, which was demonstrated.

**One new gate clause fell out of this**: with `impvalues`, the `over(cluster)` sites need
`cluster` in the key. With a cluster *finer* than `i` the collapse is **BROKEN** (sign flip,
max absolute 2.2e+02) because the collapsed rows carry an arbitrary cluster value. wooldid's
own design makes coarse clusters the normal case - the cluster must contain `i`, and `_i` is
in the key - so a conforming cluster is automatically a function of the key. Two earlier
attempts to break this were *vacuous* and were caught by the subpop metric.

---

### Does a caller's outstanding `preserve` break this? No - measured

The author asked what happens if the user has issued `preserve` before calling wooldid.
Tested directly (`scratchpad/pres/pres.do`); `preserve` is **scoped per program
invocation** and nests cleanly.

| test | result |
|---|---|
| user `preserve` outstanding, program does `preserve`/`restore` | **works** - the program restores to the user's state (N=8), the user's own `restore` still returns N=10 |
| two levels of nested programs, each preserving | **works** - inner restores to 6, outer to 10 |
| program preserves and exits WITHOUT restoring | **auto-restored** on program exit (N=10) |
| program preserves then ERRORS | caller gets its data back (N=8), and the user's own `preserve` **survives** - their later `restore` returns rc 0, N=10 |
| two `preserve`s in the SAME scope | **rc 621** ("already preserved") |

So a program-internal `preserve` is safe alongside a user's, and safe under abort. The only
constraint is the last row: **wooldid must not have two outstanding `preserve`s in the same
program scope.** Today that is free - wooldid currently issues no `preserve` anywhere - but
it is a real constraint on any future second use inside `wooldid_worker`, and worth a guard
rather than an assumption.

One cost the speed numbers do not show: `preserve` was measured to be an in-memory copy
(~37 GB/s, no temp file appears), so a user's outstanding snapshot plus wooldid's own means
**two or three copies of the dataset resident at once**. Irrelevant at 320 MB; worth
considering for a user already near their memory ceiling, and an argument for holding the
collapsed window as short as the collapse-once design already does.

---

## 3. The break-even rule

**A single `att` is already enough - by 3.7x.** Break-even is **0.27 margins calls** at
29k and 0.17 at 119k. Everything above one call is comfortably profitable.

**Margins call counts, counted by instrumentation rather than guessed:**

| option set | calls |
|---|---|
| `att` | 1 |
| `att` + event study | 3 |
| + `jointtests` | 4 |
| + `subgroups` (2 cases) | 9 |
| all four estimands | 12 |

**On horizons, the premise is confirmed and two of my own assumptions were wrong.** Extra
horizons never add calls (3 calls at every window length). But `espostlength` raises
*neither* K *nor* p - K is a property of the panel, not the reporting window - while
`esprelength` does raise p. Net: **more horizons make the collapse better**, with break-even
falling from 0.23 to 0.17 as post-length goes 0 to 6.

**Where it is not worth it.** With `R = N/K` the rows-per-cell ratio: R >= 10 collapse;
4 <= R < 10 is grey; R < 4 do not. At R=1 the collapse is a measured **net loss**. But those
cutoffs are not universal - the crossover is set by `K*p*(R-1)`, not R alone. Since the
decision is made after the regression, the code knows N, K and p exactly and needs no
threshold:

```
collapse when   ncalls * 1e-07 * p * (N - K)  >  0.008 + 2e-07 * N
```

---

## 4. What the adversarial pass broke

It did not break the method. It broke the instrument and the formula.

**The residue bound was misstated, including by me.** The perturbation is
`ULP(|xb|) ~= |xb| * 2^-52`, *not* a fixed `2^-49` - that is merely `dta_s_admin`'s number.
An amplification ladder adding a constant to the outcome (absorbed by the fixed effects, so
every treatment coefficient is algebraically unchanged) tracked the predicted binade
exactly across three decades: `2^-49` at |xb|~12, `2^-43` at ~1012, `2^-39` at ~10012.

**Consequently "clean at rtol 1e-13" is a property of the test data, not of the method.**
It holds only while |xb| <~ 4500 (compare.py's absolute floor) and |xb|/|estimand| <~ 450
(the pure relative test). At the amplified rung the collapse produced 11 difference groups
at 1e-13.

**But the correctness conclusion survives, restated scale-free.** At the matched rung, the
row-order control - pristine wooldid with rows merely permuted, no collapse anywhere -
moves `b` by the *same* absolute quantum (1.819e-12) and a *larger* relative amount
(2.2e-10 vs 4.0e-12). The right statement is: **the collapse perturbs the answer no more
than permuting the input rows does.** Standard errors never amplify at all (V and se stay
at 1e-16 across three decades of |xb|), which is the zero-Jacobian argument holding
empirically.

**The break-even formula is wrong on an unfitted case: it has no fastnose term.** Measured
flip at N=5045 with `fastnose` - same data, same N/K/p, only the option differs - where the
formula says collapse and the truth is a **net loss of 0.003s**. It is 22-23x off at high-K
under `nose`. One-token fix: scale the saving term by ~0.3 when `nose` is in effect. Stakes
are small (~1% of runtime); the defect is in the claimed precision.

**Attacks that failed, reported as failures:** V and se amplification (SEs never scale);
plain factor controls `controls(i.region)`, never previously tested as a *plain* control
(clean, 0 groups at 1e-13). A full vacuity audit of the whole corpus traced every subpop
ratio to primary logs and found **no undisclosed vacuity** - the only degenerate run is one
the break-even pass itself disclosed, and it is timing-only.

---

## 5. Open items

- **A separate, pre-existing wooldid fragility surfaced**: with a large-mean outcome
  (constant +1e5 added to y), wooldid dies `rc 198 "invalid syntax"` after its own
  all-zero-VCE warning. Attribution was checked - no instrumented block ran - so this is
  wooldid's own path, not a collapse artifact. One data point; deserves its own look.
- **`cfx` sites** (`over(cfxlatticeplot)`) carry the same structural hazard as the
  `impvalues` cluster sites and were not exercised. Treat as gated until measured.
- **`ccc_absorb()`** is targetflag-free by inspection but was not run.
- Poisson remains unmeasured, as in the first report.
