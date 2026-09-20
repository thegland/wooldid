# Fix proposals — mechanism `LAT` (`lattice()` / `contreat()`)

Finding covered: **`fp-contreat-B2`** [high] — "on the same lattice fixture the ame point
estimate is 45x the truth, and no path gives the reader a usable warning: the default path
has no standard error at all and the one path that produces one reports se 0.44 on a 45x
miss."

## The proposals at a glance

| | proposal | recommend? | cost | moves a reported number? | needs a ruling? |
|---|---|---|---|---|---|
| **L1** | the all-zero-covariance cautions stop blaming clustering on runs that clustered nothing (4 edits) | **yes — ship** | zero | no (display text only) | no |
| **L2** | free colname predicate: partial-drop tier on R6's counters | **no — I built it and it cried wolf** (78 of 312 on a healthy design) | zero | no | — |
| **L2c** | free colname predicate: excess drops over the structural baseline | **no — measured, 204 vs 144, does not separate** | zero | no | — |
| **L3** | lattice support screen: treated (cell x bin) pieces with no more observations than parameters — **implemented and run: fires on the broken fixture, silent on a healthy one, on a non-lattice run, and under `safety(off)`** | **yes — ship with L1** | one `gegen` pass; **5.05s at 5,000,000 rows**, only on `contreat()`+`lattice()` runs | no | no |
| **L3b** | wire that screen to R1's existing ame scrub | **hold for the author** | same as L3 | **YES** — blanks ame `se/t/p/ci` and `e(V_es_ame*)` on affected lattice runs | **yes** |
| **L4** | two help edits (`sthlp:1059`; the lattice caveat paragraph) | yes | zero | no | no |
| **L5** | do **not** re-add the removed `regress` refit; stop steering readers to `customvce(robust)` | yes — a "do not fix" with new evidence | zero | no | no |

**The one-line story.** A lattice gives every treated (cell x bin) piece its own intercept and
its own dose slope; a piece with no more observations than parameters is fitted exactly, so
its residuals are identically zero and every residual-based variance estimator is blind in
exactly the directions the design's information matrix says carry the most variance. On the
finding's fixture that is 101 of 273 pieces holding 157 of 847 treated observations, all with
residuals below 5e-9 against an rmse of 1.40 — which is why `reghdfe` posts an all-zero
covariance, why `customvce(robust)` posts a standard error thirty times too small, and why
neither of the free colname-based screens I tried can tell it apart from the collinearity a
healthy lattice creates as a matter of algebra.

Everything below was **measured on the pinned build** and, where it is a code change,
**patched into a private copy and run**:

| artefact | path |
|---|---|
| private patched build (L1's 4 edits + L3, both tested) | `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-LAT/patched/wooldid.ado` |
| debug build (dump after the reghdfe call) | `.../probes/prop-LAT/ado/wooldid.ado` |
| probes + logs | `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-LAT/` |
| working notes | `.../probes/prop-LAT/NOTES.md` |

The fixture throughout is the finding's own: 12 states, 5-12 counties each (110 units), 20
years, staggered onsets 8/12/15, a **unit-level lognormal dose** with 10% time noise,
`i` = **state**, true dose effect exactly +0.4 on every treated post observation.
`wooldid y st year ttre, att contreat(dose) contreatcontrols(basic) lattice(all 3) robust`.
It is a shape the package's own SI-8 caution steers users into (">100 distinct `i` values"
-> use a coarser `i`).

---

## 0. What is actually broken — five measurements, three of them new

### 0.1 The missing standard errors are Stata's variance screen, and the specification is near-singular, not the code

A debug dump inserted immediately after the `reghdfe` call (`p01_dbg.log`) shows the
covariance is **already all zero when reghdfe returns**: `max|e(V)| = 0`, `K = 750`,
`rank = 501`, `df_a = 36`, `df_r = 1663`, `N = 2200`. Re-running the identical command
outside `quietly` (`p02_reghdfe.log`) shows what the quiet call swallows:

```
warning: variance matrix is nonsymmetric or highly singular.
```

Arms on that same command, same data (`p02_reghdfe.log`):

| arm | VCE | `max|e(V)|` | rank | note |
|---|---|---|---|---|
| A | `vce(robust)` | **0** | 501 | the screen fires |
| B | `vce(unadjusted)` | **3.18e+07** | 501 | no screen; se ~5600 on some cell-bin slope |
| C | `vce(cluster unit)` | **0** | 501 | the screen fires |
| D | lattice removed from the treatment term | 21.97 | 188 | healthy |
| F | A with `tol(1e-14)` | **0** | 501 | the partialling tolerance is **not** the cause |
| G | `regress` + explicit dummies + `robust` (the `customvce` path) | 421.2 | 381 | survives the screen |

So the `ado:12851` comment is right about the mechanism, and the design really is
near-singular: under iid the *largest coefficient variance in the fit is 3.18e7*. Three
candidate causes are ruled out by measurement: the partialling tolerance (arm F),
the absorbed dose control (`ccnoabsorb` unchanged in the finding's own repro), and lattice
*construction* (a user-supplied `lattice(latu)` variable does the same).

### 0.2 Both cautions blame clustering on a run that clustered nothing — verbatim

`p01_dbg.log`, the run above, `robust`, no `cluster()` anywhere:

> Caution: estimation returned an all-zero variance-covariance matrix. The issue is likely
> related to continuous treatment effect estimation coupled with use of **analytic clustered
> SEs, with clusters occuring at a high level of aggregation**. …

> Caution: the estimation command returned an all-zero variance-covariance matrix … This is
> a property of the specification, not a numerical accident: **see the cell-containment
> caution above** for the mechanism and for the alternatives (impvalues dispersion inference,
> **a finer-than-cell clustering**, or **conditional robust standard errors**).

There is no cell-containment caution above: R1 is gated on a cluster-type VCE
(`ado:1534-1544`), so on a bare-`robust` run it never ran and `amescrub` is 0. Yet
`ado:1567` keys its pointer on `"\`contreat'" != ""`, not on whether R1 fired. And the two
remedies offered are "a finer-than-cell clustering" (nothing was clustered) and "conditional
robust standard errors" (the user typed `robust`). `e(wdid_vcvstate)` does read `absent`,
so the machine-readable side is honest; only the prose is wrong.

This is the concrete instance of `known_bydesign` §5's "an all-zero `e(V)` under `robust` …
is **not** documented and **is** a finding".

### 0.3 NEW — the one path that reports a standard error reports one 30x too small

The finding and its refuter disagreed about whether the huge ame swings are bias or a
heavy-tailed sampling distribution. Both were arguing about draws in which **the design
changes**. I froze the design instead.

`p05_condmc3.do`: dose, lattice bins, onsets, units and unit effects are **fixed**; only the
iid N(0,1) idiosyncratic error is redrawn. The estimator is then an exactly linear function
of `y` with fixed weights, so the MC sd across replications **is** its conditional sampling
sd, and a heteroskedasticity-robust standard error is consistent for exactly that. No
oracle is needed.

**24 replications** (`p05_condmc3.log`, parsed by `parse_mc.py`, summary in
`mc_summary.txt`; the run was stopped at 24 of a planned 50 — the ratios had been flat since
rep 10):

| arm | n | mean b | **MC sd(b)** | **mean reported se** | **se/sd** | 95% CI coverage of 0.4 |
|---|---|---|---|---|---|---|
| L  `lattice(all 3) customvce(robust)` | 24 | 1.243 | **17.392** | **0.4418** | **0.0254** | 0.042 |
| U  `lattice(all 3) unconditionalse robust` | 24 | 1.243 | 17.392 | 8.8352 | 0.5080 | 0.792 |
| N  no lattice, `customvce(robust)` | 24 | 0.583 | 0.0803 | 0.1181 | **1.4703** | 0.708 |

**Read the control arm first.** Arm N's `se/sd` is 1.47, not 1.00, and that is the MC design,
not a defect: the frozen-design MC redraws only `eps`, while the residuals the robust
estimator squares also contain the *fixed* county effects that `i` = state cannot absorb. So
every arm's reported `se` is inflated by the same factor relative to this benchmark, and arm
N measures that factor at 1.47. The calibrated statement is therefore

> **arm L's standard error is `0.0254 / 1.4703` = 0.017 of the estimator's conditional
> sampling sd — a factor of about 58 — on a path where the same construction on the same
> data, one option away, is correct to within the calibration constant.**

Coverage is the wrong metric on arms N and U, because the unabsorbed county effects also put
a *fixed* bias in `b` (arm N's mean is 0.583 against a truth of 0.4). `se/sd` is not affected
by a fixed bias, which is why it is the statistic quoted. Arm L's coverage of 0.042 is quoted
only because it is so far from 0.95 that no calibration argument reaches it.

**The same MC on a design with support says the remedy works.** `p07b_healthy.log`, the
identical frozen-design construction on the healthy fixture (40 units per state, 15
replications):

| arm | n | mean b | MC sd(b) | mean se | se/sd |
|---|---|---|---|---|---|
| lattice(all 3) `customvce(robust)` | 15 | 0.674 | 0.1012 | 0.1684 | **1.66** |
| no lattice, `customvce(robust)` | 15 | 0.266 | 0.0171 | 0.0433 | **2.53** |

Relative to its own control the lattice arm is at 0.66 there, against **0.017** on the broken
fixture — a factor of 38 between two designs that differ in one parameter, units per state.
So "use a coarser lattice or a finer `i`" is not folklore: it is the fix, and it is
measurable. (Both arms sit above 1.0 in absolute terms for the reason given above — the
residuals carry the frozen county effects — which is exactly why the ratio to the control,
not the raw ratio, is the statistic.)

Note also what the MC **confirms of the refuter's correction**: mean `b` on arm L is 1.24
against a truth of 0.4 with an MC sd of 17.4, i.e. **no detectable bias**. The prober's "45x
the truth" really was one draw of a very heavy-tailed distribution. What the refuter could
not see, because their draws moved the design, is that the distribution is heavy-tailed
*conditional on the design too* — the tails are sampling noise the reported standard error
does not know about.

The point estimate is identical on arms L and U to every digit — only the variance differs.
**`customvce(robust)`, the only path that prints a standard error on this design, prints one
roughly thirty times smaller than the estimator's own conditional sampling sd, on data where
nothing but the noise moved.** `unconditionalse` is about twice too small. The control arm
with no lattice is in the right ballpark.

Consequences that follow directly:

* The finding's "se 0.44 on a 45x miss" is not a *conditional-variance statement about an
  ill-conditioned design*, which is the possibility the prober honestly flagged as the reason
  they marked it uncertain. It is simply wrong, and the frozen-design MC is what separates
  the two.
* The **default engine's behaviour is the safe one**: reporting nothing is correct here.
  The shipped caution steers the reader off it, towards the worst number in the run.
* The att-rescue refit removed by author ruling on 2026-09-08 (`ado:1551-1563`) would have
  manufactured exactly that number. The ruling looks better in hindsight than it did at the
  time; see L5.

### 0.4 NEW — the one screen positioned to see this is switched off by `lattice()`, and switching it back on is the wrong fix (measured in L2)

`wooldid_ccscreen` already scans `e(b)`'s column names for treated-cell terms dropped for
collinearity. Its interaction branch is `ado:11925`:

```stata
          if `wordcountproblem' > 3  & ("`lattice'" == "" | "`verbose'" == "verbose"){
```

— and inside it, `ado:11930`:

```stata
            if `catastrophe' == 0 & "`lattice'" != "" dis as text  "If using a lattice, this is likely innocuous."
```

So with a lattice active and `verbose` off the branch is skipped entirely ("lattice drops
there are expected", `ado:11889-11891`).

**That suppression turns out to be right, and I proved it by breaking it** (L2): a lattice
drops one term per treated cell as a matter of algebra, so the per-term caution really would
fire on healthy designs — measured, 78 of 312 dropped terms on a fixture whose covariance is
perfectly fine. What is missing is not the per-term caution. It is a statistic that separates
the algebra from the pathology, and `wooldid` computes nothing of the sort today. L2c and L3
are two candidates, one free and one costing a pass.

The sibling screen R6 (`ado:11958-11976`) already counts the same terms into `amelive` and
`amedrop`, free, in a colname loop with no data pass — but it only acts when `amelive == 0`
(*every* dose slope dropped). The partial case, which is this case, falls through silently.

### 0.5 NEW — the degeneracy is visible in the data alone, before any estimation

`p03_degen.log`, pure data work on the realized estimation sample, treated rows only:

| | combos | with < 2 distinct dose values | with exactly 2 obs | median within-combo sd(dose) |
|---|---|---|---|---|
| `lattice(all 3)`: (cell x bin) | 273 | **45 (16.5%)** | 56 | 0.140 |
| no lattice: cell | 94 | **0** | 0 | 0.663 |

(sample sd(dose) = 0.623.) A piece with fewer than two distinct dose values has no slope of
its own; the engine drops it and the reported ame averages only what survived. The statistic
separates the broken arm from the healthy arm perfectly, and costs one `gegen … by()` pass.

**One number, two samples — read this once and the rest of the document is consistent.** The
counts here (273 pieces, 45 + 56 = **101** saturated, **157** treated observations in them)
are taken on the user `if` sample, which is what I had offline. The screen as implemented runs
on `esamp`, the **realized** estimation sample, and reports the same quantity restricted to it:
**228 pieces, 83 saturated, 130 observations** (`p10_l3test.log`). `esamp` is the correct
sample — it is what the engine actually fitted — so 83/228 is the shipped number and 101/273
is its offline shadow. Nothing else differs between them.

### 0.6 NEW — the mechanism, closed: **saturated cell-by-bin pieces have identically zero residuals**

This is the piece the finding, the refuter and the third lens all missed, and it makes every
other measurement above fall into one story.

A lattice gives each treated (cell x bin) piece **its own intercept and its own dose slope**:
two free parameters with `contreatpoly(1)`. A piece holding **two or fewer observations is
therefore fitted exactly.** Its residuals are not small — they are zero.

Measured on the broken fixture (`p08_predicate.log`, `regress` with explicit dummies,
`vce(robust)`, model rmse 1.398):

| treated rows | n | mean abs residual | max abs residual |
|---|---|---|---|
| in pieces with <= 2 distinct dose values | **157** | **6.57e-11** | **4.96e-09** |
| in pieces with >= 3 | 690 | 0.741 | 4.18 |
| untreated / reference rows | 1353 | 1.110 | 4.55 |

157 is not an approximation: `p03_degen.log`'s size distribution has 45 pieces of one
observation and 56 of two, and (45 pieces of one) + (56 pieces of two) = **157** exactly. **18.5% of the treated
observations contribute a residual of zero.**

Everything else follows:

* **Every residual-based variance estimator** — HC-robust, cluster-robust, the bootstrap's
  own analytic leg — assigns **zero** variance to precisely the directions in which the iid
  bread says the variance is up to **3.18e7** (`p02_reghdfe.log` arm B). A matrix holding
  both is what Stata calls "nonsymmetric or highly singular", so `reghdfe`'s post zeroes it
  wholesale (arm A) while `regress` posts a small rank-deficient one (arm G, `e(rank)` 381 of
  537) — the 0.44 of §0.3.
* **This is R1's mechanism, one geometry finer.** `ado:12840-12847` already states it for
  cell containment: "the plug-in CRVE cannot see within-cell dose variation … the ame
  standard errors it would deliver collapse toward the variance of the absorbed baseline dose
  control." R1 fires on a cell wholly inside one *cluster*. Here the cell is split by the
  *lattice* into pieces that are wholly determined by their own parameters. Same disease, and
  the cure the package already ships for it is the R1 scrub — which is why L3b is an
  extension of an accepted ruling rather than a new judgment call.
* **It gives the exact predicate**, with no threshold to argue about:
  `count(treated observations in a (cell, bin) piece) <= 1 + contreatpoly`.
  On the broken fixture that is 101 of 273 pieces; on the healthy one, 0 of 234
  (`p07b_healthy.log`: minimum 7 distinct dose values per piece, median 13).

---

# L1 — Make the all-zero-covariance cautions state the cause that actually applied

**RECOMMENDED. Free. Moves no reported number. Changes displayed text.**

Four edits, all in the same causal chain: the R4 caution and its helper have no idea whether
anything was clustered, so they say it was.

### Edit 1 of 4 — hoist the cluster extraction out of the R1 gate (`ado:1533-1541`)

```stata
  local amescrub 0
  if "`contreat'" != "" & "`impvalues'" == "" & "`unconditionalse'" == "" & "`fastnose'" == "" {
    local amescrubclus
    if "`cluster'" != "" local amescrubclus `cluster'
    local amecvw = word(`"`customvce'"', 1)
    if "`amescrubclus'" == "" & strlen("`amecvw'") >= 2 & strpos("cluster", "`amecvw'") == 1 {
      local amescrubclus : subinstr local customvce "`amecvw'" "", word
    }
    if "`amescrubclus'" != "" {
```

becomes

```stata
  * the cluster variable(s) actually in play, whichever syntax asked for them. HOISTED
  * out of the R1 gate below because R4's all-zero-covariance caution needs to know
  * whether anything was clustered at all: it used to blame "analytic clustered SEs ...
  * at a high level of aggregation" on runs that clustered nothing (contreat() +
  * lattice() + robust zeroes the covariance with no cluster in the specification).
  * R1's own gate, call and behaviour are unchanged.
  local amescrubclus
  if "`cluster'" != "" local amescrubclus `cluster'
  local amecvw = word(`"`customvce'"', 1)
  if "`amescrubclus'" == "" & strlen("`amecvw'") >= 2 & strpos("cluster", "`amecvw'") == 1 {
    local amescrubclus : subinstr local customvce "`amecvw'" "", word
  }
  local amescrub 0
  if "`contreat'" != "" & "`impvalues'" == "" & "`unconditionalse'" == "" & "`fastnose'" == "" {
    if "`amescrubclus'" != "" {
```

Pure motion: the same three statements, the same gate around the same call. `amescrubclus`
is now set on every run, which is what edits 2-4 need. Nothing else in the worker reads that
name (checked: two occurrences, both in this block).

### Edit 2 of 4 — tell the screen what the VCE was (`ado:1565`)

```stata
    wooldid_amescreen , mode(vzero) r1fired(`amescrub') contreat(`contreat')
```
becomes
```stata
    wooldid_amescreen , mode(vzero) r1fired(`amescrub') contreat(`contreat') clusvars(`amescrubclus') latflag(`=cond("`latticevar'" != "", "1", "0")')
```

`latticevar` is the worker's own local (`ado:1023`, `ado:1084`): the tempvar for a generated
lattice, or the user's variable name for `lattice(myvar)`. Empty when no lattice.

### Edit 3 of 4 — the vzero caution branches on whether anything was clustered (`ado:12982-12985`, plus `latflag` on the syntax line `ado:12863`)

`ado:12863` gains one optional token:
```stata
  syntax , mode(string) [icell(string) … r1fired(string) latflag(string)]
```

and the body becomes:

```stata
    if `vzero' & "`r1fired'" != "1" {
      if "`contreat'" != "" & "`clusvars'" != "" dis as text "Caution: estimation returned an all-zero variance-covariance matrix. The issue is likely related to continuous treatment effect estimation coupled with use of analytic clustered SEs, with clusters occuring at a high level of aggregation. Consider switching to use of bootstrap SEs, IM p-values, or randomization inference if finer clustering levels are non-viable."
      else if "`clusvars'" != "" dis as text "Caution: estimation returned an all-zero variance-covariance matrix. This usually reflects a degenerate clustering (very few clusters, or clusters at a high level of aggregation relative to the design). Consider switching to bootstrap SEs, IM p-values, or randomization inference."
      else dis as text "Caution: estimation returned an all-zero variance-covariance matrix. Nothing in this run was clustered, so this is not a clustering problem: Stata replaces a variance matrix it judges nonsymmetric or highly singular with zeros, which is what a near-singular design does to the sandwich estimator.`=cond("`latflag'" == "1", " With lattice() the usual cause is treated cohort-by-period cells that hold too few distinct values of the continuous treatment once the lattice splits them: a cell-by-bin piece with no more observations than it has parameters is fitted exactly, so its residuals are identically zero and no residual-based variance estimator can see it.", "")' Bootstrap SEs and randomization inference (permute()) do not use this matrix; impvalues does not either, but it needs a cluster variable.`=cond("`latflag'" == "1", " A coarser lattice (fewer bins), a base-period lattice type, or a coarser i also restores it.", "")'"
    }
```

The first two branches are the shipped text, now correctly gated: the clustered-contreat
sentence only fires when something really was clustered. The third is new and is the case
that previously got branch one.

### Edit 4 of 4 — the R4 caution stops pointing at a caution that never printed (`ado:1566-1568`)

```stata
    if r(vzero) == 1 {
      dis as text "Caution: the estimation command returned an all-zero variance-covariance matrix, … This is a property of the specification, not a numerical accident: `=cond("`contreat'" != "", "see the cell-containment caution above for the mechanism and for the alternatives", "the alternatives are")' (impvalues dispersion inference, a finer-than-cell clustering, or conditional robust standard errors)."
    }
```
becomes
```stata
    if r(vzero) == 1 {
      * the "see the caution above" pointer keys on whether R1 ACTUALLY FIRED, not on
      * whether contreat() was specified: R1 needs a cluster-type VCE, so on a bare-robust
      * contreat run the old text sent the reader to a caution that was never printed, then
      * offered "a finer-than-cell clustering" to a user who had clustered nothing and
      * "conditional robust standard errors" to a user who had just asked for them.
      * impvalues is offered only where it is legal: it refuses without cluster() (ado:4642)
      local wdr4alt "bootstrap() or permute()"
      if "`amescrubclus'" != "" local wdr4alt "impvalues dispersion inference, `wdr4alt', or a finer-than-cell clustering"
      if "`amescrubclus'" == "" & "`robust'" == "" & "`customvce'" == "" & "`unconditionalse'" == "" local wdr4alt "`wdr4alt', or conditional robust standard errors"
      dis as text "Caution: the estimation command returned an all-zero variance-covariance matrix, so no standard error in this run is identified. The reported estimates stand; every standard error, test statistic, p-value and interval is missing. This is a property of the specification, not a numerical accident: `=cond(`amescrub' == 1, "see the cell-containment caution above for the mechanism and for the alternatives", "the alternatives are")' (`wdr4alt')."
    }
```

Note the deliberate change of substance, not only of grammar: **`bootstrap()`/`permute()`
replace "conditional robust standard errors" at the head of the remedy list, and the
conditional-robust remedy is offered only to a run that has not already asked for a
conditional VCE.** §0.3 is the reason — on the design that produces this caution, the
conditional robust standard error is the worst number available.

### Cost

Zero. Two `strpos`/`word` calls that already ran, now running unconditionally instead of
under a gate; three locals; no data pass, no matrix work. Unmeasurable on any panel size.

### Does it move a reported number?

**No.** No `e()` result, no matrix, no estimate, no standard error changes. It changes
**displayed caution text**, so any baseline that captures the log verbatim needs re-blessing
— text only, and only for runs that hit the all-zero-covariance path.

### Blast radius

Three display branches and one macro hoist. R1's own predicate, gate and call are untouched,
which I checked by running the clustered arm (below). The only behaviour any consumer can see
is which caution string appears.

### Tested — before / after

Patched build `.../probes/prop-LAT/patched/wooldid.ado`, `p06_patchtest.log`:

| arm | expectation |
|---|---|
| B `lattice(all 3) robust` | the new non-cluster text, naming the near-singular design and the lattice; no mention of clustering |
| C no lattice, `robust` | unchanged, silent |
| D `lattice(all 3) cluster(st)` | the **original** clustered text survives verbatim (no regression) |
| E healthy lattice, 40 units per state | silent |

Results are quoted in §L6.

### Alternatives rejected

* **Leave the text and document it in the help.** The text is not merely incomplete, it is
  false on the run that prints it ("coupled with use of analytic clustered SEs" when no
  cluster was requested) and it points at a caution that is not on the screen. A help entry
  does not stop a user acting on a false sentence in their own log.
* **Detect the near-singularity before estimating and refuse.** No cheap predicate for
  "Stata's variance screen will fire" exists; arm B shows the same design has a perfectly
  finite iid covariance, so it is not even a property of `X` alone. Refusing would also
  destroy runs whose point estimates the user legitimately wants.
* **Suppress the second caution when the first one already said it.** They say different
  things (one names the mechanism, one names the consequence for every downstream object)
  and the project already has a defer-to rule between them (`r1fired`). Extending the same
  rule is cheaper than merging them.

---

# L2 — Let the partial-drop case out of the R6 counter — **I BUILT IT, RAN IT, AND IT CRIED WOLF**

**Status: rejected in its simple form, on my own measurement. Keep reading — the site is
right, the predicate is not, and the salvage is in L2c.**

### What I built

`wooldid_ccscreen` already scans `e(b)`'s colnames for treated-cell terms dropped for
collinearity, and R6 (`ado:11958-11976`) already counts them into `amelive`/`amedrop` with no
data pass at all — but acts only when `amelive == 0`. The partial case falls through. The
obvious free fix is an `else` branch. I added exactly that, gated on `lattice()` so it lands
precisely in the hole `ado:11925` leaves:

```stata
    else if `amedrop' > 0 & "`lattice'" != "" {
      local amedtot = `amelive' + `amedrop'
      dis as text "Caution: `amedrop' of `amedtot' treated cohort-period dose-slope terms for `contreat' were dropped for multicollinearity …"
    }
```

### What it did — the number that kills it

`p06_patchtest.log`:

| arm | the new caution says |
|---|---|
| B — the broken fixture (5-12 units per state) | "**114 of 306** … were dropped" |
| E — a healthy fixture, 40 units per state, `N` = 9600 | "**78 of 312** … were dropped" |

Arm E is *not* broken: its `reghdfe` covariance is finite (`max|e(V)|` = 91.5, rank 468 of
640 — `p07_healthy.log`), it prints a complete table with standard errors, and no all-zero
caution appears anywhere in it. Yet the caution fires, and loudly.

The reason is structural and the shipped code already knew it: `ado:11930` says *"If using a
lattice, this is likely innocuous."* With `L` bins the factorial gives each treated cell
`L + 1` dose terms of which exactly one is redundant, so a lattice run drops **one term per
treated cell as a matter of algebra** — 78 cells, 78 drops, 312 = 78 x 4 terms. The broken
fixture's 114 is that same structural baseline (~76) **plus ~38 real ones**. 25.0% against
37.3% is not a margin a default-on caution can live inside.

**So `ado:11925`'s lattice suppression is vindicated, not a bug.** That is worth saying
plainly, because on the reading of §0.4 alone it looks like the defect. It is not. The defect
is that nothing replaced it.

### Two things survive the failure

* The **site** is right: this is the one screen positioned to see the geometry, and it is free.
* The **statistic** is wrong: "terms the engine dropped" mixes algebra with pathology. What
  discriminates is *support* — see L3 (measured: 45 vs 0) and L2c (free, measured in §L6).

---

# L2c — a second free predicate (excess drops over the structural baseline) — **ALSO MEASURED, ALSO REJECTED**

Having lost the raw drop count to algebra, the obvious repair is to subtract the algebra: the
structural baseline is one dropped term per treated cell, and the number of treated cells is
already in the same colname loop (the `cell#tf#c.contreat` term carrying no lattice value).
Predicate: `amedrop > ncellmain`, cost one `strpos` per already-visited token, no data pass.

I implemented the counter and ran it on both designs (`p08_predicate.do`, `r6count`):

| design | live | dropped | total | cell-main terms | **excess** | drop share |
|---|---|---|---|---|---|---|
| broken (5-12 units/state) | 228 | 348 | 576 | 144 | **204** | 60.4% |
| healthy (40 units/state) | 234 | 270 | 504 | 126 | **144** | 53.6% |

**It does not discriminate.** 204 against 144, and the share 60% against 54%. Nor does the
raw share in `wooldid`'s own parameterization (37.3% against 25.0%, §L2). A default-on caution
cannot live in that margin.

*(These token counts are larger than the ones `wooldid` prints because I called `reghdfe`
directly, with Stata's default factor-variable bases rather than the bases `wooldid` sets —
576 tokens against 306. The comparison between the two designs is like-for-like, which is what
the test needed.)*

**So both free predicates are rejected, on measurement.** What the colnames record is what
the engine dropped, and the engine drops for algebraic reasons and pathological ones with the
same notation. The statistic that separates them is not in `e(b)`; it is in the data, and it
is L3.

# L3 — **RECOMMENDED**: the lattice support screen, one data pass, measured price and measured discrimination

This is the proposal I would ship alongside L1. It is the only predicate I tested that
separates the broken design from the healthy one, and §0.6 gives it an exact form with no
threshold to argue about.

### The edit — a new mode in `wooldid_amescreen`, called from the worker after `ado:1545`

Call site, immediately after the R1 block, gated so it runs only where it can matter:

```stata
  * LATTICE SUPPORT (see wooldid_amescreen mode latsupport). Only with contreat() +
  * lattice(), and never under safety(off) - one gegen pass over the treated rows.
  if "`contreat'" != "" & "`latticevar'" != "" & "`safety'" != "off" & "`wdrepmode'" == "" {
    wooldid_amescreen , mode(latsupport) icell(`_i') tcell(`_t') tflag(`targetflag') esamp(`esamp') latvar(`latticevar') contreat(`contreat') poly(`contreatpoly')
  }
```

and the mode itself:

```stata
  if "`mode'" == "latsupport" {
    * A lattice gives every treated (i, relative-period) cell x lattice-bin piece its OWN
    * intercept and its OWN dose slope: 1 + contreatpoly free parameters. A piece with no
    * more observations than parameters is fitted EXACTLY - its residuals are identically
    * zero - so every residual-based variance estimator (robust, cluster, the bootstrap's
    * analytic leg) sees no variation there at all, while the design's information matrix
    * says those directions carry enormous variance. That combination is what Stata's
    * ereturn-post screen calls "nonsymmetric or highly singular", and it is R1's
    * cell-containment mechanism one geometry finer. MEASURED on a 110-unit x 20-year panel
    * with a unit-level dose and lattice(all 3): 101 of 273 pieces are saturated, they hold
    * 157 of the 847 treated observations, and every one of those 157 residuals is < 5e-9
    * against a model rmse of 1.40.
    local lspar = 1 + `poly'
    tempvar lsn lstag
    qui gegen double `lsn' = count(`contreat') if `esamp' & `tflag' == 1, by(`icell' `tcell' `latvar')
    qui gegen byte `lstag' = tag(`icell' `tcell' `latvar') if `esamp' & `tflag' == 1
    qui count if `lstag' == 1
    local lsall = r(N)
    qui count if `lstag' == 1 & `lsn' <= `lspar'
    local lsdead = r(N)
    qui count if `esamp' & `tflag' == 1 & `lsn' <= `lspar'
    local lsobs = r(N)
    if `lsdead' > 0 {
      dis as text "Caution: lattice() splits the treated cohort-by-period cells into `lsall' cell-by-bin pieces, and `lsdead' of them hold `lspar' or fewer observations - no more than the intercept and slope that piece is given, so those pieces are fitted exactly and the `lsobs' treated observation(s) in them have residuals of zero. No residual-based standard error (robust, cluster, or the bootstrap's analytic leg) can see any variation there, so the analytic standard errors on the average marginal effects can be far smaller than the estimates' sampling variability, and the estimation command may return an all-zero covariance instead. Use fewer lattice bins, a coarser `i', or a base-period lattice type; or take inference from bootstrap() or permute()."
    }
    return scalar ndead = `lsdead'
    return scalar nobs = `lsobs'
    return scalar ncombo = `lsall'
    exit
  }
```

### Tested in the build, not only on paper

I implemented exactly the code above in the patched copy and ran it (`p10_l3test.log`):

| arm | result |
|---|---|
| broken fixture, `lattice(all 3) robust` | **fires**: "…into **228** cell-by-bin pieces, and **83** of them hold 2 or fewer observations … the **130** treated observation(s) in them have residuals of zero…" — run time **4.52s** |
| healthy fixture (40 units/state, N = 9600), same options | **silent** — run time 10.59s |
| broken fixture, no `lattice()` | **silent** |
| broken fixture, `lattice(all 3) robust safety(off)` | **silent** (the `safety(off)` gate works) |

The in-run counts (228 / 83 / 130) are smaller than the offline ones below (273 / 101 / 157)
for a reason worth keeping: the screen runs on `esamp`, the **realized** estimation sample,
while my offline count used the `if` sample. `esamp` is the right one — it is what the engine
actually fitted, after `reghdfe` dropped what it drops — and it is why the call site is after
the fit rather than before it.

Run time on the broken fixture is 4.52s against 4.78s for the same build without the screen
(§L6): the added pass is inside the noise on a 2,200-row panel, which is what the 5M-row
price below is for.

### Measured discrimination

| design | pieces | **saturated pieces** | treated obs in them | reghdfe `max|e(V)|` | reported ame se / conditional MC sd |
|---|---|---|---|---|---|
| broken, 5-12 units per state | 273 | **101** | 157 | **0** (screen fired) | 0.0254 |
| healthy, 40 units per state | 234 | **0** | 0 | 91.5 (fine) | 1.66 |

The healthy fixture's smallest piece holds 7 distinct dose values, median 13
(`p07b_healthy.log`). The screen is silent there, and silent on every non-lattice run by
construction. That is the discrimination the two free predicates could not manage.

### Cost — priced, not asserted

`p06_patchtest.log` arm T, a synthetic **5,000,000-row** frame with 60,000 (cell x bin)
groups: two `gegen … by()` calls plus three `count`s take **5.05 seconds** (`timer list 2`).
The version above is cheaper than the one timed — `count()` instead of `nunique()`, which
does not have to hash the values — so 5 seconds is an upper bound at 5M rows, and well under
a second on the panels most users bring. Three things bound it further:

* it runs **only** when both `contreat()` and `lattice()` are specified;
* it is skipped under `safety(off)`, like the screens around it;
* the `"\`wdrepmode'" == ""` gate keeps it out of every bootstrap and permutation replication,
  so it is once per run, not once per draw.

For scale: the R1 `detect` screen already standing at `ado:1542` does one `gegen nunique` plus
one `gegen tag` **per cluster variable** on every clustered `contreat` run. This is the same
cost class the package has already accepted, on a strictly narrower set of runs.

### Does it move a reported number?

**No** as written — it prints and returns. Wiring it to the scrub is L3b, and that does.

### Blast radius

`contreat()` + `lattice()` runs with at least one saturated piece. Nothing else in the package
reads `r(ndead)`/`r(nobs)`/`r(ncombo)` unless L3b is taken.

### Alternatives rejected

* **The two free colname predicates** (L2, L2c): built, run, measured, do not discriminate.
* **A relative-sd threshold (`sd < 1e-3 * sample sd`), copying the `bycoh` screen.** It asks
  the wrong question. An sd threshold asks *is this piece's slope well conditioned?*; §0.6
  says the damage is done by *saturation*, and a two-observation piece with a perfectly
  healthy within-piece sd is still fitted exactly and still contributes residuals of zero.
  Measured on the broken fixture, the threshold catches 46 of the 273 pieces — essentially the
  singletons only — and misses all 56 two-observation pieces, whose within-piece sd is
  0.02-0.10 against a sample sd of 0.62 and which hold 112 of the 157 zero-residual
  observations. The count predicate asks *is this piece saturated?*, which is exactly what no
  residual-based variance estimator can see.
* **Refusing (`error 459`) instead of cautioning.** The point estimates remain the estimates
  the user asked for, and §0.3 shows they are unbiased; it is the standard errors that are
  worthless. Blanking those is L3b; refusing the run is disproportionate and would break every
  existing lattice do-file that has one thin cell.
* **Doing it before estimation instead of after.** It must run on the *realized* estimation
  sample: `reghdfe` drops singleton groups and the screen's whole point is what the engine
  actually fitted. `esamp` exists at `ado:1505`, which is why the call site is where it is.

# L3b — (needs a ruling) withhold the ame standard errors when the support screen fires

**One extra line on L2c or L3. MOVES REPORTED NUMBERS. Do not ship without the author.**

Add to whichever screen the author picks:

```stata
      c_local amescrub 1
```

The R1 machinery then nulls `se/t/p/ci` on every ame-family capture, withholds
`e(V_es_ame*)` while `e(b_es_ame*)` still posts, and skips the ame-side `V_es` tests — what
it already does for cell containment (R1) and what R6 does one branch above for the
total-drop case.

**§0.6 changes what this is.** Before the residual measurement, wiring the screen to the
scrub looked like a new judgment call about how bad is bad enough. It is not: a saturated
cell-by-bin piece has identically zero residuals, so the plug-in variance estimator is blind
in exactly the way `ado:12840-12847` already describes for cell containment — the mechanism
R1 was built for and the author already ruled on. L3b is R1 applied to the geometry the
lattice creates, with the same consequences (`se/t/p/ci` nulled, `e(V_es_ame*)` withheld,
`e(b_es_ame*)` kept) and the same carve-outs.

**Why it deserves a ruling rather than a caution.** §0.3 measured the number this would
blank: a reported `se` of 0.4418 against a conditional sampling sd of 17.39, i.e. a
t-statistic of forty on an estimate whose own sampling distribution has a standard deviation
of seventeen. That is a *silently wrong number on realistic input*, the top of the
finding-value scale, and a caution the user may scroll past leaves it on the screen. The
package's own precedent is to blank rather than caution when the predicate is structural.

**What it costs.** Every `lattice()` run whose screen fires loses its ame standard errors, t,
p, intervals and `e(V_es_ame*)` — a behaviour change needing a re-bless of any lattice
baseline, and it will bite runs where the support failure is marginal and the standard error
is probably fine. Three dials the author would need to set:

* whether to fire at one saturated piece or above a share of them (I have no measurement
  that says where the harm becomes material; the sweep in T2 is what would produce one);
* whether `unconditionalse` is carved out, as R1 carves it out today. Measured, it is the
  least-bad reporting path on this design (`se/sd` 0.51 against arm L's 0.025) but it is not
  right either, so carving it out preserves a number that is still twice too small.

**My own recommendation** is to ship L1 + (L2c or L3) as cautions first, and hold L3b for a
ruling, because the scrub's blast radius over existing lattice do-files is unknown to me: I
have measured two designs, not the population of real ones.

---

# L4 — Help edits

Two places in `wooldid.sthlp` say something this design contradicts.

1. **`sthlp:1059`** (inside `timetrends`) ends: *"The same specification under `robust` is
   unaffected."* True of the cohort-trend case it describes, but it is the only place in the
   help where an all-zero covariance is discussed, and a reader takes away "`robust` is safe".
   Add, in the `lattice()` discussion rather than here, one sentence:

   > A `lattice()` on a continuous treatment can also leave the covariance non-identified,
   > and there it does **not** depend on `cluster()`: where the lattice splits treated
   > cohort-period cells into pieces that hold too few distinct values of the continuous
   > treatment, the design is near-singular and the estimation command may return an
   > all-zero covariance under `robust` as readily as under `cluster()`.

2. **The lattice caveat paragraph** (`sthlp:~1006-1010`, "When a specified lattice has many
   groups and almost all of the effect of Contreat_jt is captured by the intercept lattice
   terms alone, the identification of marginal effects … can be threatened", with the
   `Y = -Contreat + 100*floor(Contreat/10)` example). That paragraph describes a
   *functional-form* threat: effect loaded onto bin intercepts. The threat measured here is a
   *support* threat and has nothing to do with the shape of f() — on the fixture the true
   dose response is exactly linear with no bin structure at all. Add:

   > A second, distinct requirement is **support**: each treated cohort-period cell must
   > still hold at least two distinct values of the continuous treatment *within each lattice
   > bin it occupies*, since that is the only variation identifying that piece's slope. When
   > the continuous treatment is close to constant within units and `i` is coarser than the
   > unit, a lattice can cut cells into pieces of one or two observations; the pieces with
   > fewer than two distinct values are dropped for collinearity, the reported average
   > marginal effect then averages only the pieces that survived, and the analytic standard
   > errors on such a design can be very much smaller than the estimates' sampling
   > variability. `wooldid` cautions when it sees this. A useful rule of thumb is that the
   > number of units contributing to a treated cell should be several times the number of
   > lattice bins.

Cost zero, moves nothing, and the help is where the "ordinary applied choice" in the finding
should have been steered.

---

# L5 — Do **not** re-add the `regress` refit, and stop steering readers to `customvce(robust)`

**This is a "do not fix" with new evidence.**

The obvious reading of the finding is "the default path has no standard error; get it one".
`ado:1551-1563` records that exactly such a rescue existed and was removed on 2026-09-08:
a refit of the same specification through `regress` with explicit dummies, removed for its
cost (a dense solve over hundreds of collinear dummy columns) and its numerical noise
(1.6e-8 relative under row reordering, against 1e-13 for the partialled fit).

§0.3 adds the argument the ruling did not have: **on this design that refit's answer is the
one measured at `se/sd ≈ 0.03`.** Arm G of §0.1 is that refit — `regress` with explicit
dummies and `robust` — and arm L of §0.3 is what wooldid reports from it. Re-adding it would
convert a run that currently says "no standard error is identified" into a run that reports a
t-statistic of forty on an estimate whose conditional sampling sd is thirty times the
reported standard error. The 2026-09-08 removal should stand, and this measurement belongs in
the record next to it.

The same measurement is why L1's edit 4 demotes "conditional robust standard errors" from the
remedy list on runs that already have a conditional VCE, and puts `bootstrap()`/`permute()`
first.

**What I did not test:** whether `bootstrap(standard)` or `permute()` actually deliver
calibrated inference on this geometry. They are the remedies both the shipped text and my
replacement offer, and neither has been measured here. `known_bydesign` §7 already records
`bootstrap(standard)` over-covering at `cluster() == i` on a different design. If the author
wants the remedy list to be evidence-backed rather than conventional, that is the run to
commission: the same frozen-design MC as `p05_condmc3.do`, with `bootstrap(standard)` and
`permute(onsets)` arms. I would not assume it passes.

---

# L6 — What the patched build actually did

Patched build: `.../probes/prop-LAT/patched/wooldid.ado` (L1's four edits plus the L2 tier
as originally drafted). Probe `p06_patchtest.do`, log `p06_patchtest.log`.

### Arm B — `lattice(all 3) robust`, the finding's own run

Before (pinned build, `p01_dbg.log`):

> Caution: estimation returned an all-zero variance-covariance matrix. The issue is likely
> related to continuous treatment effect estimation coupled with use of **analytic clustered
> SEs, with clusters occuring at a high level of aggregation**…
>
> Caution: … not a numerical accident: **see the cell-containment caution above** … (impvalues
> dispersion inference, **a finer-than-cell clustering**, or **conditional robust standard
> errors**).

After (patched):

> Caution: estimation returned an all-zero variance-covariance matrix. **Nothing in this run
> was clustered, so this is not a clustering problem**: Stata replaces a variance matrix it
> judges nonsymmetric or highly singular with zeros, which is what a near-singular design
> does to the sandwich estimator. With lattice() the usual cause is treated cohort-by-period
> cells that hold too few distinct values of the continuous treatment once the lattice splits
> them… Bootstrap SEs and randomization inference (permute()) do not use this matrix; impvalues does not either, but it needs a cluster variable. A
> coarser lattice (fewer bins), a base-period lattice type, or a coarser i also restores it.
>
> Caution: … not a numerical accident: **the alternatives are** (impvalues dispersion
> inference, **bootstrap() or permute()**).

`rc = 0`, the table, `e()` and every number identical. Run time 4.78s (`timer list 1`) — the
same run on the pinned build is the same to within noise, as it must be: nothing added
touches data.

### Arm D — `lattice(all 3) cluster(st)`: no regression

The original cell-containment caution still fires verbatim ("with standard errors clustered
on st, 78 treated i-t cell(s) lie completely within a single cluster…"), `r1fired` still
suppresses the vzero caution, and the R4 caution still points the reader at it. The hoist in
edit 1 changed nothing R1 can see.

### Arm C — no lattice, `robust`: silent, unchanged.

### Arm E — the healthy lattice: the L2 tier fired anyway

"Caution: **78 of 312** treated cohort-period dose-slope terms … were dropped", on a run whose
covariance is fine and whose table is complete. This is the measurement that rejects L2; it is
written up there.

### Final text, re-verified after the `impvalues` correction

`p09_finaltext.log`, patched build, arm B, verbatim:

> Caution: estimation returned an all-zero variance-covariance matrix. Nothing in this run was
> clustered, so this is not a clustering problem: Stata replaces a variance matrix it judges
> nonsymmetric or highly singular with zeros, which is what a near-singular design does to the
> sandwich estimator. With lattice() the usual cause is treated cohort-by-period cells that
> hold too few distinct values of the continuous treatment once the lattice splits them: a
> cell-by-bin piece with no more observations than it has parameters is fitted exactly, so its
> residuals are identically zero and no residual-based variance estimator can see it.
> Bootstrap SEs and randomization inference (permute()) do not use this matrix; impvalues does
> not either, but it needs a cluster variable. A coarser lattice (fewer bins), a base-period
> lattice type, or a coarser i also restores it.
>
> Caution: the estimation command returned an all-zero variance-covariance matrix, so no
> standard error in this run is identified. … the alternatives are (bootstrap() or permute()).

and arm D on the same build still opens with "Caution: with standard errors clustered on st,
78 treated i-t cell(s) lie completely within a single cluster …", unchanged.

`impvalues` is named only in the clustered branch because it refuses without `cluster()`
(`ado:4642`, `error 198`). The shipped text offered it to unclustered runs; that is a second,
smaller instance of the same defect and it is fixed by the same edit.

### Runtime

Arm B: 4.78s on the patched build (`timer list 1`), against the pinned build's same run —
no data work was added by L1, so the two are identical to noise. L3's screen was priced
separately at 5.05s per 5,000,000 rows (arm T), and is not in the patched build.


---

# L7 — Tests that would prove these fixes correct

None of them needs an oracle.

### T1 — the caution branches (L1). Log-grep, no numbers.
Four runs on one fixture, greping the log:

| run | must contain | must NOT contain |
|---|---|---|
| `lattice(all 3) robust` | "Nothing in this run was clustered" | "clustered SEs", "cell-containment caution above" |
| `lattice(all 3) cluster(st)` | "lie completely within a single cluster" | "Nothing in this run was clustered" |
| `timetrends cluster(i)` (the `sthlp:1059` case, no `contreat`) | the shipped degenerate-clustering text | the new text |
| any run with a covariance | none of the above | all of the above |

Plus the invariant that matters: **`e(b)`, `e(wdidmainresults)` and every posted matrix are
byte-identical before and after.** `mat list, format(%21x)` on both builds.

### T2 — the screen discriminates (L3). Two fixtures, one statistic.
The two fixtures already exist and differ in exactly one parameter, units per state:
`fixB.dta` (5-12) and `healthy.dta` (40). The predicate under test is
`count(treated observations in a (cell, bin) piece) <= 1 + contreatpoly`, on `esamp`. The
screen must fire on the first and be silent on the second — both already measured — and the
sweep in between (units per state = 5, 10, 20, 40, 80) must be monotone in the saturated-piece
count.
That sweep is the test I would ask for before shipping L3b, because it is what shows the
predicate tracks the geometry rather than a threshold.

### T3 — the standard error itself (the real invariant). Frozen-design MC.
`p05_condmc3.do` is the pattern and needs no truth: freeze the design, redraw only the iid
error, and require

```
mean(reported se) / MCsd(b)   ==   the same ratio on a control arm that is known good
```

within MC error. It is cheap (24 replications was enough to see 0.025 against 1.47), it
generalizes to every estimand and every VCE path in the package, and it is the only test in
this document that would have caught the defect *before* someone noticed the missing standard
errors. **I would suggest adopting it as a general assessment shape, not only for `lattice()`.**

### T4 — cross-engine identity, which already holds and should be kept holding.
`b` on the default (reghdfe) path and on `customvce(robust)` agree to 7-8 digits on this
design (18.09100052 vs 18.09100039). None of the edits here touches an estimate, so the
identity is a regression guard for all of them.

---

# L8 — What I did not test, and the honest residue

1. **Whether `bootstrap()` / `permute()` are calibrated on this geometry.** Both the shipped
   text and my replacement send the user there. Unmeasured. `known_bydesign` §7 already has
   `bootstrap(standard)` over-covering on a different design. The frozen-design MC of T3 with
   a `bootstrap(standard)` arm would settle it, and it is the single most useful follow-up
   here, because a remedy nobody has measured is a recommendation, not a fix.
2. **~~Why `customvce(robust)`'s delta-method variance is so far off~~ — CLOSED by §0.6.**
   Saturated cell-by-bin pieces have identically zero residuals (157 treated observations,
   max |residual| 4.96e-09 against rmse 1.40), so every residual-based estimator is blind in
   exactly the directions the information matrix says are enormous. What I still have *not*
   done is the last arithmetic step — recovering `margins`' weight vector and showing how much
   of its mass sits in the null space of the posted `e(V)`. The measurement does not need it;
   a referee might want it.
3. **The population of real `lattice()` do-files.** Every blast-radius statement about L2c,
   L3 and L3b rests on two fixtures I built. The project's own materials would say more.
4. **`contreatpoly(p > 1)`, user-supplied lattices with empty bins, and
   `contreatcontrols(byCohort)`** — L2c's "one redundant term per treated cell" identity is
   measured only at `contreatpoly(1)` with a generated `lattice(all L)`. L3's predicate does
   not depend on it.
5. **The `poisson`/`ppmlhdfe` twin.** Everything here is the linear path. The same screen
   would have to be sanity-checked there before it goes in unconditionally.
6. **I did not run the project's test suite** (forbidden by the brief), so I cannot say which
   blessed baselines the caution-text change moves. It moves text, not numbers; someone with
   the suite should confirm the scope.
