# SKEPTIC's report on the hv-bootsample findings

Build asserted in every log: `which wooldid` resolves to the pinned
`scratchpad/bluesky/build/wooldid.ado` (r1_guards.log:61). Every probe generates its own
data. No repo data, no repo generator, no test suite. I read the original probe's
do-files but ran none of them; every number below comes from my own fixtures in this
directory.

## Headline

The mechanism the original agent describes is real and I reproduce it independently on
every surface they named. **All eleven findings survive in substance; none is refuted
outright.** What moves is the scoping and the severity, in both directions:

* **HV-BOOT-06 is better supported than the probe could show.** I found a direct read-out
  of the regime bit: the A10 caution at `ado:406-407` prints *iff* `wdb_nest == 1`, so
  naming `i` inside `fe()` prints the regime. It flips exactly as claimed - observed, not
  inferred from a standard error.
* **HV-BOOT-07's numerical worry is refuted.** The acceleration is a pure function of the
  jackknife, so it compares Monte-Carlo-free. It moves **0.26%**.
* **HV-PERM-02's load-bearing claim is refuted.** The "no Monte-Carlo error" sentence is
  printed alongside two cautions and a footer line that say the reference set was reduced
  and name the element; and the reduced-set p is an **exact conditional** randomization p.
* **HV-BOOT-02's severity drops to low.** `bootalpha(equal)` - the default - is
  **bit-identical** with and without the dead rows; this is confined to an option the help
  already calls "not a sound default for panel data" with "badly understated" SEs.
* **New, and the most important thing in this report for realism:** the phantom cluster
  needs neither missing data nor a singleton. **A cohort treated before the estimation
  window opens** mints one - with every variable present on every row. So does a **zero
  survey weight**. And the obvious user-side `if` fixes the missing-data trigger but
  **not** the always-treated one.

---

## Verdicts

### HV-BOOT-01 - CONFIRMED (high), three corrections
`r1_guards.do` CASE A: 20 states x 2 units x 12 years, states 1-4 (all treated) report no
outcome.

| | e(N) | header "N clusters" | e(N_clust) | e(boot_df) = e(df_r) | "<6 treated" caution |
|---|---|---|---|---|---|
| full | 384 | 16 | 16 | **19** | no |
| dead rows deleted | 384 | 16 | 16 | **15** | yes |

att = .58092551 both sides. Same under `bootstrap(bayes)` (19 vs 15). The most damning
artefact is internal: the *same output* prints `N clusters = 16` (r1_guards.log:139) and
`referred to the t(19) distribution` (:178).

**Does it move a number a user reads?** `r9_refusal.do` Q2, 20 frame clusters / 15 live,
reps 400. `lincom main_att` - which the run's own footer tells the user reproduces its
Wald figures, and which `esttab`/`coefplot`/`honestdid` take - gives
CI `[.4964515397, .6754785120]`. Forcing the correct df by hand (`bootdf(14)`, same draws,
identical `se` = .0427675382) gives `[.4942377793, .6776922724]`. Half-widths .08951349 vs
.09172725, ratio .9759, exactly `t(19).975 / t(14).975` = 2.0930/2.1448 = .9758. So the
downstream interval is **2.4% too narrow** here, and 5.8x too narrow in the CASE C limit.

**Correction (a) - the r(459) refusal leak is inconsequential.** `r9_refusal.do` Q1,
12 frame clusters and exactly 1 live: the run is stopped anyway by the *estimation-path*
check (`ado:1521`) - "cluster(state) leaves 1 cluster(s) in the estimation sample;
cluster-robust inference needs at least 2. Use a finer cluster variable, robust, or
fastnose." - which is the more informative of the two messages. The bootstrap's own
refusal is what fires once the dead rows are deleted. Drop that leg of the claim.

**Correction (b) - the doc-vs-code leg is weaker than stated.** The help for `bootdf()`
says *G* is "the number of `cluster()` groups", **not** "in the estimation sample"; and
sthlp:532 explicitly discloses that the two-cluster refusal "fires when the options are
screened, **before the estimation runs**". The claim stands on substance and on the
internal contradiction, not on a broken promise.

**Correction (c) - the standard error is not affected, and the probe correctly did not
claim it was.** `r5_size.do`: under `bootstrap(bayes) bootalpha(equal)` the SE is
`+1.efeac239772daX-004` on *both* datasets - bit-identical. Under `bootstrap(standard)`
the gap is 4% with opposite signs across fixtures, i.e. draw noise.

**User-side remedy (verified, `r11_remedy.do` A):** `if y < .` restores the correct
universe without deleting rows - `e(N)` = 384 either way, `e(boot_df)` 19 -> 15. The
bootstrap frame honours the user `if` (that is what `smp` is built from), so this works
today.

repro: `r1_guards.do`, `r9_refusal.do`, `r11_remedy.do`

### HV-BOOT-05 - CONFIRMED (high), slightly worse than described
`r1_guards.do` CASE C: 12 states, only state 3 (treated) and state 12 (never) report an
outcome. e(N) = 40 and att = .57394774 on both sides.

* full: header `N clusters = 2`, footer `t(11)`, **no caution of any kind**, and the table
  prints `se 4.018e-16  t 1.4e+15  p 0.000`.
* oracle: footer `t(1)`, **both** cautions print ("only 2 clusters are available", "only
  1 treated cluster(s)"), and the row's se/t/p/CI are left blank.

The extra finding: on the full data wooldid *prints a spuriously significant number* where
the oracle blanks the row and warns.

### HV-BOOT-04 - CONFIRMED (high; it is HV-BOOT-01 with a second trigger)
`r1_guards.do` CASE B: 10 states x 2 units x 8 years plus one extra state holding one unit
observed in a single period, **every variable present**. e(N) = 160 and e(N_clust) = 10
with and without that row; e(boot_df) = **10** with it, **9** without. Same severity as
HV-BOOT-01 because it is the same defect; my N1/N2 widen the trigger set further.

### HV-BOOT-06 - CONFIRMED (high), with a stronger instrument than the probe had
`r2b_nest.do`. 16 states x 3 units x 10 years, `i` = unit, `cluster(state)` (coarser),
`att_i`, `bootstrap(standard)`. One stray record codes unit 19 to state 8 for one year
with the outcome missing.

**Direct regime read-out.** `ado:406-407` prints the A10 caution iff `wdb_reid` is on
**and** `wdb_nest == 1`. Adding `fe(unit)` therefore prints the regime bit:

| run | A10 caution | regime |
|---|---|---|
| oracle, no stray row | **printed** (r2b_nest.log:187) | nesting |
| full, stray row with y **missing** | not printed | keep-i |
| control, stray row with y **present** | not printed | keep-i (correct there) |

The predicate at `ado:8987-8993` does not read `fe()`, and the read-out runs share the
data with the SE runs, so the regime shown transfers.

**Draw-stream isolation.** Deleting that row does not change G, so with a fixed `rseed()`
the `bsample` stream is unchanged: `bootsave()`'s per-cluster multiplicity columns give
`max|w_full - w_oracle| = 0.0000000000` over 401 saved rows, in **every** cell.

**Standard errors** (att_i; b bit-identical full vs oracle in every cell):

| data seed / rseed | FULL | ORACLE | CONTROL | full/oracle |
|---|---|---|---|---|
| 101 / 5551212 | .0978610043 | .1156219131 | .0977106359 | .846 |
| 101 / 8675309 | .0954973057 | .1138693879 | .0953086558 | .839 |
| 202 / 5551212 | .0695897059 | .0814603740 | .0702489184 | .854 |
| 202 / 8675309 | .0722971475 | .0820128301 | .0727006948 | .881 |
| 303 / 5551212 | .1153971163 | .1363411770 | .1142732772 | .846 |

The CONTROL - the same stray row carrying a real outcome, where the straddle is genuine
inside `e(sample)` and keep-i is the **correct** regime - lands on the full-data SE to
three decimals in every cell and nowhere near the oracle. That rules out the extra row
itself and pins the gap to the regime.

**Correction to the framing.** `ado:8975-8986` argues that computing the predicate "on the
raw frame instead" would misclassify. `smp` **is** the raw frame minus the user `if` and
missing clusters. The comment describes the defect it was written to prevent; the intent
is right and the sample concept it was implemented on is the wrong one.

### HV-BOOT-02 - CONFIRMED mechanically, severity revised DOWN to LOW (probable)
`r5_size.do`: 14 states x 2 units x 12 years, y missing at random on half the rows, no
cluster wholly dead. Frame rows 336, live rows 160, e(N) = 160 both sides, att
bit-identical (`+1.603ab04740d61X-001`).

`bootalpha(size)` se: **.0270485748** (full) vs **.0400620885** (dead rows deleted),
ratio **.6752** against the predicted sqrt(160/336) = .6901. The arithmetic holds.

**Why low.**
1. **Confined to `bootalpha(size)`.** `bootalpha(equal)`, the default, gives
   `+1.efeac239772daX-004` on *both* datasets - bit-identical, as Dirichlet aggregation
   requires.
2. **The help's own magnitude already describes frame counting.** sthlp:548 says the
   shrinkage is "roughly the square root of the mean cluster size". Frame n_c = 24 gives
   sqrt(24) = 4.9; the live mean is 11.4, giving 3.4. Measured in this fixture:
   `standard` .1428631773, `equal` .1210734927, `size` .0270485748 - a 5.3x understatement
   that matches the *frame* figure. So the help already documents the magnitude the code
   produces; what is imprecise is the word "observation count".
3. **Parity probably favours the code.** I could not obtain `bayesboot` (not installed
   here; not at the usual SSC path). But `bayesboot` draws its Dirichlet shape from the
   rows an `id()` holds *in memory* - it cannot know what the estimation command will
   drop - so frame counting is most likely the advertised parity. Marked **probable**.

### HV-BOOT-03 - CONFIRMED (severity revised DOWN to low)
Same `r5_size.do` run: the saved file's `_rep == 0` row holds `w_1 = 0.0029761905`,
exactly 1/336 (the frame row count), while e(N) = 160; the oracle file holds
0.00625 = 1/160. sthlp:557 says the row holds "1/*N*, *N* the number of observations in
the estimation sample". Divergence confirmed verbatim.

**Correction.** "the observed row cannot be used to recentre the draws as documented"
overreaches: recentring the coefficient draws uses the *coefficient* columns of the
`_rep == 0` row, not the `w_` columns. Only the weight-scale identity moves, and
`ado:8949` shows 1/N_frame is the scale on which the help's own companion identity
(sum_c n_c w_c = 1) actually holds. This is a help-text fix on a diagnostic file.

### HV-PERM-01 - PARTIAL (medium), and I measured what they disclaimed
`r3_perm.do` A, `r3b_perm.do` P. 20 states, 8 treated, states 1-4 (treated) report no
outcome.

Confirmed: `permsave()` writes **20** onset columns (onset_1..onset_20) on the full data
and **16** (onset_5..onset_20) on the oracle; reported possible assignments
**317,444,400** vs **43,680**; the caution "only 4 of the resampling clusters are ever
treated" fires on the **oracle only**. e(N) = 320 both.

**New: the reference distribution itself moves** (`r3b_perm.do`, weak effect so p is
interior, permreps 600, same rseed):

| data seed | full null SD | oracle null SD | ratio | full p | oracle p |
|---|---|---|---|---|---|
| 11 | .11977383 | .15011282 | .798 | .46589018 | .50582363 |
| 22 | .11848298 | .13573510 | .873 | **.04991681** | **.07986689** |
| 33 | .11484640 | .13482312 | .852 | .27287854 | .37271215 |

The robust statistic is the **null SD, 13-20% tighter on the full data against ~3% Monte
Carlo noise at 600 reps, in the same direction on all three seeds**; the p is smaller in
all three and at seed 22 the 5% verdict flips (that one gap is ~3 MC SEs, so it is an
illustration rather than the evidence). The mechanism is countable: the 20-cluster design
gives a replication 8 x (16/20) = 6.4 live treated clusters against the observed run's 4,
so its null is narrower than the observed statistic's.

**Correction - "misreported by four orders of magnitude" is not right.** 317,444,400 is
the arithmetically correct count for the design wooldid ran, and both designs are valid
randomization models (the full design's exclusion of non-estimable assignments is
conditioning on a deterministic function of the assignment with the observed assignment
retained - a valid conditional test). The defensible criticisms are: the
few-treated-cluster guard is silenced by clusters carrying no data; nothing in the output
says which design ran; and the count calibrates the "permreps() exceeds 10% of the
possible assignments" advice.

### HV-PERM-02 - PARTIAL; the load-bearing claim is REFUTED (low)
`r3_perm.do` B, my own 7-cluster fixture (unit 1 treated 2003, unit 2 treated 2004, units
3-7 never; units 1 and 7 report no outcome). e(N) = 50 and att = .42990721 both sides.

* full: `npossible = 42`, `completed = 40`, `e(perm_exact) = 1`,
  `e(perm_incomplete) = main_att`, **p = 0.2000**
* oracle: `npossible = 5`, `completed = 5`, `e(perm_incomplete)` empty, **p = 0.4000**

The design-size and p movements reproduce. But *"'exact shares with no Monte-Carlo error'
is false as printed"* does not survive: the **same output** prints, in order,

1. `Caution: 2 of 42 enumerated assignments failed to estimate, so the enumeration is
   incomplete. The p-values stay exact shares ... but over the 40 assignments that
   completed ... with each element's n in e(perm_n).`
2. `Caution: under permexact the following elements were not estimable in every enumerated
   assignment ... a reduced reference set whose size is in e(perm_n): main_att`
3. footer: `Elements not estimable in every enumerated assignment (exact shares over a
   reduced reference set; their n is in e(perm_n)): main_att`

Three disclosures naming the element. And the statistics are sound: conditioning a
randomization reference set on estimability - a deterministic function of the assignment,
with the observed assignment retained - is an **exact conditional** randomization test.
What is left is HV-PERM-01's point: two valid designs, chosen by whether the analyst
deleted unusable rows, with nothing saying which.

### HV-PERM-03 - CONFIRMED (medium)
`r3_perm.do` C and `r3b_perm.do` C, now with working plain-run controls:

| run | rc | e(N) |
|---|---|---|
| C0 ordinary run on the stray-onset file | 0 | 160 |
| C0o ordinary run, stray row deleted | 0 | 160 |
| C1 `permute()` on the stray file | **r(459)** cluster-onset nesting | - |
| C1s same + `safety(off)` | **r(459)** (not lifted) | - |
| C2ord ordinary run on the `.a`-onset file | 0 | 160 |
| C2 / C2s `permute()` on it | **r(459)** / **r(459)** | - |
| C3clean `permonsetcohorts`, clean fixture | 0 | 216 |
| C3ord ordinary run on the stray-onset file | 0 | 216 |
| C3 `permonsetcohorts` on it | **r(459)** | - |
| C3s same + `safety(off)` | **r(459)** (sthlp:520 says so) | - |
| C3o stray row deleted | 0 | 216 |

The estimation sample is identical in every successful run, so a single record wooldid
would itself drop makes the whole permutation engine unavailable, and the message calls
that record "in-sample".

**Note on the fix.** The refusal is internally *consistent*: `wooldid_permrun` builds its
design on the same unscreened `smp`, so a screen on the screened sample would let through
a design the engine then cannot run. The fix is the shared sample concept, not the screen
- and my N6 shows a markout-aware screen of exactly this shape already exists in the
worker.

### HV-VAL-01 - CONFIRMED (low), one correction
`r4_val.do` A: 20 unusable records carrying `unit == -9` with no outcome, appended to a
12-state panel.

| run | rc | e(N) |
|---|---|---|
| no `if` | **r(411)** "unit must consist of positive, non-zero integers, but does not" | - |
| `if junk == 1` | 0 | 240 |
| `safety(off)` | 0 | 240 |
| sentinel rows deleted | 0 | 240 |

**Correction the probe missed: `safety(off)` lifts it** (`ado:4959`). A remedy exists;
the defect is that the message names *neither* `safety(off)` *nor* `if()`, which is what
makes the failure read as a broken id variable. Message defect, not a missing escape
hatch - low rather than medium.

### HV-BOOT-07 - CONFIRMED mechanism (low); the numerical worry is REFUTED
`r8_fixups.do` PART 1. 18 states, 5 report no outcome. e(N) = 260 and e(N_clust) = 13 on
both sides.

* `e(boot_jk_completed) = 18` against `e(N_clust) = 13`; the footer says "delete-one-cluster
  jackknife (**18 fits**)". The structural claim is in `e()`, not just in prose.
* **The isolation the probe skipped.** `a = d3/(6 d2^1.5)` over the jackknife estimates
  alone (`ado:15786-15790`), so it does not depend on the draws. Two runs on the same data
  with *different* `rseed()`s give `e(boot_accel)` bit-identical
  (`+1.d4f20b77b80d3X-005`) - that *proves* draw-independence rather than assuming it.
  The oracle gives `+1.d3b9024548fa4X-005`.
* Decimal: **0.05724432** (18 fits) vs **0.05709506** (13 fits) - **0.26%**.

The probe's p 0.007-vs-0.021 and CI gap were draw-stream confounding, exactly as they
flagged. Real cost: five wasted refits and a footer that overstates the jackknife.
Severity stays **low**; confidence rises from uncertain to **certain**.

---

## New findings (mine)

### N1 - the phantom bootstrap cluster needs no missing data and no singleton: a cohort treated before the window opens mints one
`r7_novel.do` A / `r10_alwaystreated.do` / `r11_remedy.do` B. 16-20 states, 2 units each,
2000-2011. Some states carry onset **1996** - treated before the estimation window, which
sthlp:73 says wooldid drops by design. **Every variable is present on every row.**

| fixture | e(N) | e(N_clust) | e(boot_df) full | e(boot_df) oracle |
|---|---|---|---|---|
| r7 A (16 states, 4 always-treated) | 288 | 12 | **15** | 11 |
| r10 A (20 states, 6 always-treated) | 336 | 14 | **19** | 13 |

The `<6 treated clusters` caution fires on the **oracle only** (r10_alwaystreated.log:181).
Under `permute(shuffle)` the reported possible assignments are **10,810,800** vs **5,940**
(r7) and **465,585,120** vs **12,012** (r10), with e(N) identical on both sides; the
permutation p moves (.0732/.0832, .3661/.3827, .0033/.0100 across three seeds) though the
null-SD direction is not consistent there, so I do not claim a size distortion for this
trigger - only for the missing-outcome one (HV-PERM-01).

**And the remedy is not obvious here** (`r11_remedy.do` B): the natural hygiene filter
`if y < .` leaves `e(boot_df)` at **15** (no change); only
`if missing(onset) | onset >= <first period>` brings it to **11**. wooldid's own advice
for pre-sample onsets is to truncate the estimation window, which does not remove the
phantom clusters.

This is what changes the realism assessment of the whole family. A staggered panel with
some units already treated when the window opens is the default shape of administrative
rollout data, not an edge case.

### N2 - a zero survey weight does the same
`r7_novel.do` B. 14 states, two carrying `w == 0` on every row (`ado:686` drops `w <= 0`).
e(N) = 288 and e(N_clust) = 12 both sides; `e(boot_df)` **13** with the zero-weight rows
and **11** without. Zero-weight records are routine in weighted administrative and survey
extracts.

### N3 - `permstrata()` is built on the same unscreened set
`r7_novel.do` E. 16 states stratified by region, states 1-4 (all region 1, all treated)
report no outcome. e(N) = 288 both sides; reported possible assignments **560** vs **12**.

### N4 (negative, valuable) - `impvalues` is CLEAN
`r7_novel.do` D. 14 states, two wholly dead. `e(N)` = 288, `e(N_clust)` = 12 and the
Ibragimov-Muller `e(df_r)` = 5 are identical with and without the dead rows. The
cross-cluster IM machinery counts live clusters correctly - the one inference engine on
this axis without the defect.

### N5 (negative) - the event-study window does not drop rows
`r7_novel.do` C. `ado:9367`'s comment names "an event-study window" as one of three
triggers. `esprelength(2) espostlength(2)` on a 20-period panel leaves e(N) = 480 and
`e(boot_df)` = 11, identical to the unwindowed run. Corroborates the original probe's
negative result: of the three the comment names, missing-y and singletons reproduce and
the window does not.

### N6 - the worker's own onset-uniqueness screen IS markout-aware, so the fix already exists in the codebase
Found while building `r2_nest.do`. A stray row for unit 19 carrying a *different* onset
trips `r(459)` "wooldid requires exactly one treatment onset per group defined by unit"
**only when that row has an outcome** (r2_nest.log:336-349); with the outcome missing the
same row passes silently. So wooldid already contains a screen of exactly this shape
evaluated on the internally screened sample, while `wooldid_permctl`'s structurally
identical nesting screen (HV-PERM-03) evaluates on the unscreened one and refuses. The two
are inconsistent with each other, and the markout-aware one is the model for the fix.
