# Fix proposals — mechanism `IMP` (`impvalues`, Ibragimov-Müller)

Pinned build `c1d8979`. Working copy of the ado patched and run; the repo was not touched.
Probes, patched ados and logs:
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-IMP/`

| # | findings | what | moves a number? | cost | tested |
|---|---|---|---|---|---|
| **P1** | `mcci-01-impvalues-se`, `skep-02-im-degeneracy-message` | one geometry screen before the fit (+ a conditional tail on the degeneracy message) | **no** | 0.13 s at 1M rows; under Variant A it also **replaces a 369 s run** in the case it fires | yes — **both variants** patched + run; `%21x` identity on 100 hex lines |
| **P2** | `skep-01-impvalues-imunit` | help only — no code change; a measured lever the help can give the user | no | zero | mechanism measured, R=150 |

---

# P1 — one geometry screen, covering `mcci-01` and `skep-02`

## 1. The defect, stated so the fix follows from it

`ado:2042` is **one** `margins … nose over(`cluster')` on **one** pooled fit; `2047-48` grab the G
reads, `2052-53` mean them, `2054` t-tests them. The G reads are therefore
`W · beta`, where `beta` is the vector of treated cohort-by-relative-period coefficients and `W` is
the G x C matrix of each cluster's estimand-weight profile over those cells.

They are not G separate estimates. They **do** carry sampling variation (the refuter's correction to
`mcci-01` is right and I have kept it): every read is a linear function of the same estimated
`beta`. What is absent is *independence and count*: `rank(W) <= min(G, C)`, so when `cluster()` does
not nest `i` the G reads have at most C independent components while the t-test divides their spread
by `sqrt(G)` and reads it off `t(G-1)`. Refining `cluster()` below `i` raises G and leaves `beta`,
`W`'s row space and the data untouched.

**The discriminating condition is exactly "does `cluster()` nest `i`".** When it does, each cluster
owns a disjoint block of cohorts, `W` is block-diagonal by cluster, and rank(W) = G. When it does
not — `cluster()` finer than `i`, or crossing it — several clusters load on the same coefficients.

### I reproduced the identity on my own DGP and seeds
`p1_identity.do` / `.log`: 60 units in 12 cohort groups, 12 periods, unit shocks sd 0.40,
idiosyncratic sd 0.60. **Same data, same fit, only the cluster labels change:**

| run | `att` b | IM se | 95% CI | `e(df_r)` |
|---|---|---|---|---|
| `i=cg, cluster(cg)` (nests) | 1.0098577 | 0.06829971 | [0.85535, 1.16436] | 9 |
| `i=cg, cluster(unit)` (finer) | **1.0098577** (bit-identical, reldiff `0.0e+00`) | **0.02927130** | [0.95103, 1.06868] | 49 |
| analytic `cluster(unit)` | 0.98621117 | 0.09785480 | [0.79040, 1.18202] | 59 |

`se(unit)/se(cg) = 0.42857142857142855`; `sqrt(9/49) = 0.42857142857142860`; deviation **5.55e-17**.
Relabelling the same 720 rows from 12 clusters to 60 shrinks the reported standard error by 57% and
the interval with it. That is the whole defect in one line, and it needs no oracle.

**Be precise about the point estimate: it is invariant to `rtol 1e-12`, not always bit-identical.**
On this run (no event-study window requested) the two `b`'s agree at `reldiff 0.0e+00`. On the same
DGP and seed *with* `esprelength(1) espostlength(1)` — the `p9a.log` fixture — they are **1 ulp
apart**: `+1.d95eca99eddc2X-001` against `+1.d95eca99eddc0X-001`, about 4e-16 relative. That is F062
(`sthlp:375`): under `impvalues` the cluster variable joins the cell key, so the two clusterings
present different collapsed frames to the fit and the floating-point accumulation differs in the last
place. It does not weaken the finding — a 4e-16 move in `b` against a **57%** move in `se` is the
point — but a verifier asserting exact equality on the `b` will get a spurious failure, so the
invariant in §9 is stated at `rtol 1e-12`.

## 2. The exact edit

### 2a. The screen — insert before `ado:1447` (`  local vcespecifier`)

This site is **after** the `jointtests` `jtcoarse` block and **before** the regression at 1476, so a
refusal costs nothing and a caution prints ahead of the table. `ifstatement`, `i` and `cluster` are
all live there, and it is before the `preserve`/cell collapse at 1875 (under `impvalues` the cluster
joins the cell key at 1892, so a post-collapse predicate would be computed on a different frame).

```stata
  * --- IM GEOMETRY SCREEN (impvalues): does cluster() nest i? ------------------
  * Ibragimov-Muller needs one estimate per cluster carrying that cluster's OWN
  * sampling variation. wooldid reads its per-cluster estimates from the SINGLE
  * pooled fit with margins, over(cluster), so when several clusters share levels
  * of i they load on the same cohort-by-relative-period coefficients and differ
  * only in each cluster's composition over those cells. Same predicate and idiom
  * as the jointtests `jtcoarse' tag-count at 1425/1429, opposite polarity.
  * Computed on `ifstatement' (the rows margins will see) and BEFORE the fit.
  local imnest 1
  if "`impvalues'" != "" & "`cluster'" != "" {
    tempvar imgt
    qui gegen byte `imgt' = tag(`i') if `ifstatement'
    qui count if `imgt' == 1
    local imni = r(N)
    cap drop `imgt'
    qui gegen byte `imgt' = tag(`i' `cluster') if `ifstatement'
    qui count if `imgt' == 1
    local imnic = r(N)
    cap drop `imgt'
    local imnest = (`imni' > 0 & `imnic' == `imni')
    if `imnest' == 0 {
      dis as text "Caution: impvalues was combined with a cluster() variable that does not nest the cohort variable `i' (`imnic' distinct (`i', `cluster') combinations for `imni' values of `i'). wooldid reads its per-cluster estimates with margins, over(`cluster'), from one pooled fit, so clusters sharing levels of `i' load on the same cohort-by-relative-period coefficients and their reads differ only in each cluster's composition over those cells. Ibragimov and Muller inference requires one estimate per cluster carrying that cluster's own sampling variation, which this geometry does not deliver: the reported standard error understates the estimator's sampling variability and the interval is narrower than its stated level (measured coverage 5 to 38 percent at a nominal 95 percent), and refining cluster() below `i' shrinks it further with no change to the data. Cluster at or above the level of `i' under impvalues, or use bootstrap() or permute(), whose inference is built from replications."
    }
  }
```

**Variant A (refusal) — what the evidence actually supports.** The `if `imnest' == 0 {` body
becomes a two-branch block: refuse under default safety, caution and proceed under `safety(off)`.
**Built and run** (`patched_gated/wooldid.ado`, `p10_gated.do`):

```stata
    if `imnest' == 0 {
      * F13 convention (the jointtests screen's own comment at ado:1399): an integrity
      * screen gates on `safetyinput' - the user's own safety() - never on `safety',
      * which contreatwithin rewrites to off internally at ado:4829.
      if "`safetyinput'" != "off" {
        dis as error "Error: impvalues requires a cluster() variable that nests the cohort variable `i' (… same body as the caution …). To proceed anyway and receive this as a caution, specify safety(off)."
        error 459
      }
      dis as text "Caution: impvalues was combined with a cluster() variable that does not nest …"
    }
```

**Gate on `safetyinput`, not `safety`** — this is a correction to the obvious form. `ado:4829` is
`if "`contreatwithin'" == "contreatwithin" & "`contreat'" != "" local safety off`, and nothing
refuses `impvalues` + `contreatwithin`, so a screen written on `safety` would silently stand down for
a user who never asked for `safety(off)`. `ado:1399`'s own comment states the rule: *"so it gates on
safetyinput, not safety (contreatwithin overwrites the latter internally)"*. `safetyinput` is live in
`wooldid_worker` well before the insertion point (716, 772, 995, 1355, 1399).

`459` matches `permute()`'s refusal of the same geometry at `ado:6078`, and the caution-under-
`safety(off)` shape is the convention at 836/908/955/1713.

**Variant B (caution only)** is the block exactly as written above, with no error branch. **Both
variants are built and measured** — B in `patched/`, A in `patched_gated/` — and §3 gives the runs.

### 2b. The degeneracy message — `ado:2018` (`skep-02`)

Replace the single `local imdegmsg "…"` with:

```stata
    * skep-02: name the condition that actually fired. Under a cluster() that does not
    * nest i the collapse is a property of the cluster GEOMETRY, not of the data's noise
    local imdegtail "A collapsed IM standard error usually reflects deterministic or noiseless data - please verify your input if this was not expected."
    if `imnest' == 0 local imdegtail "Under a cluster() that does not nest `i' - the geometry the caution above names - this is a property of the cluster structure rather than of the data: the per-cluster reads are a deterministic reweighting of one pooled fit, so clusters with the same composition over the `i'-by-relative-period cells read identically however noisy the outcome is. Change the cluster variable, not the data."
    local imdegmsg "Caution: for one or more Ibragimov-Muller estimates the per-cluster estimates were numerically identical, so the cross-cluster (IM) standard error is not identified; wooldid has set the standard error, test statistic, p-value, and 95% confidence interval to missing for the affected estimate(s) while retaining the point estimate. `imdegtail'"
```

One local feeds all six IM t-test sites, so this is one edit for every `impvalues` row.
**This does not touch the guard predicate** (`se < 1e-10*|b|`, relative-only) — CORE-06 is the
protected class and stays exactly as it is. Only the sentence changes.

Under **Variant A** the geometry branch of `imdegtail` becomes unreachable on the default path
(the run is refused earlier), but it is still reached under `safety(off)`, so keep it either way.

## 3. Measured before/after (this is the part that is not an argument)

`p2a.do` (pinned build) and `p2b.do` (patched), byte-identical bodies, six cases:

| case | spec | screen fires? | `e(b)`, `e(V)`, `e(df_r)`, `e(N)` in `%21x` |
|---|---|---|---|
| A | `i=cg, cluster(cg)` | no | identical |
| B | `i=cg, cluster(unit)` (finer) | **yes** | identical |
| C | `i=cg, cluster(reg)` (crossing) | **yes** | identical |
| D | `i=unit, cluster(unit)` | no | identical |
| E | `i=cg, cluster(unit)`, **no** `impvalues` | no (screen skipped) | identical |
| F | 25 balanced regions crossing `i`, exact degeneracy | **yes** + new tail | identical |

`diff` of the hex lines from the two logs is **empty**. The caution fires in exactly 3 of 6 cases —
the 3 that should — and in none of the 3 that should not.

**Widened past those four scalars.** `p9a.do` (pinned) / `p9b.do` (patched) are the same six cases
with `mat list e(wdidmainresults), format(%21x)` and `mat list e(b), format(%21x)` added, so the
comparison covers **every reported row's b, se, t, p, lb, ub and relyr**, not just the `att` scalar.
Bodies byte-identical apart from the `adopath` line. Extracted hex lines: **100 from each log,
`diff` empty.** Caution count: 3 in the patched log, 0 in the pinned one. Case F's degeneracy tail is
the only other difference, and it is prose.

Case F before / after, on a panel whose `sd(y)` is 1.659689:

> *before:* "… A collapsed IM standard error usually reflects **deterministic or noiseless data** -
> please verify your input if this was not expected."
> *after:* "… **Under a cluster() that does not nest cg** — the geometry the caution above names —
> this is a property of the cluster structure rather than of the data … **Change the cluster
> variable, not the data.**"

**Variant A**, `p7_refuse.do` on a 24,000-row panel (2,000 units, 50 cohorts):
`cluster(unit) impvalues` → `r(459)` in **0.209 s**; `cluster(cg) impvalues` → `rc 0`, unchanged.
That form appended a bare `error 459`.

**The gated form of Variant A is now run too** (`p10_gated.do` against `patched_gated/`, the
`safetyinput` gate above, the `p9a.log` fixture — §1's DGP and seed plus `esprelength(1)
`espostlength(1)`, so its `att` differs from §1's — `%21x` against the *pinned* build's
answers in `p9a.log`):

| case | spec | rc | seconds | `b` / `se` / `df_r` / `N` vs pinned |
|---|---|---|---|---|
| G1 | `i=cg, cluster(cg)` (nests) | **0** | 0.729 | **identical** to pinned A |
| G2 | `i=cg, cluster(unit)` (finer) | **459** | **0.042** | — (refused before the fit) |
| G3 | G2 + `safety(off)` | **0** | 1.251 | **identical** to pinned B |
| G4 | `i=cg, cluster(reg)` (crossing) | **459** | — | — |
| G5 | `i=unit, cluster(unit)` (`cluster`==`i`) | **0** | — | **identical** to pinned D |
| G6 | `cluster(unit)`, **no** `impvalues` | **0** | — | **identical** to pinned E |

G3 is the load-bearing row: the `safety(off)` bypass returns the pinned build's answer **bit for
bit** (`b = +1.d95eca99eddc0X-001`, `se = +1.198c7e40347f6X-005`, `df_r = 49`, `N = 720`), so the
escape hatch is a genuine "proceed unchanged" and not a different code path. The refusal costs
**0.042 s against the 1.251 s** the run itself takes on this small fixture, and 0.209 s against
369 s on the 24,000-row one.

## 4. Cost — priced, because the author has said he does not want expensive guards

`p3_cost.do`, the screen in isolation (two `gegen tag` + two `count`, gtools):

| rows | seconds |
|---|---|
| 200,000 | 0.078 |
| 1,000,000 | 0.128 |
| 5,000,000 | 0.629 |

It runs **once per run**, only under `impvalues`, and it is the same idiom `jointtests` already pays
unconditionally at `ado:1425/1429`.

**It is cheaper than free.** Same file, same 24,000-row panel, full `impvalues` run:

| geometry | seconds |
|---|---|
| `cluster(cg)` — nests, 50 clusters | **12.885** |
| `cluster(unit)` — finer, 2,000 clusters | **369.013** |

The broken geometry is also **28.6x slower**. (The 28.6x is measured; its attribution is not
isolated — `margins, over(cluster)` is O(G) subpopulation reads, but `reghdfe`'s clustered VCE at
G = 2,000 against G = 50 is in the same run. Either way the cost is real and it is on the wrong side.)
**This saving is Variant A only.** Under Variant B the 369 s run still happens and the caution prints
on top of it; B buys the warning, not the time. Under A the user is told in 0.209 s instead of waiting
six minutes for a standard error 3.3x too small.

## 5. Does it move a reported number?

**No — neither variant.** Where the predicate passes, nothing executes but two counts; `p2a` vs `p2b`
is bit-identical on `e(b)`, `e(V)`, `e(df_r)` and `e(N)` across all six cases including the three
where the caution prints, and `p9a` vs `p9b` widens that to the **whole of `e(wdidmainresults)`**
(every row's b, se, t, p and interval) and `e(b)` — 100 hex lines each, `diff` empty. Under Variant A
the three cases that still run (`p10_gated` G1, G5, G6) are likewise `%21x`-identical to the pinned
build, and so is the `safety(off)` bypass (G3). Under Variant A the failing geometry produces an error instead of a number,
which is a behaviour change, not a moved number. **No re-bless of any passing config is required.**

## 6. Blast radius

Read-only greps of `test/suite/configs.csv`, `configs_ext.csv` and the two generators (the suite was
**not run**). Every `impvalues` config, checked against how its data is built:

| config | `i` / `cluster()` | screen |
|---|---|---|
| `impv` | `i` / `cluster(i)` | passes |
| `trig_imdeg_part` | `i` / `cluster(i)` | passes |
| `ext_bin_ext_impvalues` | `cohort` / `cluster(cohort)` | passes |
| `ext_cov10_cont_impvalues` | `cohort` / `cluster(cohort)` | passes |
| `real_impv` | `unit` / `cluster(st)` (coarser) | passes |
| `trig_imdegen` | `i` / `clus = 1+mod(floor((i-1)/3),12)`, a function of `i` → coarser | passes |
| **`adm_impv`** | **`cohid` / `cluster(state)`** | **FIRES** |
| the 6 `*_val_imp_*` / `scr_impsg` screens | refused before estimation today | unaffected |

`adm_impv` is the one hit, and it is the finding's best witness, not a false positive.
`mkdata_admin.do:54` is `gen cohid = cond(ttre == ., 1, ttre - 2000)` — **four national onset
cohorts** — run against **48 states**. A cohort spans states, so `cohid` does not nest in `state`:
this is the crossing geometry, on the project's own "realistic administrative panel". Its config
comment reads *"IM inference at 48 state clusters on realistic data"*, which is how the geometry got
past review. Under Variant B it gains one caution line (check `degeneracy_expected.csv`, which pins
caution counts). Under Variant A it must be re-specified — the natural repair, `cluster(state)` with
a state-level `i`, or `i = cohid` with `cluster(cohid)` — and its baseline re-blessed.

Elsewhere: nothing. `impvalues` is already refused with `permute()`, `bootstrap()`, `fastnose`,
`subgroups()`, `jointtests`, `adversarial`, `makecfxplots`, `pretrendpower()`, `customvce()` and
`unconditionalse`, so there is no replication path and no nested-call concern. `imnest` is a local in
`wooldid_worker` consumed 570 lines later in the same program.

**`contreatwithin` is the one place the gate choice matters.** `ado:4829` rewrites `safety` to
`off` internally for `contreatwithin` + `contreat()`, and nothing in `wooldid_validate` refuses
`impvalues` + `contreatwithin` (unlike `permute()`, which refuses all three at `ado:5802`). A
Variant A written on `safety` would therefore stand down on exactly that combination without the
user asking for it; written on `safetyinput` it does not. I did not run an
`impvalues` + `contreatwithin` case — the point is established from the two source sites and from
`ado:1399`'s own comment, which is the project's ruling on the same question for `jointtests`.

**One honest limitation.** The predicate is computed on `ifstatement`, not on the realized
`e(sample)`, so a cluster that straddles two `i` values only through rows `reghdfe` later drops as
singletons would be classified non-nesting. That is the conservative direction (it warns where it
need not) and it is the same sample concept `jointtests` uses at 1425. Computing it on `e(sample)`
would mean running it after the fit, which forfeits the 369 s saving and the pre-fit refusal.

## 7. T4 and PS-64 — why this is not the thing that was already declined

**T4 (2026-09-15) is a different defect.** It is about the *estimand*: IM reports the
cluster-equally-weighted mean whatever estimand you asked for, so non-`att` estimands under
cross-cutting clusters are the wrong quantity. The author declined both a refusal and a caution and
took the documentation instead, and the help now leads with the equal weighting. **Everything that
ruling says is still true and none of it is contradicted here.** This finding fires on `att` alone
and is about the *standard error*, which T4 never examined: the point estimate is bit-identical
across the two clusterings in my table above — it is only the interval that moves, by 57%, on
identical data.

**PS-64 (2026-09-08) is the closest precedent and it is a caution.** The joint pre-trend test at or
above `i` is oversized 11-36% at nominal 5%, and the ruling was "compute the test with a named
caution instead of refusing … a hard error sent users into a documented dead end … The caution states
what is wrong with the covariance; the number is the user's to weigh." That is Variant B, verbatim.

**The case for Variant A anyway, which the author should hear before ruling:**
1. PS-64's stated reason not to refuse was that *wooldid's own help twice recommends the clustering
   that trips the condition*, so refusal was a dead end. Here the opposite holds: the remedy —
   cluster at or above `i` — is a geometry `impvalues` fully supports and where the review measured
   coverage 0.810 at `cluster() == i` and 0.935 above it, against 0.048-0.380 below and across.
   There is no dead end to send anyone into. (P2 below shows the residual 0.810 has its own cause
   and its own lever.)
2. `permute()` already refuses this exact geometry with `r(459)` at `ado:6078`. The reason differs
   (a straddling unit would carry two onsets) but the user-visible rule would be the same across two
   of the three inference engines, and the message can borrow its shape.
3. The magnitudes are not PS-64's. A 95% interval covering 4.8% is not a number to weigh.
4. The refusal is the cheap branch: 0.209 s instead of 369 s.

## 8. Alternatives I rejected

- **G genuine per-cluster refits.** This is what IM actually is. It cannot be done here: `wooldid y
  cg t ttre, att cluster(unit)` restricted to one unit's rows exits `r(459)` — one unit has no
  control group. Where it could run it would cost G fits. Out of scope; a feature, not a fix.
- **Aggregate the reads to the join of the `i` and `cluster` partitions** (connected components of
  the bipartite graph), which gives the `cluster(cg)` answer in the finer case and one block — hence
  no inference — in the crossing case. Principled, but it is an unbounded-iteration label
  propagation over the full frame (the author's explicit no on expensive guards), and it *moves every
  reported `se`, `t`, `p` and interval* on the affected geometry, so it needs a re-bless on top of a
  ruling. The screen gets the same protection for 0.13 s.
- **Reweight the per-cluster estimates to the requested estimand's weights.** Killed by T4's own
  reasoning: it buys the estimand by spending the `t(q-1)` justification.
- **Deduplicate the reads and use the distinct count as `q`.** Works on the finer-than-`i` case,
  where the reads are exact copies to 5.6e-16, and fails silently on the crossing case, where the
  reads are distinct but still rank-deficient. A detector that catches half the problem is worse than
  the geometric one that catches both.
- **Silently fall back to the analytic path.** Changes the inference without saying so.
- **Tighten the degeneracy guard with an absolute floor** so near-degenerate reads are scrubbed.
  CORE-06 fixed that tolerance by measurement and protects the miss; and it would still not catch
  `cluster(unit)` above, where 60 reads take 12 distinct values and the se is 3.3x too small but
  nowhere near zero.

## 9. The test that proves it — no oracle needed

1. **Relabelling invariance (the primary one).** On any panel, `wooldid y i t ttre, att cluster(i)
   impvalues` and `wooldid y i t ttre, att cluster(unit) impvalues` where `unit` is any strict
   refinement of `i`: `e(b)` must agree to `rtol 1e-12` (bit-identical with no event-study
   window; 1 ulp apart with one, per F062 — see §1) while `sqrt(e(V))` differs by **exactly**
   `sqrt(df_coarse/df_fine)`. Pinned build: se-ratio deviation 5.55e-17, `b` reldiff <= 4e-16. A correct implementation cannot
   satisfy that identity, because a real per-cluster estimate would change when the clusters change.
   Under the patch the fine run is refused (A) or cautioned (B).
2. **Screen exactness.** A 6-cell fixture (`p2b.do`, cases A-F) asserting the caution fires on
   exactly {finer, crossing} and on neither {cluster==i at cohort level, cluster==i at unit level,
   no-impvalues}.
3. **No-number-moved.** `%21x` dumps of `e(wdidmainresults)` (every row's b, se, t, p and interval),
   `e(b)`, `e(V)`, `e(df_r)` and `e(N)` on the passing cases, patched vs pinned, `diff` empty.
   Already run (`p9a` vs `p9b`): 100 hex lines from each log, `diff` empty.
4. **The `safety(off)` bypass is a true no-op.** Under Variant A, the refused geometry re-run with
   `safety(off)` must return the pinned build's answer bit for bit. Already run (`p10_gated` G3):
   `b`, `se`, `df_r`, `N` all identical to the pinned `p9a` case B.
5. **Cost ceiling.** Assert the screen is under 1 s at 5M rows (measured 0.629), and that the
   refusal branch returns in well under the time the run it replaces takes (measured 0.042 s
   against 1.251 s on the small fixture, 0.209 s against 369 s on the 24,000-row one).

---

# P2 — `skep-01`: help only. There is no code fix that is still Ibragimov-Müller.

## 1. What is actually happening, measured

`skep-01` reports `se/sd ≈ 0.66-0.70` and coverage 0.82-0.86 at `i = unit, cluster(unit)` — the
geometry P1's screen **passes** and where IM is textbook-legitimate. I found the driver and it is a
lever the help can hand the user.

Every per-cluster read is *(that cluster's post-vs-pre)* minus *the shared never-treated trend*. The
second term is **common to all G reads**, so the cross-cluster dispersion — which is all IM looks at —
cannot see its sampling error. IM's `se` is therefore blind to how much control data exists, while
the estimator's actual sampling sd is not.

`p6_lever.do`, R=150 per arm, 20 treated units with a homogeneous effect of exactly 1.00, 10 periods,
`i = unit`, `cluster(unit)`, varying only the never-treated pool `Nc`:

| `Nc` | `Nc/G_t` | sd(b) | **mean(se)** | se/sd | 95% coverage | MC SE |
|---|---|---|---|---|---|---|
| 5 | 0.25 | 0.1656 | **0.0938** | 0.566 | **0.780** | 0.034 |
| 20 | 1 | 0.1153 | **0.0910** | 0.789 | **0.920** | 0.022 |
| 80 | 4 | 0.0916 | **0.0915** | 0.999 | **0.973** | 0.013 |

**`mean(se)` is flat to 3%** while `sd(b)` falls by 45%. That is the mechanism, not an inference from
it. Coverage at `Nc/G_t = 0.25` is 5 MC SE below nominal; at `Nc/G_t = 4` it is at or slightly above
nominal (the `t(G-1)` reference being conservative once the shared term is negligible). The
`sqrt(Nc/(Nc+G_t))` benchmark (0.447 / 0.707 / 0.894) tracks the measured `se/sd` in order and
roughly in size. An earlier R=40 pass (`p4`, `p5`) put the `Nc=20` ratio at 0.618 and 0.990 on two
seed blocks — **R=40 is far too noisy to quote here and I have not quoted it.**

## 2. `sthlp:583`'s first conjunct is inoperative in wooldid; the second I could not test

`sthlp:583` names a **conjunction**: IM's independence assumption "seems most likely to be satisfied
… in a bespoke diff-in-diff set up where clusters contain both treated and control units **and where
time fixed effects are not shared across clusters**." I tested both halves.

**Conjunct 1 does nothing here, structurally and measurably.** `p5_remedy.do` builds exactly it:
20 clusters, each holding one treated and one never-treated unit, against `cluster(unit)` on the same
draws. The two arms return **bit-identical** MC summaries (`mean(b) 0.99522`, `sd(b) 0.09376`,
`mean(se) 0.09283`, coverage 0.950). The reason is the IM `margins` subpop at `ado:2042`,
`subpop(if … & targetflag==1 & t >= ttre)`: a cluster's **never-treated members contribute nothing to
its read**, so no `cluster()` can give a cluster its own control arm. Putting controls inside clusters
and leaving them in their own clusters are the same computation.

**Conjunct 2 is untested and I will not claim it.** `p8b_conjunct.do` builds the design it describes —
10 regions x (2 treated + 4 never-treated) units, `i = unit`, `cluster(reg)` (which P1's screen
passes), arm 1 with wooldid's shared time effects and arm 2 adding `fe(regtime)` so time effects are
region-specific — R=100 each:

| arm | sd(b) | mean(se) | se/sd | coverage | MC SE |
|---|---|---|---|---|---|
| shared time FE | 0.1314 | 0.1248 | 0.950 | **0.960** | 0.020 |
| `fe(regtime)` | 0.1346 | 0.1444 | 1.073 | **0.960** | 0.020 |

The arm shows nothing because **the control arm was already calibrated**: 40 never-treated units
against 10 contributing clusters is `Nc/G_t = 4`, which §1's table says already covers ~0.97. This is
a null from an instrument that could not have detected the effect, and I am reporting it as one. A
design where the shared arm under-covers *and* region-by-time effects are estimable is what would
settle conjunct 2; at `i = unit` with one treated unit per cluster the region-by-time terms exhaust
the degrees of freedom, which is why I did not get one.

**What p8b does establish** is §1's lever on a second, independent geometry — `cluster()` coarser
than `i`, clusters holding both treated and control units, `Nc/G_t = 4`, coverage 0.960 against the
`Nc/G_t = 4` row's 0.973 at `cluster() == i`. Two different geometries, same lever.

## 3. Why no code change

The missing variance is the sampling error of the shared control trend. Supplying it means computing
an analytic covariance — which is the object `impvalues` exists to avoid — and the cross-cluster
t-test then no longer earns its `t(q-1)` reference. This is the same conclusion `sthlp:251` already
reached for the analytic path at or above `i`: *"There is no plug-in correction … the analytic
covariance is the wrong object rather than a mis-scaled one."* The project has ruled that shape three
times. **I recommend the same disposition here, and no code change.**

## 4. The help edits

**(a) `sthlp:1096-1097` is the one flatly false sentence and should change regardless of P1.**
It currently reads: "Valid ame inference under such designs: … **impvalues** — cross-cluster
dispersion inference, **which remains valid**". Proposed: "… `impvalues` — cross-cluster dispersion
inference, which is not disabled by the cell-containment screen; note that it carries its own
under-coverage, described under *Discussion of Ibragimov and Muller-Style Inference*, and requires a
`cluster()` that nests `i`."

**(b) `sthlp:584-586`** — the paragraph already says a shared control group "likely causes mischief
with the independence assumption". Give it the magnitude and the lever, in the style `sthlp:251`
already uses: *"How much mischief is measurable and it scales with the never-treated pool. The
cross-cluster dispersion cannot see the sampling error of the shared control trend, which every
per-cluster estimate subtracts, so the reported standard error does not respond to the size of that
pool while the estimator's own sampling variation does. Measured on a homogeneous-effect design with
20 treated clusters (150 draws): with a never-treated pool one quarter the number of treated clusters
the standard error is 0.57 of the estimator's standard deviation and a nominal 95% interval covers
0.78; at parity, 0.79 and 0.92; at four times, 1.00 and 0.97. There is no plug-in correction — the
missing term is the analytic covariance the method exists to avoid — so the remedy is a larger
never-treated pool, or `bootstrap()`/`permute()`."*

**(c) `sthlp:583`** — qualify the first conjunct. "Clusters contain both treated and control units"
is not something `cluster()` can deliver in wooldid: the per-cluster read is taken over treated
post-period observations only (`ado:2042`), so a cluster's control members contribute nothing to it
and the arm is bit-identical to leaving them in their own clusters (`p5_remedy.do`). Proposed: keep
the sentence but say what it requires — *"in wooldid the per-cluster estimate is read over that
cluster's treated post-treatment cells alone, so putting control units inside a cluster does not give
it its own control arm; what the condition asks for is a specification in which each cluster's
counterfactual is estimated from that cluster, e.g. cluster-by-time effects in `fe()`."* I did **not**
verify that the `fe()` route restores coverage — see §2 — so the sentence should describe the
requirement, not promise the remedy.

## 5. Optional, and I do not recommend it as primary

A caution when `Nc < G_t`. It is nearly free — the screen in P1 already tags `i`; two more counts
give the never-treated cohort count and the treated cluster count, ~0.1 s at 1M rows. But it is
exactly the shape T4 declined, the threshold is a judgement call, and (b) above delivers the same
information where a user reading about `impvalues` will meet it. Offered so the option is costed, not
pressed.

## 6. What I did not test

- Coverage under P1's patch on the *failing* geometry — there is nothing to measure once the run is
  refused, and under Variant B the numbers are unchanged by construction (proved in §3 of P1).
- `impvalues` + `contreat()` (the `ame` sites at 2320/2375/2715). The screen is upstream of all of
  them and geometry-only, so it applies unchanged, but I ran no continuous-treatment case.
- The `poisson` path.
- `sthlp:583`'s second conjunct (time effects not shared across clusters). `p8b_conjunct.do` ran but
  on a design that was already calibrated, so the arm had no room to move. Stated as untested.
- The behaviour of `adm_impv` itself — that would require running the suite, which was off-limits.
  Its geometry is established from the generator source, not from a run.
