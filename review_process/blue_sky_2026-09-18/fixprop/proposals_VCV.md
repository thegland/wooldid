# Fix proposals — mechanism `VCV` (e(V) construction / rank / labelling)

Pinned build under review: `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build/wooldid.ado` (sha c1d8979, 17021 lines).
Patched copy actually run: `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-VCV/build2/wooldid.ado`.
Probes (all run, all logs kept):
- `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-VCV/t1_label.do` / `.log`  — SKEP-01, patched
- `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-VCV/t2_boot.do` / `.log`   — skep-dfmismatch-03, patched, with an in-run before/after oracle
- `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-VCV/t3_boot.do` / `.log`   — t2 re-run with the final caution wording
- "before" numbers come from the finding's own pinned-build log
  `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/r4_collapse.log`

Three proposals, split by whether they move a number rather than by which finding
they came from:

| # | what | finding | moves a number? | ship |
|---|------|---------|-----------------|------|
| **VCV-1** | rescale the listwise bootstrap covariance to the per-column bootstrap SDs | skep-dfmismatch-03 | **yes** — `e(V)`, `e(V_boot)` (nothing wooldid prints) | needs a ruling |
| **VCV-2** | say what `e(boot_N_V)` was; caution when it collapses | skep-dfmismatch-03 | no | shippable now |
| **VCV-3** | test `unconditionalse` first in the `e(vce)` ladder | SKEP-01 | no (two macros change value, one is added) | shippable now |

---

## A contradiction found while writing these (report it as a finding in its own right)

The A11 ruling that `sqrt(diag(e(V_boot)))` need not equal the printed `se` column
lives **only in a code comment** (`ado:14503-14506`, `ado:9640-9644`: "A11 ruled this a
documented consequence rather than a defect"). **The help file promises the opposite,
in two places**:

- `sthlp:1410`, `e(wdid_vcvstate)`: *"`bootstrap` (the full replication covariance, with
  `e(df_r)` aligned so a hand-run `test` reproduces the `bootdisplay(wald)` table)"*
- the bootstrap `e()` block: *"`e(df_r)` is aligned with the reference distribution
  actually reported, so a `test` or `lincom` run by hand after a bootstrap fit
  **reproduces the `bootdisplay(wald)` p-values under every `bootdf()` setting**"*

On the pinned build that is false on any run that loses a single replication cell.
Measured (r4_collapse.log, the finding's own run): printed `se` for `main_att`
`.07589904`; `lincom main_att` returns `.07393592`, and across all 107 reported
estimates `sqrt(diag(e(V)))/se` ranges **0.0638 to 1.6496, mean 0.7113, with 82 of 107
columns below 0.90** (t2_boot.log, `T2_BEFORE_RATIO`). A ratio below 1 means the exported
standard error is **smaller** than the one on the screen, i.e. anti-conservative, and that
is the direction for 82 of the 107 reported estimates. So the promise is not merely
approximate — the exported standard error is on average 29% below the printed one and on
one column 94% below it. Per the coordinator note, a measured contradiction of an
existing ruling is a first-class finding; VCV-1 is the remedy that makes the sentence
true rather than the remedy that deletes it.

The finding itself offered three routes: *"Either `e(V)` should rest on the replications
the printed SEs rest on, or `e(df_r)` should describe the covariance, or the run should
say what `e(boot_N_V)` was."* VCV-1 is route 1 (as exactly as a PSD matrix permits),
VCV-2 is route 3. Route 2 is rejected below, with reasons.

---

# VCV-1 — make `sqrt(diag(e(V)))` be the bootstrap standard error the table reports

**Findings:** skep-dfmismatch-03.

### The mechanism, restated in one paragraph

`wooldid_bootfin` (`ado:14507-14580`) runs **two moment passes on different row sets**.
`se[k]` is the standard deviation of column *k* over the replications in which column *k*
alone is available (`cnt[k]`, posted as `e(boot_navail)`). `V` is the **listwise
complete-case** covariance over the rows complete across *every* non-excluded column
(`nlw`, posted as `e(boot_N_V)`). The printed table gets `se`; `e(V)` and `e(V_boot)` get
`V`. With `subgroups()` and a long event window the intersection collapses: 6 subgroups x
16 horizons = 107 reported estimates, 60 completed replications, **`e(boot_N_V)` = 4**.
Every column is individually available in nearly every replication — no `wdb_lowav`
caution fires — yet only 4 rows are complete across all 107.

### The edit

`wooldid.ado`, inside `void wooldid_bootfin()`, **immediately before the sentinel loop at
`ado:14553** (`for (k = 1; k <= K; k++) {` whose body writes `-1e100`):

```mata
  // VCV-1: rescale the listwise covariance so its DIAGONAL is the per-column bootstrap
  // standard error the table prints. The two moment passes above run on different row
  // sets, so sqrt(diag(V)) was not the reported se column and the help's promise that a
  // hand test or lincom reproduces the bootdisplay(wald) numbers was false whenever any
  // replication cell was lost. Rescaling by a positive diagonal is a congruence
  // transform: it preserves every listwise correlation and preserves positive
  // semi-definiteness. Columns with no listwise variance keep V[k,k] == 0 and fall into
  // the sentinel branch below unchanged.
  real rowvector vsc
  if (nlw >= 2) {
    vsc = J(1, K, 1)
    for (k = 1; k <= K; k++) {
      if (excl[k] == 0 & V[k,k] > 0 & se[k] < . & se[k] > 0) vsc[k] = se[k]/sqrt(V[k,k])
    }
    V = V :* (vsc' * vsc)
  }
```

`V :* (vsc' * vsc)` is `diag(vsc) V diag(vsc)`. If the K x K temporary the outer product
allocates is worth avoiding, `V = (V :* vsc') :* vsc` broadcasts instead and is
arithmetically identical; at K in the hundreds neither is measurable. Overflow is not a
concern either way: Cauchy-Schwarz bounds every rescaled off-diagonal by `se_k * se_j`.
Nothing else changes: `se`, `mn`, `cnt`,
`nlw`, `wdb_lowav`, `wdb_exclnames` are all untouched, and the sentinel loop that follows
runs on exactly the columns it ran on before (a column with `V[k,k] == 0` is left at 0 by
the guard and still gets `-1e100`). `bV` is read at only two places (`ado:10042`,
`ado:10063`), so the rescaled matrix reaches `e(V)` and `e(V_boot)` and nothing else.

### Measured before / after (t2_boot.log, patched build, r4_collapse design verbatim)

The run posts `e(bootdist)` — the raw 60 x 107 replication matrix — so the **unpatched
covariance can be rebuilt inside the same run**. That gives an exact before/after from one
launch and needs no second build.

| check | result |
|---|---|
| reconstruction is faithful (`max abs(e(wdidmainresults)[,2] - se_rebuilt)`) | **exactly 0** |
| patched `e(V)` equals `diag(sc) V_listwise diag(sc)` (max abs rel diff) | **exactly 0** |
| **invariant: `max abs(sqrt(diag(e(V))) - se)`** | **5.6e-17** (1 ulp) |
| correlations preserved (`max abs(R_old - R_new)`) | **4.4e-16** |
| PSD preserved: min/max eigenvalue | before −4.0e-16 / 0.929; after −8.3e-16 / 1.908 |
| rank | 3 before, 3 after (congruence preserves rank — VCV-1 does **not** fix rank) |
| **printed main table** | **byte-identical to the pinned run**, every row |
| `lincom main_att` se | before `.07393592`; after `.07589904` = the table's se **bitwise** |
| `lincom main_att` CI | before `[.41470132, .71489712]`; after `[.4107160, .7188825]` (hex `+1.a492bb2324da0X-002` / `+1.70115c8380b9aX-001`), which is exactly `b ± invt(35,.975)·se` = `.56479922 ± .15408325` — the documented Wald family. (The printed table's `lb_95ci` is `.44331301`, the PERCENTILE interval: `bootdisplay()` defaults to percentile, and the footer already says `e(V)` carries the standard-deviation family, so the interval a hand `lincom` must reproduce is `e(boot_wald_mainresults)`'s, not the printed one. The `se` column is the same in all three families, which is why it matches bitwise.) |
| `test ES_att_pre2 ES_att_pre3` | before F(2,35) p = **.00136**; after p = **.40448** |

That last row is worth dwelling on: the design has **no pre-trend by construction**, and
the pinned build's hand-run pre-trend test rejects at p = .0014. After the rescale it does
not. That is not a proof the test is now calibrated (see the honest caveat below) — it is
one draw — but it is the direction the mechanism predicts.

### Does it move a reported number?

**No number wooldid prints moves.** Verified byte-identical on the whole main table, and
structurally: nothing inside wooldid reads the rewired `e(V)`. The only five readers of
`e(V)` in the file are `wooldid_vsplice` (`ado:11047`, analytic path, skipped in
`wdrepmode`), `wooldid_amescreen` (`ado:12976`, analytic path), `wooldid_bootgrab`
(`ado:14433`, reads the *nested replication's own* `e(V)`), `wooldid_permcse`
(`ado:16584`, permute path) and `wooldid_bootwire`'s own preservation of the analytic
matrix (`ado:10039`). wooldid's own bootstrap joint and equality tests are computed by
`wooldid_bootqf` (`ado:6845`), which builds **its own** listwise covariance over only the
tested columns straight from the accumulator; `e(V_es_*)` is rebuilt by
`wooldid_bootves` the same way; `pretrendpower()` reads `e(V_es_*)`, not `e(V)`. All
confirmed by reading and then empirically by the identical table.

`e(V_boot)` likewise: `grep -n "V_boot" wooldid.ado` returns the posting at `ado:10066`,
the "could not be posted" caution at `ado:10056`, the footer at `ado:11787` and four
comments. **Nothing in the package reads it.**

`e(boot_wald_mainresults)` — which the regression test below keys on — is `wmsd`, and
`wmsd` is snapshotted from `wm` at `ado:9744`, immediately after the loop at
`ado:9705-9719` that writes `bse` (the bootstrap SD) into column 2, and **before** the
percentile/bc/bca substitutions at `ado:9820-9822` and the `boottstat` substitution at
`ado:9841-9880`. So its column 2 is the bootstrap standard deviation in every display
mode, `boottstat` included. That is read, not run: I did not execute a `boottstat` fit.

**What does move: `e(V)` and `e(V_boot)`, and therefore everything downstream of them** —
`esttab`/`estout` standard errors, `coefplot` bars, a hand `lincom`/`test`, and any
`honestdid` call that is pointed at `e(V)` rather than `e(V_es_*)`. On my run those move
by a factor between 1/1.6496 = 0.61 and 1/0.0638 = **15.7**, mean factor 1.41, always
toward the printed table. The author must rule, and any baseline that blesses `e(V)` or
`e(V_boot)` needs a re-bless.

### Honest caveat — what VCV-1 does *not* fix

The **correlations** still rest on `e(boot_N_V)` rows, so `rank(e(V))` is still at most
`e(boot_N_V) - 1` (3 on my run) and a hand `test` of q > 1 restrictions is still resting
on 4 draws. VCV-1 makes every **q = 1** statement exact — which is what the help
actually promises, since `e(boot_wald_mainresults)` is a per-coefficient table — and
leaves q > 1 to VCV-2's disclosure. It is not a substitute for VCV-2 and should ship with
it.

### Statement of the invariant (the correct one)

`sqrt(diag(e(V))) == e(boot_wald_mainresults)[,2]`, the **bootstrap standard deviation**.
That is the printed `se` column on every display *except* `boottstat`, where the printed
column is the **analytic** se (`ado:9871`, labelled `se_analytic`). My t2 check compares
against `e(wdidmainresults)[,2]` and is valid only because that run is not `boottstat`; a
regression test should key on `e(boot_wald_mainresults)`, which is stored whichever family
is displayed.

### Alternatives rejected

1. **Change the printed `se` to the listwise one.** Rejected: it throws away 56 of 60
   replications for every reported standard error, and it *does* move a printed number.
2. **Pairwise-complete covariance.** Rejected by A11 already ("measurably indefinite"), and
   VCV-1 gets the agreeing diagonal without reintroducing indefiniteness.
3. **Blank the off-diagonals when `e(boot_N_V)` is small**, leaving `diag(se^2)`. Rejected:
   `sthlp:1410` advertises the `bootstrap` vcvstate as *the full replication covariance* —
   that is a documented feature over the analytic `analytic-diagonal` state, and silently
   demoting it would be a second doc-vs-code divergence.
4. **Leave the code and fix the help instead** (delete the "reproduces the
   `bootdisplay(wald)` p-values" sentence). Defensible and cheap — but it documents a
   surprise rather than removing it, and it leaves `esttab` quietly publishing a standard
   error 29% away from the one on the screen with the same name. If the author prefers
   this, the help edits are in §Help below and VCV-2 still ships.

---

# VCV-2 — disclose `e(boot_N_V)`; caution when the covariance rests on under half the replications

**Findings:** skep-dfmismatch-03. **Moves no number.** Two edits.

### Edit 2a (one line): the footer names the macro but never its value

`ado:11786` is one `dis as result` with the whole sentence inside one pair of quotes, so
it prints the literal text `e(boot_N_V)`:

```stata
* NOW
   if e(boot_v_posted) == 1 dis as result "  e(V) holds the bootstrap covariance (complete-replication listwise, e(boot_N_V) replications; the analytic covariance, when one existed, remains in the hidden e(V_analytic))."
* BECOMES
   if e(boot_v_posted) == 1 dis as result "  e(V) holds the bootstrap covariance, built on the " e(boot_N_V) " replication(s) complete across every reported estimate, of " e(boot_completed) " completed (e(boot_N_V); its rank is at most one less). The analytic covariance, when one existed, remains in the hidden e(V_analytic)."
```

Measured, patched run: `e(V) holds the bootstrap covariance, built on the 4
replication(s) complete across every reported estimate, of 60 completed (e(boot_N_V); its
rank is at most one less).` Compare the pinned run, which printed `(complete-replication
listwise, e(boot_N_V) replications` — the name, never the number. This is the single
highest value-per-character edit in the whole set: the run already computes and posts the
number, it just never shows it. Line 11784 in the same block shows the macro is normally
concatenated outside the quotes (`... " completed replications, ..."`), so this looks like
a plain slip rather than a decision.

### Edit 2b: a caution and a replay-surviving flag when it collapses

After `ado:10248-10249` (`ereturn scalar boot_N_V` / `boot_v_posted`), in
`wooldid_bootwire`:

```stata
  * VCV-2: the listwise rule intersects over EVERY non-excluded column, so a run with many
  * columns (subgroups() x a long event window) can reach e(V) on a small fraction of its
  * completed replications while NO per-column low-availability caution fires - each column
  * can be available in 58 of 60 and still leave an intersection of 4. Disclose the count,
  * and leave a flag that survives into a replayed estimate after the caution has scrolled
  * away. THRESHOLD: the package's own low-availability rule (fewer than half, ado:14528),
  * so an ordinary run that loses a cell or two says nothing extra - Edit 2a already prints
  * the count on every run - and only an e(V) resting on under half its replications speaks.
  local wdb_done0 = cond("`done'" == "", 0, real("`done'"))
  if `wdb_vposted' & `wdb_nlw' < 0.5*`wdb_done0' {
    ereturn scalar boot_V_shortfall = 1
    dis as text "Caution: e(V) rests on `wdb_nlw' of the `wdb_done0' completed replications - the rows complete across all `K' reported estimates (e(boot_N_V)) - so its rank is at most `=max(`wdb_nlw'-1, 0)' and a hand-run test of more restrictions than that is not identified. `=cond(`wdbdf' < ., "e(df_r) = `wdbdf' describes the reference distribution of the printed table (the resampling-cluster count under bootdf(clusters)), not this covariance, so a test of fewer restrictions than that rank is anti-conservative relative to its nominal size.", "The reported statistics are normal/chi-squared here, which takes no account of how few replications the covariance rests on, so such a test is anti-conservative relative to its nominal size.")' The printed standard errors are NOT affected: each rests on its own e(boot_navail) replications."
  }
```

`K`, `wdbdf`, `wdb_nlw`, `wdb_vposted` and `done` are all already in scope in
`wooldid_bootwire`. The `wdb_done0` local is defensive only (`done()` is declared
`string`). The `cond()` on `wdbdf` matters: under `bootdf(normal)` `e(df_r)` is
deliberately blanked (`ado:10245`), so the "e(df_r) describes..." clause would be wrong.

This **adds one `e()` entry** (`e(boot_V_shortfall)`). No number moves, but any
ereturn-contract listing the project maintains needs the new name added.

Measured, patched run: fires, with `e(boot_V_shortfall) = 1`, text
`Caution: e(V) rests on 4 of the 60 completed replications - the rows complete across all
107 reported estimates (e(boot_N_V)) - so its rank is at most 3 ...`.

### Cost

Edit 2a: nothing — it is a line that already runs. Edit 2b: two scalar comparisons and, on
the rare firing run, one `display`. **No pass over the data, no matrix work, nothing that
scales with N, T or reps.** The author's "no guards that are expensive on large datasets"
constraint is not engaged at all.

### The size of the thing being disclosed

The finding's companion Monte Carlo (`r6_df.do`, pure simulation, no wooldid) prices what
the disclosure is about: a covariance from `N_V` draws referred to F/t with `df_r`
rejects at nominal 5% as follows — (4, 35): q=1 13.7%, q=2 31.2%, CI coverage .8647;
(6, 35): q=2 18.7%, q=4 53.8%; (14, 29): q=6 29.4%; (30, 35): q=6 10.8%; (60, 35): q=6
6.1%. A user with `e(boot_N_V) = 4` and `e(df_r) = 35` is running a nominal-5% test at a
true 14-31%.

---

# VCV-3 — test `unconditionalse` first in the `e(vce)` / `e(vcetype)` ladder

**Findings:** SKEP-01. **Moves no reported number** (verified bitwise).

### The edit

`ado:4286-4302`. The ladder currently tests `cluster` → `robust` → `customvce` →
`unconditionalse` → default. `ado:4676` refuses `unconditionalse` unless the user supplies
`robust`, `cluster()` or `customvce()` — *unless* the weight is a pweight. So the
`unconditionalse` arm at `ado:4300` is reachable by exactly one call shape, and every other
`unconditionalse` user is handed the stage-1 regression's label for a
`margins, vce(unconditional)` matrix. `ado:4283`'s own comment states the contract being
broken: *"e(vce) and e(vcetype) name what e(V) IS"*.

```stata
* NOW (ado:4286-4302)
  if `wdhasv' {
    if "`cluster'" != "" {
      ereturn local vce cluster
      ereturn local vcetype Robust
    }
    else if "`robust'" != "" {
      ...
    else if "`unconditionalse'" != "" ereturn local vce unconditional
    else ereturn local vce = cond("`poisson'" == "poisson", "oim", "ols")
  }

* BECOMES
  if `wdhasv' {
    * SKEP-01 fix: unconditionalse is tested FIRST. On that path e(V) is margins'
    * vce(unconditional) matrix whatever vce the stage-1 regress carried, and the gate
    * above forces every unconditionalse user but a bare-pweight one to supply
    * robust/cluster()/customvce(), so testing those first labelled the unconditional
    * covariance with the stage-1 vce on essentially every run. The stage-1 spelling is
    * kept as e(vce_stage1) and the design in e(clustvar): nothing readable before is lost.
    if "`unconditionalse'" != "" {
      ereturn local vce unconditional
      ereturn local vcetype Unconditional
      local wdvst1 = cond("`cluster'" != "", "cluster", cond("`robust'" != "", "robust", cond(`"`customvce'"' != "", word(`"`customvce'"', 1), cond("`weight'" == "pweight", "robust", ""))))
      if "`wdvst1'" != "" ereturn local vce_stage1 `wdvst1'
    }
    else if "`cluster'" != "" {
      ereturn local vce cluster
      ereturn local vcetype Robust
    }
    else if "`robust'" != "" {
      ... (unchanged)
    else ereturn local vce = cond("`poisson'" == "poisson", "oim", "ols")   /* the dead
        `else if "`unconditionalse'" != ""' arm is removed: it is now unreachable */
  }
```

### `e(vcetype)` = `Unconditional` is Stata's own word, not an invention

t1_label.log runs `margins, dydx(year) vce(unconditional)` directly and Stata prints

```
             |            Unconditional
             |      dy/dx   std. err.      t    P>|t|     [95% conf. interval]
```

so `esttab` will head the column "Unconditional" exactly as `margins` does. The pinned
build gave `Robust`, which is the header of a matrix this is not.

### Measured (t1_label.log, patched build; the "before" is the finding's own log)

| case | before `e(vce)`/`e(vcetype)` | after | `sqrt(e(V)[1,1])` |
|---|---|---|---|
| (i) `[pw=w], unconditionalse robust` | robust / Robust | **unconditional / Unconditional**, `e(vce_stage1)=robust` | `+1.246441c7f884cX-002` — **bitwise unchanged** |
| (iii) `[pw=w], unconditionalse` (bare) | unconditional / . | unconditional / **Unconditional**, stage1 `robust` | `+1.246441c7f884cX-002` — the same matrix as (i), now with the same name |
| (iv) `[aw=w], unconditionalse robust` | robust / Robust | unconditional / Unconditional | `+1.2565132e60078X-002` unchanged |
| **(vi) `unconditionalse cluster(state)`, no weight — the ordinary user** | cluster / Robust | **unconditional / Unconditional**, stage1 `cluster`, `e(clustvar)=state` | `+1.3225d959a9944X-002` unchanged |
| (ii)/(v) plain `robust` | robust / Robust | **unchanged** | unchanged |
| (vii) plain `cluster(state)` | cluster / Robust | **unchanged** | unchanged |
| (viii) default | ols / . | **unchanged** | — |
| (ix) `customvce(robust)` | robust / Robust | **unchanged** | — |

So after the fix: one covariance, one name. The `robust` label no longer sits on two
matrices 17.2% apart, and the clustering information a consumer might have keyed on is not
lost — `e(clustvar)` was already posted unconditionally at `ado:4285`, and `e(vce_stage1)`
is new.

### Cost, blast radius, alternatives

Cost: three `cond()` evaluations on macros, once per run. Zero data cost.

Blast radius: `e(vce)` and `e(vcetype)` on `unconditionalse` runs only, plus one new macro
(`e(vce_stage1)`) — **an added `e()` entry, so an ereturn-contract listing needs updating;
no number moves.**
An existing script that branches on `e(vce)=="cluster"` to decide "this run was clustered"
will now take the other branch on an `unconditionalse cluster()` run — it should read
`e(clustvar)`, which is what that macro is for and which is unchanged. An `esttab` table
header flips from `Robust` to `Unconditional`. Nothing numeric moves.

Rejected alternatives: (a) leave the code and document the precedence in `sthlp:1411` —
cheaper still, but it documents a broken contract that `ado:4283` itself states, and the
help currently lists `unconditional` as a value `e(vce)` takes without saying it is nearly
unreachable; (b) post a composite `e(vce) = "unconditional cluster"` — no consumer parses
two words; (c) leave `e(vcetype)` empty on the unconditional arm as the pinned build does
in case (iii) — rejected because `margins` itself has a word for it and an empty
`e(vcetype)` makes `esttab` print the default OLS header.

### Note on the brief's reviewer verdicts for SKEP-01

The verdict block under SKEP-01 in `mech/VCV.md` carries five verdicts, three of which are
about other findings entirely (cluster-at-`i` interval coverage; `permute()` +
`contreat()` pre-period dose-response rows). The two that are about SKEP-01 are the two
`sev low` ones, and both say the same thing: the label is wrong and no reported number
changes. That is what VCV-3 fixes, and its severity is correspondingly low — it is here
because it is a two-minute edit with a measured zero-risk profile, not because it is
urgent.

---

# Not proposed — and why

### Changing `e(df_r)` to describe the covariance (the finding's route 2). **Rejected.**

Setting `e(df_r) = min(wdbdf, e(boot_N_V) - 1)` is tempting and wrong:

1. It **breaks the documented alignment it is meant to restore.** `ado:10234-10236` sets
   `e(df_r) = wdbdf` precisely so a hand `test`/`lincom` reproduces the printed
   `bootdisplay(wald)` p. With VCV-1 shipped, `se` *and* `df` then both match and the
   q = 1 case reproduces exactly (measured: `lincom` CI = `b ± t(35,.975)·se`). Lowering
   `df_r` would re-break it in the opposite direction.
2. `e(df_r)` **cannot encode what is actually needed.** The correct reference for a Wald
   statistic on a covariance estimated from *n* draws with *q* restrictions is Hotelling's
   `(n-q)/((n-1)q) · T² ~ F(q, n-q)` — it depends on *q*, and `test` reads a single fixed
   `e(df_r)` whatever *q* is. There is no value of `e(df_r)` that is right for q = 1 and
   q = 6 simultaneously.
3. It would move numbers on **every** bootstrap run, including well-behaved ones, for a
   pathology that VCV-2's caution names precisely.

Disclosure is the only honest q > 1 remedy. That is VCV-2.

### Applying the same rescale to `e(V_es_*)`. **Named, not built — needs a ruling.**

`wooldid_bootves` has exactly the same two-row-set problem: its blocks are listwise over
each family's own columns, while the printed ES rows carry per-column SDs. Its listwise
count is *larger* than the global one (it intersects over ~16 columns, not 107), so the
damage is much smaller — but it is the same defect, and `e(V_es_*)` is what `honestdid`
and `pretrendpower()` consume. Rescaling it **would move a reported number**
(`pretrendpower()`'s power, detectable trend and slope figures print in the table), which
is why I have neither built nor measured it. It should be a separate ruled item.

### Widening the listwise base adaptively ("greedy column drop"). **One paragraph, not built.**

If the rank of `e(V)` matters to the author as much as its scale, the root cause is that
the intersection is taken over *all* 107 columns at once. A greedy rule — drop from the
intersection the column whose absence recovers the most rows, repeat until
`nlw >= f·done`, and give the dropped columns a diagonal-only entry (zero off-diagonals,
which is already the package's "not estimated" convention) — would take `e(boot_N_V)` from
4 to something near `done` while keeping the matrix PSD (block-PSD plus a positive
diagonal). Cost: a handful of passes over a `reps x K` missingness matrix, microseconds.
I did not build it because it needs a ruling on two things it cannot decide for itself:
the target fraction `f`, and a tie-break that makes the column order **deterministic**
across runs and seeds (a nondeterministic `e(V)` would be worse than a rank-3 one).

---

# Help edits that go with these

Each is required by a code edit above; none stands alone.

1. `sthlp:1410`, `e(wdid_vcvstate)`, the `bootstrap` clause. **Now:** *"`bootstrap` (the
   full replication covariance, with `e(df_r)` aligned so a hand-run `test` reproduces the
   `bootdisplay(wald)` table)"*. **Becomes:** *"`bootstrap` (the replication covariance:
   listwise-complete correlations over the `e(boot_N_V)` replications complete across every
   reported estimate, scaled so that each diagonal is that estimate's own bootstrap
   standard error — the one the `bootdisplay(wald)` table reports — with `e(df_r)` aligned
   so a hand-run single-coefficient `test` or `lincom` reproduces that table exactly. Its
   rank is at most `e(boot_N_V) - 1`, which can be far below the number of reported
   estimates when `subgroups()` and a long event window are combined; a joint hand-run
   `test` of more restrictions than that rank is not identified, and the run cautions when
   `e(boot_N_V)` falls below half the completed replications)"*.
1b. `sthlp:1409` carries the same phrase and must move with it: *"the blocks are in the
   hidden `e(V_analytic)` while `e(V)` holds the **full replication covariance**"* →
   *"... while `e(V)` holds the replication covariance described under
   `e(wdid_vcvstate)`"*. Fixing 1410 and leaving 1409 would leave the twin standing.
2. The bootstrap `e()` block sentence *"so a `test` or `lincom` run by hand after a
   bootstrap fit reproduces the `bootdisplay(wald)` p-values under every `bootdf()`
   setting"* — add **"single-coefficient"** before `test`, and a following sentence naming
   `e(boot_N_V)` as the rank bound for joint tests.
3. Same block, `e(boot_N_V)` "(replications behind the covariance)" — add "and therefore
   its rank, one less; the per-estimate counts behind the reported standard errors are
   `e(boot_navail)`, which are generally larger".
4. New `e()` entries to document: scalar `e(boot_V_shortfall)` (VCV-2), macro
   `e(vce_stage1)` (VCV-3).
5. `sthlp:1411`: *"`e(vcetype)` is `Robust` or `Bootstrap` where one of those fits"* →
   *"`Robust`, `Bootstrap` or `Unconditional`"*, and add a sentence: *"`unconditionalse`
   is read first, because on that path `e(V)` is `margins`' `vce(unconditional)` matrix
   whatever variance the stage-1 regression carried; the stage-1 spelling is in
   `e(vce_stage1)` and the clustering design in `e(clustvar)`."*

---

# Tests that would prove these correct (oracle-free where possible)

- **VCV-1, invariant, no oracle:** on any `bootstrap()` fit, stack column 2 of
  `e(boot_wald_mainresults)` and then of `e(boot_wald_sgresults_<sg>)` for each `sg` in
  `e(wdid_sglist)` order — that stacking is exactly how `e(b)` is built (`ado:9595-9605`; the per-subgroup Wald tables are posted at `ado:9991`)
  — and require
  `max(abs(sqrt(diagonal(e(V)))' - that vector))` to be 0 to within an ulp for every
  non-sentinel column. Fails on the pinned build on any run that loses a cell; passes on
  the patched one (measured 5.6e-17, i.e. 2^-54). Key it on `e(boot_wald_*)`, not the
  printed tables, so it also holds under `boottstat`, where the printed column is the
  analytic se. **Scope of what I actually ran:** t2/t3 compared the reconstruction to the
  printed table for the `nmain` main rows only (exact, 0), and validated the remaining ~90
  subgroup columns against the `e(bootdist)` reconstruction rather than against the printed
  subgroup tables. The regression test should cover both.
- **VCV-1, in-run oracle:** rebuild the unpatched covariance from `e(bootdist)` and check
  (a) `e(V)` equals `diag(sc)·V_listwise·diag(sc)` exactly, (b) the correlation matrices
  agree to machine precision, (c) `min eigenvalue >= -tol·max eigenvalue`. This is the
  t2/t3 probe and it is self-contained — it needs no second build and no external package.
- **VCV-1, metamorphic, no oracle:** run the same specification with and without a
  `subgroups()` variable that costs replication cells. The diagonal identity must hold in
  both; on the pinned build it holds only in the run that loses no cell.
- **VCV-1, cross-build identity:** the printed tables must be byte-identical before and
  after. Measured: they are.
- **VCV-2:** `e(boot_V_shortfall) == 1` iff `e(boot_N_V) < 0.5*e(boot_completed)`, and the
  footer's printed count equals `e(boot_N_V)`. Both checkable on one run.
- **VCV-3:** for every `unconditionalse` call shape, `e(vce) == "unconditional"`; and
  `e(V)` is **bitwise identical** across adding or removing the redundant `robust` under a
  pweight (`%21x` comparison). Measured on six shapes plus three control shapes whose
  labels must not move.

# What I could not test

- The `boottstat` interaction with VCV-1. The code path is the same (`bse` is the
  bootstrap SD in both modes and the rescale runs before any display branch), and the
  invariant is stated against `e(boot_wald_mainresults)` for exactly this reason, but I
  did not run a `boottstat` fit — it is materially slower, since every replication must
  estimate its own variance.
- Coverage of the *corrected* intervals. VCV-1 changes what `esttab`/`lincom` report; I
  measured that they now agree with the printed table, not that the printed table is
  itself well calibrated. The `bootstrap(standard)` over-coverage at `cluster() == i`
  (known_bydesign §7) is a separate, unruled matter.
- Whether the project has a blessed baseline containing `e(V)` or `e(V_boot)` for a
  bootstrap run. If it does, VCV-1 needs a re-bless; I am forbidden the test suite and
  could not check.
