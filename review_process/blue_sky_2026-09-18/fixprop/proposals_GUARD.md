# Fix proposals — mechanism `GUARD` (safety guards and preflights)

Pinned build `.../bluesky/build/wooldid.ado` (sha c1d8979, md5 444892481be8d50ec962a102d6abef51).
All line numbers are that build's.

**Every proposal below was patched into a private copy of the ado and RUN.** Three full builds
were exercised on identical do-file bodies, plus three targeted runs:

| build | contents | log |
|---|---|---|
| `base` | pinned, unmodified | `probes/prop-GUARD/run_base.log` |
| `p1b` | **G1a + G2 + G3a + G4 — the shipping candidate** | `probes/prop-GUARD/run_p1b.log`, `run_ctrl0b.log` |
| `p1` | an earlier, REJECTED form of G3a (treated-only assert) | `probes/prop-GUARD/run_p1.log`, `run_ctrl0.log` |
| `p2` | G1a + G2 + G3b + G4 | `probes/prop-GUARD/run_p2.log` |
| `p3` | p1 + G1b | `probes/prop-GUARD/run_p3.log` |
| — | cost of the pass G1a removes, 2,000,000 rows | `probes/prop-GUARD/run_cost.log` |
| — | true variable-allocation multiplier | `probes/prop-GUARD/run_bound.log`, `run_bound2.log` |
| — | a never-treated unit labelled 0, on the pinned build | `probes/prop-GUARD/run_ctrl0base.log` |

`which wooldid` is printed at the head of every log and resolves to the intended directory in
every one (`run_base.log` -> `build/`, `run_p1.log` -> `p1/`, `run_p2.log` -> `p2/`,
`run_p3.log` -> `p3/`, `run_p1b.log` -> `p1b/`).

Patcher (line-anchored, asserts the original text before replacing): `probes/prop-GUARD/patch.py`.
Probe body: `probes/prop-GUARD/body.do`; drivers `run_base.do` / `run_p1.do` / `run_p2.do`.

**Headline no-change control — the rendered tables, not just the numbers.** An ordinary
120-unit x 9-period panel, four specifications (`att att_i att_it att_itime` + `summarystats`;
`customvce(cluster id)`; `contreatwithin` + `contreat()`; `safety(off)`), plus the two clean legs
of the izero panel. For each, the ENTIRE log region wooldid prints — headers, cautions, the
`matlist` tables including the `relyr` column and its missing values, the footers — was diffed
line by line:

```
P4A   base vs p1b / p2 : RENDERED-IDENTICAL  226 lines
P4B   base vs p1b / p2 : RENDERED-IDENTICAL  144 lines
P4C   base vs p1b / p2 : RENDERED-IDENTICAL  232 lines
P4D   base vs p1b / p2 : RENDERED-IDENTICAL  124 lines
P3good base vs p1b / p2: RENDERED-IDENTICAL  234 lines
P2_B_BIN  base vs p1b  : RENDERED-IDENTICAL  149 lines
P2_B_CTW  base vs p1b  : RENDERED-IDENTICAL  273 lines
```
1,382 lines of rendered output, byte for byte, on seven specifications across three builds.
Nothing in this document changes a character on data that does not contain the defect.

---

# G1 — the maxvar preflight over-charges for absorbed fixed effects
*covers:* **SAFETY-01** and **sk-maxvar-fp** (one finding seen twice; one fix site)

### Confirmed on the pinned build first

`run_base.log`: 5,200 units x 8 years (41,600 rows), 20 ever-treated, **Stata's untouched
default `c(maxvar) = 5000`**, `wooldid y unit year ttre, att`:

```
Error: this specification requires at least 5458 variables (about 50 treated
cohort-by-period interaction terms) but maxvar is currently 5000, so estimation
would fail with r(907).                                          -> rc 907
```
The same command with `safety(off)`: **rc 0**, N = 41600, att `+1.7263d0a337fdcX+000`.
`5458 = 50 (cells) + 5200 (si8ilev) + 8 (si8tlev) + 200`. The deciding term is the unit count.

### G1a — engine-conditional demand  (PRIMARY, tested)

**File / site:** `wooldid.ado:1365-1369`, inside the SI-8 preflight block.

**Now:**
```stata
    qui gegen byte `si8tag' = tag(`t') if `ifstatement'
    qui count if `si8tag' == 1
    local si8tlev = r(N)
    cap drop `si8tag'
    local si8demand = `si8cells'*(1 + ("`contreat'" != "")) + `si8ilev' + `si8tlev' + 200
```

**Becomes:**
```stata
    * ENGINE-CONDITIONAL (guard false-positive fix): the levels of `i' and `t' are real
    * design-matrix columns ONLY on the reg/poisson path (:1485, which writes i.`i' i.`t'
    * out as regressors). On the DEFAULT path they are handed to reghdfe/ppmlhdfe's
    * absorb() (:1476/:1462) and never materialise, so charging the run for them is not a
    * lower bound: it over-predicts by roughly the unit count - one for every NEVER-TREATED
    * unit, the dimension this block's own caution calls 'comparatively cheap' - and refuses
    * specifications that estimate in seconds at that very maxvar. The predicate below is
    * the same one that selects the engine at :1457; the two must stay in step. Counting the
    * t levels is also skipped on the absorbing path, so this is one gegen pass CHEAPER than
    * before on the path essentially every run takes.
    local si8levels 0
    if "`unconditionalse'" != "" | "`customvce'" != "" {
      qui gegen byte `si8tag' = tag(`t') if `ifstatement'
      qui count if `si8tag' == 1
      local si8tlev = r(N)
      cap drop `si8tag'
      local si8levels = `si8ilev' + `si8tlev'
    }
    local si8demand = `si8cells'*(1 + ("`contreat'" != "")) + `si8levels' + 200
```
`si8ilev` is still computed unconditionally — the `> 100` granular-`i` caution at :1370 needs it,
and it still printed on every patched run (`run_p1.log:160`).

**Measured before/after**

| leg | base | p1 | p2 |
|---|---|---|---|
| 5,200-unit panel, `att`, default safety, maxvar 5000 | **r(907)**, "at least 5458" | **rc 0, 1.20 s** | **rc 0, 1.09 s** |
| same, `safety(off)` | rc 0 | rc 0 | rc 0 |
| same, `customvce(cluster unit)`, default safety | r(907), "at least 5458" | **r(907), "at least 5458"** | r(907), same |

and the estimates: **p1's and p2's default-safety run is BIT-IDENTICAL to base's `safety(off)`
run** (all 7 columns of `e(wdidmainresults)` in `%21x`). The refusal was purely the variable
count; the run it refused takes 1.2 seconds.

**Cost: NEGATIVE.** On the default path the patch *removes* one `gegen byte tag()` + `count`
pass. Priced on a 2,000,000-row panel (`run_cost.log`): `gegen tag(unit)` + `count` = **0.121 s**,
`gegen tag(year)` + `count` = **0.105 s**. So G1a saves ~0.1 s per 2M rows on the path
essentially every run takes, and costs exactly the pinned build's time on the reg/poisson path.

**Does it move a reported number?** **No.** `si8demand` is only ever compared with `c(maxvar)`;
it appears in no `e()` and no table. The new value is `<=` the old on every input, so the change
is monotone: it can convert a refusal into a run, never a run into a refusal. Proven by the
no-change control above (four specs, three builds, bit-identical) and by the `customvce` leg,
where the quoted demand string is character-for-character unchanged.

**Blast radius.** (i) The error *message* quotes `si8demand`, so a baselined log that pins the
number 5458 on an absorbing-path run would move. I cannot check the project baselines — the test
suite is off-limits to me — **flag this for a baseline grep on `requires at least`.**
(ii) The false-negative band widens by `si8ilev + si8tlev`: for cell counts in
`[c(maxvar) - si8ilev - si8tlev - 200, c(maxvar) - 200)` a run that used to be refused early now
reaches Stata's own late r(907). That is the direction the design already declares acceptable
(`ado:1342-1348`, `open_defects.md:74-79`), and G1b below closes it properly. I measured the band
is already leaky in the pinned build: at `maxvar 2048`, a 100-treated-cohort panel
(`si8cells = 1000`, pinned demand `1000+140+20+200 = 1360 < 2048`) gets **no** early refusal on
either build and dies on the raw r(907) (`run_bound.log`).

**Alternatives rejected.**
- *Drop `si8ilev + si8tlev` unconditionally.* Rejected: on the reg/poisson path those ARE columns
  (`:1485` writes `i.`i' i.`t'` out). Measured: the `customvce(cluster unit)` leg genuinely needs
  them and still refuses correctly after G1a.
- *Turn the refusal into a caution.* Rejected: it would remove a fast, correct refusal in the
  cases where it is right, and the user still loses the run.
- *Count `c(k)` (the dataset's own variables) into the demand.* It would be a legitimate
  tightening but it is orthogonal to this defect and adds a second reason for the quoted number
  to move. Left out deliberately.

**Test that would prove it correct (no oracle needed).**
1. *Cross-engine identity.* Same data, same `maxvar`: default safety must reproduce `safety(off)`
   bit-for-bit whenever `safety(off)` succeeds. Ran it; identical.
2. *Monotone-refusal invariant.* For a grid of panels, `refused(patched) => refused(pinned)`, and
   `succeeded(pinned) => succeeded(patched)`. Ran the two endpoints.
3. *Engine split.* The identical spec with `customvce()` must keep the identical refusal and the
   identical quoted demand. Ran it; identical string.

---

### G1b — name the late r(907) in the shell epilogue  (RECOMMENDED COMPANION, tested)

The preflight can only ever fail *fast when certain*. Everything else — and after G1a that
includes the band above, but also cases the pinned build already misses — surfaces as Stata's
bare `maxvar too small`, which says nothing about what in *this* command consumed the variables.

**File / site:** `wooldid.ado:100`, the thin shell, immediately before
`  if `wdsh_rc' != 0 exit `wdsh_rc'`.

**Insert:**
```stata
  * MAXVAR DIAGNOSIS (companion to the :1355 preflight): the preflight fails fast only
  * when the demand is CERTAIN, so every other route to Stata's variable limit - including
  * the sparse staggered grids the cell count under-predicts - arrives here as a bare
  * r(907) whose stock message says nothing about what in THIS command consumed the
  * variables. rc-gated: a successful run never reaches it, so it costs nothing.
  if `wdsh_rc' == 907 dis as error "Note: r(907) is Stata's variable limit; maxvar is currently `=c(maxvar)'. wooldid builds one interaction term per treated cohort-by-period cell and Stata allocates design-matrix columns for all of them, so a granular cohort variable is the usual cause. Either raise the limit with -set maxvar- (fixed at 2048 in Stata/BE and not raisable there), or specify unit-level fixed effects via fe() and swap the cohort variable to a coarser level of aggregation, such as cohorts or treatment-onset groups."
```

**Measured** (`run_p3.log`): a 100-cohort panel at `maxvar 2048` that dies late now prints the note
above Stata's own message and still exits 907; a 60-cohort panel on the same build succeeds with
**no output change at all**.

**Cost: zero.** The line is inside the existing `if `wdsh_rc' != 0` epilogue region and is gated on
`== 907`; a successful run never evaluates it. **Moves no number.** Blast radius: one extra line of
error text on r(907) failures, which is the entire point. Alternative rejected: wrapping the
`reghdfe` call at :1476 in `cap` and re-raising — that swallows reghdfe's own messages on every
other failure code and has a far wider blast radius for no extra benefit.

---

### G1c — tighten the bound to the measured allocation  (AUTHOR'S CALL, NOT in p1/p2)

While pricing G1a I measured what Stata actually allocates on the absorbing path. With
`set emptycells drop` (wooldid sets it at `:602`), at `maxvar 2048`, all cohorts sharing one onset
so `si8cells = 10K` exactly (`run_bound.log`, `run_bound2.log`):

| si8cells | 600 | 800 | 900 | **950** | 980 | 1000 | 1500 | 1800 | 2000 |
|---|---|---|---|---|---|---|---|---|---|
| outcome | rc 0 | rc 0 | rc 0 | **r(907)** | r(907) | r(907) | r(907) | r(907) | early refusal |

The break is between 900 and 950 at `maxvar 2048`, i.e. **the true demand is ≈ `2*si8cells + ~150`**,
not `1*si8cells`. `2*si8cells + 200` would have refused the 950…1800 cases *early and correctly*
(the pinned build misses all of them too) while still letting the G1 headline case through
(`2*50 + 200 = 300 << 5000`). So it is strictly better in both directions — on the binary reghdfe
path I measured.

**I am not proposing to ship this yet**, and I say so honestly: I measured the multiplier on ONE
configuration (binary, `reghdfe`, no `contreat`, no `lattice()`, no `fe()`, no controls). If any
engine allocates fewer than 2 columns per cell, a `2*` factor re-creates exactly the
false-positive class this whole proposal exists to remove. The measurement that would settle it is
cheap and mechanical: for each of {reghdfe, ppmlhdfe, reg, poisson} x {binary, contreat,
contreat+lattice}, bisect `si8cells` against a fixed `maxvar` and record the break point. If the
minimum observed multiplier is `m >= 2`, ship `floor(m)*si8cells`.

---

# G2 — `contreatwithin` + `contreat()` silently disables two input-domain screens
*covers:* **hv-ctw-ipos** (and its named sibling, the BS-34 float caution, which is the identical
one-word bug on the next block down)

`:4828` preserves the user's value as `safetyuser`; `:4829` rewrites the working `safety` to `off`
under `contreatwithin` + `contreat()`; `:4934` re-publishes the preserved copy as `safetyinput`.
Two checks below that line still gate on the rewritten working value.

**File / site 1:** `wooldid.ado:4959` — `  if "`safety'" != "off" {`
**becomes** — `  if "`safetyinput'" != "off" {`

**File / site 2:** `wooldid.ado:4973` (the float/2^24 caution) — identical one-word change.

**Rationale, in the codebase's own words.** `:1345-1350` already argues exactly this for the
neighbouring preflight: *"The relaxation contreatwithin needs is of the estimability screens, not
of a variable count."* An identifier-domain assert and a storage-type caution are no more
estimability screens than a variable count is. And `sthlp:1197-1199` already publishes the claim
that *"every data-integrity assert"* gates on the user's own value — so this edit makes the code
match a contract the help already states.

**Measured before/after** (one treated cohort's `id` made bad, everything else identical;
`contreatwithin` + `contreatcontrols(byCohort)` + `contreat(dose)`, **no `safety()` typed**):

| bad shape | base | p1 / p2 |
|---|---|---|
| `id = 1.5` | `r(452)` `__000005: factor variables may not contain noninteger values` | **`r(411)` "Error: id must consist of positive, non-zero integers, but does not"** |
| `id = -1` | `r(452)` `... may not contain negative values` | **`r(411)`, same named message** |
| `id = 0` | **rc 0, no message at all** | **`r(411)`, same named message** |
| good ids (control) | rc 0 | **rc 0, 112 hex lines BIT-IDENTICAL to base** |

**Cost: zero** on any run that already passes the screen, and **negative** on any run that fails
it (the assert now short-circuits before estimation instead of after). The assert itself is one
`assert` pass that default-safety runs already pay; the change only decides *whether* the branch
is entered.

**Does it move a reported number?** **No** for every currently-valid run — proven by the
`contreatwithin` control legs (`P3_good`, 112 lines, and `P4C`, 96 lines) being bit-identical
across base/p1/p2. It DOES change the outcome of runs that currently complete with an out-of-domain
identifier, from "a wrong or unaudited answer" to a named refusal — see G3, which is that case.

**Blast radius.** A user who today runs `contreatwithin` + `contreat()` with a non-positive or
non-integer cohort id gets an r(411) where they used to get an r(452) (no behavioural loss: both
abort) or, for `id == 0`, where they used to get silently wrong pretreatment rows. There is no
path where this turns a good run bad.

**Alternatives rejected.** *Move the check above `:4829`.* It was there before `e26b3dd` and the
CS-98 ruling deliberately hoisted the rewrite to `:4829` so that messages keyed on `safety` stop
contradicting themselves; re-splitting the block would reopen that. *Leave it and document it.*
Rejected: `sthlp:1195`'s own enumeration does not list this site, `sthlp:1209`'s does not either,
so the documentation would have to gain a site rather than lose a bug.

**Test.** Nine legs, no oracle: {no `contreatwithin`, `contreatwithin`, explicit `safety(off)`} x
{`id` non-integer, negative, zero}. The first two columns must be identical after the fix; the
third (explicit `safety(off)`) must stay as documented. Plus the bit-identity control above.

---

# G3 — a treated cohort whose identifier is literally `0` aliases wooldid's own sentinel
*covers:* **sk-izero-preweight**

### Confirmed on the pinned build first

`run_base.log`, 20 cohorts x 3 units x 8 periods, a pre-trend planted in cohort 1 only. Leg A
labels cohort 1 as `0`; leg B adds 1 to every identifier — **the same panel and the same
partition**, so every number must agree. Binary + `safety(off)`, `%21x` on all 7 columns:

```
BASE  legA vs legB : DIFFER on 6 of 72 hex lines - and they are exactly
      pretreatment:att_i      A +1.4b27a4e155ec0X-003   vs  B +1.9545f8e903fa0X-003
      pretreatment:att_itime  (same pair)
      ... with se / t / p / lb / ub all moving too.
```
Under `contreatwithin` + `contreat()` with **no `safety()` typed**, the same panel differs on 12
lines: `pretreatment:att_i`, `pretreatment:att_itime`, `pretreatment:ame_i`,
`pretreatment:ame_itime`. Every other row — `main:*`, `pretreatment:att`, `pretreatment:att_it`,
every `ES_*` — is bit-identical, which is exactly what the mechanism predicts.

**Mechanism.** `:928-929` (and its twin at `:12178-12179`) set `_i = i` on estimable treated rows
and `0` everywhere else, so `0` is wooldid's internal *"not an estimable treated cell"* code —
and `:1009` defines `targetflag = (_t > 0 & _t < .)`, i.e. the same sentinel. A cohort labelled
`0` therefore has `_i == 0` on every row; `:1602`'s classifier
`replace posttreatment = 2 if t<ttre & ttre!=. & _i!=0 & _t!=0` never promotes its pre-window rows
out of category 0; `:1620`'s `gegen w_per_i = total(w*esamp), by(posttreatment i)` then pools its
pre-window rows with its reference rows into one denominator, so `w_att_i = w/w_per_i` sums to
less than one over the pre-window and the cohort is silently downweighted. Post rows are
`posttreatment == 1` for everybody, which is why `main:*` is clean.

### G3a — refuse a cohort identifier of 0, unconditionally  (RECOMMENDED, tested in p1b)

**File / site:** `wooldid.ado`, immediately after the existing identifier block closes at `:4965`.

**Insert (a new, separate `if` — not an `else`, because Stata will not accept comments between
`}` and `else`):**
```stata

  * THE VALUE 0 IS NOT A SAFETY SCREEN. The positive-integer test above is one, and
  * safety(off) bypasses it. 0 is different: it is the code wooldid writes into its own
  * internal cohort variable (_i at :928 and :12178) to mark a row that is NOT part of
  * an estimable treated cohort-period cell. A cohort labelled 0 therefore ALIASES that
  * sentinel: the classifier at :1602 never promotes its pre-window rows, :1620 pools
  * them with the reference rows in one w_per_i denominator, and the cohort is silently
  * downweighted in pretreatment att_i / att_itime (measured 2.2x too small, with no
  * caution). There is no correct number to report, so this refuses under safety(off) too.
  * BLANKET, not restricted to units that are treated IN THE OBSERVED DATA: permute() and
  * bootstrap() re-enter this validation on a REASSIGNED onset, so a treated-only form
  * passes rep 0 and then fails inside the replication loop on the reps that hand unit 0
  * an onset - measured, 6 of 20 reps lost and the run dead at r(459). One assert pass, on
  * the branch that tested nothing, and it is the same domain default safety enforces.
  if "`safetyinput'" == "off" {
    cap assert `i' != 0 `ifassert'
    if _rc != 0 {
      dis as error "Error: `i' takes the value 0, which wooldid reserves internally to mark observations outside any estimable treated cohort-period cell. safety(off) does not bypass this: with a cohort labelled 0 the pretreatment att_i and att_itime estimates are silently wrong, and under permute() or bootstrap() a reassigned onset can make any unit that cohort. Recode the identifier so no value is 0, e.g. replace `i' = `i' + 1, or egen newi = group(`i')."
      error 411
    }
  }
```
(`i` and `ifassert` are in scope — `:4875` and `:4955-4956`.)

**Measured** (`run_p1b.log`, `run_ctrl0b.log`):

| leg | base | p1b |
|---|---|---|
| treated cohort labelled 0, binary `safety(off)` | rc 0, pretreatment att_i **2.2x too small** | **r(411), named message** |
| treated cohort labelled 0, `contreatwithin` (no `safety()` typed) | rc 0, no message at all | **r(411)** (via G2's screen) |
| leg B (ids +1) and all five ordinary specs | — | **rendered byte-identical to base**, 1,382 lines |

**Why BLANKET and not "only units treated in the observed data" — this is measured, not argued.**
My first cut asserted `i != 0 | ttre >= .`, which lets a never-treated unit numbered 0 through.
That form is not invariant under the replication engines, which re-enter `wooldid_validate` on a
*reassigned* onset (`ado:102-105`):

| leg (40 units, 20 treated, one never-treated unit labelled 0, `safety(off)`) | base | p1 (treated-only) | p1b (blanket) |
|---|---|---|---|
| plain run, leg A (control id 0) vs leg B (ids +1) | **rc 0 both, 64/64 hex lines IDENTICAL** | rc 0 | **r(411) at parse** |
| `cluster(id) permute(shuffle) permreps(20)` | **rc 0, 20/20 reps** | **6 of 20 reps lost -> `r(459)` "only 14 of 20 permutation replications estimated successfully"** | **r(411) at parse** |

So the treated-only form buys a run that works and pays for it with an opaque mid-loop failure
whose message blames the data. The blanket form refuses at parse with a one-line remedy, is
invariant under permutation and bootstrap reassignment, is one term simpler, and enforces exactly
the identifier domain default safety already enforces (`i > 0` excludes 0 for every unit).

**Cost.** ONE `assert` pass over the rows the `if` qualifier keeps, **only on the `safety(off)`
branch**, which today performs zero passes there. No sort, no `gegen`, no `by`. Order of magnitude:
the `gegen tag()` + `count` pair measured above is 0.12 s per 2,000,000 rows and a bare `assert`
is cheaper than that. Default-safety runs — essentially all runs — pay **nothing**: the branch is
not entered, and the existing `assert i > 0 & int(i) == i` already covers the same ground.

**Does it move a reported number?** **Yes, and the author must rule on it.** Two classes change:
(a) a treated cohort labelled 0 under `safety(off)` — a currently-printed `pretreatment:att_i` /
`att_itime` (and under `contreat`, `ame_i` / `ame_itime`) row disappears in favour of an error.
Those numbers are provably wrong: the panel is balanced, so `pretreatment:att_i` must equal
`pretreatment:att`; leg B satisfies that to the last hex digit and leg A does not — 0.01993
reported where 0.04404 is right. (b) a NEVER-treated unit labelled 0 under `safety(off)` — a
currently-CORRECT run (measured: 64/64 hex lines identical to the +1 relabelling) becomes an
r(411). That is the real price of the blanket form, stated plainly; the remedy is in the message
and is one line of the user's code.

**Blast radius.** Confined to `safety(off)` (and, before G2, to `contreatwithin` + `contreat()`;
after G2 that route refuses at the positive-integer screen instead). Default-safety runs are
untouched by construction — they already refuse 0.

**Why refuse rather than repair.** `_i == 0` is read as *"not a treated cell"* in at least six
other places: `:1361` (the si8 cell count, which undercounts), and the four `wooldid_hist` margins
subpops `:12360 / :12374 / :12405 / :12412`, which under `over(`_i')` would silently drop the
`i == 0` cohort from the cohort-effect histogram. G3b below repairs the *table*; it does not repair
the histogram (fixing those four needs `_t` plumbed into `wooldid_hist`, which today receives only
the user's `t` — `:1944`), and I did not audit the `wooldid_gtx` Mata name parse or
`wooldid_ccscreen` at that depth. Promising "safety(off) proceeds *correctly* with a cohort
labelled 0" is a promise I cannot underwrite; refusing is.

### G3b — make the classifier stop reading the sentinel  (tested in p2; free hardening)

**File / site:** `wooldid.ado:1602`.

**Now:**
```stata
    qui replace `posttreatment' = 2 if `t' < `ttre' & `ttre' != . & `_i' != 0 & `_t' != 0 & `esamp' == 1
```
**Becomes (a deletion):**
```stata
    qui replace `posttreatment' = 2 if `t' < `ttre' & `ttre' != . & `_t' != 0 & `esamp' == 1
```

**Why the deleted term is redundant for every well-formed identifier.** Within `esamp == 1` the
row is inside `ifstatement`, and `:923` (fixed-base regime) / `:12166` (the `wooldid_reltime`
twin) have already sent every unflagged or reference `_t` to `0`, so `_t` is never missing there;
hence `_t != 0` implies `_t > 0`, which by `:928` / `:12178` implies `_i == i`. For any `i != 0`
the two conditions coincide exactly; they differ only at `i == 0`, where the existing term is
wrong. The measurement below covers BOTH regimes: the no-change control runs use the default
(non-fixed-base) reference, i.e. the `wooldid_reltime` path.

**Measured — a metamorphic invariant, no oracle** (`run_p2.log`):

| comparison | base | p2 |
|---|---|---|
| leg A vs leg B, binary `safety(off)` | **DIFFER, 6/72 hex lines** | **IDENTICAL, 72/72** |
| p2 leg A vs base leg B | — | **IDENTICAL, 72/72** |
| base leg B vs p2 leg B | — | **IDENTICAL, 72/72** |
| ordinary panels (`P4A/B/C/D`, `P3_good`) base vs p2 | — | **IDENTICAL, 80/48/96/48/112** |

p2's leg A lands on leg B's numbers, and leg B is the leg the balanced-panel identity
(`pretreatment:att_i == pretreatment:att`) certifies. So the edit both restores relabelling
invariance and lands on the demonstrably correct side.

**Cost: zero** — one fewer term in an existing `replace`. **Moves no reported number** on any panel
without a treated cohort labelled 0, bit-identical over 336 hex lines across four specifications.

**Recommendation: ship G3a AND G3b.** G3a is the contract (there is no promise that the rest of the
`_i == 0` surface is clean); G3b is a free deletion that makes the internal invariant true anyway,
so the next relaxation of the screen cannot resurrect the defect. If the author prefers to keep
`safety(off)` meaning "proceeds", ship G3b alone and add a help sentence saying the histogram omits
a cohort labelled 0 — but I would not call that closed.

**Alternatives rejected.**
- *Recode `_i` internally (offset it away from 0).* Rejected: `_i`'s VALUES are user-visible —
  `:1766` `mkmat `_i' ...` puts them in `e(gteffects)`, and they appear in coefficient names that
  `:1671`, `:11855`, `wooldid_gtx` and `wooldid_sglbl` parse. An offset would relabel output.
- *Classify on the user's `i` plus a treated flag.* Equivalent to G3b but adds a term instead of
  removing one; `_t != 0` already IS the treated-and-flagged flag (`:1009` defines `targetflag`
  from precisely that expression).
- *A treated-only assert in G3a.* Rejected on measurement: not invariant under `permute()` /
  `bootstrap()` onset reassignment — see the table in G3a.

---

# G4 — the `relyr` column rounds fractional horizons onto the same integer
*covers:* **SKEP-02**

**File / site:** `wooldid.ado:10747`, inside `wooldid_mlfmt` (the PS-3 per-column format builder).

**Now:**
```stata
    if "`wdmc1'" == "relyr" | substr("`wdmc1'", 1, 6) == "numpre" local wdmf `wdmf' %7.0f &
```
**Becomes:**
```stata
    if "`wdmc1'" == "relyr" | substr("`wdmc1'", 1, 6) == "numpre" local wdmf `wdmf' %7.0g &
```

**Measured** (`run_base.log`, PART 0 — the same `%7.0f` / `%7.0g` pair on every value):

| value | `%7.0f` | `%7.0g` |
|---|---|---|
| -3 -2 -1 0 1 2 3 10 -10 | `-3 -2 -1 0 1 2 3 10 -10` | **byte-identical, same 7-wide field** |
| -1.5 | `-2` | `-1.5` |
| -0.5 | `-0` | `-.5` |
| 0.5 | `0` | `.5` |
| 1.5 | `2` | `1.5` |
| 2.5 | `2` | `2.5` |

`%7.0f` is round-half-to-even, which is why -0.5 and +0.5 both collapse and one prints `-0`.

**Cost: zero** (a format character). **Moves no reported number, and moves no printed character on
any integer-onset run** — proven twice: by the identity row above (every integer renders
identically in the same 7-wide field), and, more directly, by the rendered-table diff at the top of
this document. Those 1,382 byte-identical lines were produced by builds that HAVE this patch, and
they include the `relyr` column on `main:*` and `pretreatment:*` rows, where the cell is Stata's
system missing — so `%7.0g` renders `.` exactly as `%7.0f` did as well. `e(wdidmainresults)` is untouched
(it already stores -1.5 exactly). PS-3's stated intent — *"the relative-year column prints as the
integer it is"* — is preserved: `%0g` prints an integer with no decimal point.

**Blast radius.** The `numpre*` columns share this branch; they are counts, so `%7.0g` renders them
identically. The only visible change anywhere is on a fractional-onset run, which is reachable only
behind `safety(off)` (the default refuses it r(459)).

**Honest verdict on severity.** Both reviewers scored this **low** and they are right: the row
stubs (`pre1.5` / `pre.5` / `post.5`) already carry the correct horizon and the matrix stores the
exact value. *"Do nothing"* is defensible. I propose it only because it is one character, costs
nothing, and is provably a no-op everywhere else.

**Not mine.** The two "CONFIRMED (sev high)" verdicts attached to SKEP-02 in the GUARD brief —
continuous treatment + `subgroups()` + event study + bootstrap/`permanalyticse` aborting on an
internal *"This is unexpected - please report it"* assertion — describe the `wooldid_sglbl`
assertion at `ado:12632`. That finding is filed under mechanism **SG** (`fixprop/mech/SG.md`) and
is not the `relyr` format issue. No proposal here; do not let it fall between the two briefs.

---

# Help-file edits these imply

1. **`sthlp:1207`** currently lists *"the cohort identifier not being a positive integer"* under
   **"refuses, and `safety(off)` proceeds"**. With G3a, carve out the zero case: keep the
   positive-integer clause where it is and add *"a cohort identifier taking the value 0"* to the
   **"What `safety(off)` does not buy"** list at `sthlp:1214` (the ~forty-five unconditional
   refusals), with one clause of why: 0 is wooldid's own internal code for a row outside any
   estimable treated cell, and `permute()` / `bootstrap()` can make any unit a treated cohort.
2. **`sthlp:1197-1199`** already claims *"every data-integrity assert … gate[s] on your own
   value"* but its parenthetical enumeration omits the identifier-domain assert and the float
   caution. After G2 the claim is true; add *"the cohort-identifier domain"* and *"the float-storage
   caution"* to the parenthetical so the published inventory names them (the finding's own point:
   `sthlp:1195`, `1208` and `1209` between them assign this site to neither family).
3. No help change is needed for G1 or G4.

---

# What I did NOT test, stated plainly

- **The project's own baselines and battery.** Off-limits to me. Two things there need a check I
  could not run: any baseline that pins the literal string `requires at least 5458` (or any other
  `si8demand` value) on an absorbing-path run, and any baseline table carrying a `relyr` column
  (G4 is a no-op on integers, but the grep is one command).
- **G1c's multiplier on engines other than binary `reghdfe`.** Measured on one configuration only;
  that is why it is offered as a question, not a patch.
- **The other `_i == 0` sentinel sites** — `wooldid_hist`'s four margins subpops, `:1361`'s cell
  count, `wooldid_gtx`'s Mata name parse, `wooldid_ccscreen`. G3a is chosen precisely because it
  does not require them to be clean; G3b does not fix them.
- **`poisson` / `ppmlhdfe` paths** under any of these edits. G1a's engine predicate is the same one
  the code uses at `:1457`, so `ppmlhdfe` follows `reghdfe` (both absorb) and `poisson` under
  `unconditionalse`/`customvce()` follows `reg` (both write the levels out) — but I ran neither.
- **`permute()` / `bootstrap()` replication paths — partially tested.** I ran `permute(shuffle)`
  with 20 reps on three builds while choosing G3a's form (`run_ctrl0*.log`), which is what exposed
  the treated-only assert's mid-loop failure. I did NOT run a replication through G1a or G1b: the
  si8 block is already skipped inside a replication (`"`wdrepmode'" != "1"`) and G1b sits in the
  shell above the loop, so I expect no interaction, but that is reasoning, not measurement.
- **`bootstrap()` at all.** Not run. G3a's blanket form is invariant under the `reid` remap by the
  same argument that makes it invariant under `permute()`, but I measured only `permute()`.
