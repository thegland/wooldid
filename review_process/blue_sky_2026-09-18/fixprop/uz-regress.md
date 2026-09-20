# uz-regress — fix proposals for CRIT-2/CRIT-3, written from the REGRESSION-RISK side

Agent angle: *find what a fix would break.* Everything below is measured on my own copies of
the ado. **The repo was not touched.** Workdir `.../scratchpad/bluesky/probes/uz-regress`
(NOTES.md, probe01..probe05, all `.log`s, all `out_*.csv`, four builds `b_stock` / `b_fixA` /
`b_fixB` / `b_fixC`).

---

## 0. BOTTOM LINE FOR THE SYNTHESIS

1. **The comment at ado:1230-1233 does not apply under `untreatedZero`, and I proved it.** At
   `contreatpoly(1)` the untreatedZero spelling emits **no pooled `c.contreat` column at all**,
   so the "base level gets partialled out" manipulation has nothing to do. Removing the ringer
   head there is not a trade-off; it is a pure deletion of a column with no job.
2. **The fix is a zero-regression change outside untreatedZero.** 75 non-untreatedZero
   configurations (9 control sets x 7 engines + 12 poisson arms) come back **byte-identical**
   in %21x hex, and the constructed estimation-command macros are byte-identical too, for all
   three patch variants.
3. **The fitted model does not change — only the basis Stata reports it in.** On noisy data the
   head-in and head-out regressions agree in `rss` to ~1e-15 relative and report the **same
   `df_m`, `df_r`, `rank`**; across the whole untreatedZero matrix `e(N)` and `e(df_r)` are
   byte-identical stock vs patched. The head adds one column and costs one column. So the fix
   cannot move a residual df, an SE scale factor or a sample size — only what `margins` can see.
4. **A minimal fix is NOT ENOUGH: there is a second, independent defect at ado:1314.** Under
   `contreatpoly(>=2)` the untreatedZero spelling `i.cell#c.x##c.x` emits a **pooled `c.x`**
   main effect. reghdfe drops it; plain `reg`/`poisson` keep it and drop a cell slope instead.
   After the ringer fix, `untreatedZero contreatpoly(2)` **still** reports
   `att = 0.351286051773` against a truth of `2.463539274929` on every explicit-regressor
   engine. If the synthesis ships only the gate edit, this survives and will read as
   "the fix didn't work".
5. **Exactly one place where a fix moves a currently-working default-path number**, and it is
   the poly edit, not the ringer edit, and only on data that violates untreatedZero's own
   precondition. Named and measured in §6.

---

## 1. THE COMPLETE CONSUMER MAP (verified by grep + trace, nothing else reads these)

`grep -n 'contreatcontrolinserter'` returns 1215-1217 (init), 1300-1302 (assignment),
1320-1324 (the `ccnoabsorb` transfer), 1461/1462, 1468, 1475/1476, 1484/1485, and the comment
at 1896-1899. That is the whole surface.

| macro | who consumes it | when |
|---|---|---|
| `contreatcontrolinserterABS` | `absorb()` of ppmlhdfe **first** call (1462) and of reghdfe (1476) | default engine only |
| `contreatcontrolinserterFE` | explicit regressor in the ppmlhdfe **rc-3200 fallback** (1468) and in plain `reg`/`poisson` (1485) | `unconditionalse`, `customvce()`, and any poisson run whose first ppmlhdfe call fails |
| `contreatcontrolinserter` | explicit regressor in **all four** calls (1462/1468/1476/1485); non-empty only under `ccnoabsorb`, which moves FE into it at 1320 | |
| the `ringer` variable itself | `wooldid_fcol`'s cell collapse, via `fccont contreat(\`contreat' \`ringer')` at ado:1899 | **unconditional** — the ado's own comment says "averaging it costs nothing where it is unused". Confirmed harmless: the poisson arms, which keep the collapse on, recover the truth exactly after the fix. |

Nothing keys a *flag*, a *screen*, a *caution* or an `e()` posting on FE/ABS being non-empty.
The only state they touch is the estimation command. This is why the fix is cheap.

**Bootstrap / permute exposure.** `permute()` refuses `contreatcontrols(untreatedZero)` outright
(ado:5842). `bootstrap()` refuses `unconditionalse` (6410) and `customvce()` (6414), and
`wdb_ro` (ado:321-...) never emits `customvce`. So the replication loops reach the defect
**only through `poisson`**, and there the observed run and every replication are broken the same
way — measured: `bootstrap(standard) reps(15) poisson` + untreatedZero returns `ame = 0.000`
against a truth of `0.25`, repaired to `0.250000000000` by the fix.

---

## 2. THE SINGLE MOST IMPORTANT QUESTION, ANSWERED (probe01_struct.log)

The comment says `i.null##c.(ringer)` "ensures that the base level of contreat gets partialled
out. if not specified ... one of the contreat control coefficients gets dropped for colinearity
instead."

Hand-built regressions on the exact geometry, 400 obs, `d == 0` off the cells:

| # | model | what came back |
|---|---|---|
| Q1a | `reg y i.g#c.x` | `0b.g#c.x` **carries a real coefficient**. Stata's `i.a#c.z` emits a slope for **every** level, base included. There is no unpartialled base level in the `#` spelling. |
| Q2 | `reg y i.g#c.d i.g` (= untreatedZero, poly 1) | **No pooled `c.d` column exists anywhere.** `0b.g#co.d` is *omitted* because that column is identically zero. Cell slopes .69677/.70405/.70258/.70031 against a truth of .70. `df_m=8 df_r=391 rank=9` |
| Q3 | `reg y i.null##c.(r) i.g#c.d i.g` (= the stock non-reghdfe path) | `r` takes **0.70030524**, `4o.g#co.d` is dropped, every survivor becomes `truth - 0.70030524`. `df_m=8 df_r=391 rank=9`, `rss` equal to Q2's to ~1e-15. |

**Verdict: the comment's justification is VACUOUS under untreatedZero at poly 1.** The
manipulation exists to kill a margins-visible pooled `c.contreat` that `##`-spelled
`contreatinj` creates; `untreatedZero`'s `#`-spelled `contreatinj` (ado:1314) never creates one.
A naive fix therefore reintroduces nothing. **This is the answer the cell was asked for.**

**But at `contreatpoly(>=2)` the answer flips** (Q4/Q5). `contreatpolypart` is
`c.x##c.x`, and `#` binds tighter than `##`, so ado:1314 builds `(i.cell#c.x)##c.x` =
cell#x **+ pooled c.x** + cell#x#x.
- Q4, no head: the pooled `d` **survives** (0.70696) and `4o.g#co.d` is dropped instead.
  Harmless *for the ame* — `d` is the real contreat variable, so `margins` counts its
  derivative; the model is merely reparameterised.
- Q5, with head: `r` takes 0.70696 and **both** `o.d` (the margins-visible pooled dose) **and**
  `4o.g#co.d` are dropped. Two visible columns zeroed, their content parked in an invisible
  ringer. So poly>=2 is *worse* under the stock head — and, per §5, still wrong for `att`
  after the ringer fix alone.

---

## 3. WHICH PATHS RELY ON THE POOLED CONTROL FOR A LEGITIMATE REASON

- **`basic`** — yes, obviously, and DF-1 (ado:1263-1299) already established that on every
  explicit-regressor path (`ccnoabsorb`, `unconditionalse`, `customvce`, the rc-3200 fallback)
  the `i.null##` spelling is load-bearing for the **collinearity drop order**. My 75-config
  matrix leaves all of that untouched.
- **`unconditionalse` / `customvce()`** — they need the head only *because* those paths carry the
  terms as explicit regressors where a margins-visible pooled `c.contreat` would otherwise
  survive. That premise is created by the `##` in the ordinary `contreatinj` (ado:1117). Under
  untreatedZero the premise does not exist (Q2).
- **`poisson`** — no independent reason. `poisson` is in the gate for the same drop-order
  reason, plus the fact that a term can get kicked out of `absorb()` by the rc-3200 fallback.
  Two extra facts worth knowing:
  - Under untreatedZero the head is exactly what **causes** the fallback: `i.null` is
    single-level, so `absorb(i.i i.t i.null##c.(ringer))` fails. Stock printed the four-line
    "initial estimation of ppmlhdfe failed" caution on 2 of 5 poisson arms; after the fix,
    **0 of 5**. That caution disappearing is a visible, intended behaviour change.
  - The other three poisson arms were broken **silently** — `customvce`/`unconditionalse` go to
    plain `poisson` (never ppmlhdfe), and `ccnoabsorb` moves the head into
    `contreatcontrolinserter`, which the ppmlhdfe *first* call already carries as a regressor.
- **`ccnoabsorb` alone** — not in the gate, so with untreatedZero it has nothing to move and is
  **already correct today** (`att 2.463539274929`, `ame 0.9325`). Any fix must leave it alone;
  all three of mine do.

---

## 4. THE PROPOSALS

### FIX-UZ-A (one token, at the gate) — ado:1234
```
  if strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
```
becomes
```
  if (strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson") & "`contreatcontrols'" != "untreatedZero" {
```
Exact because `wooldid_validate` (ado:4797) guarantees `untreatedZero` is the sole word.

### FIX-UZ-B (three lines, inside the branch that already knows) — ado:1313-1315 — **PRIMARY**
```
    if "`contreatcontrols'" == "untreatedZero" {
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
      local contreatcontrolinserterFE
      local contreatcontrolinserterABS
    }
```
(`ccringerhead`/`ccringerheadabs` may be cleared too, defensively; they are not read again.)

**Why B is primary, from the regression-risk seat.** A and B are *numerically indistinguishable*
— `out_uz_fixA.csv` and `out_uz_fixB.csv` are byte-identical in hex on all 28 arms — so the
choice is maintenance, and B wins on three counts:
- It sits **inside** the `== "untreatedZero"` branch, so "is this an untreatedZero run?" keeps a
  single source of truth. A introduces a second, `!=`-spelled copy of that predicate 80 lines
  away that a future edit can desynchronise.
- It is robust to a future engine being added to the 1234 gate. A silently re-breaks then.
- It runs **before** the `ccnoabsorb` transfer at 1319-1325, so `ccnoabsorb` is covered for free.
- It puts the clearing where the comment at 1311-1312 already promises it
  ("just not having contreatcontrols of any kind"), making the code match its own comment.
Shipping **both** is defensible and costs nothing; if only one, ship B.

### FIX-UZ-C (the residual the ringer fix does not reach) — ado:1314 — **SHIP WITH B**
Replace the single `contreatinj` line so the polynomial is `#`-chained instead of `##`-chained,
which is what the help file already documents (sthlp:1006 defines
`f(C) = beta_it1*C^1 + ... + beta_itp*C^p` **per i-t cell**, with no pooled term):
```
      local uzdose `latticepart1hash'#c.`contreat'
      if `contreatpoly' > 1 {
        forval uzpk = 2/`contreatpoly' {
          local uzrep
          forval uzq = 1/`uzpk' {
            local uzrep `uzrep'#c.`contreat'
          }
          local uzdose `uzdose' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'`uzrep'
        }
      }
      local contreatinj   `uzdose' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
```
At `contreatpoly(1)` this is **character-for-character the current expansion** (verified in the
trace: `INJ={#c.dose i.TV1#i.TV2#i.TV3}` under both), so the edit is inert unless the user asked
for a polynomial. `margins` reads `c.x#c.x` for the derivative exactly as it reads `c.x##c.x`
(confirmed: fixC recovers `ame = 0.932500000000` and the poisson `0.250000000000` exactly).

### FIX-UZ-D (companion; the author's "very cheap" bar) — the violated-assumption screen
The help already says (sthlp:995) "If this option is specified when the always 0 condition is not
met, wooldid will return poorly identified estimates" — and today nothing checks it, and nothing
says it happened. One pass over the estimation sample, inside `wooldid_ccscreen` or beside the
other `safety`-gated screens:
```
qui count if `contreat' != 0 & `contreat' < . & (`ttre' == . | `t' < `ttre') & `ifstatement'
if r(N) > 0 dis as text "Caution: contreatcontrols(untreatedZero) asserts that `contreat' is zero outside treated cohort-periods, but `r(N)' estimation-sample observation(s) outside them have a nonzero `contreat'. No baseline dose control is fitted under untreatedZero, so the reported marginal effects absorb that variation. Use contreatcontrols(basic) (or byTime/byCohort/...) if the dose is defined outside treatment."
```
**Do NOT write this as `\`targetflag\' == 0`.** `targetflag` (ado:1007) is `_t > 0 & _t < .`
on the SHIFTED relative-period code (wooldid_reltime does `_t = _t + min_t + 1` before the
missing-to-0 fill), so under `esprelength() > 0` it is **1 on pre-treatment horizons of treated
units too** — precisely the rows untreatedZero requires to be zero. The correct predicate is the
onset test the worker already uses at ado:1602, `\`t\' < \`ttre\' & \`ttre\' != .`, widened to the
never-treated; `\`ttre\'` there is already the anticipation-adjusted onset. The version I timed is
this one (`!post` in the probe), not the `targetflag` one.

Caution, not refusal: the help documents the case as user-owned. Cost is one `count` — O(N),
one pass, no sort, no temporary variable; it found all 48 planted rows and timed at 0.000s.
Gate it on `safety != off` like its neighbours. Post `e(wdid_uzviol)` if a machine-readable
marker is wanted.

### FIX-UZ-E (documentation, mandatory with B)
The help's *model* statement is the strongest argument that the current behaviour is a bug, so
the help needs no new promise — it needs a changelog line. sthlp:1033-1034 already says
untreatedZero "specif[ies] g() = 0" and works "By not including a g() function at all". The
pooled `i.null##c.(ringer)` **is** a g(). State in the release notes that
`contreatcontrols(untreatedZero)` combined with `unconditionalse`, `customvce()` or `poisson`
previously fitted the `basic` baseline dose control and that every reported quantity moves.

---

## 5. WHAT THE FIX REPAIRS (probe03, probe05; data_uz.dta, exact-fit, truth known)

Truth: `att = 2.463539274929`, `ame = 0.932500000000`, poisson semi-elasticity `= 0.250000000000`.

| arm | stock att / ame | fixA = fixB | fixC |
|---|---|---|---|
| `robust` (default) | 2.463539274929 / .9325 | same | same |
| `cluster(unit)` | 2.463539274929 / .9325 | same | same |
| `ccnoabsorb robust` | 2.463539274929 / .9325 | same | same |
| `unconditionalse robust` | **0.351286051758 / -0.1075** | 2.4635 / .9325 | 2.4635 / .9325 |
| `customvce(robust)` | **0.351286051758 / -0.1075** | 2.4635 / .9325 | 2.4635 / .9325 |
| `customvce(cluster unit)` | **0.351286051758 / -0.1075** | 2.4635 / .9325 | 2.4635 / .9325 |
| `ccnoabsorb unconditionalse` | **0.351286051758 / -0.1075** | 2.4635 / .9325 | 2.4635 / .9325 |
| + `lattice(all 3)`, same 4 arms | **-0.308713948243 / -0.1075** | 2.4635 / .9325 | 2.4635 / .9325 |
| + `contreatpoly(2)`, same 4 arms | **0.351286051783 / -0.1075** | **0.351286051773 / .9325** | **2.4635 / .9325** |
| poisson x {robust, cluster, customvce, ccnoabsorb, uncond} | **0.10 / -0.000000000000** | 0.6078 / **0.250000000000** | same |
| ES arm, `esprelength(2) espostlength(3)` | every horizon moves (ES_ame_post0/1/2: .135/.005/-.125 vs .. | 1.175/1.045/.915 | same |

Two things worth putting in front of the author:
- **The poisson number is not just wrong, it is annihilated.** The design has a *constant* true
  semi-elasticity of 0.25 in every cell, so subtracting the dropped cell's slope leaves exactly
  `-0.000000000000` on the **default** engine. `att` survives (0.10, the true intercept shift),
  so the table looks healthy.
- **CRIT-3's "the broken number is exactly `basic`'s answer" is conditional.** It holds while
  untreatedZero's own precondition holds (then the ringer is exactly the sum of the cell dose
  columns). On data that violates it, the broken value is `basic`'s answer *minus the dropped
  cell slope* — measured 0.9325 - 1.04 = -0.1075 while `basic` itself returns 0.9325. Worth
  correcting in the write-up, because "it equals basic" is the sentence a reviewer will test.

**Curvature closure (probe06, a design with a REAL quadratic, `+0.30*dose^2` on treated rows,
noise sd 0.05; truth ame 2.183060605282).** The exact-fit designs above have a dose-linear
response, so every `c.x#c.x` coefficient there is 0 and they test structure, not differentiation.
On the curved design, fixC's poly-2 arms return the DEFAULT engine's answer:
`att 3.930967579219 / ame 2.195186690884` on `robust`, and the same to ~1e-13 relative on
`unconditionalse robust` and `customvce(robust)` (hex `+1.f729f21530b3c` vs `+1.f729f21530b55`
— the ordinary reghdfe-vs-explicit-dummies gap, not a spelling gap). So `margins` differentiates
the `#`-chained quadratic exactly as it did the `##`-chained one. Stock on the same data:
`0.329136978179 / -0.131262907883`. Poly-1 on the same data: stock non-default
`-0.705442913818 / -0.046152099304`, fixC `3.930967579219 / 2.178327293436` = the default.

**Subgroups closure (probe06).** With `sg = mod(i,2)` (constant within cohort, hence within the
(i,t) cells, so not refused): stock `unconditionalse` + `subgroups()` returns
`0.351286051758 / -0.1075`, fixC returns `2.463539274929 / 0.9325` = the default path. The
subgroup margins calls carry the defect and the fix repairs them.

Also repaired, measured: `contreatwithin` + `customvce(robust)` (-0.1075 -> 0.9325);
`semielasticity` + `unconditionalse` (-0.1075 -> 0.9325); `impvalues` + `poisson`
(0 -> 0.250000000000); `bootstrap(standard) reps(15)` + `poisson` (0 -> 0.250000000000).
`impvalues` on the default engine is **unchanged** (2.564349433284 / 0.98).

---

## 6. WHAT COULD REGRESS — the whole answer

### 6.1 Outside untreatedZero: nothing. Proven twice over.
`probe02_nonuz.do`, 75 configurations on `data_full.dta` (dose varies everywhere), all `rc=0`:
control sets `basic`, `byTime`, `basic byTime`, `byCohort`, `byTreated`, `byTtre`,
`basic`+`lattice(all 3)`, `byTime`+`lattice(all 3)`, `basic`+`contreatpoly(2)`; engines
`robust`, `cluster(unit)`, `unconditionalse robust`, `customvce(robust)`,
`customvce(cluster unit)`, `ccnoabsorb robust`, `ccnoabsorb unconditionalse robust`; plus 12
poisson arms.

- `out_nonuz_stock.csv` vs `..._fixA.csv` vs `..._fixB.csv` vs `..._fixC.csv`:
  **`diff` is empty**, in %21x hex on `att`, `ame`, `e(b)[1]`, `e(N)`, `e(df_r)`, `e(rss)`, `e(r2)`.
- `trace_nonuz_*.txt` (the constructed `FE`/`ABS`/`INS`/`INJ` macros, dumped from inside the ado
  by a trace block present in **all four** builds): **`diff` is empty**.
  Identical estimation command + identical data leaves no room for a numeric difference, which
  also sidesteps the known 1-2 ulp non-reproducibility of the clustered plain-`reg` paths.
- **The 75 configs confirm; the guard proves.** The edits touch four macros
  (`contreatcontrolinserterFE`/`ABS`, `ccringerhead`/`headabs`) and `contreatinj`, and every one
  of them is gated on the literal predicate `contreatcontrols == "untreatedZero"`. Whenever that
  is false, all five are built by byte-identical unchanged code, so the estimation command is
  identical for **every** option combination — including `fe()`, `controls()`,
  `interactivecontrols()`, `ccc_absorb()`, weights, `if`, `subgroups()` and the permute/bootstrap
  shells that I did not run. The matrix is corroboration of a structural argument, not the whole
  argument.
- This is the direct test of the DF-1 ruling's "ENGINE-CONDITIONED, and this is the part to
  keep": the `i.null##` spelling remains exactly as shipped on `ccnoabsorb`, `unconditionalse`,
  `customvce` and the rc-3200 fallback for every non-untreatedZero control set.

### 6.2 Inside untreatedZero: numbers move a lot, and that is the point
`att` moved -86% on my design; every event-study horizon moves; `ame` flips sign under poisson.
Anyone with a published `untreatedZero` + `unconditionalse`/`customvce()`/`poisson` spec has a
different number after this ships. They have a *wrong* number today — the default engine on the
same data already disagrees with them — so this is a correction, not a break. It still needs a
changelog line (FIX-UZ-E) because a user re-running an old do-file will see movement.

### 6.3 What CANNOT move, and why that matters
`e(N)` and `e(df_r)` are byte-identical stock vs patched on **every** untreatedZero arm,
including the poisson arms whose engine path changes (stock: ppmlhdfe first call fails r(3200)
-> fallback; patched: first call succeeds). Q2 vs Q3 on noisy data: same `df_m`, `df_r`, `rank`,
`rss` to 1e-15. The head adds one column and costs one column, so residual df, the SE scale
factor and the sample are untouched. Any reviewer who sees an SE move by a factor of
`sqrt(df_old/df_new)` after this fix has found a *different* bug.

### 6.4 The one currently-working number a fix moves — FIX-UZ-C only, and only on bad data
On `data_viol.dta` (untreatedZero **claimed but false**: 48 untreated rows carry a real dose):

| run | stock | fixC |
|---|---|---|
| `untreatedZero contreatpoly(2) robust` (**default engine**) | att **2.444765194869** | att **4.298562435389** |
| `untreatedZero robust` (poly 1, default engine) | att 4.298562435389 | 4.298562435389 (unmoved) |

Mechanism: the pooled `c.x` that `##` emits was acting as an accidental `basic`-style baseline
dose control, which is *right* for violated data and *not what untreatedZero means*. Removing it
makes poly-2 agree with poly-1 (both 4.2986) and with the documented `g() = 0`.

**And the number it moves is itself incoherent today.** On the same stock poly-2 default run the
two headline quantities come from two different implicit models: `att = 2.444765194869` is the
as-if-`basic` answer (the pooled `c.x` absorbing the baseline dose response, and it happens to
equal the planted `basic`-controlled truth 2.444765194869 exactly), while
`ame = 1.832500000000` is the as-if-`untreatedZero` answer (margins sums the pooled term back
into the derivative, giving nuisance 0.9 + treatment 0.9325). One table, two models. fixC makes
them consistent (4.298562435389 / 1.832500000000, both as-if-`untreatedZero`). Assessment:
**a real movement on a default path, but toward internal consistency, and it must be declared,
not discovered.**

Caveat on the mechanism: "reghdfe drops the pooled `c.x`, plain `reg` keeps it and drops a cell
slope instead" is **inferred** from the `att` outcomes plus the hand-reg column dump in probe04
(`SHIPPED col: dose = 1.04` under explicit dummies), not observed inside reghdfe. The
consequence either way is that default-path poly-2 correctness on VALID data currently rests on
a solver's collinearity drop choice — which is exactly what FIX-UZ-C removes. It is invisible on data that satisfies untreatedZero (probe03: fixC vs fixA
differ on no default-path arm). FIX-UZ-D's caution is what makes the affected user see it.
If the author would rather move nothing on the default engine, FIX-UZ-C can be narrowed to the
explicit-regressor paths only — but then `att` stays engine-dependent under poly>=2, which I do
not recommend.

### 6.5 Three regression hypotheses I tested and DISMISSED
- *"Clearing ABS changes the reghdfe `df_a` the way DF-1 did."* No: under untreatedZero the gate
  never fires on the default reghdfe path in the first place (`FE={} ABS={}` in the stock
  trace for `robust`, `cluster()` and `ccnoabsorb`), so there is nothing to clear there.
- *"The cell collapse keys on the ringer being in the model."* No: `fccont` at ado:1899 is
  unconditional by design, and the ado says so. The poisson arms keep the collapse on and
  recover the truth exactly after the fix.
- *"Some screen/caution/`e()` posting keys on FE/ABS being non-empty."* No: grep shows the only
  consumers are the four estimation calls and the `ccnoabsorb` transfer.

---

## 7. ALTERNATIVES CONSIDERED AND REJECTED

- **Refuse `untreatedZero` + `unconditionalse`/`customvce()`/`poisson` at parse time.** Rejected:
  it converts a silently wrong number into a dead end for a combination that, once fixed, is
  perfectly well defined and already works on the default engine. It would also refuse
  `poisson` + `untreatedZero` outright, a combination nothing in the help forbids and which works correctly once the head is gone.
- **Keep the head but spell it so it is not collinear.** Impossible by construction: under
  untreatedZero the dose is zero off the treated cells, so *any* pooled ringer column is exactly
  the sum of the cell dose columns. There is no non-degenerate spelling. This is the strongest
  argument that omission is the only correct action.
- **Restrict the ringer to the untreated rows.** Under a true untreatedZero the dose is
  identically zero there, so the column is identically zero and is dropped — same as omitting
  it, with extra machinery. Under a false untreatedZero it would silently become a `basic`
  control, which is the behaviour the option exists to refuse.
- **Fix it by editing `contreatpolypart` (ado:1114-1119) instead of ado:1314.** Rejected:
  `contreatpolypart` is shared with the ordinary `contreatinj` at 1117, where the `##` is
  correct and load-bearing. The untreatedZero branch is the only safe site.

---

## 8. HOW I WOULD PROVE A CANDIDATE FIX CORRECT (the test a reviewer should demand)

1. **Macro-identity leg.** Re-run the 75-config non-untreatedZero matrix with the trace block on
   stock and on the candidate; `diff` the trace and the hex CSV. Both must be empty. This is
   stronger than a numeric tolerance and immune to the clustered-`reg` ulp noise.
2. **Convergence leg.** Exact-fit continuous-treatment panel, `untreatedZero`, every engine
   (`robust`, `cluster()`, `ccnoabsorb`, `unconditionalse`, `customvce(robust)`,
   `customvce(cluster ...)`, `ccnoabsorb unconditionalse`) x `{poly 1, poly 2, lattice(all 3)}`
   x `{OLS, poisson}`: every arm must return the **default engine's** answer, and that answer
   must equal the planted truth. 26 arms; mine all pass under fixC. **Add a curved arm** — a
   design whose dose response is genuinely quadratic — or the poly legs test structure only;
   mine is probe06 (fixC poly-2 non-default agrees with poly-2 default to 1e-13 on
   `att 3.930967579219 / ame 2.195186690884`). **Add `contreatpoly(3)`.**
3. **df leg.** `e(N)` and `e(df_r)` byte-identical stock vs patched on every untreatedZero arm.
4. **Fallback leg.** Count "initial estimation of ppmlhdfe failed" in the log: stock 2, patched 0.
5. **Violation leg.** Re-run on data that violates the precondition and record the movement
   named in §6.4, so it is a declared consequence rather than a surprise.

---

## 9. UNTESTED / RESIDUAL RISK (honest list)

- `gteffects` + `unconditionalse` is refused r(198), so that pairing is untested; `gteffects` on
  the default engine under untreatedZero was not run.
- `contreatpoly(3)` and higher: FIX-UZ-C's loop is written for arbitrary `p` and is inert at
  `p = 1`, but only `p = 2` was run. A `p = 3` arm should be added to leg 2.
- `subgroups()` x untreatedZero is measured only for ONE subgroup variable on ONE design
  (probe06); the per-subgroup contrast rows and the across-subgroup joint test were not read.
- FIX-UZ-C was measured on poly 2 with and without curvature, but the `lattice()` x poly>=2
  combination under untreatedZero was not run (the loop carries `latticepart1hash` into every
  term, so it should be covered, but "should be" is not a measurement).
- FIX-UZ-D was TIMED but not IMPLEMENTED in a build; the predicate above is proposed, not run
  inside wooldid.
- Weighted runs (`aw`/`fw`/`pw`) were not exercised; the mechanism is a column-space identity
  and weights do not change collinearity, but I did not measure it.
- `honestdid`/`pretrendpower`/`adversarial` downstream of a repaired untreatedZero fit.
- Everything here is one panel shape (4 cohort groups, 10 units each, 8 periods, 320 obs) plus
  three variants of it. The mechanism is exact and engine-level, and the numbers reproduce the
  refuter's independent design, so I do not think shape matters — but I did not vary it.
- FIX-UZ-C was built on top of FIX-UZ-A, not FIX-UZ-B. They are numerically identical, so the
  composition is safe, but the exact patch the author ships should be re-run through leg 1.
