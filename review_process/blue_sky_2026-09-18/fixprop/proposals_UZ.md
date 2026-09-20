# Fix proposals — mechanism `UZ`: `contreatcontrols(untreatedZero)`

All four findings in the brief (`fp-contreat-E`, `sk-untreatedzero-regresspath`,
`FP-CTZERO-ENGINE`, `SK-UZ-POISSON`) are **one defect at one site** and are closed by **UZ-P1**,
a one-line change to the gate at `ado:1234`. Testing that patch turned up a **second, independent
defect** in the same feature that P1 does not reach — `contreatpoly(2+)` under untreatedZero —
which is **UZ-P2**. UZ-P3 records what I think should *not* be changed.

Everything below is measured. Patched builds, probe do-files and logs:
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-UZ/`
(`NOTES.md` there is the index; `build_patched/` = P1, `build_patched2/` = P1+P2).

---

## UZ-P1 — untreatedZero must not take the pooled baseline dose control on any engine

**Closes:** `fp-contreat-E` (critical), `sk-untreatedzero-regresspath` (critical),
`FP-CTZERO-ENGINE` (high, with the refuter's inversion — the *explicit* arm is the wrong one),
`SK-UZ-POISSON` (high).

### The edit

`wooldidinstall/wooldid.ado`, **line 1234**, the gate that builds `ccringerhead` /
`ccringerheadabs`.

Now:
```stata
    if strpos("`contreatcontrols'",  "basic") != 0 | "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" {
```
Becomes:
```stata
    if strpos("`contreatcontrols'",  "basic") != 0 | ( ( "`unconditionalse'" == "unconditionalse" | "`customvce'" != "" | "`poisson'" == "poisson" ) & "`contreatcontrols'" != "untreatedZero" ) {
```

Nothing else changes. `tempvar null` / `gen null = 1` live inside the gate, so under untreatedZero
they are no longer created either; `ccringertail` is already empty under untreatedZero (validate at
`ado:4797` refuses pairing it with any byX), so `contreatcontrolinserterFE` and
`contreatcontrolinserterABS` both come out empty — which is exactly what `ado:1311`'s own comment
says untreatedZero means ("the only thing to do is not control for anything!"). The gate is
upstream of the `ccnoabsorb` block at `ado:1319`, so `contreatcontrolinserter` inherits the
emptiness too.

Recommended companion (comment only, no behaviour): append one sentence to the DF-1 block that
states the rule, e.g. after `ado:1233`'s "…engine-conditioned" paragraph —
`* untreatedZero is excluded from this trigger on purpose: under it the margins-visible block`
`* carries NO pooled c.contreat (ado:1314 spells the injection with # not ##), so the pooled`
`* ringer is not a numerical twin there - it is a real column, exactly collinear with the`
`* cell#dose set whenever the dose is 0 off treated cells, i.e. always when untreatedZero is`
`* true. Rank falls, a cell slope is dropped in its place, and margins never sees the pooled`
`* slope, so att and every ES_att row lose it.`

### Why this is the right site

The gate's own comment (`ado:1225-1233`) gives the reason it includes the three engine flags: on
paths that pass these terms as **explicit regressors**, spelling `i.null##c.(ringer)` forces the
*base dose level* to be partialled out rather than letting reg drop one of the treated
dose-slope coefficients instead. That manipulation presupposes there is a pooled `c.contreat`
column in the margins-visible block to partial out. Under untreatedZero there is not:
`ado:1314` rewrites `contreatinj` with a single `#`, so the model is cell dummies + cell#dose and
nothing pooled. The ringer is then not redundant — with the dose exactly 0 off treated
cohort-periods, the pooled ringer column **is the sum of the cell#dose columns**, so it is
perfectly collinear with them, rank does not rise, and the engine drops one cell slope to pay for
it. `margins` switches the cell terms on and off but never the ringer, so the reported `att` and
every post-treatment `ES_att` row silently lose the pooled slope, and the reported `ame` becomes
`mean(beta) − beta(dropped cell)`. On the absorbing ppmlhdfe path the ringer sits inside
`absorb()`, the pooled slope is partialled out first, and under a homogeneous effect the visible
cell slopes — hence the reported `ame` — go to **exactly zero**.

### Measured, before → after

Design: exact-fit linear panel, dose ≡ 0 on every untreated observation,
`y = a_i + g_t + kappa_it + beta_it*dose`, `e(r2) == 1` asserted on every arm.
TRUE ame = 0.932500000000, TRUE att = 2.404948264178.
(`probes/prop-UZ/body.do`, PART A; `run_before.log` vs `run_after.log`.)

| arm | ame before | ame after | att before | att after |
|---|---|---|---|---|
| default (reghdfe) | 0.932500000000 | *(unchanged)* | 2.404948264178 | *(unchanged)* |
| `robust` | 0.932500000000 | unchanged | 2.404948264178 | unchanged |
| `ccnoabsorb robust` | 0.932500000000 | unchanged | 2.404948264178 | unchanged |
| `customvce(robust)` | **−0.107500000000** | **0.932500000000** | **0.375760748988** | **2.404948264178** |
| `customvce(cluster unit)` | **−0.107500000000** | **0.932500000000** | **0.375760748988** | **2.404948264178** |
| `unconditionalse robust` | **−0.107500000000** | **0.932500000000** | **0.375760748988** | **2.404948264178** |
| `lattice(all 3) customvce(robust)` | −0.107500000000 | 0.932500000000 | **−0.284239251013** | **2.404948264178** |

After the fix the six linear arms agree to ~2 ulp (`+1.33d5583f9bfa2X+001` default vs
`+1.33d5583f9bfa0X+001` explicit; ame `…a405` vs `…a3f2`, rel. 2e-15) — solver-level agreement,
far inside the project's 1e-10 indifference threshold.

Event study (`esprelength(2) espostlength(2)`, PART D). Before, `customvce(robust)` returned
`ES_att` post0/post1 = 0.756760458677 / 0.557341816386 and `ES_ame` = 0.135 / 0.005 against the
default arm's 2.874090298082 / 2.497600827286 and 1.175 / 1.045. After, the two arms return the
**same** numbers, and 1.175 / 1.045 are the closed-form truths (mean of `beta_it` at relative
year 0 and 1). This confirms the refuter's upward scope correction: the defect is not confined to
the `ame` row.

Poisson, exact-fit log-linear twin (`mu = exp(1.0 + 0.10t + 0.40*d2)`, `y = mu`; TRUE xb-scale
ame = 0.400000000000, att = 0.605020196668), PART C:

| arm | before (ame / att) | after (ame / att) |
|---|---|---|
| `poisson` (default ppmlhdfe) | −0.000000000000 / −0.000000000000 | **0.400000000000 / 0.605020196668** |
| `poisson robust` | −0.000000000000 / −0.000000000000 | 0.400000000000 / 0.605020196668 |
| `poisson customvce(robust)` | −0.000000000001 / −0.000000000002 | 0.400000000000 / 0.605020196669 |
| `poisson ccnoabsorb robust` | −0.000000000000 / −0.000000000000 | 0.400000000000 / 0.605020196668 |

On the noisy-count twin of the same DGP: ame 0.104457495291 → **0.444649524272** (truth 0.40),
att 0.131857990503 → **0.646415611201** (truth 0.605). The "exactly zero under a homogeneous true
effect" prediction in `SK-UZ-POISSON` reproduces to 1e-32 on the noiseless design.

**Two side effects, both improvements.** (1) The spurious *"a treatment period coefficient
interacted with something … was dropped for multicollinearity"* caution fired on all three
explicit UZ arms and all four poisson UZ arms before the fix; after the fix it fires on **none** of
them. That closes `FP-CTZERO-ENGINE`'s "only the explicit arm warns" asymmetry without touching the
screen — the two arms now agree *and* neither warns, because there is nothing left to drop.
(2) The poisson UZ arms were taking the documented **rc-3200 fallback** ("Caution: initial
estimation of ppmlhdfe failed, possibly due to inclusion of an absorbed categorical fixed effect
with only 1 level" — that level-1 effect is the `i.null` ringer head) and refitting with everything
pulled out of `absorb()`. After the fix that caution and that second, slower fit are gone.

**One correction to file against the brief:** `SK-UZ-POISSON` states "No multicollinearity caution
and no other warning fires on this arm." On my build of its own DGP shape the caution *does* fire
on the poisson UZ arms (`run_before.log`, arms C1–C4). The finding's substance is unaffected; the
"silent" characterisation is design-dependent and should not be leaned on.

### Does it move a reported number?

**Yes — but only for specifications the project has never blessed.** It moves `main_att`,
`main_ame`, every post-treatment `ES_att`/`ES_ame` horizon, and the corresponding se/t/p/CI, for
runs that combine `contreatcontrols(untreatedZero)` with `customvce()`, `unconditionalse` or
`poisson`. Magnitudes above: 111% and a sign flip on the linear arms; the entire effect on the
poisson arms. Pre-treatment rows do not move (they were identical across paths before, and still
are). `e(N)`, `e(r2)`, `e(df_r)`, `e(rmse)` and `e(rss)` do **not** move anywhere — the ringer column
and the dropped cell slope are a rank-neutral swap (k columns with one dropped becomes k-1
columns), which is a useful patch self-check: if the residual df moved, the patch did something
other than what is intended. Measured explicitly in `probes/prop-UZ/supp3_*.log`: on the linear
exact-fit design `e(df_r) == 293` on the default, `customvce(robust)`, `unconditionalse robust`
and `ccnoabsorb robust` arms, identical across the pinned build, P1 and P1+P2.

**Re-bless cost: zero.** Every untreatedZero call site in the repo runs on the default engine —
`test/suite/configs.csv:15,16` (`cluster(i)`), `test/suite/configs_ext.csv:270` (`robust`) and
`:298` (an r(198) pairing-refusal test), `test/runbattery.do:147,152`,
`test/suite/newfeat/fastnose_grid.do:277`, `test/test_amescreen.do:305` (`robust`),
`test/test_wave3fixes.do:65` (`robust`), `test/test_permute2.do:193` (refused at parse). None of
them carries `customvce()`, `unconditionalse` or `poisson`. I verified by reading the files, not by
running the suite. The battery baselines `test/baseline/cont_uz.csv` and `cont_uz_ss.csv` are
default-engine and should come back bit-identical.

### Runtime cost

**Zero, and slightly negative.** The edit is two macro comparisons evaluated once at
command-construction time, before any data touch. No `gen`, no `count`, no sort, no extra pass.
On the affected runs it *removes* one regressor (plus, on the ppmlhdfe path, one absorb level and
the entire rc-3200 refit), so large-data runs get faster, not slower. This is not a guard: nothing
is added to the hot path.

### Blast radius

Exactly `contreatcontrols(untreatedZero)` on the three non-default engines. `basic`, `byTreated`,
`byTtre`, `byCohort`, `byTime`, every combination with `lattice()`, and binary mode see an
identical gate — binary mode never enters the `contreatcontrols` block at all. Verified by
measurement, not by inspection: PART B of `body.do` runs nine `basic`/`byTime`/`byTreated` arms,
three `basic + lattice(all 3)` arms and two `poisson + basic` arms across default / `robust` /
`customvce` / `unconditionalse` / `ccnoabsorb`, and dumps `e(b)` in `%21x`. **Every one is
bit-identical before and after.** The only tags that changed in the whole diff are UZ tags.

### Alternatives rejected

- **Caution instead of a fix.** A caution on a point estimate that is wrong by 111% and
  sign-flipped is not a remedy, and there is a right answer available — the default path already
  computes it. Rejected.
- **Refuse `untreatedZero` + `customvce()`/`unconditionalse`/`poisson`.** Removes a legitimate and
  documented capability (`customvce()` is documented at `sthlp:270-273` as a variance switch;
  poisson + a dose that only exists once treated is precisely the case untreatedZero is *for*) to
  avoid fixing three tokens' worth of gate. Rejected.
- **Make the default reghdfe arm warn** (`FP-CTZERO-ENGINE` as originally written). The refuter
  inverted the culprit: the default arm is the correct one. Warning on it would be wrong. Also
  worth recording *why* the default arm was silent: the ringer went into `absorb()`, where reghdfe
  resolves collinearity without marking `o.` coefficients, so the screen at `ado:11890` cannot see
  it by construction. After this fix nothing absorbed is rank-deficient under UZ, so the blind spot
  is unreachable there.
- **Clear `contreatcontrolinserterFE`/`ABS` inside the untreatedZero branch at `ado:1313`.**
  Numerically equivalent, but it leaves the gate's stated rule wrong at the place a future reader
  will read it, and it leaves the `null` tempvar generated for nothing. Rejected on maintenance
  grounds only.

### The test that proves it

An **oracle-free cross-engine identity**, which is the invariant the defect violates and the one I
would put in the suite:

> For any `contreatcontrols(untreatedZero)` specification, `main_att`, `main_ame` and every
> event-study horizon must agree across `{default, ccnoabsorb, customvce(robust),
> customvce(cluster <v>), unconditionalse robust}` to rtol 1e-10 (1e-12 on the reg+clustered pair,
> per the project's non-reproducibility rule), and under `poisson` across
> `{default, customvce(robust), ccnoabsorb}`.

It needs no truth and would have failed loudly on the pinned build (ratio −0.115, and 0 vs 0.44
under poisson). Backed by the closed-form leg: the exact-fit design above makes
`ame = mean(beta_it)` and `att = mean(kappa_it + beta_it*dose)` computable in closed form with
`e(r2) == 1`, and the post-fix run hits both to 5e-16. Runner:
`probes/prop-UZ/run_before.do` / `run_after.do` (they `do` the shared `body.do`; ~3 minutes each).

---

## UZ-P2 — `contreatpoly(2+)` under untreatedZero re-injects a pooled dose term, and `att` is still engine-dependent

**New defect, found while verifying UZ-P1.** Not in the brief. UZ-P1 does **not** fix it.

### Claim

`ado:1314`'s untreatedZero rewrite reuses `contreatpolypart`, which `ado:1110` builds with the
factorial operator: `c.dose##c.dose` at `contreatpoly(2)`. Attached to the cell dummies the
spelling becomes `i._i#i._t#i.TF#c.dose##c.dose`, which Stata expands to
`i._i#i._t#i.TF#c.dose` **+ `c.dose`** + `i._i#i._t#i.TF#c.dose#c.dose`. That bare `c.dose` is a
pooled baseline dose term — exactly the `g(f(Contreat))` that `sthlp:1033` says untreatedZero sets
to zero — and it is invisible to the `att` contrast. So at `contreatpoly > 1` the untreatedZero
model has a pooled dose control after all, put there by its own spelling rather than by the
`ado:1234` gate. Verbose output confirming the spelling: `probes/prop-UZ/poly_diag.log:190`.

### Measured

Quadratic exact-fit design (`y = a_i + g_t + kappa_it + beta_it*dose + gamma_it*dose^2`, dose ≡ 0
off treated, `e(r2) == 1`), TRUE ame = 0.907598611901, TRUE att = 2.645132209592
(`probes/prop-UZ/supp2_body.do`; `supp2_p1.log` = UZ-P1 only, `supp2_p2.log` = UZ-P1 + UZ-P2):

| arm | att with P1 only | att with P1+P2 | ame (both) |
|---|---|---|---|
| `contreatpoly(2)` default / `robust` / `ccnoabsorb` | 2.645132209592 | 2.645132209592 | 0.907598611901 |
| `contreatpoly(2) customvce(robust)` | **0.264358054282** | **2.645132209592** | 0.907598611901 |
| `contreatpoly(2) unconditionalse robust` | **0.264358054282** | **2.645132209592** | 0.907598611901 |
| `contreatpoly(2) lattice(all 2) customvce(robust)` | **0.264358054282** | **2.645132209592** | 0.907598611901 |
| `contreatpoly(3) customvce(robust)` | **0.446396498494** | **2.827170654624** (= default) | 1.073333333333 |

The reported `att` is wrong by a factor of ten, `e(r2)` is 1.0000, and — unlike the UZ-P1 case —
**no caution of any kind fires** (`poly_diag.log:186-232` shows the whole run: `att .44885224`,
`se 2.058e-11`, `p 0.000`, no message). The `ame` row is *not* affected, because each cell's total
dose slope is (cell coefficient + pooled coefficient) and `margins`' derivative picks up both; it
is the `att` contrast, which switches the cell terms off but not the pooled one, that loses it.

### The edit

`ado:1314`. Now:
```stata
      local contreatinj   `latticepart1hash'#`contreatpolypart' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
```
Becomes:
```stata
      local contreatinj   `latticepart1hash'#c.`contreat'
      forvalues uzq = 2/`contreatpoly' {
        local uzterm c.`contreat'
        forvalues uzr = 2/`uzq' {
          local uzterm `uzterm'#c.`contreat'
        }
        local contreatinj `contreatinj' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'#`uzterm'
      }
      local contreatinj `contreatinj' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash'
```
Each polynomial degree becomes its own `#`-joined, cell-interacted term, so no pooled main effect
is ever created. `c.x#c.x` remains legible to `margins` for `dydx` (the `ame` numbers are
unchanged and still exact), which is why this keeps `ado:1100`'s stated requirement that "the
polynomial terms are legible to margins".

At `contreatpoly(1)` the loop body never runs and the macro reduces to
`` `latticepart1hash'#c.`contreat' i.`_i'#i.`_t'#i.`targetflag'`latticepart1hash' `` — **textually
the same string as today**, so poly-1 untreatedZero cannot move. Measured: the poly-1 arms
(`supp2_body.do` tags G1/G2) are identical across `build_patched` and `build_patched2`, and so are
the `contreatcontrols(basic) contreatpoly(2)` arms (tags Bg1/Bg2) — the edit is inside the
untreatedZero branch (ado:1313) and cannot reach any other control type.

### Moves numbers / cost / blast radius

Moves `main_att` and every `ES_att` row for `untreatedZero` + `contreatpoly(>1)` **on the explicit
engines** (`customvce`, `unconditionalse`, and — untested, see below — `poisson`). Does not move
`ame`. Does not move `contreatpoly(1)`, does not move any other `contreatcontrols` type, does not
move a single suite config (no untreatedZero config in the repo uses `contreatpoly`). Runtime: the
loop runs `contreatpoly` times over macros at construction time; the fitted model has **one fewer
column** than today. Nothing added to the hot path.

### Honest limits

- `untreatedZero + contreatpoly(2) + poisson` **is** measured (`probes/prop-UZ/supp3_*.log`, tags
  Z1-Z3; exact-fit log-linear quadratic DGP, TRUE xb-scale ame 0.357887018016). Pinned build:
  default `att -0.133854721893 / ame -0.108345031041`, i.e. wrong on the default engine too,
  because `poisson` is in the `ado:1234` gate. With UZ-P1 only: default `att 0.514894660179 /
  ame 0.357887018016` (fixed) but `customvce(robust)` `att -0.145966985190` while its ame is
  right - the UZ-P2 defect, now visible on the poisson path. With UZ-P1 + UZ-P2: every arm
  `att 0.514894660179 / ame 0.357887018016` (`ccnoabsorb` too). So UZ-P2 is needed on the poisson
  path as well, and the two patches together restore cross-engine identity there.
- I did **not** establish *why* the default reghdfe arm is right today at poly ≥ 2 — presumably
  reghdfe drops the pooled column rather than a cell slope. It is right on both seeds I ran, and it
  stays right after the edit, which is what matters for the patch; but "the default arm is right by
  luck of drop order" is a fair reading and is an argument for making the column not exist.
- `fe()` combined with untreatedZero on an explicit engine is untested under either patch.

### Alternatives rejected

- **Leave it and document it.** The failure is a silently wrong headline `att` with `R2 = 1.0000`,
  `p = 0.000` and no message. A help sentence is not proportionate to that.
- **Reuse `contreatpolypart` but strip the pooled term afterwards** (e.g. add `o.` or a separate
  omission). Fragile against Stata's expansion order and unreadable. Rejected.
- **Route poly ≥ 2 through the ringer instead** (put the pooled dose in `absorb()` so margins never
  sees it). That is the mechanism UZ-P1 just removed, and under untreatedZero it is not
  rank-neutral. Rejected.

### The test that proves it

Same cross-engine identity as UZ-P1, run at `contreatpoly(2)` and `contreatpoly(3)` and with
`lattice()`, plus the quadratic exact-fit closed form:
`ame = mean(beta_it + 2*gamma_it*dose)`, `att = mean(kappa_it + beta_it*dose + gamma_it*dose^2)`,
`e(r2) == 1`. Runner: `probes/prop-UZ/supp2_p1.do` vs `supp2_p2.do` (~4 minutes each).

---

## UZ-P3 — what should *not* change (no code)

1. **No new caution, no new refusal, no runtime screen for the untreatedZero condition.** It is
   tempting to add "dose is not 0 on every untreated observation" as a check, since after UZ-P1
   there is no pooled control left to absorb the violation. I recommend against it: it costs a full
   pass over the estimation sample on every contreat run, `sthlp:995` already states the
   consequence ("wooldid will return poorly identified estimates"), and the behaviour under
   violation is now *consistent* across engines, which it was not before. Measured on a deliberate
   misuse arm (dose non-zero off treated cells, `supp_body.do` tags S7/S8/S9): before the fix
   default 1.323333333333 vs explicit −0.086666666667; after the fix all three engines return
   1.323333333333. Consistently-mis-specified beats inconsistently-mis-specified, and the user was
   warned.
2. **No change to the multicollinearity screen at `ado:11890`.** Its blind spot is structural —
   terms placed inside `absorb()` are resolved by reghdfe/ppmlhdfe without an `o.` mark, so the
   screen can never see a rank problem created there. Closing that would mean reading engine
   internals on every run. After UZ-P1 the only thing wooldid puts inside `absorb()` that could
   create one is the DF-1 ringer set, which DF-1 measured to be rank-neutral. Not worth the cost.
3. **One optional help clause.** `sthlp:1033` ends "specifying untreatedZero makes other technical
   changes to the estimation process". Once UZ-P1 ships, the honest and useful addition is that the
   `g() = 0` model is fitted identically on every engine — reghdfe, ppmlhdfe, `reg` and `poisson` —
   so `customvce()` and `unconditionalse` remain what `sthlp:270-273` promises they are: variance
   switches. Pure documentation; it describes the post-fix code and would be false today.

---

## Untested / open

- `untreatedZero + fe()` on an explicit engine — unmeasured under either patch.
- `bootstrap()` under untreatedZero + `customvce()`: the replications re-enter the same code path,
  so the fix should carry, but I did not run one. (`permute()` is refused with untreatedZero at
  `ado:5842`, so there is no permutation leg.)
- Adjacent observation, **not** a proposal and not in scope: `contreatcontrols(basic)` on data where
  the dose *is* zero off treated cells (i.e. where untreatedZero is the right option and the user
  picked the wrong one) is degenerate in the same way — `poisson + basic` on my exact-fit dose-zero
  design returns `att`/`ame` of −1e-32 against a truth of 0.605/0.400, with only the generic
  multicollinearity caution (`run_before.log` / `run_after.log`, tags CB1/CB2 — unchanged by both
  patches). That is arguably user error the help already covers, but it is the mirror image of the
  defect above and someone may want to look at it.
