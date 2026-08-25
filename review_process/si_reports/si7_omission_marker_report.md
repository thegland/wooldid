# SI-7 Diagnosis: `ext_contract10_rcs_factor_fe` — r198 "5o. invalid name"

Date: 2026-08-25. Machine: Stata 19 MP (`/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp`).
Builds: NEW = working tree `wooldidinstall/wooldid.ado` (HEAD d0f2625, branch bugfixes); OLD = `main:wooldidinstall/wooldid.ado` (301b706), staged to `old_wooldid.ado` in this directory.
Config (configs_ext.csv line 354, prep `ext_contract10_rcs_factor_fe` = `replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .`):

```
use test/suite/data/survey_rcs.dta
wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) fe(i.group##i.t)
    esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) safety(off)
```

All probe .do/.log files referenced below live in `/Users/hegland/.claude/jobs/640b9294/tmp/si7/`.

---

## (a) Exact failing line, construct, and offending name (NEW build)

Reproduced: rc 198, error text `5o. invalid name` (probe1_new.log, RUN_A_RC=198).

Failing site (trace: probe1_new.log lines 16880–16970): **wooldid_worker, "effects by horizon for sgs" block, NEW wooldid.ado line 2546**:

```stata
local relyrinj = `relyrsg' -`min_t' -1        // line 2546 (NEW); line 2304 (OLD)
```

which expands to the illegal expression `local relyrinj = 5o. -6 -1`. The token `5o.` reaches the
expression evaluator because the upstream name-sanitizing chain (NEW 2526–2539) extracts the
subgroup event-study relative periods from the **colnames of `margins` r(table)** and strips only
`bn`, the `_t` tempvar name, and `1.targetflag`:

```stata
local names: colnames `nameget'                       // -> "5bno.__000003 7.__000003"
local names = subinstr("`names'", "bn","",.)          // -> "5o.__000003 7.__000003"
local names = subinstr("`names'", "`_t'","",.)        // -> "5o. 7."
local names = subinstr("`names'", "1.`targetflag'","",.)
...
foreach relyrsg in `names' {                          // "5o." and "7."
  local relyrinj = `relyrsg' -`min_t' -1              // "7." evaluates; "5o." -> r198
```

There is **no handling for the `o` (omitted / not-estimable) marker**. The offending column is the
subgroup-1 pre-period cell of the per-subgroup horizon margins call (NEW line 2500):

```stata
qui margins [pweight=`w_att'], dydx(`targetflag') subpop(if `ifstatement' & `targetflag'==1
    & `_t' >= 5 & `_t' <= 7 & group == 1) noestimcheck over(`_t') df(23)
```

whose r(table) on NEW is (probe3_instr.log 9275–9291):

```
            0b.tf:5o.   0b.tf:7o.   1.tf:5o.    1.tf:7.
     b          0           0           0       .69760818120352
    se          .           .           .           .
```

i.e. margins marked the relyr −2 (event-time 5) margin for subgroup 1 **omitted** (`5o.`), and the
sg-relyr parser crashed on the marker. This is *not* a dropped regression coefficient: the
underlying reghdfe e(b) contains all 253 cohort×period×targetflag dummies with **no** `o.` flags
(only `o.` tokens in the entire trace are at the margins site).

## (b) OLD-build behavior — rc 0, but the output is not trustworthy

OLD on this machine: rc 0, with or without `safety(off)` (probe2_old.log RUN_A_RC=0, RUN_B_RC=0).
At the **same site** OLD's identical margins call returns the same 4 columns, but the subgroup-1
pre-period cell is *not* flagged: colname `1.tf:5.` with **b = 8.881784197e-16 and se = .**
(probe3_instr.log 13863–13878). "5." happens to evaluate as a number, so the identical parse chain
survives by accident.

The rc-0 results are quietly degenerate on this dataset:

- **Every SE, t, p, and CI in every returned matrix is missing** — e(wdidmainresults),
  e(wdidsgresults_1), e(wdidsgresults_2), e(pwsgs) all have se=t=p=CI=`.` (probe2_old.log
  55863–55985). Root cause: an all-zero clustered VCV (bit-exact sum(e(V)) = +0.0 in both builds,
  probe3 diagnostics; every margins SE missing; and the NEW build's own screen reports
  "estimation returned an all-zero variance-covariance matrix" on the adjacent no-fe variant).
  OLD has no such screen at all — the string is absent from its source — so OLD predates the
  detection and reports the degenerate table without comment.
- The pre-trend "estimates" are machine zeros (8.88e-16), i.e. the same non-identified contrast NEW
  flags as omitted, computed as numerical noise.
- e(wdidsgresults_1) has 3 rows vs subgroup 2's 4: subgroup 1's `pretreatment` row is silently
  absent (the `truerow`/`sghaspre` logic decided sg1 "has no pre"), an asymmetry a user gets no
  message about.

So the pinned old-side "success" is a results table with no usable inference and artifact zeros for
pre-trends. The pin keeps the suite green but the behavior it preserves is not meaningfully correct.

Coefficient-name diff at the failure site, old vs new: the *regressions* hand back byte-identical
name lists (254 colnames, identical modulo tempvar labels — bvars_old.txt vs bvars_new.txt) and
bit-identical sum(e(b)). The only name difference anywhere is the margins r(table) colname
`5.`/`5bn.` (OLD) vs `5o.`/`5bno.` (NEW) for the subgroup-1 pre-period cell.

## (c) Minimal failing ingredient set and concrete user cases

Bisection on NEW (probe5a/probe5b), plus OLD cross-checks (probe6):

| variant | NEW | OLD |
|---|---|---|
| pinned config | 198 | 0 |
| − subgroups(group) | **0** | – |
| − fe(i.group##i.t) | 198 | 0 |
| fe(i.group) only | 198 | – |
| fe(i.t) only | 198 | – |
| − [pw=survey_weight] | 198 | **198** |
| − esfixedbaseperiod | 198 | – |
| safety ON (no safety(off)) | 198 | 0 |
| esprelength(0) espostlength(1) | **0** | – |

Minimal ingredient set: **`att` + `subgroups()` + an event-study window containing a horizon whose
subgroup-level margin `margins` deems not estimable** (here the pre period; esprelength(1) suffices).
fe(), pweights, esfixedbaseperiod, cluster() and the safety setting are all irrelevant — including
`safety` ON, which cannot catch it because the deliberate screen reads e(b) (clean here), not
margins output.

What a user must type to be at risk:

```stata
wooldid y i t ttre, att subgroups(g) esprelength(1) espostlength(1) ...
```

i.e. **any subgroup event study with at least one pre/post horizon** — bread-and-butter usage of two
headline options. The crash fires only on datasets/specs where some subgroup×horizon margin comes
back omission-marked. Concrete, reasonable specs that produce that:

1. `subgroups(g)` combined with `fe(i.g##i.t)` — the user's own FE absorb the subgroup×time
   variation, so subgroup pre-trends are non-identified *by construction* (this config);
2. a subgroup with no (or fully degenerate) observations at some requested horizon;
3. numerically knife-edge designs (this suite dataset), where the determination flips with
   epsilon-level input changes.

Spec 1 is a plausible real-world mistake (group-specific time effects plus effects-by-group); the
failure mode is a hard r198 with the cryptic message `5o. invalid name` after the main estimation
already succeeded — no results, no hint.

## (d) Regression vs luck: verdict — shared latent fragility; the flip is numerical luck

Evidence:

1. **The fragile code is byte-identical in both builds** (NEW 2526–2546 = OLD 2284–2304; same
   subinstr chain, same `local relyrinj =` evaluation). None of the new-build screens
   (safety/truerow/dupcol, the `o.`-scan at NEW 1339–1427) are in the failing path; safety was off,
   and with safety on the result is unchanged (198).
2. **The regressions are identical**: same reghdfe command, same 254 coefficient names, bit-identical
   sum(e(b)), no omitted coefficients, in both builds (probe3 diagnostics).
3. **The single input difference is w_att storage precision.** OLD generates the estimand weight as
   float (`qui gen `w_att' = ...`, OLD line 1100); NEW generates it double (NEW line 1304 — commit
   3cb0c41 "B21: gen double for ... estimand-weight copies", 2026-07-21). survey_weight is stored
   double; OLD's float copy truncates it (weight sums differ at the 8th significant digit:
   `+1.07851eb851eb8X+007` vs `+1.07851ef000000X+007`).
4. **Causal one-line toggle, both directions** (probe4): NEW with `gen `w_att'` (float) → rc 0;
   OLD with `gen double `w_att'` → rc 198 "5o. invalid name". Nothing else changed.
5. **OLD crashes on its own** on a *simpler* spec: the unweighted variant (no [pw]) gives r198 on
   the untouched OLD build (probe6). And the colleague's Stata-17 run of OLD on the pinned config
   also gave r198 — margins-internals version luck, same class.

Verdict: the NEW-side r198 is **not** a behavioral regression in name/omission handling. B21's
float→double weight fix (a legitimate correctness change) perturbed a knife-edge margins
estimability determination on this adversarial dataset, flipping a latent, build-independent crash
from "hidden" to "expressed". The old build's rc 0 is omission-pattern luck — it fails identically
on Stata 17, on the unweighted variant, and the moment its weights become exact.

## (e) Deliberate omission-pattern exploitation vs accidental fragility — the boundary

Deliberate (intended, robust to `o.`):

- **Pooled ES relyr extraction** (NEW 1318–1335): relative periods are derived from **e(b)
  colnames** with `regexm("`bvars'","[0-9].`_i'#`horiz'.`_t'#1.`targetflag'")` — written precisely
  so `o.`/`b.`-marked terms do not match ("anynumber not followed by a o or a b" per its comment).
  This is the author's omission-pattern reading working as designed; it is why the pooled ES path
  survives this config in both builds.
- **The safety screen** (NEW 1339–1427 / OLD 1206–1245): scans e(b) for `o.`-marked treatment
  dummies, distinguishes duplicates, empty cells, and real drops, and errors 322 with an
  explanatory message. Deliberate, but scoped to the **regression's** e(b) only.
- Base-period omission by construction (targetflag=0 at the reference period, so no `_t==6`
  dummies), and display-side skipping of `o.` columns (NEW ~4563).

Accidental (the bug): the **subgroup**-horizon path does not derive relyrs from e(b); it re-derives
them by string-stripping **margins r(table) colnames** (NEW 2526–2539) and then **arithmetically
evaluates the residue** (line 2546). That chain anticipates the `bn` base marker but not the `o`
omission marker, and its only fallback (`glevelsof` when the residue is empty) does not trigger
because `"5o. 7."` has wordcount 2. The precise boundary: **wooldid deliberately reads omission
patterns in the regression's e(b); it never anticipated omission markers in margins-over() output,
and the subgroup ES path is the one place where such a marker flows into an expression evaluator
instead of being filtered.** (For reference: the sg pairwise path did process an o-marked name
`1bno.group` without crashing in the OLD trace — probe2_old.log 40886–40889 — though its mechanism
was not traced further here.) The oversight conclusion rests on the pooled-path-vs-sg-path
contrast: the pooled ES path filters omission markers by design; the sg-horizon path evaluates
them as numbers.

## Recommended dispositions (options only — nothing implemented)

1. **Harden the sg-horizon parser** (minimal fix): strip omission/base markers robustly — e.g.
   replace the residue extraction with a regex/`_ms_parse_parts` pull of the leading digits, or add
   `local names = subinstr("`names'","o.",".",.)` *after* the `bn` strip plus a per-token
   `if regexm("`relyrsg'","o")` guard — and **drop o-marked horizons from that subgroup's matrix**
   (mirroring what the pooled path already does by regex). Optionally print the analogue of the
   pooled caution: "the relyr −2 effect for subgroup 1 was not estimable and is omitted."
2. **Align the sg path with the pooled path**: derive per-subgroup relyrs from e(b) names (or the
   existing glevelsof fallback) instead of margins colnames, removing the fragile chain entirely.
3. **Safety-screen coverage note**: the deliberate screen cannot see margins-stage omissions; if the
   author wants subgroup non-estimability to be a hard, *explained* error under safety-on, add a
   check of margins colnames for `o` markers at the sg site (an error 322 with subgroup/horizon
   named beats r198 "5o. invalid name").
4. **Spec-level warning (optional)**: detect `fe()` terms that interact the `subgroups()` variable
   with the time variable and warn that subgroup pre-trends/ES horizons are then not identified.
5. **Suite pin disposition**: re-label the pin — old rc 0 is not certified-good behavior (V≡0, no
   SEs anywhere, artifact-zero pre-trends, silently missing sg1 pretreatment row); the pin records a
   shared latent crash whose expression flips with weight precision and Stata version. After fix
   (1)/(2), the expected outcome for this config is rc 0 on NEW with the sg1 pre-period row
   reported as omitted (or a clean, explained error), not restoration of the old output.

## Evidence inventory (this directory)

- probe1_new.do/.log — NEW repro + trace (RUN_A rc198; trace to failing line, 16880–16970)
- probe2_old.do/.log — OLD traced full run (rc0), full e() dump, safety-on rc0
- probe3_instr.do/.log, new_instr.ado, old_instr.ado — instrumented site diagnostics (subpop counts,
  bit-level b/V/weight sums, r(table) both builds)
- probe4_toggle.do/.log, new_floatw.ado, old_doublew.ado — causal precision toggle (NEW+float→0,
  OLD+double→198)
- probe5a_bisect.do/.log, probe5b_bisect.do/.log — option bisection on NEW
- probe6_oldvariants.do/.log — OLD no-pw → 198, OLD no-fe → 0
- bvars_old.txt, bvars_new.txt — e(b) colname lists (identical modulo tempvar labels)
- old_wooldid.ado — staged OLD build (main)
