<!-- PS-46 research, commissioned 2026-09-10 to answer the author's conditional ruling.
Method: 7 agents, ~34 measured probes run in Stata batch mode against the shipped build.
No repository file was modified; the two remedy states were synthesized in-session.
Fixture test/data/dta_binary.dta, permute(shuffle) permreps(20) rseed(4).
Every row in the evidence table is an observed run, not an inference. -->

# Decision memo: does the permutation fast-mode default break esttab-class table tools?

## 1. Answer in one line

**Yes** — on 19 estout-family probes and 15 built-in probes, all at rc 0 for the estout family: `esttab` never refuses and instead prints a **blank standard-error row under a footer that asserts "Standard errors in parentheses" and a full star legend**, while Stata's own built-in `estimates table S, b se` **refuses loudly** with rc 321. So "breaks" resolves into two different failures — silent misreport for the estout family (including `.tex` and `.csv` files that outlive the session), loud stop for the built-ins — and the ruling's condition is met.

## 2. The evidence table

Two states side by side: the shipped fast-mode default (no `e(V)`), and the ruling's proposed remedy (all-zero 6×6 `e(V)`, synthesized outside wooldid). Analytic control shown where the probe captured it. `estimates table S, b se star` is excluded: it fails rc 198 "option star may not be combined with se, t, or p" in **both** states and on the `permanalyticse` control — a Stata syntax rule, not a VCV finding.

| Tool × command | No `e(V)` (shipped) | Verbatim observed | All-zero `e(V)` (remedy) | Verbatim observed | Analytic control |
|---|---|---|---|---|---|
| `eststo fast` | CLEAN | rc 0; stores | CLEAN | rc 0 | rc 0, identical |
| `esttab fast` | **SILENT_WRONG** | `main_att            0.824` then a row of pure whitespace — not `(.)`, not `(0.00)`. Zero stars. Footer `t statistics in parentheses` + `* p<0.05, ** p<0.01, *** p<0.001`. No warning line anywhere (verified by grep over the block). | **SILENT_WRONG** | `(.)` replaces the blank; same false footer, same legend | `0.824***` over `(13.08)` |
| `esttab fast, se` | **SILENT_WRONG** | blank SE row; no stars; footer `Standard errors in parentheses` | **SILENT_WRONG** | `main_att            0.824` over `                      (.)`; same footer | `0.824***` over `(0.0630)` |
| `esttab fast, b(3) se(3) star` | **SILENT_WRONG** | coefficients honour `b(3)`; SE row blank despite explicit `se(3)`; zero stars despite explicit `star`; legend printed anyway | **SILENT_WRONG** | `(.)`, still no stars | `0.824***` over `(0.063)`; `ES_att_post0/1` starred `***` |
| **`esttab fast, ci`** | **SILENT_WRONG** | blank bracket row under footer `95% confidence intervals in brackets` | **SILENT_WRONG — strictly worse** | `[0.824,0.824]`, `[-0.0457,-0.0457]`, `[-0.0129,-0.0129]`, `[-0.0785,-0.0785]`, `[0.842,0.842]`, `[0.781,0.781]` — fabricated zero-width intervals rendering as real, infinitely precise numbers | `[0.698,0.949]`, `[-0.166,0.0746]`, … |
| `esttab fast, p` | **SILENT_WRONG** | blank under footer `p-values in parentheses`; `e(perm_p)` holds `.0952381` for `main_att` the whole time and is never surfaced | **SILENT_WRONG** | `(.)`; same | — |
| `esttab fast using fast.tex, replace` | **SILENT_WRONG** | `main\_att    &       0.824         \\` then `            &                     \\`; no `\sym{}` anywhere; footnote `\textit{t} statistics in parentheses` + full legend written into the file | **SILENT_WRONG** | same shape with `(.)` | `main\_att    &       0.824\sym{***}\\` / `(13.08)` |
| `esttab fast using fast.csv, replace` | **SILENT_WRONG** | `="main_att",="0.824"` then `="",=""`; footer row `="t statistics in parentheses"` written as data | **SILENT_WRONG** | same with `.` | `="main_att",="0.824***"` / `="",="(13.08)"` |
| `estout fast` (default `cells(b)`) | CLEAN | `main_att         .8235285` … byte-identical to the analytic run | CLEAN | identical | identical |
| `estout ., cells(b(star) se)` | **SILENT_WRONG** | header `b/se`; se line blank | **SILENT_WRONG** | `main_att         .8235285` over `                        .` | — |
| `matrix list e(V)` | HARD_ERROR | `matrix e(V) not found  [r(111)]` | **CLEAN — zeros visible** | `symmetric e(V)[6,6]`, all 21 lower-triangle cells `0` | 6×6, diag `.0039645057 …` |
| `estat vce` | HARD_ERROR | `e(properties) does not describe the status of e(V)  [r(322)]` | **CLEAN — zeros visible** | `Covariance matrix of coefficients of wooldid model` + 6×6 grid of `0` | — |
| `di _se[main_att]` | HARD_ERROR (accurate refusal) | `_se not allowed when e(V) is not present  [r(111)]` | **SILENT_WRONG** | `0` at rc 0 — Stata's own accessor now asserts the ATT's SE is exactly zero | `.0629643` |
| `di _b[main_att]/_se[main_att]` | HARD_ERROR | same r(111) | DEGRADED_LOUD | `.` (0/0 → missing marker) | `13.08` |
| `ereturn display` | DEGRADED_LOUD (honest) | header is just `y \| Coefficient` — SE/t/p/CI columns structurally absent | **SILENT_WRONG** | full header `Coefficient  Std. err.      t    P>\|t\|     [95% conf. interval]` with `.` in all 30 cells, under the banner `(Std. err. adjusted for 80 clusters in i)`; nothing says why | populated |
| `matrix list r(table)` | DEGRADED_LOUD (honest) | 1×6, `b` row only | **SILENT_WRONG** | 9×6 with `se/t/pvalue/ll/ul` all `.` **plus `df 79` and `crit 1.9904502` populated** — downstream code reading `r(table)` gets a full-shaped matrix asserting a df and a critical value for statistics that do not exist | 9×6 real |
| **`estimates table S, b se`** | **HARD_ERROR** | `SZ estimation results do not have e(b) and e(V)  [r(321)]` | **SILENT_WRONG — worst cell** | rc 0; `main_att \|  .82352847` over `         \|          0` under `Legend: b/se` — a fabricated numeric zero, not a dot | real SE |
| `estimates table S, b t p` | HARD_ERROR | r(321) | **SILENT_WRONG** | rc 0; `t` and `p` rows both `.` | real |
| `esttab . using f.tex, se star replace` | SILENT_WRONG | rc 0, falsehood persisted to file | SILENT_WRONG | rc 0, falsehood persisted to file | — |
| **Measured aux route** (`e(se)`+`e(p)` posted, **no** `e(V)`) | — | — | **CORRECT** | `esttab ., se star` → `(0.130)`, `(0.0767)`, … ; `esttab ., p` → `(0.095)`, `(0.857)`, … ; rc 0 throughout | matches by construction |

The probes do not average out. `estout`/`eststo` are genuinely clean; `esttab` is uniformly silent-wrong; the Stata built-ins are honest without `e(V)` and dishonest with a zero one. That divergence is the finding, not noise.

## 3. The mechanism

**`e(properties)` is not the branch point.** estout never reads it. The missing `e(V)` is.

`estout.ado` line 3050 fills the variance column of its internal `bc` matrix with

```
capt mat `tmp' = vecdiag(e(V))
if _rc {
    capt confirm mat e(se)
    ...
}
```

With `e(V)` absent, the `capt` swallows rc 111; the fallback `confirm mat e(se)` also fails because wooldid posts no `e(se)`; the column is left at estout's leave-empty marker `.y` (set at 3045) and renders as whitespace. `ComputeCoefs_se/_t/_p/_ci_l/_ci_u` propagate the missing, and stars derive from `p`, so none print.

The decisive part is the footer. `esttab.ado` 578-599 builds the note **purely from the requested aux keyword** — `se` → `Standard errors`, `ci` → `95% confidence intervals`, then `' in parentheses'` (599) or `' in brackets'` (597) — with no check that any value was obtained. **The footer is a claim about the option, not about the data.** Same for the star legend, printed unconditionally.

The same code gives the remedy that works: `estout.ado` 3156 defines `local bVs2 "se var t z p ci_l ci_u"`, and 3165-3172 does `capt confirm matrix e(`v')` for each, falling through to `ComputeCoefs_`v'` (the `e(V)` path) **only on failure**. A posted `e(se)`/`e(p)` matrix is therefore used verbatim, ahead of anything derived from a covariance.

The built-ins fail differently: `estimates table` and `estat vce` read `e(properties)`, which under fast mode is `b`, hence r(321)/r(322) — honest refusals that a zero `e(V)` would convert into fabrications.

## 4. What the ruling implies

It breaks them, so the conditional fires. Two code sites, and they are not interchangeable.

**Route A, the ruling as literally written (all-zero `e(V)`).** It must be supplied at `ereturn post`, not bolted on: evidence item 5 shows `ereturn repost V =` on a b-only estimate refuses rc 503 "conformability error". Fast mode currently lands at **`wooldidinstall/wooldid.ado:3878`**, `cap ereturn post \`eb' , esample(...) depname(...) obs(...) dof(...)`, reached because `ev = diag(hadamard(se,se))` is all-missing under `fastnose` so the paired post at 3847 fails and `evpostgood == 0`. A zero matrix would go in at 3847/3878. **Gating is a real requirement, not a detail:** the ruling is about the *permutation* fast-mode default, not user-requested `fastnose` and not the replications. The parent sets `local fastnose fastnose` at **5196** under the `fnperm` predicate at **5144** (`permute()` and not `permtstat` and not `permanalyticse`); the worker would need that fact plus `wdrepmode` empty to distinguish the observed permutation run from a user `fastnose` run and from a nested rep. I am naming the requirement, not proposing a solution to it.

**Route B, the measured route (post `e(se)`/`e(p)` and no `e(V)`).** These are ordinary `ereturn matrix` calls that preserve existing `e()`, so they belong in `wooldid_permwire` beside the stripe-aligned block at **wooldidinstall/wooldid.ado:7467-7471** (`ereturn matrix perm_p = \`pv''`). Two properties matter. First, it renders correctly with no covariance posted at all. Second, it is **estout-specific**: `ereturn display`, `r(table)` and `_se[]` do not read `e(se)`, so Stata's honest refusals survive untouched — the aux route fixes exactly the tools that were silently lying and touches nothing that was already telling the truth.

**But Route B needs a ruling, because it moves a printed number.** `e(p)` has an obvious source: `e(perm_p)` is already posted (`.0952381 .85714286 .85714286 .47619048 .0952381 .0952381`) and is the number the tables print. `e(se)` has **no** designated source. `wooldid_permwire`'s own header comment states the SE column is blank under the fast-mode default *by design* — "the standard-error column is left as the worker posted it … blank under the fast-mode default." The probe used `e(perm_se_shift)` as an arbitrary stand-in. So the open choices are: (i) post `e(p)` only and leave `e(se)` unposted — `esttab ., p` becomes correct, `esttab ., se` stays blank-with-false-footer; (ii) post `e(p)` and `e(se)` from `e(perm_se_shift)`; (iii) post `e(p)`, `e(se)`, and `e(ci_l)`/`e(ci_u)` from `e(perm_ci)` (source-only observation — the CI path was not run, and it is the one path where the zero-V remedy actively fabricates, so it is the highest-value one to fix). Which SE, if any, is yours.

**A help correction is required whichever way you go.** `wooldidinstall/wooldid.sthlp` lines 271-272 currently read: "Under `fastnose`: `e(V)` is not posted (`e(properties)` is `b`), so `test`, `lincom`, **`esttab` standard errors**, and anything else needing a variance-covariance matrix **will refuse**". Measured, that is false: `esttab` does not refuse, it prints a blank row under a footer asserting standard errors are present. If you leave the code as is, that sentence must become something like: "`test` and `lincom` refuse; `esttab` and `estout` do **not** refuse — they print an empty standard-error row beneath their standard footer, so a table built from a fast-mode run will assert standard errors it does not contain." The parallel sentence at 427-428 ("`e(V)`, `e(b_es_*)` and `e(V_es_*)` are then not posted") is accurate as far as it goes but should carry the same warning.

## 5. The case against the remedy, honestly

The ruling set its own test: the artificial covariance must be **"transparently false to anyone who inspects it."** On the evidence it fails that test on the surfaces the ruling was trying to protect.

Count where the zeros are actually visible: **two**, `matrix list e(V)` and `estat vce`. Both are inspection commands a suspicious analyst runs deliberately. Everywhere else the zeros are invisible or worse. On the esttab surfaces the remedy targets, it converts a blank cell into `(.)` — marginally more visible, but the false footer and the star legend survive unchanged, so the grade does not improve. On `esttab ., ci` it is **strictly worse**: `[0.824,0.824]` is not a transparently false number, it is a plausible-looking one, and it goes into `.tex` and `.csv` files that outlive the session. On `estimates table S, b se` it converts an honest rc 321 refusal into a literal `0` printed in an SE row. On `di _se[main_att]` it converts an accurate `_se not allowed when e(V) is not present` into `0` at rc 0. On `r(table)` it manufactures a 9×6 matrix carrying `df 79` and `crit 1.9904502` for statistics that do not exist — which downstream scripts will read.

So the framing in the prompt is correct and the evidence supports it: a fabricated number in a table someone pastes into a paper is worse than a loud stop. **The all-zero matrix, as specified, does not meet the ruling's own condition.** It trades one silent failure for a new fabrication and buys visibility on exactly two inspection commands.

The reconciliation is that the ruling's *goal* — let the tools run, without misleading anyone — is achievable, just not by that instrument. The aux route was not available to weigh when the ruling was made because it had not been measured. It makes `esttab` print the real permutation numbers, leaves `e(V)` absent so `test`, `lincom`, `_se[]`, `estat vce` and `estimates table` keep refusing honestly, and fabricates nothing. That is the substitution I would put to you: **honour the ruling's intent via posted `e(se)`/`e(p)` (and possibly `e(ci_l)`/`e(ci_u)`), and drop the zero matrix.**

## 6. Residual uncertainty

The zero-V state and the posted-aux state were both **synthesized from outside wooldid** — constructed via `ereturn post e(b) J(6,6,0)` and by adding aux matrices to a posted result; the repo was not modified. The fresh `ereturn post` for the zero-V block drops `e(df_r)` and wooldid's other `e()`. That does not move any result here (`se=.` gives `t=.` regardless of df; `se=0` gives `[b,b]` regardless of critical value), but it means the zero-V rows are a faithful reconstruction rather than a shipped-code observation.

**Star rendering from `e(p)` is untested.** Every `e(perm_p)` in this fixture is ≥ `.0952`, so zero stars is *correct* in the aux probes and the star path was never exercised. That is the single most important untested thing about Route B.

**Which SE to post is unmeasured and unruled** — `e(perm_se_shift)` was an arbitrary probe choice. **`e(ci_l)`/`e(ci_u)` from `e(perm_ci)` is a source-only inference**, never run. **Non-estout writers were not probed at all**: `outreg2`, `coefplot`, `regsave`, `estimates table` under `esttab`-style wrappers, and `putexcel`-driven pipelines are unknown territory; `coefplot` in particular reads `e(V)` and would plot zero-width whiskers under Route A.

What would settle it: implement Route B on a branch, and re-run the estout battery on a fixture with at least one `e(perm_p) < .05` so the star path is exercised, plus a `coefplot` probe under both routes.
