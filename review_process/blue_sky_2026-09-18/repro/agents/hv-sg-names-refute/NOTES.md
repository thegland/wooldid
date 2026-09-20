# hv-sg-names-refute — SKEPTIC notes

Role: refute hv-sgname-01 and hv-sgname-02 from an independent repro.
All evidence below is from MY OWN do-files on MY OWN data. I read the
original prober's p1_names.do only to know what the claim was; I never
ran it.

## Files
- `r1_sgnames.do` / `.log` — subgroup e(b)/e(V) names, help-name vs true-name
  addressing, the V splice, adversarial scalar names.  DONE.
- `r1b_consumers.do` / `.log` — coefplot counts done right, e(wdidsgresults_1)
  rownames, bootstrap table names.  RUNNING.
- `r2_esnames.do` / `.log` — att_i / att_itime event-study names.  RUNNING.

## Build assertion (hazard b)
`which wooldid` in r1 log line 67 →
`/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`, `version 1.4 2026-08-19`.
Pinned build confirmed, not the PLUS v1.3.

## r1 results (my own panel: 90 units x 10 yrs, onsets 2004/2006/2008,
## 36 never-treated, sg in {1,2}, N=900)

### e(b) column names — 48 columns, ALL FLAT (coleq is `_` for every one)
```
main_att  main_att_it  main_att_i  main_att_itime
pretreatment_att ... pretreatment_att_itime
ES_att_pre2 pre1 post0 post1 | ES_att_it_pre2 ... post1
main_SG1_att  main_SG1_att_it  main_SG1_att_i  main_SG1_att_itime
pretreatment_SG1_att ... pretreatment_SG1_att_itime
ES_att_SG1_pre2 ... post1 | ES_att_it_SG1_pre2 ... post1
  (then the SG2 block, identically)
```
So: the `_SG<value>_` INFIX is real, there is NO `_s#` suffix anywhere in
e(b), and — the one refutation angle that could have bitten — there is no
equation structure either, so `[main_s1]att` / `[main]att_s1` cannot rescue
the help's name.  **hv-sgname-01's core factual claim reproduces exactly.**

### Addressing, help name vs true name (r1 log 525-648)
| reference                       | rc |
|---------------------------------|----|
| `_b[main_att_s1]`               | 111 |
| `_b[main_SG1_att]`              | 0 |
| `_b[ES_att_post0_s1]`           | 111 |
| `_b[ES_att_s1_post0]`           | 111 |
| `_b[ES_att_SG1_post0]`          | 0 |
| `test main_att_s1 = main_att_s2`| 111 |
| `test main_SG1_att = main_SG2_att` | 0 |
| `lincom main_att_s1`            | 111 |
| `lincom main_SG1_att`           | 0 |
| `esttab ., keep(*_s1*)`         | 111 |
| `esttab ., keep(*SG1*)`         | 0 |
| `esttab ., drop(*_s1*)`         | 111 |

Note: `esttab, drop(*_s1*)` ERRORS 111 here, so the failure is LOUD on
every consumer the help names except coefplot.  `coefplot, keep(*_s1*)`
is not silent either — it prints
`(.: no coefficients found, all dropped, or none kept)` and `(nothing to plot)`.
Only `coefplot, drop(*_s1*)` returns rc 0 with no message (counts pending
in r1b).  Severity-relevant: the prober's "one silent path" is one path,
and it is coefplot's own semantics for a non-matching drop pattern, not
something wooldid does.

### Three coexisting subgroup key regimes — CONFIRMED (r1 log 764-771)
| object | key | rc |
|---|---|---|
| `e(b)` / `e(V)` colnames | `_SG1_` infix | works |
| `e(b_es_att_s1)` | `_s1` suffix | rc 0 (`_SG1` → 111) |
| `e(V_es_att_s1)`, `e(escomp_s1)` | `_s1` suffix | rc 0 |
| `e(wdidsgresults_1)` | bare `_1` | rc 0 (`_s1` and `_SG1` → 111) |
| `e(adversarial_pval_att_s1)` | `_s1` suffix | present |
sthlp:1315 documents the `_s{sgcase}` form, sthlp:1318 the bare
`{sgcase}` form, sthlp:1398 the `_s[SG]` form — all three correct.
Only the e(b)/e(V) infix is undocumented. `grep -c '_SG' wooldid.sthlp` == 0
(bare "SG" appears 11 times, all as the prose abbreviation for "subgroup").

### The SG event-study covariance SPLICE lands correctly (my check, not the
### prober's).  e(V)["ES_att_SG1_pre2","ES_att_SG1_post0"] ==
### e(V_es_att_s1)["pre2","post0"] to the bit on 10 of 16 cells and to
### 1 ulp on the other 6 (e(V_es_att_s1) is itself asymmetric at the ulp
### level; ereturn post symmetrizes).  Well inside rtol 1e-10. NOT a finding.

### The printed table never shows an e(b) name
`Results for Subgroup sg == 2:` with equation-style rows `main | att`,
`ES_att | post0`.  So the help really is the user's only stated route to
the name — the prober's "aggravating" note is accurate.

## r1b / r2 / r3 / r4 / r5 results — see FINDINGS.md for the writeup
- r1b: coefplot counts 36 / 36(drop _s1, silent no-op) / 24(drop SG1) / 12(keep SG1).
  boot_wald_sgresults_1 rc 0, _s1 rc 111 → bare-value key on bootstrap sg tables.
- r2: att_i alone → ES_att_it_*; ES_att_i_post0 r(111).  Printed table header reads
  `ES_att_it` even though only att_i was requested.  att att_i together:
  ES_att_post0 = 0.5787 vs ES_att_it_post0 = 1.2262 (truths 0.5909 / 1.24).
- r3: ES_att_it_post0 bit-identical (+1.39e70e8795764X+000) under att_it/att_i/att_itime.
- r4: contreat + att_i → main_ame_i exists but ES block is ES_ame_it_*, so
  ES_ame_i_post0 r(111).  sthlp:239 makes the same false promise as sthlp:1306.
  bootdist / V_boot colnames carry the _SG infix (consistent with e(b)).
- r5: **NEW FINDING** prefix-nested names + Stata `_b[]` unique-prefix resolution:
  on an att_i-only fit `_b[main_att]` silently returns main_att_i (2.55x wrong);
  under contreat, `_b[main_ame]` silently returns main_ame_i (SIGN FLIP:
  +0.077 → -0.030).  lincom/test echo the expanded name; `di _b[]` does not;
  esttab errors.  Fires only when exactly one family member is posted.

## STATE: COMPLETE.  Verdicts: hv-sgname-01 CONFIRMED, hv-sgname-02 CONFIRMED.
## New: sk-prefix-01 (medium), sk-amename-02 (low-medium).
- r6: `set varabbrev off` (Stata default is ON) makes `_b[main_att]` and
  `lincom main_att` error r(111) on an att_i-only fit, while `_b[main_att_i]`
  still works.  No dataset variable named main_att exists (`confirm variable`
  r(111)).  So sk-prefix-01's mechanism is Stata's varabbrev applied to
  coefficient names, and there is a user-side remedy off the default path.
