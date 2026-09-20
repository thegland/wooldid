# hv-sg-names — FINDINGS

Verifier: independent reproduction of a recon claim about e(b) coefficient names.
Build asserted in every log: `.../scratchpad/bluesky/build/wooldid.ado` (line 57-59 of each .log).
All data generated from scratch inside each probe; no repo data, no repo generator, no test suite.

---

## HV-SGNAME-01 — CONFIRMED (certain). Subgroup e(b) names are an INFIX `_SG<value>_`, not the help's suffix `_s#`; and the true convention appears nowhere in the help.

**Repro:** `/private/tmp/.../probes/hv-sg-names/p1_names.do`
(`stata-mp -e do p1_names.do`; see `p1_names.log`)

**sthlp:1303-1307 says**

> The coefficient vector e(b) and its covariance e(V) carry these column names, which are what
> esttab, coefplot, lincom, test and the permtstat/bootsave() tables address a coefficient by:
> main_att ..., pretreatment_att ..., ES_att_pre# and ES_att_post# ..., and for the other
> estimand types the same names with att replaced by att_it, att_i, att_itime or ame.
> **Under subgroups() each subgroup's rows take the suffix _s#.**

**Measured** (RUN C_sgA, `subgroups(sgA)` with levels 1 and 2, four estimands, esprelength(2) espostlength(2)):

| help says | actually posted |
|---|---|
| `main_att_s1` | `main_SG1_att` |
| `pretreatment_att_itime_s1` | `pretreatment_SG1_att_itime` |
| `ES_att_pre2_s1` | `ES_att_SG1_pre2` |
| `ES_att_it_post0_s2` | `ES_att_it_SG2_post0` |

The subgroup token is an **infix placed after the family stem and before the estimand /
horizon**, and it carries the subgroup **value** (RUN D_sgB with levels 3 and 17 gives
`main_SG3_att` ... `ES_att_it_SG17_post1`; RUN L with level 0 gives `main_SG0_att`).

**`grep -c '_SG' wooldid.sthlp` == 0.** The convention actually used by e(b) is documented
nowhere in the 326 KB help file. sthlp:1303-1307 is the only passage that documents the
e(b) name scheme at all (the only other mention of any e(b) name is the `event_plot` recipe
at sthlp:1452).

### Three different subgroup keys coexist in one `e()`
Measured on the same fit (RUN C2 / RUN N):

| object | key for subgroup level 3 |
|---|---|
| `e(b)` / `e(V)` column names, `e(bootdist)`, `e(boot_tstud)`, `bootsave()` .dta varnames, `e(perm_n)`, `e(perm_p_lo)`, `e(perm_comparator_b)` | **infix `_SG3_`** (`main_SG3_att`) |
| `e(b_es_*)`, `e(V_es_*)`, `e(escomp_*)` | **suffix `_s3`** (`e(b_es_att_s3)`) — matches sthlp:1315 |
| `e(wdidsgresults_*)` | **bare `_3`** (`e(wdidsgresults_3)`) — matches sthlp:1318 |

`e(b_es_att_SG3)` -> r(111); `e(wdidsgresults_s3)` -> r(111); `e(wdidsgresults_SG3)` -> r(111).
So the help is right about the *matrices* and wrong only about *e(b)/e(V)* — which is the one
surface the help tells the user to address coefficients through.

### Consequence surface (measured, RUN C2 / p4 RUN Q)
| what a user following sthlp:1307 writes | result |
|---|---|
| `test main_att_s1 = main_att_s2` | **r(111)** (loud) |
| `lincom main_att_s1` | **r(111)** (loud) |
| `di _b[main_att_s1]` | **r(111)** (loud) |
| `esttab ., keep(*_s1*)` | **r(111)** (loud) |
| `esttab ., drop(*_s1*)` | **r(111)** (loud) |
| `coefplot, keep(*_s1*)` | **rc 0**, prints only `(nothing to plot)` — an EMPTY graph |
| `coefplot, drop(*_s1*)` | **rc 0, NO message at all, and MEASURABLY a no-op** — see `p6b_coefplot.do`: e(b) has 12 coefficients; `drop(*_s1*)` plots **12**, identical to plain `coefplot`, while the true `drop(*SG1*)` plots **8**. The user gets a plot silently containing every subgroup row they meant to remove. |
| `event_plot, stub_lag(ES_att_post#_s1) stub_lead(ES_att_pre#_s1)` | r(498) "No estimates found" |
| (the true names, as controls) `test main_SG1_att = main_SG2_att`, `esttab ., keep(*SG1*)`, `event_plot, stub_lag(ES_att_SG1_post#)` | rc 0 |

**Aggravating:** wooldid's own printed table never shows an e(b) name. The subgroup block is
headed `Results for Subgroup sgA == 1:` with rows `main | att`, so the printed output does not
teach the name either; and `esttab` abbreviates (`main_SG1_a~t`). The user's routes to the true
name are `mat list e(b)` / `: colnames e(b)` — neither of which the help points at.

**Realistic:** yes. `subgroups()` is a headline feature and the help repeatedly instructs users
to address coefficients by e(b) name (sthlp:1303, sthlp:1452).

**Severity: medium.** Tier-5 doc-vs-code. Mostly loud failures; the one silent path is
`coefplot ..., drop(*_s1*)`.

**Fix guidance.** The divergence is confined to sthlp:1305-1307. The neighbouring passages are
CORRECT and must not be changed: sthlp:1310 (`e(wdidmainresults)` ES rows "sit under ES_att,
ES_att_it, ES_ame and ES_ame_it"), sthlp:1315 (`e(b_es_{type}_s{sgcase})`) and sthlp:1318
(`e(wdidsgresults_{sgcase})`) all match the code. And the fix belongs in the HELP, not the code:
the `_SG<v>_` infix is used consistently across e(b), e(V), `e(bootdist)`, `e(boot_tstud)`,
`e(boot_ptails)`, `e(boot_wald_mainresults)`, the `bootsave()` .dta variable names (and their
`se_` twins), `e(perm_n)`, `e(perm_p_lo)` and `e(perm_comparator_b)` — eleven surfaces measured
here — so renaming to `_s#` in code would break all of them plus every existing user do-file.

**Note on the task's sthlp:484 citation.** Lines 478-492 are the `permtstat`/`permtci`/`permsave`
option prose and contain no instruction to address a coefficient by its e(b) name. The passage
that does is **sthlp:1303-1307**, which names `permtstat`/`bootsave()` among the consumers. The
substance of the claim is unaffected.

---

## HV-SGNAME-02 — CONFIRMED (certain). `ES_att_i_*` and `ES_att_itime_*` do not exist; sthlp:1305-1306 contradicts sthlp:1315.

**Repro:** `p1_names.do` RUN A_att_i / RUN A_att_itime; `p5_semantics_and_overflow.do` RUN R2.

sthlp:1305-1306: "for the other estimand types **the same names** with att replaced by
att_it, **att_i**, **att_itime** or ame" — which promises `ES_att_i_pre#` / `ES_att_itime_pre#`.

**Measured.** `wooldid y unit year ttre, att_i esprelength(2) espostlength(2) robust` posts:
```
main_att_i  pretreatment_att_i  ES_att_it_pre2  ES_att_it_pre1  ES_att_it_post0  ES_att_it_post1
```
`lincom ES_att_i_post0` -> **r(111)**; `lincom ES_att_itime_post0` -> **r(111)**;
`lincom ES_att_it_post0` -> rc 0. Same for `att_itime` alone. The LEVEL rows do keep the
estimands distinct (`main_att_i`, `main_att_itime`); only the event-study rows collapse.

This is correct behaviour and the author documents it twice — ado:6828-6835 ("the event-study
row pool gives a SINGLE equation name, ES_att_it, to att_it, att_i and att_itime alike") and
sthlp:1315 ("att_it / ame_it for the att_it, att_i, and att_itime estimands (which share a
single event-study horizon internally)"). It is mathematically right: at a fixed relative
period each cohort contributes exactly one cell, so the three equal-weighted estimands coincide
there. **The defect is that sthlp:1305-1306 contradicts sthlp:1315 in the same help file.**

**Consequence, quantified** (`p6_sharpen.do` RUN U, `p6_sharpen.log`). i = a cohort-level
group id; onset-2003 arm = 3 groups totalling 420 units, true effect 0.4; onset-2006 arm =
2 groups of 30, true effect 2.5. Requesting `att att_i` together:

| e(b) name | value | what it is |
|---|---|---|
| `ES_att_post0` | **0.65035128** | cell-weighted; design truth (420*0.4+60*2.5)/480 = 0.6625 |
| `ES_att_it_post0` | **1.23916505** | equal-weighted; design truth (3*0.4+2*2.5)/5 = 1.24 |
| `main_att_i` | 1.23351570 | (equal-weighted, for reference) |

The two differ by a factor of 1.9. A user who asks for `att_i`, looks for the
`ES_att_i_post0` the help promises, does not find it, and takes the `ES_att_*` rows as
"the event study" reads 0.650 where the att_i answer is 1.239 — a 47% error. The correct
row is there under `ES_att_it_post0`, a name the help never associates with `att_i`.
(That last inference step is the user's, so this is a plausible rather than a forced
wrong-number path — hence medium, not high.)

**Realistic:** yes. **Severity: medium** (doc self-contradiction on a headline estimand).
**Confidence: certain.**

---

## Probed and CLEAN (negative results)
- **Downstream tables DO honour the e(b) names.** `e(bootdist)`, `e(boot_tstud)`,
  `e(boot_ptails)`, `e(boot_wald_mainresults)`, the `bootsave()` .dta varnames (incl. the
  `se_` twins), `e(perm_n)`, `e(perm_p_lo)`, `e(perm_comparator_b)` all carry
  `main_SG1_att` / `ES_att_SG1_pre2`. No fourth naming scheme. (p2, p3)
- **Subgroup value is carried verbatim, not re-indexed.** Levels {3,17} give `_SG3_`/`_SG17_`;
  level 0 gives `_SG0_`. (p1 RUN D, p2 RUN L)
- **ame family follows the same rule** on an identified continuous treatment:
  `main_SG1_ame`, `pretreatment_SG1_ame_it`, `ES_ame_SG1_post0`, `ES_ame_it_SG2_pre2`, with
  `e(b_es_ame_s1)` alongside. (p4 RUN P)
- **Negative / non-integer subgroup values** are refused loudly at parse: `sg = -1` -> r(411),
  `sg = 1.5` -> r(411). No invalid-name crash. (p2 RUN J/K)
- **Negative / non-integer subgroup values are refused with wooldid-specific messages**, not a
  bare Stata error: "the subgroups() variable contains negative values on in-sample rows of
  treated units ... must be non-negative integers (or missing) there" and "the subgroup variable
  contains non-integer values within the treated estimation sample". (p2 RUN J/K)
- **The 33-character `pretreatment_SG<8 digits>_att_itime` overflow behaves exactly as
  ado:12283 documents** (`p7_overflow.do`): the length caution prints, then the run dies
  `pretreatment_SG10000001_att_itime invalid name` r(198). No silent truncation, **no collision**
  between two 8-digit values that differ only in the last digit. The 7-digit case (name length
  exactly 32, at Stata's cap) runs clean with correct names. NOTE: p5 RUN S/T did NOT reach this
  path — with no `esprelength()`/`espostlength()` no `pretreatment_*` row is built, so the
  longest name there was only 25 characters. p7 is the run that tests it.
