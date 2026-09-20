# hv-noparent — FINDINGS

Verifier for the recon claim "Eight options are silent no-ops without their parent option".
Pinned build asserted in every probe (`which wooldid` prints
`.../scratchpad/bluesky/build/wooldid.ado`, `*! version 1.4 2026-08-19`).
No repo data, no test-suite file, no `test/**` path was read or run.

---

## VERDICT ON THE SEEDED CLAIM: **TRUE, but the list is wrong in both directions.**

**Confirmed silent no-ops** (base run and orphan run produce a *byte-identical* full `e()` —
scalars, macros, matrices, all in exact hex — and a *byte-identical* isolated log):

| option | parent it needs | evidence |
|---|---|---|
| `saveplots()` | `makeplots` / `makecfxplots` / `histogramcohorteffects` | case 1, and no file written |
| `saveplots()` **with** `makeplots` but no `esprelength()`/`espostlength()` | a horizon request | case 2, and no file written |
| `ccnoabsorb` | `contreat()` | case 3 |
| `contreatpoly()` | `contreat()` | case 4 |
| `suppresssgprimaryeffects` | `subgroups()` | case 5 (only `e(wdid_suppresssgpe)`, the option echo, moves) |
| `cfxplottypes()` | `makecfxplots` | case 6 |
| `cfxplotnose` | `makecfxplots` | case 8 |
| `makecfxplotsbysg` | `makecfxplots` | case 9 |
| `latticeignoreweights` | `lattice(`*type k*`)` — the 2-token form | cases 10 (pre-made lattice var) and 11 (no `lattice()`) |

**Two names on the recon's list are REFUTED as unscreened:**
* `contreatcontrols()` without `contreat()` **is** refused: *"Error: you cannot request continuous
  treatment controls without specifying a continuous treatment variable."* r(198), ado:~5456 family.
* `lattice()` without `contreat()` is refused (ado:4818), `makecfxplots` without `contreat()` is
  refused (ado:5268) — so the cfx orphans are only reachable on a run that already has `contreat()`.

**Two are documented as inert, so the recon's "user gets no message and no documentation" framing
does not hold for them:** `makecfxplotsbysg` (sthlp:976 "silently inert without either: no error and
no caution") and `cfxplotnose` (sthlp:978 "Without makecfxplots the option is inert: nothing reads
it, and no message is printed"). The *behaviour* claim is right; the *undocumented* claim is not.
(Small nit: sthlp:976's "no caution" is itself slightly wrong — `makecfxplotsbysg` tightens the
subgroup-value-length caution from 8 digits to 7 at ado:12294 even with no `makecfxplots`.)

**Named members the recon did not have:** `suppresssgprimaryeffects`, and
`contreatelasticitytype()` — which turns out not to be inert at all (finding 2).

---

## F1 (medium, certain) — the help states a policy the program does not keep

`wooldid.sthlp:1181`, in the `{marker refusepolicy}` section, states:

> "**No accepted option is a silent no-op**, and no check downgrades itself: *safety(off)* is the
> single documented way to move a refusal to a caution or to skip a check, and it moves nothing
> else. It has never suppressed a disclosure about what was estimated - if the run changed your
> specification, it says so whatever *safety()* is set to."

Nine options above are exactly that. The help contradicts itself twice in its own option list
(sthlp:976, 978), and the other seven are documented nowhere. The author has clearly treated this
as a defect class before — there are five screens of exactly this shape already in the code
(`lattice()`→`contreat()` error ado:4818; `makecfxplots`→`contreat()` error ado:5268;
`cfxlattice()`→`makecfxplots` **caution** ado:5307; `esrelativeto()`→`esfixedbaseperiod` **caution**
ado:5371; `poisexpresults`→`poisson` error ado:5479, whose comment reads *"SI-9 (author-ruled
2026-08-25): formerly accepted silently with zero effect"*). The nine above were simply missed.

**Realism ranking, highest first.** `saveplots(`*pfx*`)` typed without `makeplots` is the one that
costs a user something concrete: *"the user believes graphs and a plot dataset were written and
nothing was."* Worse, and more likely, is `makeplots saveplots(`*pfx*`)` **with no
`esprelength()`/`espostlength()`** — the user did remember `makeplots`, and sthlp:807 says no plot
is drawn without a horizon request, but nothing is said at run time and no file appears. Next is
`cfxplottypes()`, whose own validator lives *inside* the `makecfxplots` block (ado:5188), so
`cfxplottypes(psot 1.5 nosuchvariable)` is accepted in silence without the parent and refused
r(198) with it — the wrong way round for a typo screen. `suppresssgprimaryeffects` without
`subgroups()` and `latticeignoreweights` with a pre-made `lattice(`*var*`)` are plausible copy-paste
leftovers. `ccnoabsorb` and `contreatpoly()` without `contreat()` are the least likely.

Repro: `repro_noparent.do` (self-contained; verdict lines start `===VERDICT`).
All ten pairs report `results-identical=1 output-identical=1`; the three controls behave
(garbage `cfxplottypes()` **with** `makecfxplots` → rc 198; `latticeignoreweights` on a genuinely
computed lattice → results move; `makeplots`+horizons+`saveplots(pfx3)` → `pfx3_ES.dta` and
`pfx3_ES_att.gph` appear, while `pfx1`, `pfx2`, `pfx4` write nothing).

---

## F2 (medium, certain) — `contreatelasticitytype()` without `contreat()` is not inert: it RELABELS the table

This is the member the recon listed as "at least one other" and it is the interesting one, because
it is worse than a no-op. On a run with **no continuous treatment at all**, every number is
identical to the last bit (`mreldif = 0.0000e+00` on `e(wdidmainresults)`), but the printed heading
changes:

```
  wooldid y unit t ttre, att robust
    Main Results (Full Estimation Sample):
  wooldid y unit t ttre, att robust contreatelasticitytype(eyex)
    Main Results (Full Estimation Sample): Marginal Treatment Effects Presented as Elasticities
  ... contreatelasticitytype(eydx)
    Main Results (Full Estimation Sample): Marginal Treatment Effects Presented as Semi-elasticities (eydx)
  ... contreatelasticitytype(dyex)
    Main Results (Full Estimation Sample): Marginal Treatment Effects Presented as Semi-elasticities (dyex)
```
and `e(wdid_cetype)` goes `dydx` → `eyex`. The reported `att` is 0.7595 in every one of the four
runs: a plain level ATT, printed under a heading that says it is an elasticity.

Mechanism: ado:11313-11317 (main header), ado:1930-1933 (`Ceffecttype`, used by
`histogramcohorteffects`) and ado:13306-13310 (event-study table and **graph titles** — leg E of the
repro shows it reaching an `esprelength(2) espostlength(2)` run) build the label from
`contreatelasticitytype` with **no `contreat()` gate**, while every computation that would honour
the request is inside `if "`contreat'" != ""`.

The code itself asserts the opposite. **wooldid.ado:5029**, in the BS-6 ordering note:

> "And gate the contreat side on contreat(): `contreatelasticitytype()` without `contreat()` has no
> requires-parent screen and **is inert**, so a caution there would name a dose that does not exist."

The first half is right, the second is false. The BS-6 poisson caution at ado:5041 was gated on
`contreat()` on the strength of that comment, so by code reading (NOT measured here — my poisson
legs died of ppmlhdfe separation on this data, e8/e9 in `p2` and s15/s16 in `p3`) the header should
still say "Elasticities" under `poisson` without `poisexpresults` and without `contreat()`, with the
caution that would explain the scale suppressed. Worth a separate check.

Repro: `repro_cet_mislabel.do`.

---

## F3 (medium, probable) — the elasticity reporting scale has no domain screen, and fails three different ways

Found while working on F2; not in the seeded claim. `d ln E[y]/dx` is undefined where the fitted
outcome is not positive. sthlp:243 warns only that the option "tends to require large samples and a
y-variable that is not sparsely distributed" — nothing about the sign of y, and nothing about the
result being meaningless rather than absent. Measured, on three panels of the same ordinary shape:

| data | weights | result |
|---|---|---|
| `ypos` (fitted values > 0) | none / `aw` | rc 0, semi-elasticity 0.1027 / 0.1049 — sensible (level ATT 0.7906 ÷ mean y ≈ 7.7) |
| `yneg` = `ypos` − 7 (fitted values cross 0) | none | **rc 0, att = 1.327828, se .090441, t 14.68, p 0.000, 95% CI [1.149, 1.506]** — full inference row, **no caution of any kind** |
| `yneg` | `aw` / `pw` | rc 0, att = 1.351847 |
| second panel, min(y) = −3.58 | none / `aw` | r(459), Stata's generic *"could not calculate numerical derivatives — discontinuous region with missing values encountered"* |
| third panel, min(y) = −3.17 | none | r(459) |
| third panel, min(y) = −3.17 | `aw` | **r(3598)**, `variable __00001O not found` / `st_data(): 3598` / `__marg_cr_eval()` — an internal tempvar name in a Mata traceback |

A constant shift of the outcome cannot change the level ATT and does not (0.790571 both times), so
the 1.33 is entirely a reporting-scale artefact of dividing by a mean near zero. Which of the three
behaviours a user gets depends on the data and on whether a weight is present. The same family
covers `contreatelasticitytype(eyex|eydx)`, which the code's own F-collapse text (ado:5074, sthlp:370)
groups with `semielasticity` as "divide by the prediction".

Realistic: outcomes that cross zero are ordinary (net income, profit, a change or a log difference,
a demeaned index). Severity is capped at medium and confidence at *probable* because I have NOT
discriminated the two mechanisms that both fit 1.33: (a) `margins` averaged `1/yhat` over every
treated cell — a complete but meaningless number; or (b) `margins` silently dropped the rows with
non-positive `yhat`, in which case 1.33 is the semi-elasticity of a **shifted subpopulation** and is
a silently wrong number for the requested estimand. Deciding that is the single check that would
move this to `high`: compare the reported figure against a hand-built `margins, eydx` on the same
fit with and without the non-positive rows, and check `e(N)` / the margins subpop count.
What IS established: a full inference row with no caution at all for a quantity undefined on part of
the treated sample, and three different behaviours for one option depending on data and weights.

Repro: `repro_semielas_domain.do`.

---

## PROBED HARD, CAME BACK CLEAN

* **`latticeignoreweights` was a false alarm and the controls caught it.** On two panels it made
  no difference *with* its parent `lattice(byTreated 4)` and aweights, even though the weighted and
  unweighted dose quartiles differed (p75 6.011 vs 7.793; a hand-built `gquantiles ... xtile
  nquantiles(4) by(evertreated)` differed on 65 of 900 rows). On a shape where the lattice actually
  bites — a repeated cross-section with four rows per unit-period cell, so the lattice subdivides
  cells instead of being absorbed by them — the option is **exactly right**: the weighted run is
  byte-identical to passing the hand-computed **weighted** lattice as a pre-made variable, and the
  `latticeignoreweights` run is byte-identical to passing the hand-computed **unweighted** one
  (`p7_latticewt_controls.do`, 0 diff lines on both). On the one-row-per-cell panels the lattice
  terms are collinear with the unit-by-period cell dummies, so nothing could have moved.
  **Do not file this.**
* `simultaneousbands` with no event-study horizons: **not** silent — a named caution fires.
* `histogramcohorteffects` with only `att` requested: works, posts `e(histogramestimates)`.
* `altmethod_adversarial` without `adversarial`, `contreatwithin` without `contreat()`,
  `poisexpresults` without `poisson`, `pretrendpower()` without `esfixedbaseperiod`, and all nine
  `boot*` and ten `perm*` companions without their engine: all properly refused r(198).
* The `e()` comparison was full-surface (every scalar, macro and matrix in `%21x`), so a difference
  anywhere in the returned results would have shown; nothing but `e(cmdline)` (and, for
  `suppresssgprimaryeffects`, its own echo macro `e(wdid_suppresssgpe)`) ever moved.
