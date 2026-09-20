# hv-noparent-refute — working notes (SKEPTIC)

Build asserted in every log: `which wooldid` resolves to
`/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`. All data generated in-probe.
The project test suite was never touched.

## Probe inventory (all in this directory, each with its .log)

| file | what it settles |
|---|---|
| `r_f2_cet.do` | F2: 4 paired runs + ES + subgroups + histogram + poisson leg |
| `r_f1_noparent.do` | F1: 9 paired orphan-option legs, file-write test, cfxplottypes validation gap |
| `r_f3_semielas.do` | F3: scale-invariance + additive-shift legs (weights mis-typed — see below) |
| `r_f3b_semantics.do` | F3: standalone `margins, eydx()` oracle; weight syntax fixed; F1 loose ends |
| `r_f3c_pole.do` | F3: the POLE sweep — the decisive run |
| `r_f3d_scope.do` | F3: event-study + att_it legs on a straddling panel |
| `r_new_validgate.do` | my own unseeded sweep of value-domain validation gated behind a parent |

**Harness note for a successor:** in `r_f3_semielas.do` I wrote `wooldid ... , att robust semielasticity [aw=w]`.
Stata weights must precede the comma; that produced spurious `r(198)` on the Q2/Q3/Q4/R2 legs.
Fixed in `r_f3b_semantics.do` (`wooldid y unit t ttre [aw=w], ...`) — all those legs are rc 0.
Also: `capture dis` SUPPRESSES the display. Use `capture noisily dis`.

## F1 — PARTIAL

Behaviour half reproduces exactly (`r_f1_noparent.log`): nine legs, all rc 0,
`mreldif(e(wdidmainresults), baseline) = 0` for every one, nothing printed beyond the
ordinary table, `saveplots(p1)`/`makeplots saveplots(p2)` wrote **no** files while the
control `makeplots esprelength(2) espostlength(2) saveplots(p3)` wrote both `p3_ES.dta`
and `p3_ES_att.gph`.

Doc half largely fails. **Five** of the nine are explicitly documented as inert and
silent, in the very words the finding says the help forbids:

- `ccnoabsorb` — sthlp:1121 "the option is inert and silent on a run with no continuous treatment"
- `contreatpoly()` — sthlp:917 "without `contreat()` the option is otherwise inert: nothing reads it and nothing is said about it"
- `cfxplotnose` — sthlp:978 "Without `makecfxplots` the option is inert: nothing reads it, and no message is printed"
- `makecfxplotsbysg` — sthlp:976 "silently inert without either: no error and no caution"
- `makeplots` with no horizon request — sthlp:807 "No plot is drawn unless at least one event-study horizon is requested"

The finding conceded only two of these. With five explicit entries in the same file,
sthlp:1181's "No accepted option is a silent no-op" cannot be read as the promise the
finding needs; its own closing clause scopes it ("if the run **changed your
specification**, it says so"), and an unreached child option changed nothing.

The one item in F1 that survives at full strength is **not a no-op claim at all**:
`cfxplottypes()`'s token validator sits inside `if "makecfxplots" != ""` (ado:5188).
Clean paired proof, same data, same options:
- `... contreat(dose) contreatcontrols(untreatedZero) cfxlattice(dose0) cfxplottypes(psot) safety(off)` → **rc 0**
- `... + makecfxplots` → **r(198)** "cfxplottypes() contains an unrecognized entry: psot"
sthlp:974 states that rejection with **no** parent qualifier, and the sibling entry
sthlp:917 says contreatpoly() "is validated whether or not contreat() is given".
ado:5188's own VAL-G13A comment records that garbage tokens used to be swallowed; the
fix was placed inside the parent block, so without the parent the pre-fix behaviour is live.

## F2 — CONFIRMED, with three surface corrections

Core measured (`r_f2_cet.log`): four runs on identical data, att hex
`+1.d169f113af0c8X-001` in **all four**, `mreldif(e(wdidmainresults))` = 0 for each,
header moves to "Marginal Treatment Effects Presented as Elasticities /
Semi-elasticities (eydx) / (dyex)", `e(wdid_cetype)` `dydx` → the request.

1. **The ES graph titles are NOT relabeled.** Both saved `.gph` files carry
   `ytitle("ATT Effects on y by Horizon")`. ado:13308-13310 sets `elasinj2`, which feeds
   only `ame_it_text`/`ame_text` — the AME graphs, which do not exist without contreat().
   The ATT graph title comes from `elasinj` at 13305-06, driven by `semielasticity` only.
   The ES *table* is relabeled only because its rows sit under the same Main Results
   header (ado:11368).
2. **`histogramcohorteffects` is NOT relabeled.** `Ceffecttype` (ado:1930-33) is consumed
   only at ado:12391/12428 and 2506-2602, every one gated on `if "contreat" != ""`.
3. **Worse than stated:** the relabel also reaches the per-subgroup headers (ado:11528)
   and the pairwise-comparison header (ado:11588) —
   `Results for Subgroup sgv == 0: Marginal Treatment Effects Presented as Elasticities`,
   `Pairwise Comparisons of Different Subgroups' Estimates: ... Elasticities`.

The finding's untested poisson prediction is **correct**: `poisson` without
`poisexpresults`, no contreat(), `contreatelasticitytype(eyex)` →
`Main Poisson Regression Results (Full Estimation Sample): Marginal Treatment Effects
Presented as Elasticities`, **no BS-6 caution**, att `.28180154` identical to the plain
poisson run, `e(wdid_cetyperaw)` not posted.

## F3 — CONFIRMED, two mechanism corrections, open question answered

`r_f3c_pole.log` is the decisive run. One panel, one model, one additive shift `c`
applied to y. Level ATT = `.7543601` at **every** `c` (invariance control):

| c | min fitted-ish y over treated-post | semielasticity result |
|---|---|---|
| 0 | +7.655 | .08617444 (se .0053) rc 0 |
| -6 | +1.655 | .27606125 (se .0171) rc 0 |
| **-8** | **-0.345** | **1.1572378 (se .0864, t 13.39, p 0.000) rc 0, NO domain caution** |
| -8.4 | -0.745 | **r(3598)** `variable __00001O not found` / `__marg_cr_eval()` |
| -8.7, -9, -10, -20 | deeper | r(3598) |

Adding `esprelength(2) espostlength(2) jointtests` to the c=-8 run turns rc 0 into
**r(459)** "could not calculate numerical derivatives — discontinuous region"
(`r_f3d_scope.log`). `att_it` on the same straddling data: rc 0, `1.2005979` (se .0817).

**Correction 1 — the r(3598) is not the weighted leg.** Every r(3598) above is
*unweighted*. With the weight syntax corrected, `[aw=w]` and `[pw=w]` legs on a
straddling panel return rc 0 (`r_f3b_semantics.log`, QB2/QB3). The discriminator is how
far the fitted values run past zero, not the presence of a weight.

**Correction 2 — it is knife-edge.** A 0.4 change in an arbitrary additive constant,
same model, same option, flips a confident complete inference row into a raw Mata
traceback.

**Open question answered — it is (a), not (b).** Standalone oracle, no wooldid
(`r_f3b_semantics.log`): hand-computed mean of `_b[x]/yhat` over the subpop =
`1.0868654`; `margins, eydx(x) subpop(...)` returned `1.0868654` to every printed digit;
the ratio-of-means alternative was `1.0039811`. So `margins` eydx is the **mean of
per-observation ratios** and drops nothing. The number wooldid prints on a straddling
panel is a complete average over terms of opposite sign, exploding near yhat = 0.

**Not established (by them or by me):** the `contreatelasticitytype(eyex|eydx)` legs of
the same family. My contreat designs died on wooldid's multicollinearity screen r(322)
(`r_f3d_scope.log`, C legs). The family-membership claim remains an assertion.

## My own unseeded sweep (`r_new_validgate.log`)

The generalizable class is "value-domain validation gated behind a parent", not "inert
options". Sweeping every other valued option:

- `esrelativeto(3)` — a **positive** value, where the help documents `(-int)` and the
  gated error message says "must be a negative integer" — is accepted **rc 0** without
  `esfixedbaseperiod` and refused **r(198)** with it. BUT wooldid prints
  "Caution: esrelativeto() has no effect without esfixedbaseperiod and is being ignored",
  so the inertness *is* announced. Low severity; it is the announcement that saves it.
- `pretrendpower()` bad values: r(198) with and without a window. Screened.
- `customweightmultiplier()` with negatives: r(411) (matches the ruled VAL-G06).
- `cfxlattice()` without `makecfxplots`: my leg died earlier on the contreatcontrols
  screen (r 198), so the ruled caution was not reached — untested here.

Net: `cfxplottypes()` is the only valued option I found that is **both** unvalidated and
unannounced when its consumer is not reached. That is the honest scope of F1's live part.

---

## CORRIGENDUM (after the corrected oracle, r_f3e_oracle.do) — read this over the F3 section above

Two claims of MY OWN in the F3 section above are wrong and are superseded:

1. **"The discriminator is how far the fitted values run past zero" is refuted by my own
   logs.** Panel Q reaches roughly -2.7 in fitted value and returns rc 0 with a number;
   pole shift c = -8.4 (min y over treated-post only -0.745) dies r(3598). Depth is not
   the rule. What the crashing configurations share is fitted values tightly CLUSTERED
   near zero rather than spread across it — a hypothesis, not a measured rule. The honest
   statement is that the crash-vs-number switch is data-dependent and uncharacterized by
   either party.

2. **The r_f3b oracle never straddled zero** (min yhat = .56): `u` was per-observation
   noise, which `reg y x i.g` averaged out within (g,x) cells. Redone in `r_f3e_oracle.do`
   with a GROUP-level `u`: min yhat -1.990, max +4.013, 70 of 200 subpop rows negative.
     hand mean(_b[x]/yhat) = -.24488522 ; ratio-of-means = 1.00068
     margins, eydx(x) subpop(if x==1) = -.24488522, rc 0    <- exact match
   So (a) is confirmed ON THE ACTUAL DOMAIN IN QUESTION: mean of per-observation ratios,
   nothing dropped. And the demonstration is sharper than the finding's: a treatment with
   a true effect of +1.0 on an outcome whose fitted values average +1.0 is reported as a
   semi-elasticity of **-0.24**. Wrong sign, meaningless magnitude, rc 0, no message.

3. **Primary evidence for the silent-number leg is panel Q (r_f3b_semantics.log, QB1:
   rc 0, 1.0437756), not pole c = -8.** At c = -8 min y over treated-post is -0.345
   against a residual sd of 0.3, so I did not establish that its fitted values cross zero
   — the same inference I flagged as unreliable when reading the prober's panel A. The
   pole run remains the evidence for the knife-edge and for the unweighted r(3598).

4. I did not independently check `e(wdid_suppresssgpe)`; I accept the prober's observation
   that it is the one echo macro that moves, and my silence claim is about printed output
   and the results matrices only.
