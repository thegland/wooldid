# The four critical findings

**Status line, dated 2026-09-18.** All four were reproduced by an agent that did not write the
original probe; two were additionally reproduced by the coordinator. Repro paths are in the review
scratchpad and are self-contained.

---

## CRIT-1 — `impvalues` is not Ibragimov-Müller inference when `cluster()` is finer than, or crosses, `i`

**Verdict chain:** prober CONFIRMED -> independent refuter CONFIRMED (own DGP, own seeds) -> third lens.

`impvalues` needs G *separate* per-cluster estimates whose dispersion carries the sampling variation.
`ado:2042` is **one** `margins ... nose over(cluster)` on **one** pooled fit; `2047-48` grab the G
reads, `2052-53` mean them, `2054-59` t-test them. They are deterministic re-weightings of a single
fit, so their spread is composition x effect heterogeneity.

R=200, all arms on identical draws, analytic arm as positive control:

| `cluster()` geometry | IM 95% coverage | analytic coverage |
|---|---|---|
| **finer than `i`** | **0.380** [.316,.449] | 0.925 |
| **crossing `i`** | **0.250** [.195,.314] | 0.960 |
| equal to `i` | 0.810 | 0.720 |
| coarser than `i` | 0.935 | 0.335 |

**The mechanism is an exact identity.** With K treated cohorts and n units each, the per-unit reads
are exact copies of the per-cohort reads, so `b` is unchanged and `se` is multiplied by
`sqrt((K-1)/(nK-1))`. Measured over 200 draws: `se(unit)/se(cg)` in
[0.431787769588372, 0.431787769588373]; `sqrt(11/59)` = 0.431787769588; max deviation **5.6e-16**.
Refining `cluster()` changes **no data, only labels**, and shrinks the reported SE by 57%.

**Impossibility check.** `wooldid y cg t ttre, att cluster(unit)` on one unit's rows exits **r(459)** —
wooldid cannot estimate from one unit. Yet at `cluster(unit)` it reports 60 "per-cluster estimates"
and t-tests them on 59 df. The honest one-unit SE is 0.2497; the whole-sample IM se is **0.069**.

**What the user sees:** nothing. The degeneracy scrub (`ado:2066-77`) is *relative only*
(`se < 1e-10*|b|`), so 60 reads taking 12 distinct values pass it. `e(df_r)` = 59 is posted for
downstream `test`/`lincom`. When the scrub does fire it misattributes the cause — "usually reflects
deterministic or noiseless data" — on data that is not deterministic.

**Why it is realistic:** sthlp:67/255 are the package's headline guidance and name this exact
geometry. sthlp:257 restricts `impvalues` only to a single cluster variable. sthlp:375 normalises the
crossed case on a measurement of read **accuracy** (1e-12), which is not IM **validity**.

**Boundary, stated fairly:** at and above `i`, `impvalues` is sound and *beats* the analytic arm
(0.935 vs 0.335). That is why the option exists. This is not the ruled T4 item, which concerns
non-`att` estimands; this fires on `att` alone.

---

## CRIT-2 — `contreatcontrols(untreatedZero)` acquires the very control it is defined to exclude

**Verdict chain:** prober CONFIRMED -> independent refuter CONFIRMED on a different design ->
coordinator reproduced on a third design.

The gate at `ado:1234` tests `strpos(cc,"basic") | unconditionalse | customvce() | poisson` with no
`untreatedZero` exclusion; the `untreatedZero` branch at `ado:1313` rewrites `contreatinj` but never
clears `contreatcontrolinserterFE/ABS`. Because the dose is exactly 0 off treated cohort-periods, the
pooled column **equals the sum of the cell#dose columns**, so the engine drops one cell slope and
every reported quantity returns as `truth - (slope of the dropped cell)`.

- Refuter's unseen design: true ame 0.907333333333, broken paths 0.187333333333, implied
  `beta_dropped` 0.7200000000 against a true 0.72 — **1.0e-13 relative**. Dropped cell identified
  three independent ways including wooldid's own `verbose`.
- **`att` and every event-study horizon move too**, not just `ame` (main att 2.118 -> 1.149, -46%).
- **Poisson breaks on its DEFAULT engine, no opt-in.** `i.null` is single-level, so ppmlhdfe's first
  attempt always fails and the rc-3200 fallback carries the term as an explicit regressor.
  Coordinator's own run, exact-fit log-linear, true semi-elasticity +0.25 everywhere, plain
  `att poisson`: reported ame **-0.0284, p = 0.820** — sign-flipped and insignificant.
- **The broken number is exactly `contreatcontrols(basic)`'s answer** — the pairing
  `wooldid_validate` (`ado:4797`) refuses when a user asks for it.

**Not a regression.** `git log -L1234,1236` touches two commits: `fc17051` (2023-08-15, the original
upload) already carried the disjunction; `1df8613` (SI-20) only removed an unrelated qualifier.

**Why it survived:** `review_process/cell_collapse_contreat_2026-09-02.md` read this exact block and
recorded the **opposite** conclusion — that `untreatedZero` "puts the ringer nowhere". The AMESE
adversarial verifier's P4 probed `untreatedZero` on the **default** engine only and marked it OK.

---

## CRIT-3 — `untreatedZero` is silently ignored on every non-reghdfe path (the same mechanism, stated as the user meets it)

**Source:** found by a refuter while verifying CRIT-2. Kept separate because it is the user-facing
statement of the defect and the cleaner test.

Any run that leaves the reghdfe path — `customvce(anything)`, `unconditionalse`, or `poisson` —
returns the **`contreatcontrols(basic)` point estimates**, with no caution, the same `e(N)`, the same
`e(r2)` and the same `e(df_r)` as the correct run, and pre-treatment rows identical to 10 decimals.

On a hospital-year panel with a hand-built saturated oracle implementing `untreatedZero`'s own
documented equation (sthlp:997): default path `main_ame` 0.9127941639 matches the oracle to 10 dp;
`customvce(ols)`, `customvce(robust)`, `customvce(cluster hid)` and `unconditionalse` all return
0.4766147418 — which is exactly what `contreatcontrols(basic)` returns on the default path. Every
post-treatment ES horizon moves with them; every pre-treatment row is byte-identical. **The only
signal that the model changed is a generic multicollinearity caution.**

---

## CRIT-4 — The ruled subgroup-degeneracy guard misses on ordinary ragged panels, and a contrast prints p=0.000 where the truth is 0.25

**Source:** found by a refuter. Independent bootstrap and robust oracles on the same data.

On a staggered panel with holes in it — some `(unit,year)` rows missing, varying units per cluster —
with `i` = state, `cluster(state)` and `subgroups()`, the analytic pairwise-contrast standard error is
**4x to 15x smaller** than the same contrast's bootstrap and robust SEs and than its own realized
sampling SD. `wooldid_sgdeg`'s ruled tolerance of `1 - 1e-2` does not fire, so the row and the
across-subgroup equality test print with **no caution**.

Measured (20 states, 2-5 units each, 10 years, 8% of rows dropped, 640 obs):

| estimator | se | p |
|---|---|---|
| analytic `cluster(state)` | 0.0259 | **0.0000000** |
| `bootstrap(standard)` reps(399), same data and cluster | 0.3105 | **0.2506** |
| `robust` | 0.2734 | 0.2400 |

Ratio **11.995**; implied rho 0.9786; `e(wdid_sgdegenact)` = NONE. The run prints
*"p-value for joint test that all SGs have the same att: 0.0000"*.

Sweep over 8 independent ragged designs (G=20/30, T=8/10, 8% or 15% of rows dropped): **8 of 16
contrast rows blanked correctly, 8 of 16 leaked**, and on every leaked row the analytic SE is
4.1x-14.5x below both oracles. Five of the eight reject at 5%; three print p <= 5e-5.

This is the most likely of the four to be *acted on* by a user, because a subgroup contrast with
p < 0.000001 is exactly the kind of result that goes into a paper's abstract.

---

## CRIT-5 — `contreatcontrols(basic)` returns a sign-flipped answer on the DEFAULT path, no opt-in

**Status: added 2026-09-19.** Found by the untreatedZero blast-radius probe; **confirmed by the
coordinator on an independent exact-fit design.** This is NOT the untreatedZero defect and is not
fixed by any change to the ado:1234 gate — it needs its own decision.

CRIT-2/CRIT-3 framed the contreat control defect as reachable only by leaving the reghdfe path.
That framing is wrong. On data where the dose has no variation outside the treated cohort-periods —
the exact shape `untreatedZero` is documented for — `contreatcontrols(basic)` is silently wrong on
**every** engine including the plain default.

Coordinator verification (`repro/c5b.do`), exact-fit design, R2 = 1.0000, true ame = 1.1233333333:

| `contreatcontrols()` | default reghdfe path, reported ame |
|---|---|
| `untreatedZero` | **1.1233333333** — exactly the truth |
| `basic` | **-0.2266666667** — wrong, and the sign is flipped |

The probe measured the same shape independently (0.373750000 against a truth of 0.973750000, on all
seven engines) and extended it: **`byTime`, `byTreated` and `byTtre` are wrong on every engine too**,
by a DIFFERENT mechanism (the `ccringertail`, which no gate guards), so they are not repaired by the
untreatedZero fix either. `byCohort` is the only member of the family already protected, by
`wooldid_amescreen` mode(bycoh), r459 — which is exactly the screen the proposed FIX-2 generalises.

**Two facts that decide how this should be written up.**
1. ~~**Every pre-treatment row and every event-study pre row is BIT-IDENTICAL between the correct and
   the broken run.**~~ **CORRECTED 2026-09-19 — this is not universal.** It holds on the designs
   measured here, but the FIX-2 fire-rate study found a counterexample: under a relabelled cohort
   identifier with `esprelength(4)`, the pre-treatment `ame` comes back **-0.575000000000 against a
   truth of 0**. So a pre-trend check *can* catch this — on some labellings, which is worse than
   never, because it means the pre-trend evidence is itself label-dependent. The original sentence
   explains why this survived a campaign with extensive pre-trend machinery; it does not license
   "a placebo test cannot detect this".
2. The mitigating reading is that a user running `basic` on dose-zero-off-treatment data is already
   choosing the wrong option — the help documents `untreatedZero` for that shape. But a silently
   sign-flipped point estimate is wrong by any standard, and nothing tells the user.

**Label dependence: my null is WITHDRAWN — the probe was right.** I recorded here that the probe's
factor-ordering claim was "NOT reproduced by the coordinator" and "should not be repeated without its
own repro". It now has one. On a family-D shape with `byTime`, `grp2 = 8 - grp` is a pure relabel —
identical partition, identical rows, identical onsets — and it moves the headline `ame` from
0.752777777778 to **0.508333333333** at `esprelength(0)`, and to **-0.155555555556** at
`esprelength(4)`. My null was on a design where the relabel did not change which cell sorted last;
that is not a refutation.

This matters beyond the record. Whether a given rank deficiency damages `main:ame`, damages the
pre-treatment rows, or damages nothing is **not a function of the data** — it is a function of
Stata's factor sort order. So "per-period false positive" is not a stable category, and a screen
cannot be narrowed to "only the deficiencies that do damage". The run-level statement is the only
defensible one, and it is clean: every run whose screen fired was wrong somewhere, and no run whose
screen stayed silent was wrong anywhere.

**Proposed remedy: FIX-2, now specified and measured — see `fixprop/FIX2_FINAL.md`.** Generalise
`wooldid_amescreen` mode(bycoh) into a reference-pool identification screen covering
basic/byTime/byTreated/byTtre. The open question was the `byTime` fire rate on ordinary staggered
panels; it has been answered, and it reframes the question. The screen fires **structurally** — 100%
of draws on the families where it fires at all — but with **0 run-level false positives and 0 false
negatives** across 48 (dose-shape x type) cells at `e(r2) == 1` plus 16 boundary designs, and the
predicate's count of dangerous groups equals wooldid's own dropped-coefficient count exactly in every
case checked. `basic`/`byTreated`/`byTtre` fire 0 times in 17 of 20 panel families. And an
unconditional `r(459)` on ordinary data is already shipped precedent: `byCohort`'s screen at
`ado:1383` has no `safety()` gate and fires today on a plain panel under `esprelength(4)`.

One change in `FIX2_FINAL.md` is a ship gate: the screens **must** fire under `contreatwithin`,
reversing the proposal. Measured, `contreatcontrols(basic) contreatwithin` on a within-treated-
variation design returns **-0.238888889 against a truth of 1.411111111 — a sign flip — with zero
cautions printed**, and `byCohort`'s shipped screen already fires there.
