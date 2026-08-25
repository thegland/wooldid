# SI-6: `ccc_absorb()` r198 — resolution report (2026-08-25)

Answering the author's challenge on standing issue SI-6 (suite config `cov_ccc_absorb`
errs bare r198 "invalid syntax" in the margins phase on both builds). All runs on the
dev build (`wooldidinstall/wooldid.ado`, loaded explicitly), StataMP 19. Probes and
logs in this directory (`p_batch.do/.log`, `p_trace.do/.log`, `p_verify.do/.log`).

## 1. The author's example: confirmed working

Reconstruction (the author gave the panel shape and dose formula, not the ttre/y DGP;
mine: 30 cohorts x 6 units x 12 periods, onsets 3..9 via `3+mod(i,7)`, one third of
cohorts never-treated, true post effect 1.5, `dose = exp(0.8*runiform()+0.2*mod(i,5))`;
seed 20260825 — `p_batch.do` P1):

- `wooldid y i t ttre, att ccc_absorb(dose)` — rc 0, ATT 1.42 (se .25) vs true 1.5. Sensible.
- `wooldid y i t ttre, att coarse(c.dose)` — rc 0, ATT 1.42 (se .25), point-identical
  to the ccc_absorb call (the absorb-side and interaction-side implementations coincide
  on this DGP — a useful cross-check that ccc_absorb estimated what it should).

The author's claim "that works just fine" is fully reproduced. Note the working example
passes **cohort-level `i`** (6 units per cohort-period cell).

## 2. The failing config: reproduced, bisected, root-caused

Verbatim `cov_ccc_absorb` (dataset `test/suite/data/dta_s_covars.dta`, varlist
`y i t ttre` where `i` is the 150-unit id):
`wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) ccc_absorb(x2)` → r198 bare. Reproduced.

Bisect (all on the dev build):

| probe | change from verbatim | rc |
|---|---|---|
| B1 | drop cluster + ES windows (`att ccc_absorb(x2)`) | 461 "no observations in subpop() subpopulation" |
| B2 | drop ES windows only | 461 |
| B3 | drop cluster only | 198 |
| B4 | drop ccc_absorb only | 0 |
| B5 | clean x2 clone (no injected missings) | 198 |
| B6 | cohort-level id (`coh = 1+mod(i,3)`), all options kept | **0**, ATT 1.30 (cf. B4's 1.28) |
| C1 | author's demo panel, unit-level id, `att ccc_absorb(dose)` | 461 |
| C2 | author's demo panel, cohort id, all failing-config options | 0 |

Minimal failing ingredient: **`ccc_absorb()` + a unit-level (or otherwise
singleton-cell) `i`**. cluster() and the x2 missings are irrelevant; the ES-window
options only select which downstream line dies (198 vs 461).

### Mechanism (verified directly, `p_verify.log`)

`ccc_absorb(x2)` becomes `absorb(... i._i#i._t#c.(x2))` (wooldid.ado L942, regression
at L1141) — a separate x2 slope per treated cohort-by-relative-period cell. With
unit-level `i`, `_i` is the unit id (L520-523), so **every treated `_i`x`_t` cell is a
single observation** (verified: 700 treated rows, max cell size 1). reghdfe drops
one-obs slope groups as singletons: "(dropped 677 singleton observations)" — **the
entire treated sample vanishes from e(sample)** (0 treated rows survive; the 23 others
fell to x2 missings). e(b) then contains only `0o._i#0o._t#0o.tf` and `_cons` — no
treated dummy at any level. With cohort-level `_i` the same command keeps all 677
treated rows and estimates cleanly.

### Exact failing line

Trace (`p_trace.log` L14446-14472): with no treated dummies in e(b), the horizon scan
at wooldid.ado **L1324-1331** matches nothing, `relyrlist` stays empty, and
**L1334 `matrix `relyrcolumn' = (`relyrlist')'` executes as `matrix __00000N = ()'`
→ bare "invalid syntax" r198.** Without ES windows the run instead reaches the first
margins call (L1597) whose `subpop(... targetflag==1 ...)` is empty → r461.

Why no existing screen catches it: the dropped-treated-dummy screen (L1360-1420, r322)
iterates over `o.`-marked coefficients **present** in e(b); here the treated levels are
absent entirely (rows dropped from the sample, not coefficients omitted), so the screen
has nothing to inspect and the wreck sails into the margins phase.

## 3. Verdict: (b) — the author is right on substance; wooldid delivers the refusal badly

The specification is genuinely unestimable **by construction**: ccc_absorb's
cell-specific control slopes require multiple observations per treated cohort-period
cell, and unit-level `i` makes every such cell a singleton. This is not an estimator
bug, and the author's "ccc_absorb isn't suitable for every specification" framing is
correct — his example works, and the same options succeed on the same data once `i` is
cohort-level (B6). But wooldid neither documents nor detects the constraint: it burns
the regression, silently loses 100% of the treated sample, then crashes with a bare
r198 at an internal matrix line (or a raw margins r461). That is an informative-refusal
gap plus one small genuine robustness defect (L1334 can emit `matrix = ()'`).

### Recommended disposition (three pieces, in priority order)

1. **Primary screen (catches the whole class):** post-estimation, right after
   `esamp`/`df_r` are grabbed (~L1170): if no treated observation survived in
   `e(sample)`, refuse informatively (r459, sufficiency-family style; house precedent:
   the bycoh screen L1101-11). Draft:
   > Error: the estimation command dropped every treated observation from the
   > estimation sample (e.g., reghdfe singleton drops), so no treatment effects are
   > estimable. [if ccc_absorb active:] ccc_absorb() estimates a separate control
   > slope for every treated `i'-by-period cell, so each cell needs multiple
   > observations; with a unit-level `i' every treated cell is a singleton and the
   > whole treated sample is dropped. Specify `i' at a coarser level (e.g., cohort or
   > onset group, adding unit fixed effects via fe()), or use coarsecohortcontrols()
   > instead of ccc_absorb().
   This also backstops other estimator-side sample losses (ppmlhdfe separation,
   SI-8-adjacent shapes).
2. **Defensive guard at L1334:** `if "`relyrlist'" == ""` → informative r459 ("none of
   the requested event-study horizons appear in the estimated model") instead of the
   bare-198 `matrix ()'`. Cheap; kills this bare-198 class for good.
3. **Optional pre-flight (best UX):** when ccc_absorb is active, a cheap treated-cell
   size check (gegen count by `_i` `_t` on the screened sample) before the regression,
   refusing/steering before the fit is burned. Help-file note on the cell-size
   requirement either way.

**Suite config:** not (c) user error in intent, but `cov_ccc_absorb` as written asks
for an inherently unidentifiable spec. After the screen lands: keep it as
screen-coverage (rc_old stays 198; rc_new becomes the informative 459) and add a
happy-path companion (prep `gen coh = 1+mod(i,3)`, varlist `y coh t ttre`, same
options, rc 0) — ccc_absorb currently has zero passing coverage anywhere.

## 4. T-M1 adjacency

No interaction. T-M1 is the `evertreated` coding inversion in the
contreatcontrols/lattice blocks plus L1014's pooled-term condition omitting
`ccnoabsorb` — a contreat-seam issue; `ccnoabsorb` (moves contreat controls out of
absorb) is unrelated to `ccc_absorb` (moves coarse cohort controls into absorb), and
SI-6 reproduces with neither contreat nor lattice active (B1/B3). The author's
"ignore" ruling on T-M1 stands unaffected. Only connection: the proposed
treated-sample-survival screen (piece 1) would incidentally backstop any future
absorb-boundary corner that strands treated rows.
