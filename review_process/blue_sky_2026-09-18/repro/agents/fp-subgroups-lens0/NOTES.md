# fpsg-01 — THIRD LENS (independent reproduction). VERDICT: CONFIRMED, high.

Built from the claim alone. Neither prior agent's script was run or read for design.

Repro: `l0_repro.do` / `l0_repro.log` (9 arms), `l0_ereturn.do` / `l0_ereturn.log` (ereturn audit).
Run: `cd <this dir> && /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do l0_repro.do`

## Design — a THREE-way separation, so the mechanism is pinned, not just the direction
45 firms x 12 years (2000-2011), 540 rows. Six treated cohorts, onsets 2004-2009, five firms
each, cohort-constant true effects 6,5,4,3,2,1. Fifteen never-treated.
Exactly ONE firm per cohort carries `sector==1` (ids 1,6,11,16,21,26), so each sector-1 firm's
own effect is separately identified and the six are all different. Verified in-log:

    id   P   Llab  eff   share      (P = post rows, Llab = post rows still carrying the label)
     1   8      1    6   .125
     6   7      2    5   .28571429
    11   6      6    4   1
    16   5      5    3   1
    21   4      4    2   1
    26   3      3    1   1

Item nonresponse: the sector field stops being recorded for firm 1 after 2004 and firm 6 after
2006. **12 blanked rows out of 540 = 2.2%.** Nothing else changes.

Three arithmetically distinct candidate answers for `att_i(sector==1)` in the blanked arm:
  equal-unit (documented meaning) = 3.5000000000
  share-of-own-post-weight (claim) = 2.7611336032
  row-weighted (innocent alternative) = 3.1428571429

## Results

| arm | data | subgroup var | att_i(sg1) | att_itime(sg1) | att / att_it (sg1) |
|---|---|---|---|---|---|
| A | noiseless | clean | **3.5000000000** | 3.5000000000 | 4.0303030303 |
| B | noiseless | 2.2% missing | **2.7611336032** | 2.7611336032 | 3.1428571429 |
| B2 | noiseless | missing, `safety(off)` | 2.7611336032 | — | — |
| C | noiseless | value CHANGES 1->0 | **r(459) refused** | — | — |
| C2 | noiseless | changes, `safety(off)` | **2.7611336032** | 2.7611336032 | — |

Deviations: B vs SHARE = **+4.44e-15**; B vs equal-unit = -7.389e-01 (-21.1%);
B vs row-weighted = -3.817e-01. A vs equal-unit = 8.88e-16.
`att`/`att_it` in arm B match row-weighted to 5.33e-15 — **correct for them**, that is the
documented "belong to no subgroup" semantic. (Does not bear on the skeptic's separate skp-01,
which needs several rows per (i,t) cell; my design has one.)

### Aggregation-only, proved
`e(N)` = 540 in both A and B. `e(N_dropped)` = 0 in both. **`mreldif(e(gteffects)_A, e(gteffects)_B)
= 0.000000000000e+00`** — the estimation stage is bit-identical. Full-sample `att_i` = 3.5 in
both. `att_i(sg0)` = 3.5 in both (no sg0 firm loses its label). The movement is entirely the
`over(sg)` renormalisation.

### The decisive new fact: missing IS the refused case, minus the guard
Arm C (label CHANGES on exactly the same rows) is refused r(459) with the text
"...att_i and att_itime enter a split unit into each subgroup weighted by its share of sample
weight there, rather than with the equal unit weighting those estimands normally imply."
Arm C2 (same, `safety(off)`) returns att_i = 2.7611336032, hex `+1.616cd36ee6177X+001`.
Arm B returns att_i = 2.7611336032, hex `+1.616cd36ee6177X+001`.
**B - C2 = 0.000000000000e+00. Bit-identical.**
So the missing case reaches the same number as the case wooldid refuses by default and warns
about under an explicit opt-in — with neither the refusal nor the caution.

### Caution / disclosure audit
The string "aution" appears **once** in the entire 1830-line log: arm C2's split-unit caution.
Arms B and D1 emit none. Header reads `N = 540 (0 obs dropped from initial sample)`.
A full uncaptured `ereturn list` after arm B (`l0_ereturn.log`) holds nothing recording the
missing-subgroup rows or the lost weight mass: `e(N_dropped)=0`, `e(N)=540`, and no subgroup
scalar or macro beyond the joint tests. `safety(off)` on arm B changes nothing and warns nothing.

### Noisy arm, independent NON-CIRCULAR oracle
Same design + N(0,0.3). Oracle built from single-firm `subgroups()` runs (E0 = label on all
rows, E1 = label on the surviving rows only). With ONE firm in a subgroup, share weighting and
equal-unit weighting coincide **identically** by construction, so these measure each firm's own
post-cell mean without assuming which theory is right and without borrowing `w_att_i` or
`e(histogramestimates)`.

    D0 clean  att_i(sg1) = 3.5547643993   vs equal-unit oracle 3.5547643993   dev  0.00e+00
    D1 missing att_i(sg1) = 2.8252546366  vs share oracle      2.8252546366   dev  8.88e-16
                                          vs equal-unit oracle               dev -7.295e-01 (-20.5%)

### Inferential consequence, truth known to be EXACTLY ZERO
Both sectors contain the same spread of cohort effects (6,5,4,3,2,1), so the true equal-unit
`att_i` contrast sector1 - sector0 is exactly 0. On identical outcome data:

    D0 (clean)   att_i contrast = -0.00908452, se .07956633, t = -0.11, p = 0.909   correct null
    D1 (missing) att_i contrast = -0.73859429, se .08912133, t = -8.29, p = 0.000
                 CI [-0.91393411, -0.56325447] — excludes the true 0 by 8.3 SE

**12 rows (2.2%) of item nonresponse in one classification variable manufacture a subgroup
difference of -0.739 at p < 0.001 where the true difference is zero, with no caution.**

## Doc-vs-code
`sthlp:887` "Observations with missing subgroup values remain allowed and belong to no subgroup"
is exactly right for `att`/`att_it` and materially incomplete for `att_i`/`att_itime`: the unit
is not merely excluded on those rows, it is down-weighted **as a unit** inside its own subgroup.
`sthlp:241` simultaneously names share weighting as the distortion that makes the split-unit
case require an explicit opt-in. The help therefore documents the same distortion as
opt-in-worthy in one place and benign in another.

## Corrections to the two prior agents
- Prober's arm-D oracle (`e(histogramestimates)`) was the wrong comparator, as the skeptic said;
  my single-subgroup oracle is independent of both and lands on share to 8.9e-16. Skeptic right.
- Prober's note that `att`/`att_it` are unaffected holds in a one-row-per-cell design and is
  reproduced here to 5.3e-15. It is a statement about that design, not a general one.
- Neither agent showed B == C2 bit-for-bit. That is the single fact that closes "same
  consequence, different treatment".

## What would have refuted this and did not
Landing on 3.1428571429 (silent row-weighting, a different and lesser defect); any caution or
`e()` record of the lost mass; a changed `e(N)`/`e(gteffects)` making it a sample effect rather
than an aggregation effect; `safety(off)` disclosing it; or B != C2, which would have meant the
missing path is something other than the refused split-unit path.
