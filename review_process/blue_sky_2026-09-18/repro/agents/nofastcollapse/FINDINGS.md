# nofastcollapse — differential-testing findings

Agent: NFC.  Method: differential testing (arm A = collapse default, arm B = `nofastcollapse`),
same pinned build, same data re-`use`d fresh for each arm, every `e()` scalar/macro/matrix dumped
in %21x hex and compared elementwise (dimensions, row/column names and presence compared
separately, so a horizon present in one arm only would show even with no numeric diff).
Harness noise floor measured WITHOUT controls (same arm twice, both arms): EXACTLY 0.
The same check WITH the control present is in pcrit.do section (1) - see NOTES.md for whether
it had run by the time this file was written.

---

## NFC-01 (found by me, NOT seeded - `controls()` is absent from the seeded cell list) — `controls()` with a continuous covariate puts the fast and slow collapse paths 4-5 orders past the project's own tolerance, and the disagreement reaches the PRINTED digits

**Surface**: `nofastcollapse` x `controls()` (binary treatment, ordinary linear path).

> **READ THE "LATE RESULTS (pcrit)" SECTION AT THE BOTTOM FIRST.**  It settles determinism
> (exactly 0), and then the CORRECTION section at the very bottom retracts an over-reach:
> sthlp:357's "a row permutation moves it ~19x more than the collapse" is about a COEFFICIENT
> and is CORRECT as written (measured 15.5).  What survives is a scope gap - coefficients are
> fine everywhere (1e-16..1e-15), variances are not (1e-7..5e-6) - plus the new arm-G result
> that the collapse makes the reported VARIANCE row-order dependent (6.6e-16 -> 8.7e-8) where
> it is not otherwise.  It also records the two caveats that keep severity at
> `medium`: no decision ever flipped, and the effect is DESIGN-dependent (it does not fire on a
> simple hand-built panel, and it does not grow with rows-per-cell).

**Claim**: On a plain staggered panel with one ordinary continuous covariate in `controls()`,
the default cell collapse and `nofastcollapse` disagree by up to **5.5e-6 relative on `e(V)`**
(a variance, i.e. ~2.7e-6 on a standard error), up to **3.5e-6 on `e(wdidmainresults)`**, and the
disagreement changes digits that wooldid **prints**.  Point estimates are unaffected (~1e-15).
The project's indifference threshold is rtol 1e-10.

**Contradicts four documented, MEASURED statements in wooldid.sthlp** (the `nofastcollapse` entry;
line numbers verified on the pinned build):
 1. `sthlp:332-333` "The difference is bounded, in **absolute** terms, by about one unit in the
    last place of the linear prediction ... and every difference measured was a small integer
    multiple of that one quantum. On a panel whose predictions are of order 10, that is about
    2e-15."  Measured here: the moves are ~1e8-1e9 ulps of the prediction, not a small integer
    multiple of one.
 2. `sthlp:338` "Standard errors and variances **do not amplify at all: they were measured at
    about 1e-16 relative** regardless of the scale of the predictions."  Measured here: **5.5e-6**
    on a diagonal of `e(V)`.
 3. `sthlp:350` "For almost all work, not at all: **no digit wooldid displays is affected**."
    (repeated at `sthlp:347` for the continuous case).  Measured here: displayed standard errors,
    displayed CI bounds and displayed event-study point estimates all change in their last
    printed digit.
 4. `sthlp:344` (continuous treatment) "Measured over 305 paired continuous-treatment runs ... **no
    cell moved by more than 1e-9** in relative terms."  Cell w09 — `contreat(dose)
    contreatcontrols(basic) controls(xdollar)` — moves **8.75e-7** on `e(wdidmainresults)`.  Every
    continuous-treatment cell I ran WITHOUT a continuous `controls()` term was clean at <=1e-13,
    so the 305-run measurement almost certainly did not include one.

**Observed** (probe p07, data/dCc.dta, 150 units in 3 treated cohorts + never-treated, T=10,
i = cohort; `att att_i att_it esprelength(3) espostlength(3) jointtests`):

| cell | control | max rel diff, `e(V)` | max rel, `e(wdidmainresults)` | rel diff `e(joint_horiz)` |
|---|---|---|---|---|
| w01_c1      | `controls(x1)`         (sd 1)      | **1.81e-07** | 1.14e-07 | - |
| w05_pct     | `controls(xpct)`       (mean 50)   | **5.49e-06** | 3.46e-06 | 8.3e-10 |
| w03_dollar  | `controls(xdollar)`    (mean 52000)| **2.01e-06** | 2.60e-06 | 3.9e-10 |
| w04_dollar2 | two dollar-scale controls          | 7.68e-08 | 5.33e-08 | 1.7e-09 |
| w10_rob     | `controls(xdollar) robust`         | 2.01e-06 | 1.33e-06 | 5.4e-10 |
| w11_clu     | `controls(xdollar) cluster(region)`| 2.56e-06 | 2.32e-06 | **2.5e-08** |
| w09_ct      | `contreat(dose) controls(xdollar)` | 2.82e-07 | 8.75e-07 | - |
| c09_ctrl    | `controls(x1 x2)`                  | 5.34e-09 | 2.82e-09 | 2.9e-10 |

The printed table, two continuous controls (p07.log, blocks DISPLAY_A / DISPLAY_B):

```
                A (collapse, the default)          B (nofastcollapse)
 main  att   | -.07134316  .11903432 ... -.30483918   |  ... -.30483917    <- CI bound
 ES_att pre1 | -.10003275  .16194921 ...    -.41771   |  -.10003275  .1619492 ... -.41770999
                          ^^^^^^^^^ printed se        |            ^^^^^^^^ printed se
 ES_att att  |  .16215286                             |   .16215285         <- printed estimate
 ES_att pre1 |   .2176445                             |   .21764449         <- printed estimate
```

**Mechanism**.  What is *established* is the three-arm isolation below.  The numerical-analysis
story underneath it (see the last paragraph) is *consistent with* the evidence but is not itself
measured here; `dx01.do` in the workdir is the diagnostic that would settle it.
| arm | result |
|---|---|
| `controls(xcohortc)` — the SAME control made constant within the (cohort, period) cell | **3.58e-14 CLEAN** |
| `fe(region)` — absorbed, so it carries no coefficient | **3.60e-14 CLEAN** |
| `controls(i.region)` — a FACTOR control | **5.04e-14 CLEAN** |
| no controls at all | 2.9e-14 CLEAN |

So the trigger is precisely: *a continuous control that carries a coefficient and varies within the
collapse's own (cohort, relative-period, targetflag) cell.* The collapsed row carries ONE row's
value of that control. The gate's whole argument is that such a control cancels between
xb(targetflag=1) and xb(targetflag=0) — which it does for the POINT estimate (1e-15 here).  It
does not cancel in `margins`' delta-method Jacobian: the Jacobian entry for a control coefficient
is analytically zero but is formed as a difference of two nearly equal predictions, so what
survives is cancellation residue that depends on WHICH row carries the control value, and
`J V J'` then squares and cross-multiplies it.

**Why this is the cell nobody measured**: the ado's gate comment (wooldid.ado:5588-5610) reasons
only about the derivative of the *prediction*, never about the Jacobian of the *variance*, and the
sthlp's continuous-treatment paragraph reports 305 paired runs but its binary paragraph reports the
1-ulp bound — the combination `binary + continuous controls()` is exactly where the two arguments
meet and neither covers it.

**Scope caveat, measured**: the magnitude is DATASET-dependent, as cancellation residue is.  The
same option on a different panel (`m10_messy`, data/dCm.dta, `controls(x1 x2)`) came back at
2.85e-13.  It is *reproducible per dataset* — the harness noise floor is exactly 0 — not run-to-run
noise.  So the claim is "this can and does happen on ordinary input", not "it always happens".

**realistic**: TRUE. `controls(x)` on a continuous covariate is the commonest option in applied DiD.
**severity**: medium  **confidence**: certain (see the four-line argument at the end of this entry)

**PRE-REGISTERED upgrade criteria** (written before the remaining probes returned).  Upgrade to
`high` if any of these lands:
 - p10 (rows-per-cell ladder, up to 512 rows per cell — the help's own administrative scale is
   ~500) shows the disagreement growing into the 4th-5th displayed digit;
 - p08 x05-x07 or p11 show a **reported set** changing rather than a low digit: a horizon dropped
   in one arm only, `e(wdid_fcguard)` firing in arm A only, a `permexact` p jumping by 1/n, or the
   Ibragimov-Muller degeneracy guard scrubbing SEs in one arm only.
Downgrade to `low` if p12's row-shuffle arm shows a plain row permutation (collapse OFF in both
runs) moves `e(V)` by the same ~1e-6 — that would make this a scope error in the help's measured
bounds rather than anything the collapse does.
**repro**: `<WORKDIR>/repro/repro_NFC01_controls_V.do` (self-contained: builds its own 1200-row
panel, prints a noise floor of exactly 0 with the control present, then the four mechanism arms,
then the two printed tables to diff).

---

## Verified CLEAN (negative results with the measured floor)
Listed in NOTES.md; every cell below matched at <= 1e-12 with the harness noise floor at exactly 0.

### Cells verified CLEAN against `nofastcollapse` (worst relative difference over EVERY posted
### e() scalar and every element of every e() matrix; harness noise floor exactly 0)

Binary linear path, i = cohort, 4 estimands + event study + jointtests (data/dC.dta):
plain 2.9e-14 | robust 3.6e-14 | cluster(region) 2.3e-14 | [aw=w] 1.5e-14 | [pw=w] 7.8e-15 |
subgroups(sgc) 1.4e-14 | esfixedbaseperiod+adversarial 2.2e-14 | esrelativeto(-2) 7.2e-15 |
timetrends 1.5e-14 | fe(unit) 7.2e-14 | anticipation(1)+oostreatedcontrols(2) 5.5e-15 |
fastnose 3.5e-14 | histogramcohorteffects+summarystats+gteffects 2.9e-14 | an `if` restriction
1.4e-14 | controls(i.region) 5.0e-14 | fe(region) 3.6e-14.

`impvalues` — the case the sthlp says left the ineligible list on 2026-09-06, and the one place the
cluster variable joins the cell key: cluster CROSSING i 3.0e-15, cluster NESTING i 4.3e-16.
CONFIRMS the sthlp's "measured at 1e-12 both for a cluster nesting i and for one crossed with it".

Continuous treatment (data/dC.dta, dose varying within the cell):
contreatcontrols(basic) 6.3e-15 | +byTime 5.6e-13 | +byTtre 5.4e-15 | ccnoabsorb 8.8e-15 |
contreatelasticitytype(dyex) 5.9e-15 | contreat+[aw=w] 2.6e-14 | contreat+customweightmultiplier
with exact zeros 8.2e-14 | contreat+subgroups 1.6e-14.

Poisson (the BS-6/F032 remap path, where the gate lifts the ban):
poisson 2.4e-14 | poisson+contreat 8.9e-15 | poisson+contreat+dyex 8.9e-15 |
poisson+semielasticity 1.2e-14 — and `e(wdid_fcol)` is 1 in all four, so the remap really does
turn the collapse back on, as sthlp:249 claims.

`customweightmultiplier()` with a multiplier that is exactly 0 on 10% of units (the case the ado's
own comment at 5613-5645 says broke an earlier form of the dose mean): 1.5e-14 binary,
8.2e-14 with contreat, 1.7e-13 on the messy panel.  The `!= 0` guard holds.

The `wooldid_fcdeg` degeneracy screen, live: `subgroups(sgc) cluster(cohort)` FIRED the screen
(`e(wdid_fcguard)` == "fallback" in arm A) and the fallback then reproduced the `nofastcollapse`
arm to **4.58e-16** across every result.  The screen and its restore-and-re-issue machinery work.

Messy administrative shape (entry, exit, interior gaps, 5% missing y, 5% missing x1):
4 estimands 2.85e-13 | contreat 2.2e-14 | cwm-with-zeros 1.7e-13.

Repeated cross-sections (4 rows per unit-period, i = unit): **EXACTLY 0** on the 4-estimand run
and on subgroups+cluster(region) — the collapse is bit-identical there.

The 2026-09-16 option consolidation is CLEAN in code: no surviving reference to `fastcollapse` or
`fullprecisionslow` anywhere in wooldid.ado outside explanatory comments (only `NOFASTCOLLAPSE` in
either `syntax` line), and nothing in wooldid.sthlp.  All eight `restore, preserve` sites in the
fcdeg fallback reset BOTH `local fcnoes` and `local fcsg` on the next two lines, so the fallback
cannot leave a later `margins` running `noesample` against full data.

No data-state leak in any successful run: `datasignature` before and after is identical in both
arms of all 41 cells, and `_N` and `c(k)` are unchanged.

### Amplification probes — does the residue ever change a DECISION?  (so far: no)
- `permute(shuffle) permexact cluster(cohort) controls(xdollar)`: the whole permutation
  distribution `e(permdist_b)` and every p match at **4.26e-14**.  A rank-based p does not move.
- `bootstrap(standard) reps(200) rseed(1) cluster(region) controls(xdollar)`: every bootstrap
  b/se/ci/p matches; only the RETAINED ANALYTIC `e(V_analytic)` moves (2.32e-07) — NFC-01 again.
- `impvalues cluster(region) controls(xdollar)`: **1.85e-15 CLEAN.**  impvalues sets `nose`, so
  there is no delta-method VCE for the residue to live in.  This is a strong corroboration that
  the residue enters through the delta-method Jacobian rather than through the estimates.
- **No caution or error text ever differed between the arms.**  Across the 5 probe logs
  (p01,p02,p03,p07,p11 - 46 paired runs), every `Caution:`/`Error:` line appears an EVEN number of
  times, i.e. in both arms.  The collapse never changed WHAT wooldid says, only low digits of
  what it computes.  Structural objects (`e(escomp)`, matrix dimensions, row/column names,
  the reported horizon set) were identical in every cell.

### What the author would want to know about NFC-01
1. The `wooldid_fcol` gate's argument (ado:5588-5610) is entirely about the *derivative of the
   prediction*: "every regressor that does not interact with targetflag contributes identically
   to both terms and cancels".  That is true, and the point estimates confirm it (1e-15).  The
   argument is never extended to `margins`' delta-method *Jacobian*, where the same cancellation
   is between two nearly-equal quantities and leaves residue instead of zero.  Every measured
   failure is in a V-derived quantity; every `nose` path (`impvalues`, `fastnose` with no ES
   variance) is clean.
2. The collapse ALREADY carries a cell mean for exactly this reason in the `contreat()` case
   (wooldid_fcol's `contreat()` option averages the dose and the ringer onto the retained row).
   Continuous `controls()` terms get no such treatment - the retained row keeps one arbitrary
   row's value.  Whether averaging them would close the gap is not obvious (the uncollapsed
   Jacobian residue is an average of per-row residues, not the residue of the average), so the
   cheap, certainly-correct remedies are either to put a continuous `controls()` term that varies
   within the cell on the ineligible list, or to widen the sthlp's measured bounds to say the
   1-ulp / 1e-16 / no-displayed-digit figures were measured WITHOUT continuous controls.
3. The measured evidence in the sthlp looks scoped exactly this way: the binary paragraph's
   1-ulp bound and the continuous paragraph's "305 paired runs, no cell above 1e-9" both hold in
   every cell I ran that had no continuous `controls()` term, and both fail as soon as one is
   added (w09 is a continuous-treatment run at 8.75e-7).

### NFC-02 (low) — sthlp: "The full data is restored before any plotting or display work"
`wooldid_hist` (`histogramcohorteffects`) is invoked at ado:1944, which is INSIDE the
`preserve` ... `restore` window opened at ado:1876 and closed at ado:3496; the cfx `marginsplot`
calls at ado:2492-2625 are inside it too.  No measured consequence: the histogram is built in its
own frame from a margins result matrix (so it does not read the collapsed rows), `m07_hist`
matched at 2.9e-14 including `e(histogramestimates)`, and the cfx sites are unreachable under the
collapse because `contreat()` + `makecfxplots` is on the ineligible list.  Cosmetic doc-vs-code.
**realistic**: true but inconsequential.  **severity**: low.  **confidence**: certain.

### NFC-01, the full cell-by-cell picture on ONE run (w05_pct: `controls(xpct)`, xpct = 50 + 8*x1,
### i.e. an ordinary percentage-scale covariate).  `e(wdidmainresults)`, every cell above rtol 1e-10.
wooldid prints estimates, standard errors and CI bounds to 8 decimal places, so these ARE the
printed digits.

```
 family        row      column    rel diff   A (collapse, default)   B (nofastcollapse)
 pretreatment  att      se        2.41e-06   0.125623208243          0.125622905190   <- 7th dp
 pretreatment  att      lb_95ci   1.99e-06  -0.298298927330         -0.298298332867   <- 6th dp
 pretreatment  att      ub_95ci   3.06e-06   0.194542173605          0.194541579143   <- 6th dp
 pretreatment  att_it   se        2.74e-06   0.125301806901          0.125302150605   <- 7th dp
 pretreatment  att_it   ub_95ci   3.46e-06   0.194664248342          0.194664922547   <- 6th dp
 pretreatment  att_i    (identical to att_it)
 ES_att        pre2     se        7.30e-08   0.159360561954          0.159360550314
 ES_att        pre2     ub_95ci   8.40e-08   0.271659161721          0.271659138888
 ES_att        post0    p         2.84e-06   1.38632744935e-24       1.38632351379e-24
 ES_att_it     post0    p         2.81e-06   7.17329739260e-25       7.17331758487e-25
```
The POOLED pre-treatment row is the worst hit (2.4e-6 to 3.5e-6); the individual event-study
horizons move at 2e-8 to 8e-8; the far-tail p-values at 2.8e-6 (sthlp:338 offers 3e-11 on a
p of 2e-07 as the far-tail figure).  The t and p columns of the main rows are printed to 2 and 3
decimals so they are unchanged on screen; the se and the two CI bounds are not.

### NFC-01 frequency, honestly stated
Every paired cell I ran that contained a continuous `controls()` term, worst relative difference
over all of `e()`:

```
 5.49e-06  w05_pct      dCc  controls(xpct)                    FAIL
 2.56e-06  w11_clu      dCc  controls(xdollar) cluster(region) FAIL
 2.01e-06  w03_dollar   dCc  controls(xdollar)                 FAIL
 2.01e-06  w10_rob      dCc  controls(xdollar) robust          FAIL
 8.75e-07  w09_ct       dCc  contreat + controls(xdollar)      FAIL
 2.32e-07  v03_boot     dCc  bootstrap + controls(xdollar)     FAIL (analytic V only)
 1.81e-07  w01_c1       dCc  controls(x1)                      FAIL
 7.68e-08  w04_dollar2  dCc  controls(xdollar xdollar2)        FAIL
 5.34e-09  c09_ctrl     dC   controls(x1 x2)                   FAIL
 2.85e-13  m10_messy    dCm  controls(x1 x2)                   pass
 2.18e-14  m11_messy_ct dCm  contreat + controls(x1)           pass
 1.85e-15  v04_imp      dCc  impvalues + controls(xdollar)     pass (impvalues runs `nose`)
```
9 of 12 exceed the project's rtol 1e-10; one of the three passes is structural (no delta-method
VCE exists under `impvalues`' `nose`), two are dataset-specific.  So: this fires on most ordinary
datasets, not all, and when it fires it is 1e-9 to 5e-6.  Every cell with NO continuous
`controls()` term (32 of them, listed above) passed at <= 5.6e-13.

### NFC-01: why this cannot be run-to-run noise (four independent lines)
1. `n0` (nofastcollapse twice) and `n1` (collapse twice), both with the full 4-estimand +
   event-study + jointtests request, are **EXACTLY 0** across every posted result.  The margins
   stage is deterministic.
2. In every FAILING cell the REGRESSION-STAGE scalars are **bit-identical** between the arms.
   For w05_pct: `F N N_dropped df_r escomp_* level ll ll_0 numpre_att numpre_att_it r2 r2_within
   r2_withinadj r2adj rmse rss wdid_vblockn` all match to the last bit; only `joint_horiz`,
   `joint_horiz_i`, `e(V)` and `e(wdidmainresults)` move.  So reghdfe ran identically and the
   divergence is confined to the `margins` stage - which is exactly the stage the collapse
   touches, and nothing else.
3. The three mechanism arms (control made cell-constant, control absorbed by `fe()`, control made
   a factor) all return to ~1e-14.  A reproducible, systematic dependence on a structural property
   of the control is not what noise looks like.
4. Magnitudes are coherent across related cells (xdollar 2.01e-06 plain and 2.01e-06 under
   `robust`, 2.56e-06 under `cluster(region)`).

**Confidence upgraded to `certain` on the strength of (1)-(4).**

### THE WHOLE SWEEP, one line per cell (worst relative difference over every posted e() object)
No cell anywhere showed a STRUCTURAL divergence: no result present in one arm only, no matrix
dimension difference, no row/column-name difference, no macro difference outside the three
expected (`wdid_fcol`, `wdid_fcskip`, `cmdline`) - checked at tolerance 1e-30 across all 51 cells.

```
 5.49e-06  w05_pct        controls(xpct)                              <- has continuous controls()
 2.60e-06  w03_dollar     controls(xdollar)                           <- "
 2.56e-06  w11_clu        controls(xdollar) cluster(region)           <- "
 2.01e-06  w10_rob        controls(xdollar) robust                    <- "
 8.75e-07  w09_ct         contreat + controls(xdollar)                <- "
 2.32e-07  v03_boot       bootstrap + controls(xdollar)               <- "
 1.81e-07  w01_c1         controls(x1)                                <- "
 7.68e-08  w04_dollar2    controls(xdollar xdollar2)                  <- "
 2.53e-08  v05_impdet     impvalues + controls, NOISELESS design      <- "
 5.34e-09  c09_ctrl       controls(x1 x2)                             <- "
 1.98e-09  v06_det        controls, NOISELESS design                  <- "
 ----------------------------- project threshold rtol 1e-10 ------------------------------
 5.64e-13  k02_ct_byT        2.85e-13 m10_messy       1.71e-13 m12_messy_cwm
 8.22e-14  k13_cwm0_ct       7.17e-14 m03_feunit      5.04e-14 w07_fvctrl(i.region)
 4.26e-14  v01_permex        3.62e-14 c02_robust      3.60e-14 w06_fe(region)
 3.58e-14  w02_cell(cell-constant control)             3.52e-14 m06_fastnose
 2.89e-14  m07_hist          2.89e-14 c01_plain       2.60e-14 k07_ct_aw
 2.41e-14  k08_pois          2.29e-14 c03_clreg       2.18e-14 m11_messy_ct
 2.16e-14  m04_esfbp         1.58e-14 k14_sg_ct       1.50e-14 m01_timetrend
 1.45e-14  k12_cwm0          1.45e-14 c04_aw          1.44e-14 c10_sg
 1.42e-14  m08_ifres         1.20e-14 k11_pois_semi   8.85e-15 k10_pois_dyex
 8.85e-15  k09_pois_ct       8.84e-15 k05_ct_ccno     7.78e-15 c06_pw
 7.24e-15  m09_esrel         6.30e-15 k01_ct_basic    5.91e-15 k06_ct_dyex
 5.51e-15  m05_antic         5.35e-15 k04_ct_byTt     2.98e-15 c07_impx
 1.85e-15  v04_imp           4.58e-16 k15_sg_clu      4.28e-16 c08_impnest
 0         n0 (nofastcollapse twice)   0  n1 (collapse twice)
 0         m13_rcs (repeated cross-sections)           0  m15_rcs_sg
```
**EVERY cell above the threshold has a continuous `controls()` term.  EVERY cell without one is
at or below 5.6e-13.**  The converse does NOT hold: the rows-per-cell ladder cs4..cs512 and the
first repro panel all carry a continuous `controls()` term and are clean, so a continuous
control is necessary but not sufficient - see LATE RESULTS (d).

### Two ways to reproduce NFC-01
1. **Self-contained, one file** (builds its own 1200-row panel, prints a noise floor, the four
   mechanism arms and the two printed tables to diff):
   `<WORKDIR>/repro/repro_NFC01_controls_V.do`
   NOTE: the machine stayed at 8-9 Stata instances (the project cap) for this whole session and
   this file was still queued when the report was written, so it has not been executed
   end-to-end.  Its ARMS are not untested - they are exactly the cells `p07.do` ran and are the
   source of every number quoted above; what is unvalidated is only the packaging.
2. **Already executed, evidence on disk**: `<WORKDIR>/mkdata.do` then `<WORKDIR>/p07.do`.
   `p07.log` contains the whole table above, the eleven `WDCELL|` sentinel lines, and the
   `DISPLAY_A` / `DISPLAY_B` printed-table blocks.  The per-cell dumps are in `<WORKDIR>/out/`
   (`w05_pct.A.txt` vs `w05_pct.B.txt` etc.) and are compared with
   `python3 <WORKDIR>/cmp.py out/w05_pct.A.txt out/w05_pct.B.txt w05_pct 1e-10`.

===========================================================================
## LATE RESULTS (pcrit) — these settle the two open questions and CHANGE the framing
===========================================================================

### (a) DETERMINISM: settled, exactly 0
On dCc, the dataset where the effect fires, WITH `controls(xdollar)` present:
```
  collapse arm run twice              worst rel = 0   (exactly)
  nofastcollapse arm run twice        worst rel = 0   (exactly)
  collapse arm, reversed session order worst rel = 0   (exactly)
```
NFC-01 is fully deterministic.  Confidence `certain`.

### (b) THE FRAMING CELL — this is the headline, and it INVERTS a documented calibration
All four runs use `controls(xdollar)` on the same data:
```
  F  row permutation, collapse OFF in both runs   2.08e-13   <- what a row shuffle costs
  C  collapse vs nofastcollapse, SAME row order   1.40e-06   <- what the collapse costs
  G  row permutation, collapse ON in both runs    8.70e-08
     (row permutation, NO controls, collapse off  1.95e-11)
```
`wooldid.sthlp` (the "When should you care?" paragraph) tells the user:

  "merely permuting the order of the rows in your dataset - leaving the data otherwise
   untouched, with the collapse off - moves the same near-zero coefficient about **19 times
   more** than the collapse does.  So if a conclusion is sensitive at this scale, it is
   already sensitive to the order your observations happen to be stored in."

**RETRACTED - see the CORRECTION section at the bottom of this file.**  That sentence is about a
COEFFICIENT, and on the actual near-zero coefficient the measured ratio is 15.5 against the
stated ~19, i.e. correct.  The numbers above are variances and do not test it.  What the
comparison DOES establish is the benign-alternative kill ("margins' delta-method V is just row-order sensitive once a
control is present") - with the collapse OFF, a full row shuffle costs only 2.08e-13.

**And the second-order consequence, arm G**: *with the collapse ON*, a row permutation moves
`e(V)` by 8.70e-08 instead of 2.08e-13 - about 400,000x more.  The collapse keeps each cell's
FIRST row, so which control value survives depends on the storage order of the user's
observations.  Re-`sort`ing a dataset therefore changes wooldid's reported standard errors in
the 7th digit, and the help's reassurance points the user in exactly the wrong direction.

### (c) SEVERITY: the pre-registered upgrade criteria did NOT fire — severity stays `medium`
- **Rows-per-cell ladder** (4 / 8 / 32 / 128 / 512 rows folded into one cell, same DGP):
  `2.45e-15, 5.07e-15, 4.62e-14, 4.27e-14, 2.85e-13` with the control; `0, 0, 0, 0, 0` without.
  **The disagreement does NOT grow with cell size.**  It is a property of the DESIGN, not of how
  much the collapse folds.
- **Rank-truncation flip** (cluster(cohort): 4 clusters against 9 horizons, so the ES covariance
  is rank deficient and the adversarial line test's 1e-12 truncation is live): x05_advc moves
  1.69e-06 WITH the control and 1.42e-14 without - but the reported SET is identical, no rank
  flip, no horizon dropped, no `e(wdid_fcguard)` asymmetry.
- **subgroups() + controls** (x07_sgc): 4.81e-06 on `e(V)`, and the subgroup joint tests move
  (`joint_horiz_s0` 1.5e-08, `joint_horiz_s1` 1.1e-08, `joint_pre_att_allsgssame` 2.5e-08) -
  again with no structural change.
- No decision anywhere flipped.  Severity **medium**.

### (d) WHY THE FIRST REPRO ATTEMPT FAILED, and what that tells us (important caveat)
My first self-contained repro built a simple hand-made panel (120 firms, y = a_i + trend +
0.5*x + 1.2*post + noise) and every arm came back at ~6e-16 - the effect did NOT fire.  The
rows-per-cell ladder uses the same simple DGP and is also clean at every rung.  So:
**the failure is design-dependent, not universal.**  It fires on the recon generator's panels
(dC / dCc: unit effects + a common period shock + horizon-varying effects) and not on a
hand-built panel with a single covariate.  It fires on 3 of the 4 base designs I built
(dC/dCc yes, dCm no, cs* no).  The shipped repro therefore regenerates the panel with the
recon generator so that it does fire.  Honest statement of scope: *this can and does happen
on ordinary applied specifications, but it does not happen on every dataset.*

### (e) two more control sets, resolved
```
  coarsecohortcontrols(c.xdollar)  1.38e-06   FAILS  (continuous, coefficient-bearing)
  ccc_absorb(xdollar)              5.76e-15   clean  (absorbed, no coefficient)
  interactivecontrols(xdollar)     EXACTLY 0  clean  (ineligible: collapse off in both arms,
                                                      which also verifies the ineligible list)
```

### (f) determinism, re-measured in an INDEPENDENT session (pnf.do) on `controls(xpct)`
```
  collapse arm run twice                        0   (exactly)
  nofastcollapse arm run twice                  0   (exactly)
  collapse arm again after an intervening -reg- 0   (exactly)
  nofastcollapse arm again later in session     0   (exactly)
  ---------------------------------------------------------
  the two arms against each other          5.49e-06   on e(V)
```
5.49e-06 is the same number p07 measured in a different session.  There is no run-to-run,
session-order or session-state component: the two paths simply compute different variances.

### (g) the framing result, apples-to-apples on ONE matrix element — `e(V)[2,2]`, `controls(xdollar)`
```
  row permutation, collapse OFF in both runs     6.63e-16    <- pure ULP.  nothing.
  collapse vs nofastcollapse, same row order     1.40e-06
  row permutation, collapse ON in both runs      8.70e-08
```
- The collapse moves this element **2.1e9 times more** than a row permutation does.
  `wooldid.sthlp` states a row permutation moves it about **19 times more than the collapse**.
  The documented calibration is inverted by ~4e10 on this element.
- Turning the collapse ON makes the same element **1.3e8 times more sensitive** to the order the
  user's observations are stored in (6.63e-16 -> 8.70e-08).  That is the direct consequence of
  keeping each cell's FIRST row: `sort`ing a dataset differently changes which control value
  survives, and therefore changes the reported standard error in the 7th digit.  The help's
  "if a conclusion is sensitive at this scale, it is already sensitive to the order your
  observations happen to be stored in" is precisely backwards - with the collapse off it is not.

===========================================================================
## CORRECTION to (b) and (g) — I over-reached, and the corrected version is cleaner
===========================================================================
`sthlp:357`'s sentence is about "the same **near-zero coefficient**".  My 2.1e9 ratio compared a
row permutation against the collapse on `e(V)[2,2]`, a **variance** — i.e. I checked a coefficient
claim against a variance measurement.  Recomputed on the actual near-zero coefficient
(`e(b_es_att)`'s smallest pre-treatment element, pre2 = -0.0417), same data, `controls(xdollar)`:

```
  quantity                      row-perm, collapse OFF | collapse vs nofc | row-perm, collapse ON
  b_es_att[pre2]  (coefficient)        2.06e-14        |     1.33e-15     |      2.13e-14
  e(V)[2,2]       (variance)           6.63e-16        |     1.40e-06     |      8.70e-08
```
Ratio for the coefficient: a row permutation moves it **15.5 times more** than the collapse does.
`sthlp:357` says "about 19 times more".  **THE DOCUMENTED CALIBRATION IS CORRECT AS WRITTEN.**
Retract the "inverted by 4e10" claim; it was an apples-to-oranges comparison.  Every b_es_att
column behaves the same way (collapse 1.9e-16 to 1.3e-15, row permutation 1.0e-15 to 5.7e-14).

**What survives, and is the accurate statement of the whole finding:**
1. On COEFFICIENTS the collapse is exactly as good as the help says, everywhere I looked —
   1e-16 to 1e-15, well under a row permutation.  All four documented bounds hold for `e(b)`.
2. On VARIANCES, with a continuous `controls()` term, it is not: `e(V)` moves 1e-7 to 5e-6,
   against `sthlp:338`'s measured "about 1e-16 relative", and that carries into the printed
   `se`, `lb_95ci` and `ub_95ci` columns (6th-7th decimal) against `sthlp:350`'s "no digit
   wooldid displays is affected".  This is a single, uniform SCOPE GAP: the documented figures
   were evidently measured on coefficients and/or without a continuous control, and the gate's
   own argument (ado:5588-5610) reasons only about the derivative of the prediction, never
   about `margins`' delta-method Jacobian.
3. **NEW, and the part with a practical edge (arm G, unaffected by the correction):** turning the
   collapse ON introduces a row-order dependence into the reported VARIANCE that does not exist
   with it off — `e(V)[2,2]` moves 6.63e-16 under a row permutation with the collapse off and
   8.70e-08 with it on, a factor of 1.3e8 — while leaving the coefficients' row-order behaviour
   untouched (2.06e-14 vs 2.13e-14).  The collapse keeps each cell's FIRST row, so `sort`ing the
   dataset differently changes which control value survives and hence the reported standard
   error in the 7th digit.  With the collapse off, re-sorting costs nothing at all.

===========================================================================
## REPRO v2 VALIDATED — `repro/repro_NFC01_controls_V.do` was executed and reproduces
===========================================================================
```
arm                                              fcol A/B | max rel e(V) | max rel e(b)
N  noise floor: collapse arm run TWICE             1/1    |  0.00e+00    |  0.00e+00
A  no controls                                     1/0    |  5.78e-16    |  1.7e-...
C  controls(xpct)  <- varies within the cell       1/0    |  5.49e-06    |  6.0e-...
D  controls(xcellconst)  <- constant in the cell   1/0    |  1.29e-15    |  2.6e-...
E  fe(region)  <- absorbed, no coefficient         1/0    |  5.47e-16    |  2.0e-...
```
Arm C's 5.49e-06 is the SAME number measured in p07 and again in pnf, in three separate Stata
sessions.  The noise floor is exactly 0.  The two mechanism arms (D, E) are clean.  The finding
reproduces from a single self-contained file.

===========================================================================
## SECOND CORRECTION — arms F and G of the validated repro RETRACT the arm-G claim
===========================================================================
The full repro run (all seven arms, `controls(xpct)`, estimands `att att_i att_it`):
```
N  noise floor: collapse arm run TWICE             1/1    |  0.00e+00
A  no controls                                     1/0    |  5.78e-16
C  controls(xpct)  <- varies within the cell       1/0    |  5.49e-06
D  controls(xcellconst)  <- constant in the cell   1/0    |  1.29e-15
E  fe(region)  <- absorbed, no coefficient         1/0    |  5.47e-16
F  ROW SHUFFLE, collapse OFF in both               0/0    |  2.92e-07
G  ROW SHUFFLE, collapse ON in both                1/1    |  2.28e-07
```
**F is 2.92e-07, not the ~1e-15 the earlier pcrit `y01` cell showed.**  The two are not
comparable: `y01` requested `att` only, this requests `att att_i att_it`.  So:

**RETRACTED: "the collapse introduces a row-order dependence into the variance that does not
exist with it off."**  It does exist with it off — F (2.92e-07, collapse off) and G (2.28e-07,
collapse on) are the same size.  A plain `sort` of the user's data moves `e(V)` by ~3e-7 whether
or not the collapse runs.

**THE SURVIVING, CORRECTED STATEMENT OF THE FINDING:**
Whenever a continuous, coefficient-bearing control varies within the collapse's cell, wooldid's
reported VARIANCE is reproducible only to about 1e-7 - 5e-6 relative:
  - between the two collapse paths                       5.49e-06   (arm C)
  - between two row orderings of the same data           2.92e-07   (arm F, collapse OFF)
Coefficients are unaffected throughout (1e-16 to 1e-15), and every other configuration I ran -
no controls, a cell-constant control, an absorbed control, a factor control - is at 1e-14 or
better.  So:
  - `sthlp:338` "Standard errors and variances do not amplify at all: they were measured at
    about 1e-16 relative regardless of the scale of the predictions" is wrong by ~10 orders in
    this configuration, and it is wrong for BOTH paths, so it is not solely a collapse issue.
  - `sthlp:350` "no digit wooldid displays is affected" is wrong: the 6th-7th printed decimal of
    `se`, `lb_95ci` and `ub_95ci` differs between the default run and `nofastcollapse` (measured
    verbatim in p07.log's DISPLAY_A / DISPLAY_B blocks).
  - `sthlp:332` the 1-ulp absolute bound holds for coefficients and fails for variances.
  - `sthlp:357`'s row-permutation calibration (~19x) is CORRECT for the near-zero coefficient
    (measured 15.5) and says nothing about variances.
The attribution is therefore: the collapse is the LARGER of two perturbations of the same
underlying instability in `margins`' delta-method variance under continuous controls, not the
sole cause.  For the `nofastcollapse` surface specifically, what matters is that the two paths
disagree by 4-5 orders more than the project's rtol 1e-10 and more than the help promises, so a
user reproducing an archived result with `nofastcollapse` does NOT get the same printed table as
the default run, contrary to `sthlp:350`.
