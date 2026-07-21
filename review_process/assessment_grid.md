# wooldid coverage grid

Every assessment and verification report must end with this grid, marking each cell:
`T` tested empirically, `S` covered by static/structural argument (say why), `U` untested
but relevant, `-` not applicable to the change under review. Untested-but-relevant cells
are acceptable to leave — hiding them is not: one verifier in every wave is assigned the
union of the assessors' `U` cells. The B15 pweight over-rejection and the B14 sd-guard
corner — the only two defects that reached committed code in the B1–B23 campaign — were
both unexamined cells of this grid.

Not every cell combination is required. The **mandatory minimum** for any fix:
- every value of each dimension the change directly touches;
- the full **weights** row and both **treatment modes** (these two dimensions produced the
  campaign's escapes);
- at least one **pathological data shape** relevant to the change;
- one clean bit-identity config per treatment mode (old vs new, %21x).

## Dimensions

**D1. Estimands**: att | att_it | att_i | att_itime | all four together (multi-estimand
loops carry cross-iteration state — B7).

**D2. Treatment mode**:
- binary (default)
- continuous: contreat() with contreatcontrols(basic) | contreatcontrols(untreatedZero)
  + safety(off) | other documented contreatcontrols values
- continuous extras: lattice(all N) | lattice(by* N) | lattice(baseTreated* N) |
  lattice(varname) passthrough | latticeignoreweights | contreatpoly(>1) | contreatwithin |
  contreatelasticitytype variants | ccnoabsorb

**D3. Weights**: none | [aw=] | [pw=] | non-integer/wide-dispersion weights |
customweightmultiplier(). (aw and pw differ in VCE mechanics — pw forces robust — and in
what summarize/egen accept; never test only one.)

**D4. Estimator path**: reghdfe (default) | ppmlhdfe (poisson) — incl. its rc-3200
fallback | plain reg (customvce) | plain poisson (customvce) | unconditionalse path.

**D5. Inference options**: cluster() | robust | unconditionalse | customvce(robust/cluster)
| customvce(non-robust) | impvalues | jointtests.

**D6. Sample restriction**: no if | comparison if (t<=8) | compound if | bare-flag if with
values {0,1,2,.} (the if-qualifier vs (if)==1 distinction — B23) | if referencing the
asserted/used variable.

**D7. Event-study structure**: none | espre/espost pooled reference | esfixedbaseperiod
with esrelativeto(-1) | esrelativeto(< -1) | timetrends | adversarial |
altmethod_adversarial.

**D8. Subgroups**: none | subgroups() with values incl. 0, gaps (3, 900), missings (. and
.a), realistic large ID-style codes | suppresssgprimaryeffects | makecfxplotsbysg.

**D9. Plot paths** (all produce e() content, not just pictures): makeplots |
histogramcohorteffects | makecfxplots with cfxplottypes(post | pre | integer relyrs |
mixes) | cfxplotnose | cfxlattice(varname) | saveplots(). Run with `set graphics off`.

**D10. Misc options**: safety(off) | oostreatedcontrols(>0) | verbose | cleanmatrices |
oldsyntax | poisexpresults | semielasticity | regtol.

**D11. Environment**: gtools active (the only supported path since `40414c7` made gtools
a hard requirement) | dependency-check negative test (a COPY with a package renamed in
the entry check's list must error r(198) informatively) | back-to-back runs in one
session (stale-state bugs — B7, B12, B19).

**D12. Data shape** — use `test/data/dta_path_*.dta` (built by `test/mkpathological.do`)
or construct equivalents:

Baseline behavior on `wooldid y i t ttre, att cluster(i) espre(2) espost(2)` as of commit
`1fd84c0`-era code (2026-07-21 smoke test) is recorded per shape so deviations are
detectable:

| shape | dataset | what it stresses | baseline rc |
|---|---|---|---|
| well-behaved balanced staggered panel | `test/data/dta_binary`, `dta_cont` | reference / bit-identity | 0 |
| tiny panel | `dta_tiny` | small-sample paths | 0 |
| no never-treated units (staggered-only identification) | `dta_path_nonever` | control-group construction, `ttre == .` logic | **459** — errors even though not-yet-treated controls exist mid-sample; see the below-the-line observation logged 2026-07-21 |
| no controls at all (every unit treated, same onset) | `dta_path_nocontrols` | the no-sample/estimability error paths — must fail informatively, not crash | 459 (informative, as intended) |
| cohort with no pre-period (treated at first sample period) | `dta_path_nopre` | reference-period selection, estimability screens | 0 |
| cohort with no post-period (onset after sample end) | `dta_path_nopost` | out-of-sample-treatment handling, oostreatedcontrols interaction | 0 |
| cohort missing specific relative years — incl. relyr −1 (the default reference) and a mid-window post gap | `dta_path_relyrgap` | horizon construction, reference-period fallbacks, ES plots, adversarial line; try esfixedbaseperiod variants | 0 (pooled ref) |
| exactly 1 observation per ttre-year cell (single-unit cohorts) | `dta_path_onepercell` | cell-weight denominators, singleton drops, clustering with few clusters | 0 |
| deeply imbalanced panel (units span 2–12 periods) + an fe() singleton-inducing variable | `dta_path_imbalanced` | e(sample) shrinkage (B10 class — `fe(extraFE)` drops 2 obs), ragged horizons | 0 |
| poisson separation (all-zero count cells/units) | `dta_path_sep` | ppmlhdfe drops, count models | 0 linear; **322** under poisson (wooldid's own dropped-dummy guard — expected) |
| extended missings (.a–.z) sprinkled per user variable | build ad hoc | input screens (B9 class) | — |
| calendar gaps (years absent from the panel entirely) | build ad hoc | `_t` construction, mint/maxt logic | — |

## Report format

End your report with the grid as a table: one row per dimension, listing the values you
tested (`T`), argued (`S` + one-line reason), left untested (`U`), or n/a (`-`). Example:

```
D3 Weights:  none T | aw T | pw U (not run — flag for verifier) | cwm -
D12 Shape:   balanced T | relyrgap T | onepercell U | others -
```
