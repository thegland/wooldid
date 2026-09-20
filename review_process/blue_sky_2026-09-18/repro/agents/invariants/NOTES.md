# Recon agent: METAMORPHIC INVARIANT CATALOG — working notes

Output file: ../../recon/invariants.md
Workdir files:
  wdgen.do        shared DGPs (wdgen_d1 unit-level, wdgen_d2 coarse-cohort) + wdsnap (%21x e() dump)
  cmp.py          compares tagged blocks in a wdsnap file
  p01_harness.do / p01.txt / p01_harness.log   relabel/roworder/dup/affine/y-shift/weights
  p02_degen.do   / p02.txt                     fcol, sg-one-level, anticipation, if-vs-drop, t*2 no-ES, expand-vs-aw, aw*1024
  p03_misc.do    / p03.txt                     determinism control, sg relabel, simultaneousbands, gteffects

## MEASURED (all on d1: 60 units, 10 periods, 3 onset arms + never-treated, cluster(clus)=3 units)
Baseline command: `wooldid y i t ttre , att att_i att_it att_itime esprelength(3) espostlength(3) cluster(clus)`

| transform | verdict |
|---|---|
| same command run twice                   | BIT-IDENTICAL (77/77). bit comparison is a valid tool on the reghdfe path |
| random row permutation                   | b moves ~2 ulp (4e-16 rel); V moves ~1e-14. OK |
| relabel i (i -> 1000-7i)                 | b moves ~2 ulp; V ~1e-13. OK |
| duplicate EVERY row                      | b identical to 2 ulp; **V shrinks by 1.72x** — consistent with reghdfe's (N-1)/(N-K) with K~355, N=600. df_r stays G-1=19 |
| t -> t + 500 (ttre too)                  | **BIT-IDENTICAL, 83/83** |
| t -> 2t (ttre too), WITH esprelength(3) espostlength(3) | **att moves 4.2%**, ES matrices lose rows |
| t -> 2t (ttre too), NO event-study window | **BIT-IDENTICAL** -> the 4.2% above is entirely the calendar-unit window semantics |
| y -> y + 100                             | b identical to 2 ulp. OK |
| fastcollapse vs nofastcollapse (+controls(x) varying in cell) | **BIT-IDENTICAL 77/77** |
| subgroups(constant) vs no subgroups()    | main rows **BIT-IDENTICAL**; adds 20 SG keys |
| subgroups level relabel {0,1}->{7,3}     | **BIT-IDENTICAL after renaming SG0->SG7, SG1->SG3** (111/111) |
| anticipation(2) vs hand-shifted ttre-2   | **BIT-IDENTICAL 77/77** |
| `if t<=2006` vs `keep if t<=2006` first  | **BIT-IDENTICAL 77/77** |
| expand 2 if i<=4  vs  [aw=ww], ww=2 for i<=4 | b agrees to 1e-14..3e-16; V differs (dof). OK |
| [aw=w] vs [aw=1024*w]                    | b bit-identical, **V moves ~3.3e-9 relative** (above the 1e-10 project threshold) |
| aw vs pw on identical weights UNDER cluster() | b and V **bit-identical** |
| simultaneousbands: add att_it to an att-only run | e(esband_crit_att) unchanged; c(rngstate) after the run identical with and without simultaneousbands |

## LEADS THAT SURVIVED (handed to the catalog)
L1  esprelength()/espostlength() are counted in CALENDAR units of t, so att itself depends on
    the units t is coded in (t=1,2,3 vs t=2000,2002,2004). Help says "relative-period slots".
L2  aw scale (exact power of two) moves V by 3e-9. Probably reghdfe tol(1e-9) but point
    estimates did NOT move, which is odd for a convergence explanation. Worth a dedicated probe.
L3  The analytic V is extremely sensitive to N through (N-1)/(N-K) because K/N ~ 0.6 when i is
    the unit id — the same regime as the confirmed undercoverage.

## DELIVERABLE WRITTEN
/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/invariants.md
13 sections, ~110 invariants, tolerance classes 0-3 measured, ranked hit list at s11.
d1.dta is saved in this directory and is the reference dataset the catalog's sketches use.

## POST-REVIEW CORRECTIONS APPLIED (do not undo)
- permute() invariance is THREE-TIER (assignment matrix Class 0 via permsave; permdist_b
  Class 1; perm_p exactly identical). It is NOT "everything Class 0" - each rep is a full
  wooldid run and row order alone moves b by 2 ulp.
- bootstrap cross-leg parity (E1) is CLASS 1, not bit-identical: the physical leg gives
  duplicated clusters separate unit/cell dummies while the count path shares them under fw.
- `att boottstat robust` vs `att boottstat cluster(clus)` is NOT a parity pair (different
  resampling units; clusv falls back to an obs-level id at ado:8918-8921). Use
  `att boottstat robust` vs `att att_i boottstat robust`.
- leap-year %td split downgraded to medium/UX (it surfaces as visible post365/post366 rows).
- new hit-list item 15: simultaneousbands can be NARROWER than the pointwise CI at few
  clusters (simulated sup-t crit vs t_{G-1}).
