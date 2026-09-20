# SKEPTIC NOTES — hv-lattice-maxvar refutation

## Orientation (done)
- brief.md, known_bydesign.md read.
- known_bydesign SECTION 5: "`si8demand` is a deliberate LOWER bound ... The design accepts
  false negatives to guarantee ZERO FALSE POSITIVES. Opened as a defect and RETRACTED the
  same day, 2026-09-16."
- Retraction read verbatim: /Users/hegland/Documents/Research/wooldid_dev/review_process/open_defects.md:63-85
  It is SILENT on the help wording. It justifies the design on "zero false positives, because
  a false positive would refuse a specification that would actually have estimated fine."
- guard_cost_2026-09-16.md:172-175, 342, 404 already records the under-prediction (>2x) on a
  200x20 panel. So hv-lat-01's MECHANISM is in the project record already.
- sthlp:1167 / sthlp:1207 confirmed to contain no lower-bound hedge (grepped "lower bound",
  "at least", "variable count" across the whole sthlp).

## BIG RESULT (s1_falsepos.do / .log) — MY OWN, contradicts known_bydesign section 5
ARM 1: 2200-unit x 12-year panel, 20 treated, i=unit, maxvar 2048, DEFAULT safety
   -> wooldid: "Error: this specification requires at least 2542 variables (about 130
      treated cohort-by-period interaction terms) but maxvar is currently 2048" -> rc 907
ARM 2: IDENTICAL data, IDENTICAL maxvar 2048, only safety(off) added
   -> rc 0, att = 1.4700209 (truth 1.5), e(N)=26400, **1.91 SECONDS**
=> FALSE POSITIVE. The guard refused a specification that estimates in under two seconds at
   that very maxvar. si8demand 2542 vs a design matrix of ~130 treated cells: overshoot ~19x.
Mechanism: ado:1476 is `reghdfe ... , absorb(i.`i' i.`t' `fe')`. i and t levels are ABSORBED,
   never design-matrix columns. si8ilev + si8tlev is pure additive slack growing one-for-one
   with the NEVER-TREATED unit count. Formula reproduced exactly: 130 + 2200 + 12 + 200 = 2542.
Also confirms hv-lat-03's counter-observation: Stata's canned r(907) text DOES print after
   wooldid's two lines when the preflight fires.

## TODO
- s2: default maxvar 5000 arm (5200 units) -- same false positive with NOTHING set.
- s3: my own lattice repro of hv-lat-01.
- s4: hv-lat-02 family check.

## s2/s3 (default maxvar arm) -- CONFIRMED
- 5200 units, untouched default maxvar 5000, DEFAULT safety:
  "requires at least 5542 variables (about 130 treated cohort-by-period interaction terms)" -> rc 907
- identical data + maxvar, safety(off): rc 0 in 3.75 s, att = 1.3745078
- identical data + DEFAULT safety at maxvar 8000: rc 0, att = 1.3745078, hex +1.5fdfbe7391aecX+000
  => safety(off) changed NOTHING about the estimate; the refusal was purely the variable count.
- hand-built wooldid-shaped design on the same data: reghdfe reports F(220, 56969), i.e. 220
  estimated columns of 240 (20 omitted). Claimed "at least 5542". Overshoot ~23x.
- 2200-unit and 5200-unit panels share si8cells=130; si8demand moves 2542 -> 5542, exactly the
  3000-unit change in si8ilev. The deciding term costs zero columns.
- MECHANISM PRECISION: the ado comment's "lower bound" attaches to si8cells*(1+contreat), the
  treated interaction-cell count. si8demand then ADDS si8ilev + si8tlev + 200, and on the
  default path (reghdfe absorb(i.`i' i.`t')) those are not columns at all. So the invariant the
  retraction relies on is a property of one TERM of si8demand, not of si8demand.
- IRONY worth reporting: at 5200 units the SAME block prints, one line above the refusal,
  "a granular unit can dominate runtime long before it approaches Stata's variable limits" --
  and then refuses on the variable limit a run that takes 3.75 seconds.
- sthlp never contains the words "lower bound" (grep count 0).

## s4 (my own lattice panel, 40 states x 10 counties x 18 yrs)
- si8demand PREDICTED 1030 (cells=386, ilev=40, tlev=18); guard silent at maxvar 5000 -- correct,
  because the actual is also below 5000 here.
- ACTUAL columns: no lattice ~775; lattice(all 4) ~3685; lattice(all 5) ~4191.
  => under-prediction factor 3.6x / 4.1x. Reproduces the prober's magnitude on my own data.
- RUN 1 (no lattice) rc 0 in 13.15 s.
- RUN 2 (lattice all 4, ~3685 columns, 7200 rows) was still running after ~4 MINUTES and I
  KILLED IT MYSELF to free a Stata slot (not a deadline kill, and never finished).
  Relevant to hv-lat-04: >=18x the 13.15 s no-lattice arm and still unfinished, under
  8-process contention. Column ratio 4.75x, 4.75^1.9 ~ 19. Weakens hv-lat-04.

## s6 (scaled own lattice panel, 50 states x 8 counties x 22 yrs) -- hv-lat-01 REPRODUCED
- si8demand 1568 (cells=648, ilev=50, tlev=22) -> silent. Actual: no lattice ~1299,
  lattice(all 5) ~6827. 1568 < 5000 < 6827. Under-prediction 4.35x.
- RUN1 rc 0 in 44.00 s; RUN2 rc 907 in **0.194 s**.
- CORRECTION: the death costs 0.19 s, not an expensive partial estimation. "mid-run" (sthlp:1167,
  and the prober's framing) overstates the harm: the loss is the MESSAGE, not compute.
- dataset survives (_N=8800, 11 vars), WOOLDID_LASTRC=907, bare replay refused r(198) with
  "the last wooldid call failed with r(907); the stored results predate that failed call."

## s5 (family, own data, maxvar 2048, NO lattice) -- hv-lat-02 REPRODUCED
- si8demand 642 on EVERY arm (cells=199, ilev=30, tlev=14); guard silent throughout.
- ARM1 baseline rc 0 | ARM2 +interactivecontrols(x1..x8) rc 907 | ARM3 +coarsecohortcontrols(i.reg)
  rc 907 | ARM4 binary+interactivecontrols rc 0 (si8demand 443).

## Mechanism precision (checked in the ado)
- ado:928-929: `_i' = i if _t>0 & _t<., else 0. ado:1009: targetflag = (_t>0). So wooldid's
  own design carries ONLY post-treatment cells as nonzero -- si8cells is the right count of THAT
  term, and my s3 hand-build (which put _i=unit on all rows of treated units, 240 cols) is a
  generous UPPER proxy. The overshoot in s2 is therefore >= 23x, reported conservatively.
- The "lower bound" in the ado comment attaches to si8cells*(1+contreat). si8demand then adds
  si8ilev + si8tlev + 200, which on the reghdfe/ppmlhdfe path are absorbed and not columns.
  On the reg/poisson path (ado:1490, `unconditionalse`/`customvce()`) they ARE columns, so the
  formula is right there and wrong on the path essentially every run takes.

## s7 (3,143 US counties x 18 yrs, 100 treated, contreat) -- PARTIAL, honestly
- ARM 1 DEFAULT safety, DEFAULT maxvar 5000: "Error: this specification requires at least 5761
  variables (about 1200 treated cohort-by-period interaction terms)" -> rc 907.
  Formula check: 1200*2 + 3143 + 18 + 200 = 5761. Exact.
- ARM 2 (same + safety(off)) did NOT complete inside the 500 s deadline and was killed.
  So at county scale WITH contreat I did NOT prove the run is feasible; the design there is
  genuinely expensive (~2,401 columns on 56,574 rows). Reported as arithmetic only.
  => s8 is the cheap realistic continuous-treatment anchor instead (5,000 ZIPs, 60 exposed).

## s8 (5,000 ZIPs x 12 yrs, 60 exposed, contreat) -- FALSE POSITIVE ON THE CONTINUOUS PATH
- ARM1 DEFAULT safety, DEFAULT maxvar 5000: "requires at least 5992 variables (about 390 treated
  cohort-by-period interaction terms)" -> rc 907.   (390*2 + 5000 + 12 + 200 = 5992, exact)
- ARM2 identical + safety(off), SAME maxvar: rc 0 in 68.18 s, estimate .63452766, N = 60000.
- True demand ~781 columns. Overshoot ~7.7x.
