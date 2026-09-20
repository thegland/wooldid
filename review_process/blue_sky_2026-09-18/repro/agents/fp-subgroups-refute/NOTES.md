# SKEPTIC notes — refuting fp-subgroups findings fpsg-01 / fpsg-02

Build asserted in every log: `which wooldid` -> .../bluesky/build/wooldid.ado (r01 log line 76).

## State so far

### fpsg-01 — CONFIRMED (independently, on my own designs)
- `r01_sgmissing.do` / `.log`. My design, not theirs: 5 units in sg1, effects 1..5,
  observed shares 0.125/0.25/0.5/1/1. equal-unit target 3.0, share-weighted target
  3.8695652174.
  - RUN A (clean sg):  att_i = att_itime = 3.0000000000 exactly (dev -7.1e-15).
  - RUN B (sg missing on some post rows): att_i = att_itime = 3.8695652174,
    dev from SHARE = -3.1e-15, dev from EQUAL = +0.8696 (+29%).
  - RUN C (same rows, sg CHANGES 1->0): r(459). Asymmetry confirmed.
  - NO caution anywhere in RUN B. Header prints "0 obs dropped".
- `r02_mechanism.do` / `.log` PART 1 — the mechanism is *identified*, not just the error.
  UNBALANCED panel pulls apart three candidate weightings:
  EQUAL 5.0 / SHARE 3.6666666667 / ROW-WEIGHTED 3.0.
  att_i and att_itime = 3.6666666667 (dev from SHARE 1.3e-15) -> SHARE, exactly as claimed.
  att and att_it = 3.0 -> row-weighted, correct for them.
  So it is NOT merely "att_i collapses onto att"; it is specifically
  w / (unit's TOTAL post weight over the whole estimation sample).
- Code corroboration (read-only): `w_per_i` at ado:1621 is `gegen total(w*esamp), by(posttreatment i)`
  — no subgroup term; `subgroups` is deliberately absent from the must-be-nonmissing
  markout list (ado:~700, comment "normally, it is ok to have missing values in the sg
  variable"); the F13 split guard at ado:12266-12280 folds missing via
  `cond(sg>=., ., sg)` with the comment "missing sg rows belong to no subgroup by design".
  Subgroup margins at ado:2819 use `over(subgroups)` on the full post subpop.
- Not in known_bydesign.md. §3 covers the WITHIN-CELL straddle (refused) only.

### fpsg-02 — CONFIRMED
- `r03_names.do` / `.log`. sg levels 1,3,7,99.
  All six help-shaped names (`main_att_s1`, `main_att_s2`, `main_att_it_s1`,
  `pretreatment_att_s1`, `ES_att_post0_s1`, `main_att_s99`) -> **r(111)** from `lincom`.
  All five infix names (`main_SG1_att`, `main_SG3_att`, `main_SG99_att`,
  `ES_att_SG1_post0`, `pretreatment_SG1_att`) -> rc 0.
  `e(b_es_att_s1)`, `e(b_es_att_s99)`, `e(escomp_s3)`, `e(numpre_att_s99)` all exist:
  the `_s<value>` token is real but belongs to the e() NAME family.
- DUPLICATE: agent `hv-sg-names` filed the identical divergence as HV-SGNAME-01,
  with more surface coverage (`grep -c '_SG' wooldid.sthlp` == 0). Verdict unchanged,
  but the coordinator should merge.
- ONE CORRECTION to the claim: "carries the subgroup VALUE, not a sequence number" is not
  itself a defect — `#` is Stata help idiom for "a number", and `_s99` satisfies `_s#`.
  The defect is POSITION (infix, not suffix) and SPELLING (`SG`, not `s`).

### NEW (mine) — the PRE-treatment analogue
`r02_mechanism.do` PART 2. sg missing on a unit's PRE rows only.
Clean: pretreatment att/att_it/att_i/att_itime all 0.7; ES_att:pre3 = pre2 = 1.4.
Missing: pretreatment rows all 0.2666666667, ES pre3=pre2=0.6.
0.2666666667 is exactly the SHARE-weighted pre figure (0.6/2.25); the equal-unit pre
figure is 0.2. MAIN att_i is UNAFFECTED (5.0 exactly).
=> the same unit is weighted 1.00 in the `main:att_i` row and 0.25 in the
`pretreatment:att_i` row of the SAME printed table, silently. This is the surface that
feeds `jointtests` pre-trend tests and the honestdid handoff.

## TODO
- r04: realistic noisy magnitude, MCAR vs effect-correlated missingness, oracle from
  e(gteffects) (NOT e(histogramestimates) — hist margins at ado:12405 has no sg term in
  its subpop, so it averages over ALL post cells, in-subgroup or not).

## FINAL STATE (all probes run, FINDINGS.md complete)
- r01_sgmissing.do/.log  — fpsg-01 independent repro. CONFIRMED.
- r02_mechanism.do/.log  — PART1 mechanism identified (share, not row). PART2 = NEW skp-02.
- r03_names.do/.log      — fpsg-02 CONFIRMED (dup of hv-sg-names HV-SGNAME-01).
- r04_realistic.do/.log  — magnitude: MCAR -0.115 SE; label-stops-being-collected -1.468 SE.
- r05_contrast.do/.log   — NEW skp-03: e(pwsgs) contrast CI stops covering (8.3 SE).
- r06_attit.do/.log      — NEW skp-01: att_it ALSO distorted (contradicts "att_it unaffected").
Verdicts: fpsg-01 CONFIRMED high; fpsg-02 CONFIRMED medium.
