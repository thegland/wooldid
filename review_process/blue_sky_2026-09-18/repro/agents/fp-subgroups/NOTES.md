# fp-subgroups — running notes

Target: subgroups(), contrasts, homogeneity tests.
Method: independent reimplementation + degenerate-limit identities.
Build asserted in every log via `which wooldid` -> the pinned build path.

## Probes run

| file | what | verdict |
|---|---|---|
| p01_onelevel.do | one-level subgroups(sgall) vs no-subgroups | CLEAN, bit-identical |
| p02_sgtruth.do  | noiseless DGP (gysd=0), exact recovery of subgroup truths | CLEAN to 1e-15 |
| p03_identities.do | aggregation identity A, contrast identity C, 2-sg joint==pairwise D, covariance use S | CLEAN |
| p04_joint4.do | FULL reimplementation of the 4-level homogeneity test | CLEAN to 1e-22 |

## Results so far

### p01 — one-level subgroups
`subgroups(sgall)` with sgall==1 constant.
`e(wdidmainresults)` is **bit-identical** (%21x) to the no-subgroups run in every
column (estimate/se/t/p/lb/ub/relyr), and so are e(b_es_att), e(V_es_att).
e(N), e(df_r) identical. The only difference is the extra 16 SG1 rows appended to
e(b)/e(V). NO LOCALISED BUG.

### p02 — exact recovery on a noiseless DGP
gysd=0, effects made subgroup-dependent. i=unit, cluster(unit).
All four estimands, pooled and per-subgroup, match the design truth (computed by
independent reweighting of `truetau`) to <= 1.8e-15 absolute.
`e(pwsgs)` estimate == SG1 - SG0 exactly. Row key is `1vs0.sg2` (VALUES, not labels).
NOTE: on this NOISELESS design every joint test came back missing (`.`) — the
covariance is numerically singular when the fit is exact. That is a pathological-data
artifact, not a defect; re-ran everything with noise in p03/p04.

### p03 — algebraic identities with noise (2 subgroups)
i = cohort x region (20 levels), fe(i.unit), cluster(unit) (240 clusters, FINER than i).

* **Identity A (aggregation).** pooled estimand == subgroup-weight-weighted mean of
  the subgroup estimands, with W = {sum w; # post cells; # units; # units} for
  {att, att_it, att_i, att_itime}.  **reldiff = 0.000e+00 for all four.**  This is
  the strongest available check that subgroups() only reweights the same cell
  effects, and that the weights are the documented ones.
* **Identity C.** e(pwsgs) estimate == b1 - b0 to 6e-16.
* **Identity S.** se_contrast is NOT sqrt(se0^2+se1^2): implied rho = 0.16-0.20.
  The off-diagonal covariance IS used (margins pwcompare supplies it).
* **Identity D.** the "all SGs the same" p == the single pairwise p to 1.1e-16.

### p04 — independent reimplementation of the homogeneity test, 4 levels (1 3 7 99)
Recovered the full 4x4 level covariance from the 4 level SEs + 6 pairwise SEs via
C_ab = (s_a^2+s_b^2-s_ab^2)/2, rebuilt W = (Rb)'(RVR')^-1(Rb) with R the 3 chain
restrictions, and compared.

* every pairwise estimate == b_A - b_B **exactly** (d = 0.000e+00, all 24 rows)
* every implied rho in [0.04, 0.17] — in range; rebuilt V positive definite
  (eigenvalues 2.4e-2 .. 3.8e-1)
* **rebuilt p = Ftail(3, e(df_r), W/3) matches e(joint_*_allsgssame) to 1e-22**
  for all four estimands.  The chi2(3) reading does NOT match (off by 4e-10), so
  the reference distribution is confirmed to be F(k-1, df_r), which is right.

=> the ANALYTIC subgroup level / contrast / homogeneity machinery is CORRECT.

### NAMING CONTRACT (verified by running, seeded lead CONFIRMED)
`sthlp:1307`: "Under subgroups() each subgroup's rows take the suffix _s#."
FALSE for e(b)/e(V) column names. Actual, enumerated from a 4-level run:
  main_SG1_att  main_SG3_att_it  main_SG99_att_itime
  pretreatment_SG7_att  pretreatment_SG99_att_itime
  ES_att_SG3_pre2  ES_att_SG99_post1  ES_att_it_SG1_post0
i.e. an **infix `_SG<value>_`**, and the token is the subgroup VALUE (99), not a
sequence number (#).
The `_s<value>` SUFFIX convention is real but applies to a different family —
e() matrix/scalar NAMES: e(b_es_att_s99), e(V_es_att_it_s3), e(escomp_s7),
e(numpre_att_s99). Those the help documents correctly elsewhere.
Consequence: a user following sthlp:1307 writes `lincom main_att_s3` (error) or
`test main_att_s1 = main_att_s2` (error) instead of `main_SG1_att`.

## Harness gotchas hit (for a successor)
* `margins ... pwcompare` names the FIRST contrast row `3vs1bn.sgq` — the `bn`
  base marker is part of the rowname. Strip non-digits when parsing.
* Mata `substr(s,i,.)` does NOT mean "to end"; use the 2-arg form.
* A `matrix X[i,j] = ` assignment inside a nested do-file foreach after certain
  local reuse threw "invalid syntax"; moved the whole reconstruction into Mata.

## Still to do
- triple-diff separability oracle (block-diagonal design: subgroup estimate MUST
  equal the subsample run)
- suppresssgprimaryeffects: contrasts must be identical to the full run
- BS-15 degeneracy screen: does it actually fire at cluster() == i?
- negative / zero / long subgroup values; missing sg; tiny subgroup
- weights (aw/pw/fw) x subgroups
- subgroup event-study rows and per-subgroup pre-trend tests

## FINAL STATE (all probes run)

| file | what | verdict |
|---|---|---|
| p01_onelevel.do | one-level subgroups vs none | CLEAN, bit-identical |
| p02_sgtruth.do | noiseless exact recovery of subgroup truths | CLEAN 1e-15 |
| p03_identities.do | aggregation A, contrast C, k=2 joint==pairwise D, covariance S | CLEAN |
| p04_joint4.do | full reimplementation of the 4-level homogeneity test | CLEAN 1e-22 |
| p05_separable.do | block-separable triple-diff: subsample identity, rho==0, suppresssgprimaryeffects | CLEAN (rho = 0.00000000) |
| p06_espretrend.do | per-subgroup pre-trend Wald, e(V) splice, ES aggregation | CLEAN 5.6e-17 / 3.5e-18 |
| p07_stress.do | degeneracy screen at cluster==i, aweights, missing sg, tiny sg | screen FIRES as documented; N7 -> fpsg-01 |
| p08_ragged.do | ragged subgroup ES horizons, exact LABEL oracle | CLEAN 1e-14, no label slip |
| p09_edgevalues.do | negative / non-integer / 8-digit / 7-digit / never-treated-only level / no-pre-period | all as documented; N7 -> fpsg-01 |
| p10_partialmissing.do | **REPRO for fpsg-01** | share weighting confirmed, 8e-14 |

FINDINGS.md holds the three findings and the twelve negative results.

Deliverables: FINDINGS.md + p01..p10 do-files and logs in this directory.
The only self-contained repro a verifier needs for the one substantive finding is
p10_partialmissing.do; p04_joint4.do carries the naming-contract evidence.
