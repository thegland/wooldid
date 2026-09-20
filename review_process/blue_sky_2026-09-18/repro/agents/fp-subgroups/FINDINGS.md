# fp-subgroups — FINDINGS

Agent: FIRST-PRINCIPLES ORACLE for `subgroups()`, contrasts, homogeneity tests.
Method: independent reimplementation + degenerate-limit identities.
Every probe is a self-contained do-file in this directory; every log opens with
`which wooldid` asserting the pinned build.

---

## fpsg-01 — a partially-MISSING subgroup variable silently turns subgroup
## `att_i` / `att_itime` into a SHARE-weighted average  [FOUND BY ME]

**claim.** When the `subgroups()` variable is missing on only SOME of a unit's
rows (item nonresponse — the unit's subgroup is known in some years and not in
others), wooldid accepts the run with no error and no caution, and the subgroup
`att_i` / `att_itime` it reports are NOT the equal-unit averages those estimands
are defined to be. The unit enters its own subgroup weighted by the share of its
post-treatment sample weight on which the subgroup happens to be observed.

wooldid REFUSES (r459, default safety) the ADJACENT case — a unit whose subgroup
value *changes* over time — and its own error text gives exactly this reason:
  "att_i and att_itime enter a split unit into each subgroup weighted by its
   share of sample weight there, rather than with the equal unit weighting
   those estimands normally imply"
The missing-value case produces the same distortion and fires no guard. The help
says only "Observations with missing subgroup values remain allowed and belong to
no subgroup", which reads as benign.

**mechanism.** `w_att_i = w / (total post weight of i)` is built once on the whole
estimation sample. `margins, over(sg)` renormalises within the over-group, so a
unit observed in the subgroup on a fraction s of its post rows contributes total
weight s instead of 1.

**repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-subgroups/p10_partialmissing.do`
run:  cd to that directory, then
      "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p10_partialmissing.do

**observed.** Noiseless design, 5 units in subgroup 1 with unit effects
(5, 1, 1, 1, 1); unit 1's subgroup value blanked in 5 of its 6 post years.
```
run A (clean sg)       att_i(sg==1) = 1.8000000000     <- correct
run B (sg partly .)    att_i(sg==1) = 1.1600000000     <- reported
     equal-unit definition          = 1.8000000000     d = -0.64
     share-weighted                 = 1.1600000000     d =  8.19e-14
run C (sg CHANGES 1->0 on the same rows)  = r(459), refused with the text above
```
No caution of any kind in run B (only the unrelated identical-SE notice that the
noiseless fit produces).
Read the run-B row carefully: `att` and `att_it` also print 1.16, and that is
CORRECT for them — att_it is 29/25 by cell count and att is the row-weighted
mean. The coincidence is an artifact of this design, not evidence that all four
are affected. Run D below is the discriminating number.
Noisy confirmation in the same file (run D, y + N(0,0.7)): wooldid reports
att_i(sg==1) = 0.9545 where the equal-unit average of the five units' own
`histogramcohorteffects` estimates is 1.7485 — a 45% departure. (The exact
share-weighted figure there is 0.9486 rather than 0.9545 because unit 1 also
contributes only its single in-subgroup CELL rather than its 6-period average;
the noiseless run B is the exact algebraic check.)

**expected.** The strongest argument is the asymmetry: wooldid's OWN r(459) text
names this exact distortion as the reason it refuses the value-changes case, and
the missing-value case reaches the same place unguarded. So: either the
equal-unit average over the units in the subgroup (the documented meaning of
`att_i` / `att_itime`), or the same refusal / caution the adjacent case gets.
The author may instead rule this a documentation gap; the number itself is still
not the estimand the table's row label names.
`att` and `att_it` are unaffected and correctly so: they are row- and
cell-weighted, so dropping rows is the right behaviour for them.

realistic: **true** (item nonresponse in a subgroup variable is ordinary in
survey and administrative panels, and passing the raw variable is the natural
user action).
severity: **high** — a silently wrong number on realistic input.
confidence: certain (exact noiseless oracle + the refused adjacent case).
surface: `subgroups()` x `att_i` / `att_itime` x missing subgroup values.

---

## fpsg-02 — sthlp:1307 documents the WRONG e(b) row names for subgroups
[SEEDED — recon reported it from reading; verified here by running, and the real
name set enumerated]

**claim.** `sthlp:1307` says "Under `subgroups()` each subgroup's rows take the
suffix `_s#`." The e(b)/e(V) column names actually use an **infix `_SG<value>_`**
and carry the subgroup VALUE, not a sequence number. That paragraph is the one
that says these names "are what `esttab`, `coefplot`, `lincom`, `test` and the
`permtstat`/`bootsave()` tables address a coefficient by", so a user following it
gets "not found".

**repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-subgroups/p04_joint4.do`
(section "=== e(b) colnames (NAMING CONTRACT) ===" of `p04_joint4.log`).

**observed.** subgroup levels 1, 3, 7, 99:
```
main_SG1_att        main_SG3_att_it      main_SG7_att_i      main_SG99_att_itime
pretreatment_SG1_att                     pretreatment_SG99_att_itime
ES_att_SG3_pre2     ES_att_SG99_post1    ES_att_it_SG1_post0
```
The `_s<value>` SUFFIX convention is real but applies to a DIFFERENT family — the
e() matrix and scalar NAMES: `e(b_es_att_s99)`, `e(V_es_att_it_s3)`,
`e(escomp_s7)`, `e(numpre_att_s99)`. The help has transplanted that rule onto the
coefficient rows.

**expected.** `main_att_s1` etc. per sthlp:1307.
realistic: true. severity: medium. confidence: certain.
surface: `subgroups()` + e(b)/e(V) naming contract.

---

# PROBED HARD, CAME BACK CLEAN  (negative results)

### N0. the `bn` base marker in the stored e(pwsgs) key is HARMLESS — retracted
I noticed that the first pairwise row of each family is STORED as `3vs1bn.sgq`
while the display and `matrix list` show `3vs1.sgq`, and drafted it as a finding.
`p11_pwkey.do` settles it: `rownumb(e(pwsgs), "att:3vs1.sgq")` returns 1 and
`e(pwsgs)["att:3vs1.sgq","estimate"]` returns the right number. Stata resolves
the operator-free key against the stored one, so a user indexing by the printed
key is fine. Only a raw `st_matrixrowstripe` string comparison sees a difference.
NOT a finding; noted so nobody re-files it.

### N1. one-level `subgroups()` reproduces the no-subgroups run EXACTLY
`p01_onelevel.do`. With a constant subgroup variable, `e(wdidmainresults)` is
**bit-identical in %21x** in all 7 columns and 16 rows, as are `e(b_es_att)`,
`e(V_es_att)`, `e(N)`, `e(df_r)`. Would have localised any subgroup code path
that perturbs the pooled margins. (On a larger design with `fe()` and partially
missing subgroup values the pooled table moves by 1-3 ulps — well inside the
project's 1e-10 indifference threshold, so not filed. `p07_stress.do` section C.)

### N2. exact recovery of every subgroup estimand on a noiseless DGP
`p02_sgtruth.do`. All four estimands, pooled and per subgroup, match the design
truth (independent reweighting of `truetau` with the four documented weights) to
<= 1.8e-15.

### N3. the aggregation identity — the sharpest structural check available
`p03_identities.do`, and again with `[aw=w]` in `p07_stress.do` section B.
Because `subgroups()` never touches the regression, the pooled estimand must be
the subgroup-weight-weighted mean of the subgroup estimands, with
W = {sum w, # post cells, # units, # units} for {att, att_it, att_i, att_itime}.
**reldiff 0.000e+00 unweighted; <= 9.1e-16 under aweights.**

### N4. the contrast uses the real off-diagonal covariance, and the zero limit
`p03`/`p04`/`p05`. `se_contrast != sqrt(se0^2+se1^2)` on ordinary designs
(implied rho 0.04-0.20, always inside [-1,1], rebuilt level covariance positive
definite). And in a deliberately BLOCK-SEPARABLE triple-diff design
(`p05_separable.do`: disjoint onset cohorts per subgroup + `fe(i.sg##i.t)` +
`cluster(unit)` nested in sg) the true covariance is exactly zero and wooldid
returns **implied rho = 0.00000000, se_c/sqrt(se0^2+se1^2) - 1 = -8.4e-16**.

### N5. full reimplementation of the homogeneity test at 4 levels
`p04_joint4.do`. Recovered the whole 4x4 level covariance from the 4 level SEs
and 6 pairwise SEs via C_ab = (s_a^2+s_b^2-s_ab^2)/2, rebuilt
W = (Rb)'(R V R')^{-1}(Rb), and got
**p = Ftail(3, e(df_r), W/3) == e(joint_*_allsgssame) to 1e-22** on all four
estimands. chi2(3) does NOT match (4e-10 off), confirming the reference
distribution is F(k-1, df_r). Every pairwise estimate == b_A - b_B exactly
(0.000e+00, 24 rows). Values 1/3/7/99 are non-contiguous and value-labelled and
the rows are keyed by VALUE, correctly. At k=2 the joint p equals the single
pairwise p to 1.1e-16 (`p03`).
Would have caught: wrong df, chi2/F confusion, a chain on the wrong columns, the
`half` column split misfiring, the rank-guard fallback reporting a residual-rank p.

### N6. the triple-diff recipe the help documents
`p05_separable.do`. Each subgroup's estimate equals the estimate from running
wooldid on that subgroup's SUBSAMPLE alone to 2.6e-15 (post) / 4.7e-13 (pre), so
the difference of subgroup estimates really is the triple-diff. (The SEs differ
by 3.0-3.5% — the cluster small-sample factor is computed on different N and G.
Not filed.)

### N7. `suppresssgprimaryeffects` leaves the contrasts alone
`p05_separable.do` S4: `e(pwsgs)` estimate, se and p are **bit-identical**, and
`e(joint_att_allsgssame)` identical to 17 digits.

### N8. per-subgroup pre-trend joint tests and the e(V) splice
`p06_espretrend.do`. Rebuilt W = b_pre' V_pre^-1 b_pre from
`e(b_es_att_s<v>)` / `e(V_es_att_s<v>)`: p = Ftail(npre, df_r, W/npre) matches
`e(joint_horiz_s<v>)` and `e(joint_horiz_i_s<v>)` to 5.6e-17 for both subgroups
and both families, and the pooled twins likewise. The SUBGROUP ES blocks spliced
into `e(V)` match `e(V_es_att_s<v>)` cell for cell (25/25, max 3.5e-18).
`e(wdidsgresults_v)` ES rows == `e(b_es_att_s<v>)` exactly.

### N9. RAGGED subgroup event-study horizons, with an exact LABEL oracle
`p08_ragged.do`. Two subgroups with disjoint onset windows so one reaches a
pre-horizon the other cannot (`numpre` pooled 4, sg0 3, sg1 4). On a noiseless
DGP whose true effect is base(sg)+slope(sg)*h, **every row labelled post<h>
comes back as base+slope*h to <= 1.0e-14 for both subgroups**, so the positional
filter behind the "not all sgs have all relyrs" code path does not slip. The
pooled ES at each horizon equals the weighted mean over the subgroups that HAVE
that horizon (d <= 2.2e-16), including the horizon only one subgroup reaches.

### N10. the documented rank-deficiency screen actually fires
`p07_stress.do` section A, with `cluster(grp)` AT the level of `i`.
`e(wdid_sgdegen)=4`, `e(wdid_sgdegenr)=1.0000`, all four `e(pwsgs)` inference
columns missing, all `allsgssame` p's missing, and the full named caution
printed. Point estimates untouched. Matches the `sgdegen` help text exactly.

### N11. subgroup value edge cases — all behave as documented
`p09_edgevalues.do`.
* a level existing only among NEVER-TREATED units is correctly excluded from
  `sglist`; contrasts among the remaining levels stay exact (N1, rc 0)
* NEGATIVE value on treated rows -> r(411) with a named, informative message
  (no mangled `main_SG-1_att`)
* NON-INTEGER value on treated rows -> r(411), named
* 8-digit values + `att_itime` -> the length caution prints FIRST, then a clean
  `pretreatment_SG48000001_att_itime invalid name` r(198). The same values with
  `att` alone run fine (rc 0). 7-digit values + `att_itime` run fine (rc 0) —
  the 32-character boundary is exactly where the help says it is
* a subgroup with no estimable pre period, and a tiny subgroup (2 of 20 i's),
  both complete without error
* name COLLISION from truncation was reasoned about and is not reachable: a
  collision needs two values agreeing past ~17 digits, which is the excluded
  giant-integer class and fails margins' factor-level cap first

### N12. two code-reading suspicions that dissolved on inspection
* the EFF3 pre-period fallback (`ado:~3175`) omits `& ttre != .` from its
  subpop — redundant: `targetflag` is `_t > 0` (ado:1009) and never-treated rows
  carry `_t == 0` (ado:923), so controls cannot enter
* the pooled pre-trend joint tests use `jointmarginsdof` while the PER-SUBGROUP
  twins (ado:3028, 3045) use `marginsdof` — the two differ only under
  `impvalues` and `fastnose`, both of which are refused alongside
  `jointtests`/`subgroups()`, so it is latent, not live

---

## Explicitly NOT filed
* joint tests returning `.` on a noiseless (gysd=0) design — the covariance is
  numerically singular when the fit is exact; pathological input
* SE differences between an in-pooled subgroup estimate and a subsample run —
  different cluster small-sample factors
* 1-3 ulp movement of the pooled table when `subgroups()` is added on a design
  with `fe()` — inside the 1e-10 threshold
* anything in `known_bydesign.md` §1-§6 I did not measure a contradiction to
