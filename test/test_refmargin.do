* test_refmargin.do - SI-20: the event-study REFERENCE CELL must carry a structurally zero
* margin. Hunts for a non-zero reference-cell margin with EXTERNAL detectors (the suite
* cannot see the margin itself) across binary/continuous treatment, weights, poisson,
* subgroups, esfixedbaseperiod, lattice types, permute()/bootstrap() and one pathological
* shape. Run from test/:  "<stata>" -e do test_refmargin.do   -> expect REFMARGIN_ALL_PASS
*
* THE TWO MECHANISMS (fix commit 1df8613, investigation wooldid_ado_review.md "Findings
* (Claude, 2026-09-05)"): (1) lattice() + contreatcontrols(basic) skipped the baseline dose
* control on the default reghdfe path, so every ame-family row carried the baseline dose
* slope; (2) interactivecontrols() gave the pooled reference cell its own slope column
* (i.targetflag instead of 1.targetflag), exactly collinear with i.t#c.x plus the treated-
* cell slopes, so the reported horizons were solver- and order-dependent and the
* reference-cell margin was non-zero wherever the covariate varies within cohort over
* time, under weights, or on an unbalanced panel.
*
* DETECTORS
* (a) JOINT-WALD IDENTITY. For every case with jointtests and >= 2 pre horizons: b = the
*     "pre*" columns of e(b_es_<k>), V = the matching block of e(V_es_<k>),
*     W = b invsym(V) b', p = Ftail(npre, e(df_r), W/npre) - or chi2tail(npre, W) when
*     e(df_r) is missing (ppmlhdfe posts no e(df_r), so the poisson joint tests are z/chi2
*     based) - compared with e(joint_horiz) [att], e(joint_horiz_i) [att_it],
*     e(joint_Choriz) [ame], e(joint_Choriz_i) [ame_it] and the subgroup twins
*     e(joint_horiz_s#)/e(joint_Choriz_s#) (e() MACROS) against e(b_es_att_s#)/e(b_es_ame_s#).
*     Pass iff |d| <= 1e-8 max(|p1|,|p2|) AND (|d| <= 1e-12 whenever min(p1,p2) < 1e-6).
*     The posted joint test is a contrast(overjoint) over a set that INCLUDES the reference
*     cell, so the two coincide iff the reference-cell margin is a structural zero - the
*     SI-20 symptom. Computed in Mata from the posted matrices, never from displayed digits.
*     The "all pre margins equal" alternative (npre-1 df) is printed as a diagnostic only.
*     CALIBRATION (C0, first case): a plain binary run with no interactivecontrols and no
*     lattice, where the anchor was ALWAYS a structural zero: (a) must pass on both builds;
*     if it fails there the detector premise is wrong. C0c (continuous, basic, no lattice)
*     and C0p (poisson, no ic) are the same parity check for the ame family and the chi2
*     convention.
* (b) ES-PLOT REFERENCE PIN. After makeplots saveplots(refm_es) replace, the saved
*     <name>_ES.dta row with rownamesvar=="reference" has estimate==0 and lb==ub==0.
*     PIN ONLY: wooldid_esplots plants that 0 unconditionally, so (b) cannot detect either
*     defect on any build; it is NOT counted toward non-vacuity.
* (c) THE UNCONDITIONALSE TWIN. Where the models coincide (no fe(), or fe(i.var) only -
*     never bare fe(var): under unconditionalse fe() is passed bare to reg, a separate
*     out-of-scope defect) the point estimates (estimate column of e(wdidmainresults),
*     e(b_es_att), e(b_es_att_it), e(b_es_ame), e(b_es_ame_it)) of the default reghdfe run
*     equal the unconditionalse (plain reg) run's, row-matched by FULL name, with Stata's
*     reldif(x,y) = |x-y|/(|y|+1) <= 1e-8 on every matched cell; the matched count is
*     printed and zero matches is a failure. TOLERANCE: wooldid's regtol(9) default runs
*     reghdfe's iterative solver at 1e-9, so a reghdfe-vs-reg identity cannot be asserted
*     tighter than ~1e-8 in general; the first draft used 1e-10 and the fixed build failed
*     it on D1 (dta_cont lattice(all 3), R2 = 1.0000) at max reldif 3.415e-10 on 2 of 12
*     cells while every other case sat at <= 4.6e-13. The SI-20 defects move these cells by
*     0.09-2 (four to nine orders above 1e-8): the default path differs from the reg path
*     (0.09-0.10 on binary ic; the base dose slope on lattice+basic), so (c) detects both
*     mechanisms on the pristine build. The B8 twin (att att_it) carries safety(off)
*     because the ado refuses unconditionalse with a non-att estimand type (rc 198)
*     otherwise; safety(off) only suppresses cautions/screens there (the _t construction is
*     identical and the e(b) omission scan it skips feeds ame scrubbing only).
* (d) DGP-TRUTH CHECKS on the investigation's fixtures, where the defect produces grossly
*     wrong printed rows - usable also under bootstrap()/permute() where (a) does not exist:
*     dta_cont lattice(all 3)/(all 2) + basic (true pre ame 0, DGP base slope 2): every
*     ES_ame pre row and pretreatment:ame within [-0.3, 0.3] (pristine prints ~2.0); the
*     trending-covariate binary fixture rebuilt inline from the investigation's gen2.do
*     (true pre-effects 0, true ATT 1.0-1.5): every ES_att pre row and pretreatment:att
*     within [-0.5, 0.5] and 0.5 < main:att < 2.5 (pristine prints pre -2.96/-2.10/-0.81,
*     att 3.92).
* (e) REPLICATION ENGINES. An inline NO-EFFECT continuous panel (mkperm0: 30 units x 8
*     periods x 10 obs per unit-period, onsets 4/never, dose ~ U(0,1) for every unit,
*     y = 2 dose + unit FE + period FE + N(0, .02), no treatment effect) under lattice(all 2)
*     + basic, cluster(i), permute(shuffle) permreps(20): rc 0 AND e(perm_offcentre) EMPTY.
*     The guard (wooldid_permtail) fires when the central 95% band of the raw permutation
*     draws excludes 0: on the pristine build every ame draw sits at the baseline dose
*     slope (~2) under any reassignment, so the guard fires; on the fixed build the draws
*     are mean-zero. Plus (d) rows on the observed run: ES_ame pre1/pre2, pretreatment:ame
*     AND main:ame within [-0.3, 0.3] (true 0; pristine ~2). dta_cont is NOT used here: its
*     DGP carries treated*intensity (a true effect of 1 x dose on top of the base slope 2),
*     so the sharp null is false under reassigned onsets and the ame permutation
*     distribution is legitimately off-centre on ANY build (the first draft flagged
*     ES_ame_pre1 on the fixed build for exactly that reason). dta_cont stays under
*     bootstrap(standard) reps(10): rc 0, finite summaries, (d). The same two engines on
*     the gen2 fixture with interactivecontrols(x): rc 0, finite summaries, (d).
* CASES: C0 C0c C0p (calibration/parity), B1-B8 (binary interactivecontrols on
* dta_s_covars: unit-level robust; the FLAGSHIP cohort-level fe(i) cluster(i); aw; pw;
* poisson; subgroups; esfixedbaseperiod; att att_it), L1-L10 (continuous lattice+basic on
* continuous_control_dense: baseTreated 2, all 3, byCohort 2, the pre-made lattice
* variable, aw, pw, subgroups, esfixedbaseperiod, poisson (parity class: (a) must hold on
* both builds because poisson always absorbed the baseline dose control), contreatpoly(2)),
* D1-D2 (the truth fixtures), E1-E4 (permute/bootstrap), P1 (pathological shape:
* dta_path_nopre with a constructed time-varying covariate under interactivecontrols():
* rc 0 and (a) holds, OR an informative error rc in {198, 322, 459, 2001}).
* DEVIATIONS FROM THE TASK BRIEF, deliberate: permute() requires cluster() (wooldid_permctl
* refuses a run without it), so E1/E3 use cluster(i) rather than robust; the poisson cases
* carry esprelength(2) because a single pre horizon produces no joint test at all; the (e)
* permute case runs on the inline no-effect fixture instead of dta_cont (see (e)); the (c)
* tolerance is 1e-8, not 1e-10 (see (c)); the B8 twin carries safety(off) (see (c)).
* NON-VACUITY: live evaluations of (a), (c) and (d) are counted and the suite fails if any
* of the three evaluated zero times.
* MEASURED 2026-09-05 (mirror layouts, this file unchanged apart from this header; both
* runs in parallel on the same machine; live evaluations on both builds: (a) 43, (c) 304,
* (d) 25 over 28 cases; rc 0 everywhere on both builds):
*   build 310fa61 (pristine): fails 52 checks in 1368 s. Failing checks (52 over 24
*     cases; only C0, C0c, C0p and the poisson-lattice parity case L9 pass, as required):
*     B1 (a)(c); B2 (a)(c) [the flagship: posted joint p .0510 vs pre=0 Wald .7014; 12 of
*     12 cells, max reldif 5.94e-2]; B3 (a)(c); B4 (a)(c); B5 (a); B6 (a) pooled, sg1, sg2;
*     B7 (a)(c); B8 (a) att, (a) att_it, (c); L1-L6 (a) ame, (c) [10 of 20 cells = the ame
*     family, max reldif .89-.95 = the base dose slope]; L7 (a) ame pooled, sg1, sg2, sg3,
*     (c); L8 (a)(c); L10 (a)(c); D1 (d) [pre ame 1.97/2.01/1.99], (c) [6 of 12, max 1.99];
*     D2 (a) [1.88e-38 vs 4.03e-38], (d) x2 [pre -2.96/-2.10/-.81, att 3.92], (c) [12 of
*     12, max 2.88]; E1 (d) [ame rows 2.006/2.017/2.012/1.964], (e) [off-centre: main_ame
*     pretreatment_ame ES_ame_pre2 ES_ame_pre1 ES_ame_post0 ES_ame_post1]; E2 (d); E3 (d)
*     x2; E4 (d) x2; P1 (a) [.0096 vs .0029].
*   build 1df8613 (fixed):    passes (the ALL_PASS sentinel) in 1426 s (max (c) reldif
*     3.415e-10 on D1, <= 4.6e-13 elsewhere; E1 e(perm_offcentre) empty, ame rows
*     .006/.017/.011/-.036).
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 20260905
set sortseed 20260905
do "../wooldidinstall/wooldid.ado"
local fails 0
local na 0
local nc 0
local nd 0
local ncase 0

* ============================== Mata kernels ==============================
mata:
// pre = 0 Wald from the posted event-study matrices of family <key>; results in
// refm_p (recomputed p), refm_np (pre count), refm_rk (rank used as df1), refm_W,
// refm_peq (diagnostic: all pre margins equal), refm_st (0 ok, 1 matrices absent,
// 2 no pre columns, 3 missing/degenerate block, 5 dimension mismatch)
void refm_wald(string scalar key, real scalar dfr)
{
  real rowvector b, bp
  real matrix V, Vp, Vi, C, Ci
  real colvector idx, D
  string matrix cs
  real scalar np, W, p, i, rk, Weq, peq
  st_numscalar("refm_p", .)
  st_numscalar("refm_np", 0)
  st_numscalar("refm_rk", .)
  st_numscalar("refm_W", .)
  st_numscalar("refm_peq", .)
  st_numscalar("refm_st", 1)
  b = st_matrix("e(b_es_" + key + ")")
  V = st_matrix("e(V_es_" + key + ")")
  if (rows(b) == 0 | rows(V) == 0) return
  if (cols(V) != cols(b) | rows(V) != cols(b)) {
    st_numscalar("refm_st", 5)
    return
  }
  cs = st_matrixcolstripe("e(b_es_" + key + ")")
  idx = J(0, 1, .)
  for (i = 1; i <= cols(b); i++) {
    if (substr(cs[i, 2], 1, 3) == "pre") idx = idx \ i
  }
  np = rows(idx)
  st_numscalar("refm_np", np)
  if (np == 0) {
    st_numscalar("refm_st", 2)
    return
  }
  bp = b[1, idx]
  Vp = V[idx, idx]
  Vp = (Vp + Vp') / 2
  if (hasmissing(bp) | hasmissing(Vp)) {
    st_numscalar("refm_st", 3)
    return
  }
  Vi = invsym(Vp)
  rk = np - sum(diagonal(Vi) :== 0)
  if (rk < 1) {
    st_numscalar("refm_st", 3)
    return
  }
  W = bp * Vi * bp'
  if (dfr < .) p = Ftail(rk, dfr, W / rk)
  else p = chi2tail(rk, W)
  st_numscalar("refm_p", p)
  st_numscalar("refm_W", W)
  st_numscalar("refm_rk", rk)
  st_numscalar("refm_st", 0)
  if (np >= 2) {
    C = (I(np - 1), J(np - 1, 1, -1))
    D = C * bp'
    Ci = invsym(C * Vp * C')
    Weq = D' * Ci * D
    if (dfr < .) peq = Ftail(np - 1, dfr, Weq / (np - 1))
    else peq = chi2tail(np - 1, Weq)
    st_numscalar("refm_peq", peq)
  }
}
// name-matched comparison of two Stata matrices: byrow=1 compares column 1 across rows
// matched by "eq:name"; byrow=0 compares row 1 across columns matched the same way.
// reldif = |x-y|/(|y|+1). results in refm_cn (matched), refm_cbad (> tol), refm_cmax
void refm_cmp(string scalar A, string scalar B, real scalar byrow, real scalar tol)
{
  real matrix MA, MB
  string matrix sa, sb
  string colvector na, nb
  real scalar i, j, n, bad, mx, x, y, d
  MA = st_matrix(A)
  MB = st_matrix(B)
  n = 0
  bad = 0
  mx = 0
  if (rows(MA) == 0 | rows(MB) == 0) {
    st_numscalar("refm_cn", 0)
    st_numscalar("refm_cbad", 0)
    st_numscalar("refm_cmax", .)
    return
  }
  if (byrow) {
    sa = st_matrixrowstripe(A)
    sb = st_matrixrowstripe(B)
  }
  else {
    sa = st_matrixcolstripe(A)
    sb = st_matrixcolstripe(B)
  }
  na = sa[., 1] :+ ":" :+ sa[., 2]
  nb = sb[., 1] :+ ":" :+ sb[., 2]
  for (i = 1; i <= rows(na); i++) {
    for (j = 1; j <= rows(nb); j++) {
      if (na[i] == nb[j]) {
        x = (byrow ? MA[i, 1] : MA[1, i])
        y = (byrow ? MB[j, 1] : MB[1, j])
        if (x < . & y < .) {
          d = abs(x - y) / (abs(y) + 1)
          n++
          if (d > tol) bad++
          if (d > mx) mx = d
        }
        break
      }
    }
  }
  st_numscalar("refm_cn", n)
  st_numscalar("refm_cbad", bad)
  st_numscalar("refm_cmax", mx)
}
// column 1 of the row of matrix M whose "eq:name" equals rname -> refm_v (. if absent)
void refm_getrow(string scalar M, string scalar rname)
{
  real matrix X
  string matrix s
  real scalar i
  st_numscalar("refm_v", .)
  X = st_matrix(M)
  if (rows(X) == 0) return
  s = st_matrixrowstripe(M)
  for (i = 1; i <= rows(s); i++) {
    if (s[i, 1] + ":" + s[i, 2] == rname) {
      st_numscalar("refm_v", X[i, 1])
      return
    }
  }
}
end

* ============================== helpers ==============================
* (a): recompute the pre = 0 Wald for family <key> and compare with the posted joint p
cap program drop refm_a
program define refm_a, rclass
  syntax , label(string) key(string) [posted(string)]
  return scalar bad = 0
  return scalar n = 0
  if "`posted'" == "" | "`posted'" == "." {
    di as text "`label' (a) key=`key': no posted joint p (not evaluated)"
    exit
  }
  local dfr = e(df_r)
  mata: refm_wald("`key'", `dfr')
  if refm_st != 0 {
    di as error "`label' (a) FAIL: key=`key' has a posted joint p (`posted') but the identity is not computable from e(b_es_`key')/e(V_es_`key') (status " refm_st ", npre " refm_np ")"
    return scalar bad = 1
    exit
  }
  scalar refm_pp = `posted'
  scalar refm_dd = abs(refm_p - refm_pp)
  local ok 1
  if refm_dd > 1e-8 * max(abs(refm_p), abs(refm_pp)) local ok 0
  if min(refm_p, refm_pp) < 1e-6 & refm_dd > 1e-12 local ok 0
  local rkflag
  if refm_rk != refm_np local rkflag " RANK-DEFICIENT pre block (df1 = rank)"
  di as text "`label' (a) key=`key': npre=" refm_np " df1=" refm_rk " df_r=`dfr' W=" %21.15g refm_W " p_recomputed=" %21.15g refm_p " p_posted=" %21.15g refm_pp " |diff|=" %10.3e refm_dd " [diag p_allequal=" %21.15g refm_peq "]`rkflag'"
  return scalar n = 1
  if !`ok' {
    di as error "`label' (a) FAIL: key=`key' posted joint p " %21.15g refm_pp " vs pre=0 Wald recomputed from e(b_es_`key')/e(V_es_`key') " %21.15g refm_p " (npre " refm_np ", W " %21.15g refm_W ", |diff| " %10.3e refm_dd ")"
    return scalar bad = 1
  }
end

* (c): compare the stored default-run matrices (refm_M0, refm_E0_*) with the current e().
* tol 1e-8: reghdfe runs at wooldid's regtol(9) = 1e-9, so the reghdfe-vs-reg identity
* holds only to ~1e-8 (fixed build: D1 max 3.415e-10, all others <= 4.6e-13; the defect
* moves cells by 0.09-2)
cap program drop refm_c
program define refm_c, rclass
  syntax , label(string) [tol(real 1e-8)]
  return scalar bad = 0
  return scalar n = 0
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    di as error "`label' (c) FAIL: the unconditionalse twin posted no e(wdidmainresults)"
    return scalar bad = 1
    exit
  }
  local n 0
  local nb 0
  local mx 0
  cap matrix drop refm_M1
  matrix refm_M1 = e(wdidmainresults)
  mata: refm_cmp("refm_M0", "refm_M1", 1, `tol')
  local n = `n' + refm_cn
  local nb = `nb' + refm_cbad
  local mx = max(`mx', refm_cmax)
  di as text "`label' (c) e(wdidmainresults) estimates: matched " refm_cn " rows, " refm_cbad " beyond reldif `tol', max reldif " %10.3e refm_cmax
  foreach k in att att_it ame ame_it {
    cap confirm matrix refm_E0_`k'
    if _rc != 0 continue
    cap confirm matrix e(b_es_`k')
    if _rc != 0 {
      di as error "`label' (c) FAIL: e(b_es_`k') is posted by the default run but absent on the unconditionalse twin"
      local ++nb
      continue
    }
    cap matrix drop refm_E1
    matrix refm_E1 = e(b_es_`k')
    mata: refm_cmp("refm_E0_`k'", "refm_E1", 0, `tol')
    local n = `n' + refm_cn
    local nb = `nb' + refm_cbad
    local mx = max(`mx', refm_cmax)
    di as text "`label' (c) e(b_es_`k'): matched " refm_cn " cols, " refm_cbad " beyond reldif `tol', max reldif " %10.3e refm_cmax
  }
  return scalar n = `n'
  if `n' == 0 {
    di as error "`label' (c) FAIL: zero cells matched between the default run and the unconditionalse twin (vacuous)"
    return scalar bad = 1
    exit
  }
  if `nb' > 0 {
    di as error "`label' (c) FAIL: `nb' of `n' name-matched cells differ between the default (reghdfe) run and the unconditionalse (reg) twin beyond reldif `tol' (max reldif " %10.3e `mx' ")"
    return scalar bad = 1
  }
end

* (d): named rows of e(wdidmainresults) within [lo, hi]; a missing required row fails
cap program drop refm_d
program define refm_d, rclass
  syntax , label(string) rows(string) lo(real) hi(real)
  return scalar bad = 0
  return scalar n = 0
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    di as error "`label' (d) FAIL: no e(wdidmainresults) to read"
    return scalar bad = 1
    exit
  }
  cap matrix drop refm_M1
  matrix refm_M1 = e(wdidmainresults)
  local n 0
  local badl
  foreach rw of local rows {
    mata: refm_getrow("refm_M1", "`rw'")
    if refm_v >= . {
      local badl `badl' `rw'=MISSING
      continue
    }
    local ++n
    di as text "`label' (d) `rw' = " %21.15g refm_v "  bounds [`lo', `hi']"
    if refm_v < `lo' | refm_v > `hi' local badl `badl' `rw'=`=string(refm_v, "%12.6g")'
  }
  return scalar n = `n'
  if "`badl'" != "" {
    di as error "`label' (d) FAIL: rows outside [`lo', `hi'] or missing: `badl'"
    return scalar bad = 1
  }
end

* (b) pin: the saved ES plot data's reference row is 0/0/0
cap program drop refm_b
program define refm_b, rclass
  syntax , label(string) file(string)
  return scalar bad = 0
  cap confirm file "`file'"
  if _rc != 0 {
    di as error "`label' (b) FAIL: `file' was not written"
    return scalar bad = 1
    exit
  }
  cap frame drop refm_fr
  frame create refm_fr
  frame refm_fr: use "`file'", clear
  frame refm_fr: qui count if rownamesvar == "reference"
  local nref = r(N)
  frame refm_fr: qui count if rownamesvar == "reference" & estimate == 0 & lb_95ci == 0 & ub_95ci == 0
  local nok = r(N)
  frame refm_fr: qui count
  local nall = r(N)
  frame drop refm_fr
  di as text "`label' (b) `file': `nall' rows, `nref' reference row(s), `nok' pinned at 0/0/0"
  if `nref' == 0 | `nok' != `nref' {
    di as error "`label' (b) FAIL: `nref' reference row(s), `nok' with estimate==lb==ub==0"
    return scalar bad = 1
  }
end

* finite summaries: every row has a non-missing estimate and a non-missing column <col>
cap program drop refm_fin
program define refm_fin, rclass
  syntax , label(string) col(integer)
  return scalar bad = 0
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    di as error "`label' finite FAIL: no e(wdidmainresults)"
    return scalar bad = 1
    exit
  }
  cap matrix drop refm_M1
  matrix refm_M1 = e(wdidmainresults)
  local R = rowsof(refm_M1)
  local nm 0
  forval r = 1/`R' {
    if refm_M1[`r', 1] >= . | refm_M1[`r', `col'] >= . local ++nm
  }
  di as text "`label' finite: `R' rows, `nm' with a missing estimate or column `col'"
  if `R' == 0 | `nm' > 0 {
    di as error "`label' finite FAIL: `nm' of `R' rows carry a missing estimate or column `col'"
    return scalar bad = 1
  }
end

* data loaders / fixtures
cap program drop mkcoh
program define mkcoh
  qui gen coh = cond(ttre == 5, 1, cond(ttre == 9, 2, 3))
end
* the investigation's trending-covariate fixture (refute/gen2.do, verbatim generator):
* 120 units, onsets 5/8/never (40 each), 10 periods, untreated-state covariate slope 0.3t,
* true ATT = 1.0 + 0.2(x - xb), true pre-effects 0
cap program drop mksim2
program define mksim2
  clear
  set seed 777
  set sortseed 777
  qui set obs 120
  qui gen i = _n
  qui gen ttre = cond(i<=40, 5, cond(i<=80, 8, .))
  qui gen coh = cond(i<=40, 1, cond(i<=80, 2, 3))
  qui gen double ai = rnormal()
  qui gen double xb = rnormal()
  qui expand 10
  bysort i: gen t = _n
  qui gen double x = xb + 0.3*t + 0.5*rnormal()
  qui gen byte post = t >= ttre & ttre < .
  qui gen double y = 1 + ai + 0.2*t + 0.3*t*x + 1.0*post + 0.2*(x - xb)*post + rnormal(0,0.5)
  qui gen double yp = exp(0.1*y)
end
* the (e) permute fixture: a NO-EFFECT continuous panel with a non-zero baseline dose
* slope. 30 units (15 treated at t = 4, 15 never), 8 periods, 10 obs per unit-period,
* dose ~ U(0,1) for every unit (never-treated included), y = 2 dose + unit FE + period FE
* + N(0, .02). The true ame is 0 at every horizon; under the SI-20 lattice+basic defect
* every ame row carries the base slope 2 on the observed run and on every permutation draw
cap program drop mkperm0
program define mkperm0
  clear
  set seed 31415
  set sortseed 31415
  qui set obs 30
  qui gen i = _n
  qui gen ttre = cond(i <= 15, 4, .)
  qui gen double ai = rnormal()
  qui expand 8
  bysort i: gen t = _n
  qui gen double pt = 0.4*t - 0.05*t^2
  qui expand 10
  qui gen double dose = runiform()
  qui gen double y = 2*dose + ai + pt + rnormal(0, 0.02)
end
cap program drop mkxtv
program define mkxtv
  set seed 4444
  set sortseed 4444
  qui gen double xtv = 0.3*t + rnormal()
end
cap program drop refm_load
program define refm_load
  args data prep
  ereturn clear
  if substr("`data'", 1, 3) == "mk:" {
    local p = substr("`data'", 4, .)
    `p'
  }
  else use "`data'", clear
  if "`prep'" != "" {
    `prep'
  }
end

* the case driver: one wooldid call, the rc check, then every detector that applies;
* returns r(bad) (failed checks), r(na)/r(nc)/r(nd) (live evaluations), r(rc), r(secs)
cap program drop refm_case
program define refm_case, rclass
  syntax , label(string) data(string) cmd(string asis) [prep(string) akeys(string) sgkeys(string) uncond twin(string asis) drows(string) dlo(real 0) dhi(real 0) drows2(string) dlo2(real 0) dhi2(real 0) fin(integer 0) offcentre pin(string) rcok(numlist) seed(integer 0)]
  local bad 0
  local na 0
  local nc 0
  local nd 0
  local otherbad 0
  local abad 0
  local cbad 0
  local dbad 0
  foreach mm in refm_M0 refm_M1 refm_E0_att refm_E0_att_it refm_E0_ame refm_E0_ame_it refm_E1 {
    cap matrix drop `mm'
  }
  refm_load "`data'" "`prep'"
  if `seed' > 0 {
    set seed `seed'
    set sortseed `seed'
  }
  local dcmd `cmd'
  if "`pin'" != "" local dcmd `dcmd' makeplots saveplots(`pin') replace
  di as text _n "==== `label': wooldid `dcmd'"
  timer clear 9
  timer on 9
  cap noi wooldid `dcmd'
  local rc = _rc
  timer off 9
  qui timer list 9
  local secs = r(t9)
  return scalar rc = `rc'
  * --- rc
  local rcpass = (`rc' == 0)
  local rcobs 0
  if "`rcok'" != "" {
    foreach k of numlist `rcok' {
      if `rc' == `k' local rcobs 1
    }
  }
  if `rcpass' di as text "`label' rc=0"
  else if `rcobs' di as text "`label' rc=`rc': an informative refusal on this pathological shape (accepted: rcok `rcok')"
  else {
    local rcalt
    if "`rcok'" != "" local rcalt " or one of `rcok'"
    di as error "`label' FAIL: rc=`rc' (want 0`rcalt')"
    local ++otherbad
  }
  if `rcpass' {
    * --- (a) on the pooled families
    foreach k of local akeys {
      local pj
      if "`k'" == "att" local pj `=e(joint_horiz)'
      if "`k'" == "att_it" local pj `=e(joint_horiz_i)'
      if "`k'" == "ame" local pj `=e(joint_Choriz)'
      if "`k'" == "ame_it" local pj `=e(joint_Choriz_i)'
      refm_a , label(`label') key(`k') posted(`pj')
      local abad = `abad' + r(bad)
      local na = `na' + r(n)
    }
    * --- (a) on the subgroup twins (e() macros)
    if "`sgkeys'" != "" {
      local sgl `e(wdid_sglist)'
      if "`sgl'" == "" {
        di as error "`label' (a) FAIL: sgkeys requested but e(wdid_sglist) is empty"
        local ++abad
      }
      foreach sg of local sgl {
        foreach k of local sgkeys {
          local pj
          if "`k'" == "att" local pj `e(joint_horiz_s`sg')'
          if "`k'" == "att_it" local pj `e(joint_horiz_i_s`sg')'
          if "`k'" == "ame" local pj `e(joint_Choriz_s`sg')'
          if "`k'" == "ame_it" local pj `e(joint_Choriz_i_s`sg')'
          refm_a , label(`label'_sg`sg') key(`k'_s`sg') posted(`pj')
          local abad = `abad' + r(bad)
          local na = `na' + r(n)
        }
      }
    }
    * --- (d)
    if "`drows'" != "" {
      refm_d , label(`label') rows(`drows') lo(`dlo') hi(`dhi')
      local dbad = `dbad' + r(bad)
      local nd = `nd' + r(n)
    }
    if "`drows2'" != "" {
      refm_d , label(`label') rows(`drows2') lo(`dlo2') hi(`dhi2')
      local dbad = `dbad' + r(bad)
      local nd = `nd' + r(n)
    }
    * --- finite summaries / off-centre guard (replication engines)
    if `fin' > 0 {
      refm_fin , label(`label') col(`fin')
      local otherbad = `otherbad' + r(bad)
    }
    if "`offcentre'" != "" {
      local oc `e(perm_offcentre)'
      di as text "`label' (e) e(perm_offcentre)=[`oc'] (want empty)"
      if "`oc'" != "" {
        di as error "`label' (e) FAIL: the permutation distribution is off-centre for: `oc' - the reference cell carries a non-zero margin under reassigned onsets"
        local ++otherbad
      }
    }
    * --- store the default run's estimates for (c), then the (b) pin
    matrix refm_M0 = e(wdidmainresults)
    foreach k in att att_it ame ame_it {
      cap confirm matrix e(b_es_`k')
      if _rc == 0 matrix refm_E0_`k' = e(b_es_`k')
    }
    if "`pin'" != "" {
      refm_b , label(`label') file(`pin'_ES.dta)
      local otherbad = `otherbad' + r(bad)
    }
    * --- (c) the unconditionalse twin
    if "`uncond'" != "" | `"`twin'"' != "" {
      if `"`twin'"' != "" local tcmd `twin'
      else {
        local tcmd : subinstr local cmd "jointtests" "", all
        local tcmd `tcmd' unconditionalse
      }
      refm_load "`data'" "`prep'"
      di as text _n "==== `label' twin: wooldid `tcmd'"
      timer clear 9
      timer on 9
      cap noi wooldid `tcmd'
      local rct = _rc
      timer off 9
      qui timer list 9
      local secs = `secs' + r(t9)
      if `rct' != 0 {
        di as error "`label' (c) FAIL: the unconditionalse twin returned rc=`rct'"
        local ++cbad
      }
      else {
        refm_c , label(`label')
        local cbad = `cbad' + r(bad)
        local nc = `nc' + r(n)
      }
    }
  }
  local bad = `otherbad' + `abad' + `cbad' + `dbad'
  di as text "REFM_CASE `label' rc=`rc' secs=" %6.1f `secs' " a_eval=`na' a_bad=`abad' c_eval=`nc' c_bad=`cbad' d_eval=`nd' d_bad=`dbad' other_bad=`otherbad'"
  return scalar bad = `bad'
  return scalar na = `na'
  return scalar nc = `nc'
  return scalar nd = `nd'
  return scalar secs = `secs'
end

* ============================== cases ==============================
local BIN data/dta_binary.dta
local CONT data/dta_cont.dta
local COV suite/data/dta_s_covars.dta
local CCD suite/data/continuous_control_dense.dta
local LOPT contreat(dose_common) contreatcontrols(basic) esprelength(2) espostlength(2) robust jointtests safety(off)

* ---- C0: CALIBRATION. Binary, no interactivecontrols, no lattice: the anchor was always a
* structural zero, so (a) and (c) must pass on BOTH builds. If (a) fails here the detector
* premise (overjoint over a set including the reference cell == pre = 0 Wald) is wrong.
refm_case , label(C0) data(`BIN') cmd(y i t ttre, att robust esprelength(3) espostlength(2) jointtests) akeys(att) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- C0c: continuous parity - basic WITHOUT a lattice absorbs the baseline dose control on
* both builds, so the ame-side identity must hold on both
refm_case , label(C0c) data(`CONT') cmd(y i t ttre, att contreat(intensity) contreatcontrols(basic) robust esprelength(2) espostlength(2) jointtests) akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- C0p: poisson parity (no ic, no lattice) - pins the chi2 convention of (a)
refm_case , label(C0p) data(`BIN') cmd(ycount i t ttre, att poisson robust esprelength(2) espostlength(1) jointtests) akeys(att)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B1: binary interactivecontrols, unit-level i, robust (+ the (b) pin)
refm_case , label(B1) data(`COV') cmd(y i t ttre, att robust interactivecontrols(x1) esprelength(2) espostlength(2) jointtests) akeys(att) uncond pin(refm_es)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B2: THE FLAGSHIP. Cohort-level i with unit fixed effects fe(i) and cluster(i):
* pristine joint p 0.0510 vs pre = 0 Wald 0.70. (c) via an explicit fe(i.i) twin (fe(i.var)
* survives the bare-fe reg handoff as a factor term).
refm_case , label(B2) data(`COV') prep(mkcoh) cmd(y coh t ttre, att cluster(i) fe(i) interactivecontrols(x1) esprelength(2) espostlength(3) jointtests) akeys(att) twin(y coh t ttre, att cluster(i) fe(i.i) interactivecontrols(x1) esprelength(2) espostlength(3) unconditionalse)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B3/B4: analytic and probability weights
refm_case , label(B3) data(`COV') cmd(y i t ttre [aw=w], att robust interactivecontrols(x1) esprelength(2) espostlength(2) jointtests) akeys(att) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase
refm_case , label(B4) data(`COV') cmd(y i t ttre [pw=w], att robust interactivecontrols(x1) esprelength(2) espostlength(2) jointtests) akeys(att) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B5: poisson (ppmlhdfe path), cohort-level with fe(i)
refm_case , label(B5) data(`COV') prep(mkcoh) cmd(ycount coh t ttre, att poisson robust fe(i) interactivecontrols(x1) esprelength(2) espostlength(1) jointtests) akeys(att)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B6: subgroups (the subgroup joint twins are e() macros); cluster(i) at the unit
* level trips the SI-A1 screen, hence safety(off)
refm_case , label(B6) data(`COV') cmd(y i t ttre, att cluster(i) subgroups(sg) interactivecontrols(x1) esprelength(2) espostlength(2) jointtests safety(off)) akeys(att) sgkeys(att)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B7: esfixedbaseperiod esrelativeto(-2)
refm_case , label(B7) data(`COV') cmd(y i t ttre, att robust esfixedbaseperiod esrelativeto(-2) interactivecontrols(x1) esprelength(2) espostlength(2) jointtests) akeys(att) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- B8: att and att_it together (e(joint_horiz_i) / e(b_es_att_it)). The twin needs
* safety(off): the ado refuses unconditionalse with att_it otherwise ("unconditional
* standard errors may be unreliable with non-att estimand types", rc 198); safety(off)
* suppresses cautions/screens only, so the point estimates stay a fair twin.
refm_case , label(B8) data(`COV') cmd(y i t ttre, att att_it robust interactivecontrols(x1) esprelength(2) espostlength(2) jointtests) akeys(att att_it) twin(y i t ttre, att att_it robust interactivecontrols(x1) esprelength(2) espostlength(2) unconditionalse safety(off))
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- L1-L4: lattice types on continuous_control_dense (L1 carries the ame (b) pin)
refm_case , label(L1) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(baseTreated 2) `LOPT') akeys(att ame) uncond pin(refm_es2)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase
refm_case , label(L2) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(all 3) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase
refm_case , label(L3) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(byCohort 2) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase
refm_case , label(L4) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(custom_lattice_common) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- L5/L6: weights
refm_case , label(L5) data(`CCD') cmd(y_cont_common cohort t ttre [aw=analytic_weight], att lattice(baseTreated 2) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase
refm_case , label(L6) data(`CCD') cmd(y_cont_common cohort t ttre [pw=probability_weight], att lattice(baseTreated 2) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- L7: subgroups (esprelength(2) so the subgroup joint twins exist; espostlength(1))
refm_case , label(L7) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(baseTreated 2) subgroups(sg_three) contreat(dose_common) contreatcontrols(basic) esprelength(2) espostlength(1) robust jointtests safety(off)) akeys(att ame) sgkeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- L8: esfixedbaseperiod
refm_case , label(L8) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(baseTreated 2) esfixedbaseperiod esrelativeto(-2) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- L9: poisson lattice - PARITY class: the guard's poisson disjunct always absorbed the
* baseline dose control, so (a) must hold on both builds
refm_case , label(L9) data(`CCD') cmd(y_count cohort t ttre, att poisson lattice(all 2) contreat(dose_common) contreatcontrols(basic) esprelength(2) espostlength(1) robust jointtests safety(off)) akeys(att ame)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- L10: contreatpoly(2)
refm_case , label(L10) data(`CCD') cmd(y_cont_common cohort t ttre, att lattice(all 2) contreatpoly(2) `LOPT') akeys(att ame) uncond
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- D1: DGP truth on dta_cont (true pre ame 0, base slope 2): pristine prints ~2.0
refm_case , label(D1) data(`CONT') cmd(y i t ttre, att contreat(intensity) contreatcontrols(basic) lattice(all 3) robust esprelength(2) espostlength(2) jointtests) akeys(att ame) uncond drows(ES_ame:pre1 ES_ame:pre2 pretreatment:ame) dlo(-0.3) dhi(0.3)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- D2: DGP truth on the trending-covariate fixture (the investigation's exact spec):
* pristine prints pre -2.96/-2.10/-0.81, att 3.92, joint p 1.9e-38
refm_case , label(D2) data(mk:mksim2) cmd(y i t ttre, att robust jointtests esfixedbaseperiod esrelativeto(-1) esprelength(3) espostlength(2) interactivecontrols(x)) akeys(att) uncond drows(ES_att:pre4 ES_att:pre3 ES_att:pre2 pretreatment:att) dlo(-0.5) dhi(0.5) drows2(main:att) dlo2(0.5) dhi2(2.5)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- E1: permute on the no-effect fixture, lattice(all 2)+basic: the off-centre guard
* must NOT fire and every ame row (pre AND main) is within 0.3 of the true 0 (pristine:
* the guard fires and the rows sit at the base slope ~2). permute() requires cluster();
* i is the unit and ttre is unit-constant. The seed is set right before the call.
refm_case , label(E1) data(mk:mkperm0) seed(20260905) cmd(y i t ttre, att contreat(dose) contreatcontrols(basic) lattice(all 2) cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20)) fin(4) offcentre drows(ES_ame:pre1 ES_ame:pre2 pretreatment:ame main:ame) dlo(-0.3) dhi(0.3)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- E2: bootstrap(standard) on the same model (att/ame only, so observation resampling
* under robust is permitted)
refm_case , label(E2) data(`CONT') cmd(y i t ttre, att contreat(intensity) contreatcontrols(basic) lattice(all 2) robust esprelength(2) espostlength(2) bootstrap(standard) reps(10) rseed(20260905)) fin(2) drows(ES_ame:pre1 ES_ame:pre2 pretreatment:ame) dlo(-0.3) dhi(0.3)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- E3/E4: the binary interactivecontrols fixture under permute and bootstrap
refm_case , label(E3) data(mk:mksim2) seed(20260905) cmd(y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(3) espostlength(2) interactivecontrols(x) permute(shuffle) permreps(20)) fin(4) drows(ES_att:pre4 ES_att:pre3 ES_att:pre2 pretreatment:att) dlo(-0.5) dhi(0.5) drows2(main:att) dlo2(0.5) dhi2(2.5)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase
refm_case , label(E4) data(mk:mksim2) cmd(y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(3) espostlength(2) interactivecontrols(x) bootstrap(standard) reps(10) rseed(20260905)) fin(2) drows(ES_att:pre4 ES_att:pre3 ES_att:pre2 pretreatment:att) dlo(-0.5) dhi(0.5) drows2(main:att) dlo2(0.5) dhi2(2.5)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ---- P1: pathological shape (dta_path_nopre: one cohort treated in the first period, so
* no baseline exists for it) with a constructed time-varying covariate under
* interactivecontrols(): correct-or-observable-failure
refm_case , label(P1) data(data/dta_path_nopre.dta) prep(mkxtv) cmd(y i t ttre, att robust interactivecontrols(xtv) esprelength(2) espostlength(1) jointtests) akeys(att) rcok(198 322 459 2001)
local fails = `fails' + r(bad)
local na = `na' + r(na)
local nc = `nc' + r(nc)
local nd = `nd' + r(nd)
local ++ncase

* ============================== non-vacuity, cleanup, sentinel ==============================
di as text _n "REFM_TOTALS cases=`ncase' a_eval=`na' c_eval=`nc' d_eval=`nd' failed_checks=`fails'"
if `na' == 0 {
  di as error "NONVACUITY FAIL: detector (a) evaluated zero times"
  local ++fails
}
if `nc' == 0 {
  di as error "NONVACUITY FAIL: detector (c) evaluated zero times"
  local ++fails
}
if `nd' == 0 {
  di as error "NONVACUITY FAIL: detector (d) evaluated zero times"
  local ++fails
}
* erase everything this suite wrote (saveplots outputs)
local fl : dir "." files "refm_*"
foreach f of local fl {
  cap erase "`f'"
}
cap frame drop refm_fr
if `fails' == 0 di as result "REFMARGIN_ALL_PASS"
else di as error "REFMARGIN_FAILS=`fails'"
