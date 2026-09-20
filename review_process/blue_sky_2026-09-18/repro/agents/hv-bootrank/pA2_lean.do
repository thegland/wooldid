*=======================================================================
* pA2_lean.do -- the MECHANISM half of the verification, on a small fast
* panel (40 counties in 20 states).  Answers:
*  (1) is rank(e(V)) a function of the REPLICATION COUNT (intrinsic to
*      any bootstrap) or of the `ereturn repost' bypass?
*  (2) would `ereturn post' have refused the matrix that repost accepts?
*  (3) on a TRUE NULL (this DGP has no pre-trends at all), what does a
*      hand `test' give at small reps?
* SELF-CONTAINED: needs only gen_panel.do next to it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

local pN 40
local pT 12
local pG 20
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank/gen_panel.do"

*---------------------------------------------------------------
* (1) rank vs reps.  K is fixed; only reps() moves.
*---------------------------------------------------------------
foreach r in 4 8 15 30 60 200 {
  dis "=== RS_`r'_BEGIN ==="
  qui cap wooldid y unit t ttre, att att_i att_it att_itime esfixedbaseperiod ///
      esprelength(3) espostlength(3) cluster(clu) bootstrap(standard) reps(`r') rseed(11)
  local rc = _rc
  dis "RS_`r'_RC=" `rc'
  if `rc' == 0 {
    dis "RS_`r'_COMPLETED=" e(boot_completed)
    dis "RS_`r'_NV=" e(boot_N_V)
    dis "RS_`r'_DFR=" e(df_r)
    mata: V = st_matrix("e(V)")
    mata: printf("RS_K=%f RANK=%f NSENT=%f\n", cols(V), rank(V), sum(diagonal(V) :< -1e99))
  }
  dis "=== RS_`r'_END ==="
}

*---------------------------------------------------------------
* (3) TRUE-NULL hand test at a range of reps.  The DGP has NO
*     pre-trends, so every ES_att_pre* is exactly 0 in truth.
*---------------------------------------------------------------
foreach r in 8 30 60 200 {
  dis "=== NULLTEST_`r'_BEGIN ==="
  qui cap wooldid y unit t ttre, att esfixedbaseperiod ///
      esprelength(5) espostlength(3) cluster(clu) bootstrap(standard) reps(`r') rseed(909)
  dis "NULLTEST_`r'_RC=" _rc
  if _rc == 0 {
    dis "NULLTEST_`r'_NV=" e(boot_N_V)
    dis "NULLTEST_`r'_DFR=" e(df_r)
    mata: V = st_matrix("e(V)")
    mata: printf("NULLTEST_K=%f RANK=%f\n", cols(V), rank(V))
    mata: st_matrixcolstripe("e(b)")
    cap noisily test ES_att_pre2 ES_att_pre3 ES_att_pre4 ES_att_pre5 ES_att_pre6
    dis "NULLTEST_`r'_TESTRC=" _rc
    dis "NULLTEST_`r'_F=" r(F)
    dis "NULLTEST_`r'_chi2=" r(chi2)
    dis "NULLTEST_`r'_df=" r(df)
    dis "NULLTEST_`r'_dfr=" r(df_r)
    dis "NULLTEST_`r'_p=" r(p)
  }
  dis "=== NULLTEST_`r'_END ==="
}

*---------------------------------------------------------------
* (2) Would ereturn post have refused what repost accepts?
*---------------------------------------------------------------
clear
set obs 300
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()
qui regress yy x1 x2
tempname bb vv
matrix `bb' = e(b)
matrix `vv' = e(V)
matrix `vv'[1,3] = 0
matrix `vv'[2,3] = 0
matrix `vv'[3,1] = 0
matrix `vv'[3,2] = 0
matrix `vv'[3,3] = 0
mata: printf("SCRATCH_RANK_A=%f of %f\n", rank(st_matrix("`vv'")), cols(st_matrix("`vv'")))
dis "=== POSTA_BEGIN ==="
cap noisily ereturn post `bb' `vv'
dis "POSTA_RC=" _rc
dis "=== POSTA_END ==="

* a genuinely rank-deficient PSD matrix built exactly the way a
* B-draw bootstrap builds one: K=5 columns from 3 draws -> rank 2
clear
set obs 3
mata: A = (1, 2, 0.5, -1, 3 \ 0.2, -1, 2, 0.7, 1 \ -1, 0.4, 1.5, 2, -2)
mata: d = A :- mean(A)
mata: VB = (d' * d) / 2
mata: st_matrix("VB", VB)
mata: printf("SCRATCH_RANK_B=%f of %f MINEIG=%g\n", rank(VB), cols(VB), min(symeigenvalues(makesymmetric(VB))))
mata: st_matrix("BB", J(1,5,0.3))
matrix colnames BB = c1 c2 c3 c4 c5
matrix colnames VB = c1 c2 c3 c4 c5
matrix rownames VB = c1 c2 c3 c4 c5
dis "=== POSTC_BEGIN ==="
cap noisily ereturn post BB VB
dis "POSTC_RC=" _rc
cap noisily test c1 c2
dis "POSTC_TEST_RC=" _rc
dis "POSTC_TEST_p=" r(p)
dis "POSTC_TEST_df=" r(df)
cap noisily test c1 c2 c3 c4 c5
dis "POSTC_TESTALL_RC=" _rc
dis "POSTC_TESTALL_p=" r(p)
dis "POSTC_TESTALL_df=" r(df)
dis "=== POSTC_END ==="

* the -1e100 sentinel on the diagonal: post vs repost
clear
set obs 300
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()
qui regress yy x1 x2
tempname b2 v2
matrix `b2' = e(b)
matrix `v2' = e(V)
matrix `v2'[1,1] = -1e100
matrix `v2'[1,2] = 0
matrix `v2'[2,1] = 0
matrix `v2'[1,3] = 0
matrix `v2'[3,1] = 0
dis "=== POSTB_BEGIN ==="
cap noisily ereturn post `b2' `v2'
dis "POSTB_RC=" _rc
dis "=== POSTB_END ==="

qui regress yy x1 x2
matrix `v2' = e(V)
matrix `v2'[1,1] = -1e100
dis "=== REPOSTB_BEGIN ==="
cap noisily ereturn repost V = `v2'
dis "REPOSTB_RC=" _rc
cap noisily test x1
dis "REPOSTB_TEST_RC=" _rc
dis "REPOSTB_TEST_F=" r(F)
dis "REPOSTB_TEST_p=" r(p)
cap noisily lincom x1
dis "REPOSTB_LINCOM_RC=" _rc
dis "=== REPOSTB_END ==="

dis "PROBE_A2_DONE"
