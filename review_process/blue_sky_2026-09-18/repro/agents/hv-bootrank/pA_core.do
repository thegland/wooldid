*=======================================================================
* pA_core.do -- VERIFY the seeded claim: reps(8) posts a rank-deficient
* e(V) over 20 coefficients and a hand test() returns a tiny p with no
* warning.  ALSO measure the mechanism (rank vs reps; ereturn post's
* screen) so the cause can be stated correctly.
* SELF-CONTAINED: run it from this directory.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

dis "WHICH_START"
which wooldid
dis "WHICH_END"

local pN 200
local pT 12
local pG 40
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank/gen_panel.do"

qui count
dis "NOBS=" r(N)
qui levelsof clu, local(cl)
dis "NCLU=" `: word count `cl''

*---------------------------------------------------------------
* 1. THE SEEDED COMMAND, verbatim shape
*---------------------------------------------------------------
dis "=== RUN8_BEGIN ==="
cap noisily wooldid y unit t ttre, att att_i att_it att_itime esfixedbaseperiod ///
    esprelength(3) espostlength(3) cluster(clu) bootstrap(standard) reps(8) rseed(11)
dis "RUN8_RC=" _rc
dis "=== RUN8_END ==="

if _rc == 0 {
  dis "K=" colsof(e(b))
  dis "PROPS=" e(properties)
  dis "VCVSTATE=" e(wdid_vcvstate)
  dis "BOOT_N_V=" e(boot_N_V)
  dis "BOOT_FAILED=" e(boot_failed)
  dis "BOOT_VPOSTED=" e(boot_v_posted)
  mat list e(b), format(%12.6f)
  mata: V = st_matrix("e(V)")
  mata: printf("RANKV=%f\n", rank(V))
  mata: printf("KV=%f\n", cols(V))
  mata: printf("MINDIAG=%21x\n", min(diagonal(V)))
  mata: printf("MAXDIAG=%21x\n", max(diagonal(V)))
  mata: printf("NSENT=%f\n", sum(diagonal(V) :< -1e99))
  mata: ev = symeigenvalues(makesymmetric(V))
  mata: ev
  mata: printf("SUMNEG=%f\n", sum(ev :< 0))
  * list the coefficient names so a skeptic can re-target the test
  mata: st_matrixcolstripe("e(V)")
}

*---------------------------------------------------------------
* 2. The hand test the claim names (3 df, POST periods)
*---------------------------------------------------------------
dis "=== TEST3_BEGIN ==="
cap noisily test ES_att_post0 ES_att_post1 ES_att_post2
dis "TEST3_RC=" _rc
dis "TEST3_chi2=" r(chi2)
dis "TEST3_F=" r(F)
dis "TEST3_df=" r(df)
dis "TEST3_dfr=" r(df_r)
dis "TEST3_p=" r(p)
dis "=== TEST3_END ==="

* the PRE periods: a TRUE null under this DGP (no pre-trends at all),
* so any rejection here is pure size distortion
dis "=== TESTPRE_BEGIN ==="
cap noisily test ES_att_pre1 ES_att_pre2 ES_att_pre3
dis "TESTPRE_RC=" _rc
dis "TESTPRE_chi2=" r(chi2)
dis "TESTPRE_df=" r(df)
dis "TESTPRE_p=" r(p)
dis "=== TESTPRE_END ==="

* a test with MORE constraints than the rank of V: does Stata drop them?
dis "=== TESTBIG_BEGIN ==="
cap noisily test ES_att_pre1 ES_att_pre2 ES_att_pre3 ES_att_post0 ES_att_post1 ///
    ES_att_post2 ES_att_post3 (att) (att_i) (att_it)
dis "TESTBIG_RC=" _rc
dis "TESTBIG_chi2=" r(chi2)
dis "TESTBIG_df=" r(df)
dis "TESTBIG_p=" r(p)
dis "=== TESTBIG_END ==="

*---------------------------------------------------------------
* 3. MECHANISM A: rank as a function of reps
*---------------------------------------------------------------
foreach r in 4 8 15 30 60 {
  dis "=== REPSCAN_`r'_BEGIN ==="
  qui cap noisily wooldid y unit t ttre, att att_i att_it att_itime esfixedbaseperiod ///
      esprelength(3) espostlength(3) cluster(clu) bootstrap(standard) reps(`r') rseed(11)
  local rc = _rc
  dis "REPSCAN_`r'_RC=" `rc'
  if `rc' == 0 {
    dis "REPSCAN_`r'_NV=" e(boot_N_V)
    dis "REPSCAN_`r'_FAILED=" e(boot_failed)
    mata: V = st_matrix("e(V)")
    mata: printf("REPSCAN_K=%f RANK=%f\n", cols(V), rank(V))
  }
  dis "=== REPSCAN_`r'_END ==="
}

*---------------------------------------------------------------
* 4. MECHANISM B: would ereturn post have caught it?
*    (a) rank-deficient but PSD V   (b) V with a -1e100 diagonal
*---------------------------------------------------------------
preserve
clear
set obs 200
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()
qui regress yy x1 x2
tempname bb vv
matrix `bb' = e(b)
matrix `vv' = e(V)

* (a) rank-deficient PSD: zero out the last row/col entirely
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

qui regress yy x1 x2
matrix `bb' = e(b)
matrix `vv' = e(V)
* (b) the -1e100 sentinel on one diagonal
matrix `vv'[1,1] = -1e100
matrix `vv'[1,2] = 0
matrix `vv'[2,1] = 0
matrix `vv'[1,3] = 0
matrix `vv'[3,1] = 0
dis "=== POSTB_BEGIN ==="
cap noisily ereturn post `bb' `vv'
dis "POSTB_RC=" _rc
dis "=== POSTB_END ==="

* and the repost route with the same sentinel
qui regress yy x1 x2
matrix `vv' = e(V)
matrix `vv'[1,1] = -1e100
dis "=== REPOSTB_BEGIN ==="
cap noisily ereturn repost V = `vv'
dis "REPOSTB_RC=" _rc
cap noisily test x1
dis "REPOSTB_TEST_RC=" _rc
dis "REPOSTB_TEST_chi2=" r(chi2)
dis "REPOSTB_TEST_p=" r(p)
dis "=== REPOSTB_END ==="
restore

dis "PROBE_A_DONE"
