*=======================================================================
* pD2_sentinel.do -- does `ereturn post' REALLY refuse / alter the
* -1e100 sentinel, as the repost-site comment (ado:10044-10053) and
* wooldid_bootfin's comment (ado:14554-14560) both assert?
* Both are run inside a proper e-class program so `ereturn repost' is
* legal, which it is not at the top level of a do-file.
* SELF-CONTAINED, runs in seconds.
*=======================================================================
clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916

clear
set obs 300
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()

cap prog drop sentprobe
prog def sentprobe, eclass
  qui regress yy x1 x2
  tempname bb vv
  matrix `bb' = e(b)
  matrix `vv' = e(V)
  matrix `vv'[1,1] = -1e100
  matrix `vv'[1,2] = 0
  matrix `vv'[2,1] = 0
  matrix `vv'[1,3] = 0
  matrix `vv'[3,1] = 0
  dis "PRE_V11=" `vv'[1,1]

  * --- (1) ereturn post with the sentinel ---
  cap noisily ereturn post `bb' `vv'
  dis "POST_SENT_RC=" _rc
  if _rc == 0 {
    matrix vgot = e(V)
    dis "POST_SENT_V11=" vgot[1,1]
    dis "POST_SENT_V22=" vgot[2,2]
    cap noisily test x1
    dis "POST_SENT_TEST_RC=" _rc
    dis "POST_SENT_TEST_p=" r(p)
    cap noisily lincom x1
    dis "POST_SENT_LINCOM_RC=" _rc
  }

  * --- (2) ereturn repost with the sentinel, from a posted state ---
  qui regress yy x1 x2
  tempname b3 v3
  matrix `b3' = e(b)
  matrix `v3' = e(V)
  ereturn post `b3' `v3'
  matrix `v3' = e(V)
  matrix `v3'[1,1] = -1e100
  cap noisily ereturn repost V = `v3'
  dis "REPOST_SENT_RC=" _rc
  if _rc == 0 {
    matrix vgot2 = e(V)
    dis "REPOST_SENT_V11=" vgot2[1,1]
    cap noisily test x1
    dis "REPOST_SENT_TEST_RC=" _rc
    dis "REPOST_SENT_TEST_p=" r(p)
    dis "REPOST_SENT_TEST_F=" r(F)
    dis "REPOST_SENT_TEST_chi2=" r(chi2)
  }

  * --- (3) ereturn post with a rank-deficient bootstrap-style V ---
  qui regress yy x1 x2
  tempname b4 v4
  matrix `b4' = e(b)
  matrix `v4' = e(V)
  * force exact rank 1: V = a a'
  mata: a = (0.05 \ 0.03 \ 0.07)
  mata: st_matrix("V4", a * a')
  matrix colnames V4 = x1 x2 _cons
  matrix rownames V4 = x1 x2 _cons
  mata: printf("RANK1_CHECK=%f of %f\n", rank(st_matrix("V4")), cols(st_matrix("V4")))
  cap noisily ereturn post `b4' V4
  dis "POST_RANK1_RC=" _rc
  if _rc == 0 {
    matrix vg3 = e(V)
    dis "POST_RANK1_V11=" vg3[1,1]
    cap noisily test x1 x2
    dis "POST_RANK1_TEST_RC=" _rc
    dis "POST_RANK1_TEST_df=" r(df)
    dis "POST_RANK1_TEST_p=" r(p)
  }
end

sentprobe
dis "PROBE_D2_DONE"
