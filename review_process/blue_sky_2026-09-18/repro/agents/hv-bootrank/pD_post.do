*=======================================================================
* pD_post.do -- MECHANISM test, no wooldid needed, runs in seconds.
* The repost-site comment (ado:10044-10053) says `ereturn repost V ='
* is used ON PURPOSE because it "does NOT re-run `ereturn post''s
* positive-definiteness screen".  Two questions:
*   (a) does `ereturn post' actually refuse a rank-deficient PSD
*       matrix?  (if not, the repost cannot be the cause of the
*       rank-7-over-20 e(V) the seeded claim reports)
*   (b) does `ereturn post' actually refuse the -1e100 sentinel?
*       (if not, the comment's stated premise is wrong)
* SELF-CONTAINED.
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

* ---- (a1) rank-deficient PSD by zeroing a row/col -------------------
qui regress yy x1 x2
tempname bb vv
matrix `bb' = e(b)
matrix `vv' = e(V)
matrix `vv'[1,3] = 0
matrix `vv'[2,3] = 0
matrix `vv'[3,1] = 0
matrix `vv'[3,2] = 0
matrix `vv'[3,3] = 0
mata: printf("RANK_A=%f of %f\n", rank(st_matrix("`vv'")), cols(st_matrix("`vv'")))
dis "=== POSTA_BEGIN ==="
cap noisily ereturn post `bb' `vv'
dis "POSTA_RC=" _rc
dis "=== POSTA_END ==="

* ---- (a2) a covariance built EXACTLY as a B<K bootstrap builds one --
* K=5 coefficients from 3 draws -> rank 2, PSD, no zero rows
clear
set obs 3
mata: A = (1, 2, 0.5, -1, 3 \ 0.2, -1, 2, 0.7, 1 \ -1, 0.4, 1.5, 2, -2)
mata: d = A :- mean(A)
mata: VB = (d' * d) / 2
mata: st_matrix("VB", VB)
mata: printf("RANK_B=%f of %f MINEIG=%g\n", rank(VB), cols(VB), min(symeigenvalues(makesymmetric(VB))))
mata: st_matrix("BB", J(1,5,0.3))
matrix colnames BB = c1 c2 c3 c4 c5
matrix colnames VB = c1 c2 c3 c4 c5
matrix rownames VB = c1 c2 c3 c4 c5
dis "=== POSTC_BEGIN ==="
cap noisily ereturn post BB VB
dis "POSTC_RC=" _rc
cap noisily test c1 c2
dis "POSTC_TEST2_RC=" _rc
dis "POSTC_TEST2_df=" r(df)
dis "POSTC_TEST2_p=" r(p)
cap noisily test c1 c2 c3 c4 c5
dis "POSTC_TEST5_RC=" _rc
dis "POSTC_TEST5_df=" r(df)
dis "POSTC_TEST5_p=" r(p)
dis "=== POSTC_END ==="

* ---- (b) the -1e100 sentinel: post vs repost ------------------------
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

dis "PROBE_D_DONE"
