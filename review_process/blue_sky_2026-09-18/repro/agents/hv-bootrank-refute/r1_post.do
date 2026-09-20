*=======================================================================
* r1_post.do -- SKEPTIC's own test of what `ereturn post' and
* `ereturn repost' actually accept.  No wooldid, no data, runs in seconds.
* Tests the premise of hv-bootrank-01 and -02 from scratch.
*=======================================================================
clear all
set more off
set processors 2
set seed 20260917
set sortseed 20260917
dis "STATA_VERSION=" c(stata_version) " FLAVOR=" c(flavor) " BORN=" c(born_date)

mata:
void shmat(string scalar tag, string scalar nm)
{
  real matrix W
  W = st_matrix(nm)
  if (rows(W) == 0) { printf("%s_ABSENT=1\n", tag); return; }
  printf("%s_ALLZERO=%f\n", tag, all(W :== 0))
  printf("%s_NMISS=%f\n", tag, sum(W :>= .))
  printf("%s_RANK=%f of %f\n", tag, rank(editmissing(W,0)), rows(W))
  printf("%s_V11=%21.13g\n", tag, W[1,1])
  if (cols(W) > 1) printf("%s_V12=%21.13g\n", tag, W[1,2])
  if (cols(W) > 2) printf("%s_V33=%21.13g\n", tag, W[3,3])
  printf("%s_SYM=%f\n", tag, all(W :== W'))
}
end

program define shwv
  args tag
  cap confirm matrix e(V)
  if _rc {
     dis "`tag'_EV_ABSENT=1"
     exit
  }
  tempname W
  mat `W' = e(V)
  mata: shmat("`tag'", "`W'")
end

qui set obs 200
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double x3 = rnormal()
gen double y  = x1 + rnormal()

* ==================== A. rank-deficient PSD, symmetric ====================
mat b = (1, 2, 3)
mat colnames b = x1 x2 x3
mat A = (0.05 \ 0.03 \ 0.07)
mat V = A * A'
mat colnames V = x1 x2 x3
mat rownames V = x1 x2 x3
mata: shmat("A_INPUT", "V")
cap noisily ereturn post b V, depname(y) obs(200)
dis "A_POST_RC=" _rc
shwv A
cap noisily test x1 x2 x3
dis "A_TEST_RC=" _rc
dis "A_TEST_df=" r(df) " chi2=" r(chi2) " p=" r(p)

* ==================== B. -1e100 sentinel diagonal ====================
mat V2 = J(3,3,0)
mat V2[1,1] = -1e100
mat V2[2,2] = 0.0025
mat V2[3,3] = 0.0049
mat colnames V2 = x1 x2 x3
mat rownames V2 = x1 x2 x3
cap noisily ereturn post b V2, depname(y) obs(200)
dis "B_POST_RC=" _rc
shwv B
cap noisily lincom x1
dis "B_LINCOM_RC=" _rc " se=" r(se)
cap noisily test x1
dis "B_TEST_RC=" _rc " df=" r(df) " p=" r(p)
cap noisily ereturn post b V, depname(y) obs(200)
cap noisily ereturn repost V = V2
dis "B_REPOST_RC=" _rc
shwv BR

* ==================== C. NONSYMMETRIC ====================
mat V3 = V
mat V3[1,2] = V3[1,2] + 0.01
mata: shmat("C_INPUT", "V3")
cap noisily ereturn post b V3, depname(y) obs(200)
dis "C_POST_RC=" _rc
shwv C
cap noisily ereturn post b V, depname(y) obs(200)
cap noisily ereturn repost V = V3
dis "C_REPOST_RC=" _rc
shwv CR

* ==================== D. MISSING VALUE in V ====================
mat V4 = V
mat V4[1,1] = .
cap noisily ereturn post b V4, depname(y) obs(200)
dis "D_POST_RC=" _rc
shwv D
cap noisily ereturn post b V, depname(y) obs(200)
cap noisily ereturn repost V = V4
dis "D_REPOST_RC=" _rc
shwv DR

* ==================== E. finite NEGATIVE diagonal ====================
mat V5 = V
mat V5[1,1] = -0.004
cap noisily ereturn post b V5, depname(y) obs(200)
dis "E_POST_RC=" _rc
shwv E
cap noisily lincom x1
dis "E_LINCOM_RC=" _rc " se=" r(se)

* ==================== F. bootstrap-shaped rank-2-of-5 ====================
mata: D = rnormal(3,5,0,1); D = D :- mean(D); VV = (D'*D)/2; st_matrix("V6", VV)
mat colnames V6 = c1 c2 c3 c4 c5
mat rownames V6 = c1 c2 c3 c4 c5
mat b6 = (1,1,1,1,1)
mat colnames b6 = c1 c2 c3 c4 c5
mata: shmat("F_INPUT", "V6")
cap noisily ereturn post b6 V6, depname(y) obs(200)
dis "F_POST_RC=" _rc
shwv F
cap noisily test c1 c2 c3 c4 c5
dis "F_TEST_RC=" _rc " df=" r(df) " chi2=" r(chi2) " p=" r(p)
cap noisily test c1 c2
dis "F_TEST2_RC=" _rc " df=" r(df) " chi2=" r(chi2) " p=" r(p)

dis "R1_DONE"
