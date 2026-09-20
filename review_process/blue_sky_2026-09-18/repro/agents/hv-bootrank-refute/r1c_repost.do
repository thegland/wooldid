*=======================================================================
* r1c_repost.do -- the REPOST half, run where repost is legal (inside an
* eclass program).  My r1b showed `ereturn repost' returns rc 152 at the
* top level of a do-file, which is a harness fact, not a Stata fact.
* Question: does `ereturn repost' accept what `ereturn post' refuses,
* and does it zero a nonsymmetric matrix the way post does?
* SELF-CONTAINED, no wooldid, seconds.
*=======================================================================
clear all
set more off
set processors 2
set seed 20260917
set sortseed 20260917
dis "STATA_VERSION=" c(stata_version) " MP=" c(MP)

mata:
void shmat2(string scalar tag, string scalar nm)
{
  real matrix W
  W = st_matrix(nm)
  if (rows(W) == 0) { printf("%s_ABSENT=1\n", tag); return; }
  printf("%s_ALLZERO=%f  NMISS=%f  RANK=%f of %f  V11=%21.13g  V12=%21.13g\n",
         tag, all(W :== 0), sum(W :>= .), rank(editmissing(W,0)), rows(W), W[1,1], W[1,2])
}
end

qui set obs 300
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()

cap prog drop rprobe
prog def rprobe, eclass
  * ---- (R1) -1e100 sentinel through REPOST
  qui regress yy x1 x2
  tempname b1 v1
  matrix `b1' = e(b)
  matrix `v1' = e(V)
  ereturn post `b1' `v1'
  tempname vs
  matrix `vs' = e(V)
  matrix `vs'[1,1] = -1e100
  matrix `vs'[1,2] = 0
  matrix `vs'[2,1] = 0
  matrix `vs'[1,3] = 0
  matrix `vs'[3,1] = 0
  cap noisily ereturn repost V = `vs'
  dis "R1_REPOST_RC=" _rc
  mata: shmat2("R1_EV", "e(V)")
  cap noisily lincom x1
  dis "R1_LINCOM_RC=" _rc " se=" r(se)

  * ---- (R2) the SAME sentinel matrix through POST
  qui regress yy x1 x2
  tempname b2 v2
  matrix `b2' = e(b)
  matrix `v2' = e(V)
  matrix `v2'[1,1] = -1e100
  matrix `v2'[1,2] = 0
  matrix `v2'[2,1] = 0
  matrix `v2'[1,3] = 0
  matrix `v2'[3,1] = 0
  cap noisily ereturn post `b2' `v2'
  dis "R2_POST_RC=" _rc
  mata: shmat2("R2_EV", "e(V)")
  cap noisily lincom x1
  dis "R2_LINCOM_RC=" _rc " se=" r(se)

  * ---- (R3) NONSYMMETRIC through REPOST
  qui regress yy x1 x2
  tempname b3 v3
  matrix `b3' = e(b)
  matrix `v3' = e(V)
  ereturn post `b3' `v3'
  tempname vn
  matrix `vn' = e(V)
  matrix `vn'[1,2] = `vn'[1,2] + 0.01
  mata: printf("R3_INPUT_SYMDIFF=%g\n", max(abs(st_matrix("`vn'") - st_matrix("`vn'")')))
  cap noisily ereturn repost V = `vn'
  dis "R3_REPOST_RC=" _rc
  mata: shmat2("R3_EV", "e(V)")

  * ---- (R4) MISSING value through REPOST (post refuses this with r504)
  qui regress yy x1 x2
  tempname b4 v4
  matrix `b4' = e(b)
  matrix `v4' = e(V)
  ereturn post `b4' `v4'
  tempname vm
  matrix `vm' = e(V)
  matrix `vm'[1,1] = .
  cap noisily ereturn repost V = `vm'
  dis "R4_REPOST_RC=" _rc
  mata: shmat2("R4_EV", "e(V)")
  cap noisily lincom x1
  dis "R4_LINCOM_RC=" _rc " se=" r(se)

  * ---- (R5) MISSING value through POST, for the comparison
  qui regress yy x1 x2
  tempname b5 v5
  matrix `b5' = e(b)
  matrix `v5' = e(V)
  matrix `v5'[1,1] = .
  cap noisily ereturn post `b5' `v5'
  dis "R5_POST_RC=" _rc
  mata: shmat2("R5_EV", "e(V)")

  * ---- (R6) rank-deficient-but-symmetric (a a') through POST inside eclass
  qui regress yy x1 x2
  tempname b6
  matrix `b6' = e(b)
  mata: aa = (0.05 \ 0.03 \ 0.07); st_matrix("VS6", aa*aa')
  matrix colnames VS6 = x1 x2 _cons
  matrix rownames VS6 = x1 x2 _cons
  cap noisily ereturn post `b6' VS6
  dis "R6_POST_RC=" _rc
  mata: shmat2("R6_EV", "e(V)")
end
rprobe
dis "R1C_DONE"
