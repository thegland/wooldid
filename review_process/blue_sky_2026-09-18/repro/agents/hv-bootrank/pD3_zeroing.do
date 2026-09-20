*=======================================================================
* pD3_zeroing.do -- the OTHER half of the PS-62 premise.
* review_process/held_for_ruling_2026-09-09.md:99 asserts
*   "Stata's `ereturn post' screen does zero a nonsymmetric or highly
*    singular VCE wholesale"
* and only `repost' was ever measured (probe_ps62.do, 2026-09-10).
* Here BOTH halves are measured on `post': highly singular (symmetric)
* and NONSYMMETRIC.  SELF-CONTAINED, seconds.
*=======================================================================
clear all
set more off
set processors 2
set seed 20260916
clear
set obs 300
gen double yy = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()

cap prog drop zprobe
prog def zprobe, eclass
  * (A) HIGHLY SINGULAR but symmetric: V = a a' (rank 1 of 3)
  qui regress yy x1 x2
  tempname bb
  matrix `bb' = e(b)
  mata: a = (0.05 \ 0.03 \ 0.07)
  mata: st_matrix("VS", a * a')
  matrix colnames VS = x1 x2 _cons
  matrix rownames VS = x1 x2 _cons
  mata: printf("A_RANK=%f of %f\n", rank(st_matrix("VS")), cols(st_matrix("VS")))
  cap noisily ereturn post `bb' VS
  dis "A_POST_RC=" _rc
  if _rc == 0 {
    matrix g = e(V)
    dis "A_V11=" g[1,1] "  A_V12=" g[1,2] "  A_V33=" g[3,3]
    mata: printf("A_ALLZERO=%f\n", all(st_matrix("e(V)") :== 0))
  }

  * (B) NONSYMMETRIC V
  qui regress yy x1 x2
  matrix `bb' = e(b)
  tempname vn
  matrix `vn' = e(V)
  matrix `vn'[1,2] = `vn'[1,2] + 0.01
  mata: VN = st_matrix("`vn'")
  mata: printf("B_SYMDIFF=%g\n", max(abs(VN - VN')))
  cap noisily ereturn post `bb' `vn'
  dis "B_POST_RC=" _rc
  if _rc == 0 {
    matrix g2 = e(V)
    dis "B_V11=" g2[1,1] "  B_V12=" g2[1,2] "  B_V21=" g2[2,1]
    mata: printf("B_ALLZERO=%f\n", all(st_matrix("e(V)") :== 0))
  }

  * (C) NONSYMMETRIC via repost, for the comparison the record never made
  qui regress yy x1 x2
  tempname b3 v3
  matrix `b3' = e(b)
  matrix `v3' = e(V)
  ereturn post `b3' `v3'
  matrix `v3' = e(V)
  matrix `v3'[1,2] = `v3'[1,2] + 0.01
  cap noisily ereturn repost V = `v3'
  dis "C_REPOST_RC=" _rc
  if _rc == 0 {
    matrix g3 = e(V)
    dis "C_V12=" g3[1,2] "  C_V21=" g3[2,1]
    mata: printf("C_ALLZERO=%f\n", all(st_matrix("e(V)") :== 0))
  }
end
zprobe
dis "PROBE_D3_DONE"
