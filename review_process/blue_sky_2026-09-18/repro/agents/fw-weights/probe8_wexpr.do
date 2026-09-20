* probe8: THE WEIGHT EXPRESSION HAZARD.
* Stata lets any command take [aw = <expression>]; every official command materialises
* the expression before using it.  wooldid instead keeps the raw expression text in the
* macro `w' and PASTES it into arithmetic:
*     ado:1619  gegen `w_per_it' = total(`esamp' * `w') , by(`i' `t')
*     ado:1630  gen double `w_att_i'     = `w' `customweighttext'/`w_per_i'
*     ado:1637  gen double `w_att_itime' = `w' `customweighttext' /(`w_per_it'*`t_per_i')
*     ado:1644  gen double `w_att_it'    = `w' `customweighttext'/`w_per_it'
*     ado:1650  gen double `w_att'       = `w' `customweighttext'
* so an expression whose top-level operator is + or - binds WRONG:
*     [aw = ka+kb]  ->  w_att_it = ka + kb/W_it   instead of   (ka+kb)/W_it.
* METAMORPHIC TEST: ka+kb is constructed to equal k EXACTLY, so [aw=ka+kb] must be
* indistinguishable from [aw=k].
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

mata:
void wd_rel(string scalar a, string scalar b, string scalar tag) {
  real matrix A, B, D
  A = st_matrix(a); B = st_matrix(b)
  if (rows(A)==0 | rows(B)==0) {
    printf("REL_%s=EMPTY\n", tag)
  }
  else if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    printf("REL_%s=SHAPEMISMATCH %g x %g vs %g x %g\n", tag, rows(A),cols(A),rows(B),cols(B))
  }
  else {
    D = abs(A:-B) :/ (abs(A):*(abs(A):>1e-12) :+ (abs(A):<=1e-12))
    printf("REL_%s=%10.3e  maxabs=%10.3e\n", tag, colmax(rowmax(D)), colmax(rowmax(abs(A:-B))))
  }
}
end

program define pairrun
  syntax , fa(string) ca(string asis) fb(string) cb(string asis) tag(string)
  tempname b1 b2
  qui use "`fa'", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 matrix `b1' = e(b)
  qui use "`fb'", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 matrix `b2' = e(b)
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    mata: wd_rel("`b1'","`b2'","`tag'")
    dis "BVEC_A_`tag':"
    matlist `b1', format(%16.10g)
    dis "BVEC_B_`tag':"
    matlist `b2', format(%16.10g)
  }
end

use "p2_A.dta", clear
gen double ka = floor(k/2)
gen double kb = k - ka
assert ka + kb == k
gen double two = 2
gen double kx2 = 2*k
save "p8.dta", replace

* control: a product expression is safe (\* binds tighter than /)
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=kx2], att att_i att_it att_itime) fb(p8.dta) cb(wooldid y i t ttre [aw=2*k], att att_i att_it att_itime) tag(CTRL_PROD)

* the hazard: a SUM expression
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=k], att att_i att_it att_itime) fb(p8.dta) cb(wooldid y i t ttre [aw=ka+kb], att att_i att_it att_itime) tag(SUMEXPR)
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=k], att) fb(p8.dta) cb(wooldid y i t ttre [aw=ka+kb], att) tag(SUMEXPR_ATTONLY)
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=k], att customweightmultiplier(two)) fb(p8.dta) cb(wooldid y i t ttre [aw=ka+kb], att customweightmultiplier(two)) tag(SUMEXPR_CWM)
pairrun , fa(p8.dta) ca(wooldid y i t ttre [pw=k], att att_i att_it att_itime) fb(p8.dta) cb(wooldid y i t ttre [pw=ka+kb], att att_i att_it att_itime) tag(SUMEXPR_PW)

* with an event study (also exercises ado:1792 `total `w'' and ado:3833)
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=k], att att_it esprelength(2) espostlength(3) summarystats) fb(p8.dta) cb(wooldid y i t ttre [aw=ka+kb], att att_it esprelength(2) espostlength(3) summarystats) tag(SUMEXPR_ES)

* a minus expression
gen double kbig2 = k + 5
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=k], att att_i att_it att_itime) fb(p8.dta) cb(wooldid y i t ttre [aw=kbig2-5], att att_i att_it att_itime) tag(MINUSEXPR)

* spaces around the operator
pairrun , fa(p8.dta) ca(wooldid y i t ttre [aw=k], att att_i att_it att_itime) fb(p8.dta) cb(wooldid y i t ttre [aw = ka + kb], att att_i att_it att_itime) tag(SUMEXPR_SPACED)

* what an official Stata command does with the same expression (sanity: reghdfe agrees)
use "p8.dta", clear
qui reghdfe y post [aw=k], absorb(i t)
scalar bk = _b[post]
qui reghdfe y post [aw=ka+kb], absorb(i t)
scalar bs = _b[post]
dis "REGHDFE_k=" %21x bk
dis "REGHDFE_sum=" %21x bs
dis "REGHDFE_SAME=" (bk == bs)

dis "PROBE8_DONE"
