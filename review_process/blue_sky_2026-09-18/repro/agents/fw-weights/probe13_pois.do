* probe13: redo the poisson leg of probe7 on a fixture with genuine Poisson noise.
* probe7's ypos was a deterministic exp() of the model, so ppmlhdfe fit it exactly and every
* standard error was machine zero (se 1.2e-11, z 1.7e10); only e(b) was meaningful there.
* poisson accepts pweights and refuses aweights, so pw is the only weighted poisson path,
* and the metamorphic claim is that [pw=k] reproduces k physically duplicated rows for the
* POINT ESTIMATES (the pweight robust meat legitimately differs).
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
    printf("REL_%s=SHAPEMISMATCH\n", tag)
  }
  else {
    D = abs(A:-B) :/ (abs(A):*(abs(A):>1e-12) :+ (abs(A):<=1e-12))
    printf("REL_%s=%10.3e  maxabs=%10.3e\n", tag, colmax(rowmax(D)), colmax(rowmax(abs(A:-B))))
  }
}
end

use "p2_A.dta", clear
set seed 31337
gen double mu = exp(0.3 + 0.2*post + 0.08*t + 0.3*ufe)
gen double ycount = rpoisson(mu)
save "p13_A.dta", replace
expand k
drop k
save "p13_B.dta", replace

program define pairrun
  syntax , fa(string) ca(string asis) fb(string) cb(string asis) tag(string)
  tempname b1 b2 m1 m2
  qui use "`fa'", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 {
    matrix `b1' = e(b)
    cap matrix `m1' = e(wdidmainresults)
  }
  qui use "`fb'", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 {
    matrix `b2' = e(b)
    cap matrix `m2' = e(wdidmainresults)
  }
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    mata: wd_rel("`b1'","`b2'","`tag'_b")
    mata: wd_rel("`m1'","`m2'","`tag'_main")
    dis "MAIN_A_`tag':"
    matlist `m1', format(%14.8g)
    dis "MAIN_B_`tag':"
    matlist `m2', format(%14.8g)
  }
end

pairrun , fa(p13_A.dta) ca(wooldid ycount i t ttre [pw=k], att poisson) fb(p13_B.dta) cb(wooldid ycount i t ttre, att poisson) tag(POIS)
pairrun , fa(p13_A.dta) ca(wooldid ycount i t ttre [pw=k], att att_i att_it att_itime poisson cluster(i)) fb(p13_B.dta) cb(wooldid ycount i t ttre, att att_i att_it att_itime poisson cluster(i)) tag(POISALL)
pairrun , fa(p13_A.dta) ca(wooldid ycount i t ttre [pw=k], att att_it poisson poisexpresults esprelength(2) espostlength(2)) fb(p13_B.dta) cb(wooldid ycount i t ttre, att att_it poisson poisexpresults esprelength(2) espostlength(2)) tag(POISES)
dis "PROBE13_DONE"
