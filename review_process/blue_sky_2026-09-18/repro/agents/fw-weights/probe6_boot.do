* probe6: the fw-identity under the RESAMPLING ENGINES.
*   B*: [aw=k] vs physical duplication under bootstrap(), cluster(i) (same 24 clusters in
*       both frames, so a common rseed() must draw the same multiplicities).  Point
*       estimates match (probe2), so every bootstrap family must match too.
*   L*: is the probe5 att_i leg gap (reid leg vs phys leg) systematic or Monte Carlo?
*       Repeated at three seeds and 150 reps.
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
    matlist `m1', format(%12.6g)
    dis "MAIN_B_`tag':"
    matlist `m2', format(%12.6g)
  }
end

*--- B: aw=k vs physical duplication under the bootstrap engines ---------------
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att cluster(i) bootstrap(standard) reps(100) rseed(555)) fb(p2_B.dta) cb(wooldid y i t ttre, att cluster(i) bootstrap(standard) reps(100) rseed(555)) tag(BCNT)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att_i cluster(i) bootstrap(standard) reps(60) rseed(555)) fb(p2_B.dta) cb(wooldid y i t ttre, att_i cluster(i) bootstrap(standard) reps(60) rseed(555)) tag(BREID)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att cluster(i) bootstrap(standard) reps(60) rseed(555) boottstat) fb(p2_B.dta) cb(wooldid y i t ttre, att cluster(i) bootstrap(standard) reps(60) rseed(555) boottstat) tag(BTST)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att cluster(i) bootstrap(bayes) reps(100) rseed(555)) fb(p2_B.dta) cb(wooldid y i t ttre, att cluster(i) bootstrap(bayes) reps(100) rseed(555)) tag(BBAYES)

*--- L: is the unclustered att_i leg gap systematic? ---------------------------
foreach s in 4242 777 31337 {
  pairrun , fa(p2_A.dta) ca(wooldid y i t ttre, att_i bootstrap(standard) reps(150) rseed(`s') boottstat) fb(p2_A.dta) cb(wooldid y i t ttre [aw=one], att_i bootstrap(standard) reps(150) rseed(`s') boottstat) tag(LEG`s')
}

dis "PROBE6_DONE"
