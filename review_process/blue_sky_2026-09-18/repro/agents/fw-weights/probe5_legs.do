* probe5: (A) [aw=one] must equal NO WEIGHT at all  (the iw-vs-aw margins plumbing);
*         (G) BOOTSTRAP LEG EQUIVALENCE: unclustered boottstat with no weight takes the
*             internal FWEIGHT count leg (ado:9048); with a user weight identically 1 it
*             takes the PHYSICAL-EXPANSION leg (ado:379/9236).  The data is the same, so
*             every studentized quantity must agree.  ado:9040 asserts the two are
*             "arithmetically IDENTICAL ... on all 400 pairs".
*         (F) permute() at a common rseed: [aw=k] vs physical duplication.
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
  tempname b1 b2 m1 m2 v1 v2
  qui use "`fa'", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 {
    matrix `b1' = e(b)
    cap matrix `m1' = e(wdidmainresults)
    cap matrix `v1' = e(V)
    local n1 = e(N)
  }
  qui use "`fb'", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 {
    matrix `b2' = e(b)
    cap matrix `m2' = e(wdidmainresults)
    cap matrix `v2' = e(V)
    local n2 = e(N)
  }
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    dis "PR_`tag'_N=" `n1' "," `n2'
    mata: wd_rel("`b1'","`b2'","`tag'_b")
    mata: wd_rel("`m1'","`m2'","`tag'_main")
    mata: wd_rel("`v1'","`v2'","`tag'_V")
  }
end

*--- (A) aw=one vs no weight ---------------------------------------------------
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=one], att att_i att_it att_itime) fb(p2_A.dta) cb(wooldid y i t ttre, att att_i att_it att_itime) tag(AWONE)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=one], att att_it esprelength(2) espostlength(3) jointtests) fb(p2_A.dta) cb(wooldid y i t ttre, att att_it esprelength(2) espostlength(3) jointtests) tag(AWONEES)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=one], att cluster(i)) fb(p2_A.dta) cb(wooldid y i t ttre, att cluster(i)) tag(AWONECL)

*--- (G) leg equivalence: unclustered boottstat --------------------------------
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre, att bootstrap(standard) reps(50) rseed(4242) boottstat) fb(p2_A.dta) cb(wooldid y i t ttre [aw=one], att bootstrap(standard) reps(50) rseed(4242) boottstat) tag(LEGATT)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre, att_i bootstrap(standard) reps(50) rseed(4242) boottstat) fb(p2_A.dta) cb(wooldid y i t ttre [aw=one], att_i bootstrap(standard) reps(50) rseed(4242) boottstat) tag(LEGATTI)

* same comparison WITHOUT boottstat (both legs are then the count leg) - control
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre, att bootstrap(standard) reps(50) rseed(4242)) fb(p2_A.dta) cb(wooldid y i t ttre [aw=one], att bootstrap(standard) reps(50) rseed(4242)) tag(NOBTSATT)

*--- (F) permute at a common seed: aw=k vs physical duplication ----------------
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att cluster(i) permute(shuffle) permreps(60) rseed(909)) fb(p2_B.dta) cb(wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(60) rseed(909)) tag(PERMATT)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att att_i cluster(i) permute(shuffle) permreps(60) rseed(909) permtstat) fb(p2_B.dta) cb(wooldid y i t ttre, att att_i cluster(i) permute(shuffle) permreps(60) rseed(909) permtstat) tag(PERMTS)

dis "PROBE5_DONE"
