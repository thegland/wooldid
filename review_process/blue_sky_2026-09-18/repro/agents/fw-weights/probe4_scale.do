* probe4: (a) per-COLUMN attribution of the escomp/gteffects divergences probe2 found;
*         (b) WEIGHT-SCALE INVARIANCE.  ado:1654 asserts "weight scale here actually
*             doesn't matter - no need to multiply weights by scalar adjustment factors".
*             A user whose survey weights are scaled to population totals instead of to
*             mean one must get identical output.  Everything must match: b, se, N, df,
*             escomp, gteffects.
*         (c) continuous treatment + lattice under aw=k vs physical duplication
*             (weighted quantiles are a classic aw-vs-fw divergence site).
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

mata:
void wd_relcol(string scalar a, string scalar b, string scalar tag) {
  real matrix A, B
  real scalar j, m
  string rowvector nm
  A = st_matrix(a); B = st_matrix(b)
  if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    printf("COL_%s=SHAPEMISMATCH %g x %g vs %g x %g\n", tag, rows(A),cols(A),rows(B),cols(B))
  }
  else {
    nm = st_matrixcolstripe(a)[.,2]'
    for (j=1; j<=cols(A); j++) {
      m = colmax(abs(A[.,j]:-B[.,j]))
      if (m > 1e-9) printf("COL_%s DIFF col=%s maxabs=%10.3e\n", tag, nm[j], m)
    }
    printf("COL_%s=SCANNED %g cols\n", tag, cols(A))
  }
}
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

*--- (a) per-column attribution ------------------------------------------------
use "p2_A.dta", clear
qui wooldid y i t ttre [aw=k], att att_it esprelength(2) espostlength(3) gteffects
matrix EA = e(escomp)
matrix GA = e(gteffects)
use "p2_B.dta", clear
qui wooldid y i t ttre , att att_it esprelength(2) espostlength(3) gteffects
matrix EB = e(escomp)
matrix GB = e(gteffects)
mata: wd_relcol("EA","EB","ESCOMP")
mata: wd_relcol("GA","GB","GTEFF")

*--- (b) scale invariance ------------------------------------------------------
program define scalecmp
  syntax , opts(string) tag(string) w1(string) w2(string)
  tempname b1 b2 v1 v2
  qui use "p2_A.dta", clear
  cap gen double kbig  = 1000*k
  cap gen double khalf = k/7
  cap noi wooldid y i t ttre `w1' , `opts'
  local r1 = _rc
  if `r1' == 0 {
    matrix `b1' = e(b)
    matrix `v1' = e(V)
    local n1 = e(N)
    local d1 = e(df_r)
    tempname M1
    cap matrix `M1' = e(wdidmainresults)
  }
  qui use "p2_A.dta", clear
  cap gen double kbig  = 1000*k
  cap gen double khalf = k/7
  cap noi wooldid y i t ttre `w2' , `opts'
  local r2 = _rc
  if `r2' == 0 {
    matrix `b2' = e(b)
    matrix `v2' = e(V)
    local n2 = e(N)
    local d2 = e(df_r)
    tempname M2
    cap matrix `M2' = e(wdidmainresults)
  }
  dis "SCALE_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    dis "SCALE_`tag'_N=" `n1' "," `n2' "  DF=" `d1' "," `d2'
    mata: wd_rel("`b1'","`b2'","`tag'_b")
    mata: wd_rel("`v1'","`v2'","`tag'_V")
    mata: wd_rel("`M1'","`M2'","`tag'_main")
  }
end

scalecmp , opts(att att_i att_it att_itime) tag(AW_X1000) w1([aw=k]) w2([aw=kbig])
scalecmp , opts(att att_i att_it att_itime) tag(AW_DIV7)  w1([aw=k]) w2([aw=khalf])
scalecmp , opts(att att_i att_it att_itime) tag(PW_X1000) w1([pw=k]) w2([pw=kbig])
scalecmp , opts(att att_it esprelength(2) espostlength(3) jointtests) tag(AWES_X1000) w1([aw=k]) w2([aw=kbig])
scalecmp , opts(att cluster(i)) tag(AWCL_X1000) w1([aw=k]) w2([aw=kbig])
scalecmp , opts(att) tag(AWONE_VS_NONE) w1([aw=one]) w2()
scalecmp , opts(att att_i att_it att_itime esprelength(2) espostlength(3)) tag(AWONE_VS_NONE_ALL) w1([aw=one]) w2()

*--- (c) continuous treatment + lattice, aw=k vs duplication --------------------
foreach o in "contreat(dose) contreatcontrols(basic)" "contreat(dose) contreatcontrols(basic) lattice(dose 4)" "contreat(dose) contreatcontrols(basic) lattice(quantile 4)" {
  tempname cA cB
  qui use "p2_A.dta", clear
  cap noi wooldid y i t ttre [aw=k], att `o'
  local ra = _rc
  if `ra' == 0 matrix `cA' = e(b)
  qui use "p2_B.dta", clear
  cap noi wooldid y i t ttre , att `o'
  local rb = _rc
  if `rb' == 0 matrix `cB' = e(b)
  dis "CONT_RC=" `ra' "," `rb' "  opts=`o'"
  if `ra' == 0 & `rb' == 0 {
    mata: wd_rel("`cA'","`cB'","CONT")
  }
}

dis "PROBE4_DONE"
