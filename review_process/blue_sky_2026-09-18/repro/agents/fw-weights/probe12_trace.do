* probe12: locate the weight-EXPRESSION failure, and check whether safety(off) turns it
* from an error into a silent wrong number.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

use "p8.dta", clear
* baseline
qui wooldid y i t ttre [aw=k], att att_i att_it att_itime
matrix BREF = e(b)
dis "BASELINE_OK"

* does safety(off) change anything?
cap noi wooldid y i t ttre [aw=ka+kb], att att_i att_it att_itime safety(off)
dis "SAFOFF_SUM_RC=" _rc
if _rc == 0 {
  matrix BSO = e(b)
  matlist BREF, format(%18.12g)
  matlist BSO, format(%18.12g)
}
cap noi wooldid y i t ttre [aw=2*k], att safety(off)
dis "SAFOFF_PROD_RC=" _rc

* other realistic expression shapes
foreach e in "k/2" "k*1" "(ka+kb)" "k+0" "int(k)" "cond(k>2,k,1)" "ka" {
  cap noi wooldid y i t ttre [aw=`e'], att
  dis "EXPR_`e'_RC=" _rc
}

* parenthesised sum: does it get further, and is it RIGHT?
cap noi wooldid y i t ttre [aw=(ka+kb)], att att_i att_it att_itime
dis "PARENS_RC=" _rc
if _rc == 0 {
  matrix BPAR = e(b)
  mata:
  A = st_matrix("BREF"); B = st_matrix("BPAR")
  nm = st_matrixcolstripe("BREF")[.,2]'
  for (j=1; j<=cols(A); j++) {
    printf("PARENS %-12s  aw=k %18.12g   aw=(ka+kb) %18.12g   reldiff %10.3e\n",
       nm[j], A[1,j], B[1,j], abs(A[1,j]-B[1,j])/abs(A[1,j]))
  }
  end
}

* where does it die?  trace the failing call at depth 1
set trace on
set tracedepth 1
cap noi wooldid y i t ttre [aw=ka+kb], att
set trace off
dis "TRACE_RC=" _rc

dis "PROBE12_DONE"
