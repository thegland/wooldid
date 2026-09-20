* probe17: the expression guard that stops [aw=ka+kb] is REGHDFE's, not wooldid's
*   (probe16 trace: reghdfe.Estimate runs `unab exp : ka+kb, min(1) max(1)' -> "+ invalid name").
* wooldid has three estimation paths that DO NOT go through reghdfe - `unconditionalse',
* `customvce()' and `ccnoabsorb', which use plain `regress'.  `regress' accepts any weight
* expression.  If the expression reaches wooldid's own estimand-weight arithmetic there,
*     ado:1619  gegen `w_per_it' = total(`esamp' * `w') , by(`i' `t')   -> total(esamp*ka+kb)
*     ado:1644  gen double `w_att_it' = `w' `customweighttext'/`w_per_it' -> ka + kb/W_it
* the cell-weighted estimands are silently wrong.  ka+kb == k exactly, so [aw=ka+kb] must
* equal [aw=k].
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

* does plain regress accept the expression? (the premise)
use "p8.dta", clear
cap noi qui reg y post i.i i.t [aw=ka+kb]
dis "REG_EXPR_RC=" _rc
cap noi qui reghdfe y post [aw=ka+kb], absorb(i t)
dis "REGHDFE_EXPR_RC=" _rc
cap noi qui poisson y post i.i i.t [pw=ka+kb]
dis "POISSON_EXPR_RC=" _rc

program define pairrun
  syntax , ca(string asis) cb(string asis) tag(string)
  tempname b1 b2
  qui use "p8.dta", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 matrix `b1' = e(b)
  qui use "p8.dta", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 matrix `b2' = e(b)
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    mata: wd_rel("`b1'","`b2'","`tag'")
    dis "BA_`tag':"
    matlist `b1', format(%18.12g)
    dis "BB_`tag':"
    matlist `b2', format(%18.12g)
  }
end

* the non-reghdfe paths
pairrun , ca(wooldid y i t ttre [aw=k], att unconditionalse robust) cb(wooldid y i t ttre [aw=ka+kb], att unconditionalse robust) tag(UNCOND_ATT)
pairrun , ca(wooldid y i t ttre [aw=k], att att_it att_i att_itime unconditionalse robust safety(off)) cb(wooldid y i t ttre [aw=ka+kb], att att_it att_i att_itime unconditionalse robust safety(off)) tag(UNCOND_ALL)
pairrun , ca(wooldid y i t ttre [aw=k], att customvce(robust)) cb(wooldid y i t ttre [aw=ka+kb], att customvce(robust)) tag(CVCE_ATT)
pairrun , ca(wooldid y i t ttre [aw=k], att att_it att_i att_itime customvce(robust)) cb(wooldid y i t ttre [aw=ka+kb], att att_it att_i att_itime customvce(robust)) tag(CVCE_ALL)
pairrun , ca(wooldid y i t ttre [aw=k], att att_it customvce(cluster i)) cb(wooldid y i t ttre [aw=ka+kb], att att_it customvce(cluster i)) tag(CVCE_CL)
* ccnoabsorb needs contreat; try it with the dose
pairrun , ca(wooldid y i t ttre [aw=k], att att_it contreat(dose) contreatcontrols(basic) ccnoabsorb) cb(wooldid y i t ttre [aw=ka+kb], att att_it contreat(dose) contreatcontrols(basic) ccnoabsorb) tag(CCNOABS)
* poisson path (ppmlhdfe)
pairrun , ca(wooldid y i t ttre [pw=k], att poisson) cb(wooldid y i t ttre [pw=ka+kb], att poisson) tag(POISEXPR)

dis "PROBE17_DONE"
