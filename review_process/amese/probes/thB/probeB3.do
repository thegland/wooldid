* isolate: does the adversarial option toggle F2 posting on this d1-style draw?
clear all
set more off
adopath + "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall"
set seed 987654321
set sortseed 987654321
set obs 180
gen i = ceil(_n/6)
gen unitid = 1000 + _n
gen dose = rnormal(1,1)
gen cohfe = rnormal(0,1)
bysort i (unitid): replace cohfe = cohfe[1]
expand 12
bysort unitid: gen t = _n
gen ttre = 6 if mod(i,2)==0
gen y = cohfe + .05*t + rnormal(0,1)
tempfile d1
save `d1'

capture program drop f2check
program define f2check
  foreach mm in b_es_att V_es_att b_es_ame V_es_ame {
    cap matrix chk = e(`mm')
    if _rc == 0 & rowsof(chk) + colsof(chk) > 2 dis as result "  e(`mm'): " rowsof(chk) " x " colsof(chk)
    else dis as result "  e(`mm'): ABSENT"
    cap matrix drop chk
  }
  dis as result "  adversarial_pval_ame = " e(adversarial_pval_ame)
end

use `d1', clear
dis as result "===== A: basic cluster(i), NO adversarial ====="
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod
f2check
tempname VA
matrix `VA' = e(V_es_ame)
mata: V = st_matrix(st_local("VA")); C = corr(V); printf("A: dim %g, min corr = %12.9g\n", rows(V), min(C))
mata: ev = symeigenvalues(V); printf("A: eig max = %10.6g, eig2 = %10.3g, eig min = %10.3g\n", max(ev), sort(ev',1)[rows(V)-1], min(ev))
mata: printf("A: diag min %10.6g max %10.6g\n", min(diagonal(V)), max(diagonal(V)))

use `d1', clear
dis as result "===== B: basic cluster(i), WITH adversarial ====="
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
f2check

dis as result "PROBE B3 COMPLETE"
