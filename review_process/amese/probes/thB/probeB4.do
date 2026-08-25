* theoryB probe 4: d3 deterministic behavior (cell-constant dose) + d1 robust/unit-cluster sanity
clear all
set more off
set maxvar 20000
qui do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"
set seed 24681357
set sortseed 24681357

* base panel
set obs 180
gen i = ceil(_n/6)
gen unitid = 1000 + _n
gen udose = rnormal(1,1)          // unit-level dose for d1
gen cohfe = rnormal(0,1)
bysort i (unitid): replace cohfe = cohfe[1]
expand 12
bysort unitid: gen t = _n
gen ttre = 6 if mod(i,2)==0
gen y = cohfe + .05*t + rnormal(0,1)
* d3a: cohort-level dose, time-invariant
gen d3adose = mod(i,7) + 0.5
* d3b: cohort-by-period dose (cell-constant, time-varying)
gen d3bdose = mod(i,7) + 0.3*mod(t,4) + 0.5
tempfile base
save `base'

capture program drop f2check
program define f2check
  foreach mm in b_es_ame V_es_ame {
    cap matrix chk = e(`mm')
    if _rc == 0 & rowsof(chk) + colsof(chk) > 2 dis as result "  e(`mm'): " rowsof(chk) " x " colsof(chk)
    else dis as result "  e(`mm'): ABSENT"
    cap matrix drop chk
  }
  cap matrix chk2 = e(V_es_ame)
  if colsof(chk2) > 1 {
    mata: V = st_matrix("chk2"); printf("  diag min %10.6g max %10.6g\n", min(diagonal(V)), max(diagonal(V)))
    mata: ev = symeigenvalues(st_matrix("chk2")); printf("  eig max %10.6g min %10.3g\n", max(ev), min(ev))
  }
  cap matrix drop chk2
  dis as result "  adversarial_pval_ame = " e(adversarial_pval_ame)
end

use `base', clear
dis as result "===== D3A basic cluster(i)  [cohort-constant dose] ====="
cap noi wooldid y i t ttre, att cluster(i) contreat(d3adose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
dis as result "RC = " _rc
if _rc == 0 f2check

use `base', clear
dis as result "===== D3A basic robust ====="
cap noi wooldid y i t ttre, att robust contreat(d3adose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
dis as result "RC = " _rc
if _rc == 0 f2check

use `base', clear
dis as result "===== D3B basic cluster(i)  [cell-constant, time-varying dose] ====="
cap noi wooldid y i t ttre, att cluster(i) contreat(d3bdose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
dis as result "RC = " _rc
if _rc == 0 f2check

use `base', clear
dis as result "===== D1 basic robust ====="
cap noi wooldid y i t ttre, att robust contreat(udose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
dis as result "RC = " _rc
if _rc == 0 f2check

use `base', clear
dis as result "===== D1 basic cluster(unitid) ====="
cap noi wooldid y i t ttre, att cluster(unitid) contreat(udose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
dis as result "RC = " _rc
if _rc == 0 f2check

dis as result "PROBE B4 COMPLETE"
