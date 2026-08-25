* probeA.do - theory-A mechanical verification probes (single dataset draws, NO Monte Carlo)
* verifies: (A) byCohort+cluster(i) rank-0 zero-V prediction + FWL identity gamma_gh = s_gh - s_g,base
*           (B) byTime+cluster(i) full-rank horizon-varying prediction
*           (C) within == basic bit-identity
*           (D) d3 (cell-constant dose) drop behavior: slopes-omitted caution vs error 322
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"

* ---- d1 dataset (grid DGP, single draw) ----
set seed 2010
qui set obs 30
gen i = _n
gen ttre = cond(mod(_n,2)==0, 6, .)
qui expand 6
bysort i: gen unit = _n
gen dose = exp(0.8*runiform() + 0.2*mod(i,5))
qui expand 12
bysort i unit: gen t = _n
gen y = 0.1*i + 0.05*t + rnormal()
gen unitid = (i-1)*6 + unit
sort i unit t
tempfile d1
save `d1'

* ================= RUN D: basic + cluster(i) (baseline) =================
dis _n as result "==== RUN D: basic + cluster(i) ===="
qui wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod
matrix BD = e(b_es_ame)
matrix VD = e(V_es_ame)
dis as text "colnames b_es_ame: " as result "`: colnames BD'"
forval j = 1/8 {
  dis as text "  h`j'  b = " %10.6f BD[1,`j'] "  se = " %14.10f sqrt(VD[`j',`j'])
}
mata: V = st_matrix("VD"); printf("maxdiag = %12.6e ; max|offdiag-diag spread| = %12.6e\n", max(diagonal(V)), max(abs(V)) - min(abs(V)))

* ================= RUN A: byCohort + cluster(i) =================
use `d1', clear
dis _n as result "==== RUN A: byCohort + cluster(i) ===="
cap noi wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(byCohort) esprelength(4) espostlength(4) esfixedbaseperiod
dis as result "rc = " _rc
matrix BA = e(b_es_ame)
matrix VA = e(V_es_ame)
forval j = 1/8 {
  dis as text "  h`j'  b = " %10.6f BA[1,`j'] "  se = " %14.10e sqrt(abs(VA[`j',`j']))
}
mata: V = st_matrix("VA"); printf("byCohort cluster(i): max|V| entry = %12.6e (basic maxdiag was above)\n", max(abs(V)))
dis as text "reported SE on main ame / pretreatment ame rows shown in table above"

* ---- FWL identity check: b_es_ame[post0] =?= mean_g( s_g,t6 - s_g,t5 ) ----
dis _n as result "---- FWL check: AME_post0 vs mean_g(s_g6 - s_g5) ----"
use `d1', clear
local acc = 0
forval g = 2(2)30 {
  qui reg y dose if i==`g' & t==6
  local s6 = _b[dose]
  qui reg y dose if i==`g' & t==5
  local s5 = _b[dose]
  local acc = `acc' + (`s6' - `s5')
}
local manual = `acc'/15
* post0 should be column 5 (pre5 pre4 pre3 pre2 post0 ...)
dis as result "manual mean_g(s_g6-s_g5) = " %12.8f `manual' "   b_es_ame col5 = " %12.8f BA[1,5] "   diff = " %12.4e (`manual' - BA[1,5])

* ================= RUN B: byTime + cluster(i) =================
use `d1', clear
dis _n as result "==== RUN B: byTime + cluster(i) ===="
cap noi wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(byTime) esprelength(4) espostlength(4) esfixedbaseperiod
dis as result "rc = " _rc
matrix BB = e(b_es_ame)
matrix VB = e(V_es_ame)
forval j = 1/8 {
  dis as text "  h`j'  b = " %10.6f BB[1,`j'] "  se = " %14.10f sqrt(VB[`j',`j'])
}
mata: V = st_matrix("VB"); ev = symeigenvalues(V); printf("byTime cluster(i): eig max = %10.6g  eig min = %10.6g (full rank if min >> 0)\n", max(ev), min(ev))
mata: C = corr(st_matrix("VB")); _diag(C, J(8,1,0)); printf("max |offdiag corr| = %10.6g\n", max(abs(C)))

* ================= RUN C: within (basic+contreatwithin+safety(off)) + cluster(i) =================
use `d1', clear
dis _n as result "==== RUN C: within + cluster(i) ===="
cap noi wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) contreatwithin safety(off) esprelength(4) espostlength(4) esfixedbaseperiod
dis as result "rc = " _rc
matrix BC = e(b_es_ame)
matrix VC = e(V_es_ame)
dis as result "within vs basic: mreldif(b) = " %12.4e mreldif(BC, BD) "   mreldif(V) = " %12.4e mreldif(VC, VD)

* ================= RUN E: byCohort + robust (structure check) =================
use `d1', clear
dis _n as result "==== RUN E: byCohort + robust ===="
qui wooldid y i t ttre, att robust contreat(dose) contreatcontrols(byCohort) esprelength(4) espostlength(4) esfixedbaseperiod
matrix VE = e(V_es_ame)
forval j = 1/8 {
  dis as text "  h`j'  se = " %12.8f sqrt(VE[`j',`j'])
}
mata: C = corr(st_matrix("VE")); _diag(C, J(8,1,0)); printf("byCohort robust: mean offdiag corr = %10.6g (predict ~ +0.4-0.6)\n", sum(C)/56)

* ================= RUN F: d3 (dose constant within cohort) + basic + cluster(i) =================
dis _n as result "==== RUN F: d3 cohort-constant dose, basic + cluster(i) ===="
clear
set seed 2010
qui set obs 30
gen i = _n
gen ttre = cond(mod(_n,2)==0, 6, .)
qui expand 6
bysort i: gen unit = _n
gen dose = exp(0.2*mod(i,5) + 0.1*i)
qui expand 12
bysort i unit: gen t = _n
gen y = 0.1*i + 0.05*t + rnormal()
sort i unit t
cap noi wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod
dis as result "d3 rc = " _rc
cap confirm matrix e(b_es_ame)
if _rc == 0 {
  matrix BF = e(b_es_ame)
  matrix VF = e(V_es_ame)
  mata: printf("d3: max|b| = %12.6e   max|V| = %12.6e\n", max(abs(st_matrix("BF"))), max(abs(st_matrix("VF"))))
}
else dis as result "d3: no e(b_es_ame) posted (run errored before F2)"

dis _n as result "PROBEA_DONE"
exit, clear
