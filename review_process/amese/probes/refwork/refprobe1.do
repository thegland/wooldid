* refprobe1.do — ADVERSARIAL REFEREE probe (1 of max 2)
* PART A: F1-predicate FALSE-NEGATIVE demo. d1 data; mixed cluster variable
*   cmix = i for cohorts 1..20, 100+unitid for cohorts 21..30. Treated cohorts are
*   even i: 10 of 15 have every cell CONTAINED in one cmix cluster, 5 of 15 have
*   cells SPLIT across 6 unit-clusters. F1's predicate as drafted ("EVERY treated
*   (i,relyr) cell inside a single cluster") does NOT fire here.
*   PREDICTION (written before running): cluster(cmix) posts a FULL-RANK,
*   horizon-varying V_es_ame with NO caution, but horizon SEs materially BELOW the
*   cluster(unitid) SEs on the same data (partial annihilation: contained cohorts'
*   s-noise estimated as 0; if the s-part dominates, ratio ~ sqrt(5/15) ~ 0.58).
*   Anchor runs on the SAME data: cluster(i) = rank-1 collapse; cluster(unitid) = full.
* PART B: byTime spoof equality robustness. DIFFERENT shape (20 cohorts x 9 units
*   x 12 periods) and fresh seed. PREDICTION: each reported horizon ame SE under
*   contreatcontrols(byTime)+cluster(i) equals SE_cl(phi_t) from a visible
*   ibn.t#c.ringer refit, to displayed digits (blind theorist saw this on 30x6x12).
clear all
set more off
set maxvar 20000
set linesize 240
do /Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado
do /Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh

* ============================ PART A ============================
set seed 99000001
amese_mkdata d1
gen double cmix = cond(i <= 20, i, 100 + unitid)
qui tab cmix
dis as result "PART A: n clusters cmix = " r(r)
save /Users/hegland/.claude/jobs/640b9294/tmp/amese/refwork/refA.dta, replace

* --- run 1: cluster(cmix) — the partial-containment case ---
dis as result "===== A1: cluster(cmix) partial containment ====="
wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) cluster(cmix)
matrix VMIX = e(V_es_ame)
mata:
V = st_matrix("VMIX")
printf("A1 cmix: ES-ame SEs (sqrt diag V_es_ame):\n")
sqrt(diagonal(V))'
printf("A1 cmix: eigenvalues:\n")
symeigenvalues(V)
st_matrix("SEMIX", sqrt(diagonal(V))')
end

* --- run 2: cluster(unitid) on the same data ---
use /Users/hegland/.claude/jobs/640b9294/tmp/amese/refwork/refA.dta, clear
dis as result "===== A2: cluster(unitid) all-split anchor ====="
wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) cluster(unitid)
matrix VCU = e(V_es_ame)
mata:
V = st_matrix("VCU")
printf("A2 cu: ES-ame SEs:\n")
sqrt(diagonal(V))'
printf("A2 cu: eigenvalues:\n")
symeigenvalues(V)
st_matrix("SECU", sqrt(diagonal(V))')
end

* --- run 3: cluster(i) on the same data — full-containment rank-1 anchor ---
use /Users/hegland/.claude/jobs/640b9294/tmp/amese/refwork/refA.dta, clear
dis as result "===== A3: cluster(i) full containment anchor ====="
wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) cluster(i)
matrix VCI = e(V_es_ame)
mata:
V = st_matrix("VCI")
printf("A3 ci: ES-ame SEs:\n")
sqrt(diagonal(V))'
printf("A3 ci: eigenvalues:\n")
symeigenvalues(V)
end

* --- ratio table ---
mata:
semix = st_matrix("SEMIX")
secu  = st_matrix("SECU")
printf("A RATIOS cmix/cu per horizon (prediction ~0.58 if s-part dominates):\n")
semix :/ secu
end

* ============================ PART B ============================
clear
set seed 99110001
qui set obs 20
gen i = _n
gen ttre = cond(mod(i,2)==0, 6, .)
gen a = rnormal()
qui expand 9
bysort i: gen unit = _n
gen unitid = (i-1)*9 + unit
gen dose = exp(0.8*runiform() + 0.2*mod(i,5))
qui expand 12
bysort i unit: gen t = _n
gen y = a + .05*t + rnormal()
drop a
sort i unitid t
save /Users/hegland/.claude/jobs/640b9294/tmp/amese/refwork/refB.dta, replace

dis as result "===== B1: byTime + cluster(i) on 20x9x12 ====="
wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(byTime) cluster(i)
matrix VBT = e(V_es_ame)
mata:
V = st_matrix("VBT")
printf("B1: reported ES-ame SEs (rows pre5 pre4 pre3 pre2 post0 post1 post2 post3):\n")
sqrt(diagonal(V))'
printf("B1: eigenvalues:\n")
symeigenvalues(V)
end

* --- visible per-period ringer refit ---
use /Users/hegland/.claude/jobs/640b9294/tmp/amese/refwork/refB.dta, clear
gen rel = t - ttre
gen _t = rel + 6
replace _t = 0 if _t == . | _t == 5
gen byte targetflag = _t > 0 & _t < .
gen double _i = i if targetflag
replace _i = 0 if _i == .
gen double ringer = dose
reghdfe y ibn.t#c.ringer i._i#i._t#i.targetflag##c.dose , absorb(i.i i.t) cluster(i) tol(1e-9)
dis as result "B2: visible phi_t SE_cl by calendar period (horizons: pre = t 1..4, post = t 6..9):"
forval tt = 1/12 {
  cap dis "t=`tt'  phi_t = " %12.8f _b[`tt'.t#c.ringer] "   SE_cl = " %12.8f _se[`tt'.t#c.ringer]
}
exit, clear
