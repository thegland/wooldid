* theoryB verification probe (single dataset, NOT a Monte Carlo)
* Purpose: verify three deterministic/algebraic structural predictions on one d1-style draw:
*   (P1) byCohort + cluster(i): visible-slope CRVE block annihilated -> ame SEs ~ machine zero
*   (P2) byTime  + cluster(i): V_es_ame full rank, horizon-varying SEs, adversarial POSTS
*   (P3) within-arm (basic+contreatwithin+safety(off)) numerically identical to basic
* plus re-confirm (P0) basic + cluster(i) rank-1 collapse with all SEs = SE(phihat).
clear all
set more off
set maxvar 20000
* IMPORTANT: a stale wooldid.ado (2024-07-31) lives in the PLUS adopath and wins over
* `adopath +` (which appends). Load the repo build directly so in-memory programs win.
qui do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"
set seed 987654321
set sortseed 987654321

* ---- d1-style DGP: 30 cohorts, 6 units/cohort, t=1..12, even i treated at t=6
set obs 180
gen i = ceil(_n/6)
gen unitid = 1000 + _n
gen dose = rnormal(1,1)          // unit-level dose, constant within unit over time (d1)
gen cohfe = rnormal(0,1)
bysort i (unitid): replace cohfe = cohfe[1]
expand 12
bysort unitid: gen t = _n
gen ttre = 6 if mod(i,2)==0
gen y = cohfe + .05*t + rnormal(0,1)   // dose NEVER enters y
tempfile d1
save `d1'

capture program drop vsummary
program define vsummary
  * summarize e(V_es_ame): diag range, eigenvalues, min offdiag corr
  cap matrix VV = e(V_es_ame)
  if _rc != 0 {
    dis as result "V_es_ame: NOT POSTED"
    exit
  }
  matrix bb = e(b_es_ame)
  local k = colsof(VV)
  dis as result "V dims: " rowsof(VV) " x " `k'
  * count missing entries; build a missing-patched copy for eigen
  local nmiss 0
  matrix VP = VV
  forval a = 1/`k' {
    forval bq = 1/`k' {
      if VV[`a',`bq'] >= . {
        local ++nmiss
        matrix VP[`a',`bq'] = 0
      }
    }
  }
  dis as result "n missing entries: " `nmiss'
  local dmin = .
  local dmax = .
  forval a = 1/`k' {
    if VV[`a',`a'] < `dmin' local dmin = VV[`a',`a']
    if (VV[`a',`a'] > `dmax' & VV[`a',`a'] < .) | `dmax' == . local dmax = VV[`a',`a']
  }
  local cmin = .
  local cmax = .
  forval a = 1/`k' {
    forval bq = 1/`k' {
      if `a' < `bq' & VV[`a',`bq'] < . & VV[`a',`a'] > 0 & VV[`bq',`bq'] > 0 {
        local cc = VV[`a',`bq']/sqrt(VV[`a',`a']*VV[`bq',`bq'])
        if `cc' < `cmin' | `cmin' == . local cmin = `cc'
        if (`cc' > `cmax' & `cc' < .) | `cmax' == . local cmax = `cc'
      }
    }
  }
  dis as result "V diag min/max: " `dmin' " / " `dmax'
  dis as result "offdiag corr min/max: " `cmin' " / " `cmax'
  cap matrix symeigen X L = VP
  if _rc == 0 {
    dis as result "eigenvalues (missing-patched V):"
    matrix list L, format(%14.0g)
  }
  else dis as result "symeigen failed rc " _rc
  dis as result "raw V_es_ame:"
  matrix list VV, format(%12.0g)
  dis as result "b_es_ame:"
  matrix list bb, format(%12.0g)
  dis as result "adversarial_pval_ame = " e(adversarial_pval_ame)
  dis as result "es_ndegen_ame = " e(es_ndegen_ame)
end

* ---------------- RUN 1: basic + cluster(i)
use `d1', clear
dis as result "===== RUN1 basic cluster(i) ====="
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
vsummary
matrix B1 = e(b_es_ame)
matrix V1 = e(V_es_ame)

* ---------------- RUN 2: byCohort + cluster(i)
use `d1', clear
dis as result "===== RUN2 byCohort cluster(i) ====="
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(byCohort) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
vsummary

* ---------------- RUN 3: byTime + cluster(i)
use `d1', clear
dis as result "===== RUN3 byTime cluster(i) ====="
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(byTime) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
vsummary

* ---------------- RUN 4: within-arm = basic + contreatwithin + safety(off) + cluster(i)
use `d1', clear
dis as result "===== RUN4 within cluster(i) ====="
wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) contreatwithin safety(off) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
vsummary
matrix B4 = e(b_es_ame)
matrix V4 = e(V_es_ame)
dis as result "within-vs-basic mreldif(b) = " mreldif(B4,B1)
dis as result "within-vs-basic mreldif(V) = " mreldif(V4,V1)

* ---------------- RUN 5: byCohort + robust (sanity: full rank, positive SEs)
use `d1', clear
dis as result "===== RUN5 byCohort robust ====="
wooldid y i t ttre, att robust contreat(dose) contreatcontrols(byCohort) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
vsummary

dis as result "PROBE COMPLETE"
