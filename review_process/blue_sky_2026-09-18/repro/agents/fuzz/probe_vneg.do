* Does sthlp:1616-1617's claim hold?  "each affected estimate carries a placeholder
* entry that displays as a missing standard error (and that test DROPS and lincom
* PROPAGATES AS MISSING). The placeholder is stored as a large negative number"
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

use "`D'/ds_dose_lattice.dta", clear
ereturn clear
qui wooldid y unit t ttre if region<=3 [pw=w], att_i cluster(region) controls(x1) interactivecontrols(x1) fe(i.region) espostlength(4) contreat(dose) contreatcontrols(untreatedZero) contreatpoly(1) lattice(byCohort 3) ccnoabsorb contreatelasticitytype(dydx) makecfxplots cfxplotnose jointtests summarystats anticipation(2) oostreatedcontrols(1) level(99)
dis "RC=" _rc
matrix VV = e(V)
matrix BB = e(b)
local K = colsof(BB)
local cn : colnames BB
local ce : coleq BB
dis "K=`K'"
local affected ""
local clean ""
forvalues k = 1/`K' {
  local nm : word `k' of `cn'
  local eq : word `k' of `ce'
  local d = VV[`k',`k']
  dis "COL `k' : `eq':`nm'  V=" %21x `d' "  = " `d'
  if `d' < -1e50 local affected "`affected' `eq':`nm'"
  else if `d' > 0 local clean "`clean' `eq':`nm'"
}
dis "AFFECTED:`affected'"
dis "CLEAN:`clean'"

local a1 : word 1 of `affected'
local c1 : word 1 of `clean'
dis "===== test on an AFFECTED coefficient (`a1') ====="
capture noisily test _b[`a1'] = 0
dis "TEST_RC=" _rc "  p=" r(p) "  F=" r(F) " chi2=" r(chi2) " df=" r(df)
dis "===== lincom on an AFFECTED coefficient ====="
capture noisily lincom _b[`a1']
dis "LINCOM_RC=" _rc "  se=" r(se) " estimate=" r(estimate) " lb=" r(lb) " ub=" r(ub)
dis "===== lincom mixing a CLEAN and an AFFECTED coefficient ====="
capture noisily lincom _b[`c1'] + _b[`a1']
dis "LINCOM2_RC=" _rc "  se=" r(se) " estimate=" r(estimate)
dis "===== test joint of clean+affected ====="
capture noisily test (_b[`c1'] = 0) (_b[`a1'] = 0)
dis "TESTJ_RC=" _rc " p=" r(p) " df=" r(df)
dis "===== nlcom on an AFFECTED coefficient ====="
capture noisily nlcom 2*_b[`a1']
dis "NLCOM_RC=" _rc
dis "===== test on a CLEAN coefficient (control) ====="
capture noisily test _b[`c1'] = 0
dis "TESTC_RC=" _rc " p=" r(p)
