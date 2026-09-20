clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/datashapes/gen"

foreach S in onepre onsetgap {
  di as result _n "@@@@ SMOKE `S'"
  local SHAPE `S'
  qui include "`GENDIR'/gen_shapes.do"
  qui count
  di as result "@@@@ BUILT `S' rows=" r(N)
  capture noisily wooldid ync unit t ttre, att att_i cluster(unit) espre(2) espost(2)
  di as result "@@@@ RESULT `S' rc=" _rc
}

* dose shapes with the required contreatcontrols()
foreach S in dose_spike dose_cellconst {
  di as result _n "@@@@ SMOKE `S'_cc"
  local SHAPE `S'
  qui include "`GENDIR'/gen_shapes.do"
  capture noisily wooldid ync unit t ttre, att contreat(dose) contreatcontrols(untreatedZero) cluster(unit)
  di as result "@@@@ RESULT `S'_cc rc=" _rc
}

* admin at i=cohort (the coarse-id route past r(907))
di as result _n "@@@@ SMOKE admin_icohort"
local SHAPE admin
qui include "`GENDIR'/gen_shapes.do"
capture noisily wooldid ync cohort t ttre, att att_i cluster(cohort) espre(2) espost(2)
di as result "@@@@ RESULT admin_icohort rc=" _rc

* poisson at i=cohort
di as result _n "@@@@ SMOKE poisson_icohort"
local SHAPE poisson
qui include "`GENDIR'/gen_shapes.do"
capture noisily wooldid ycount cohort t ttre, att poisson cluster(cohort)
di as result "@@@@ RESULT poisson_icohort rc=" _rc

* the exactness harness at i=cohort and i=unit on the SAME zero-noise data
di as result _n "@@@@ SMOKE exact_axis"
local SHAPE canon_exact
qui include "`GENDIR'/gen_shapes.do"
capture noisily wooldid ync unit t ttre, att att_i att_it att_itime cluster(unit)
di as result "@@@@ EXACT unit rc=" _rc
if _rc == 0 {
  matrix bU = e(b)
  matrix list bU, format(%21x)
}
capture noisily wooldid ync cohort t ttre, att att_i att_it att_itime cluster(cohort)
di as result "@@@@ EXACT cohort rc=" _rc
if _rc == 0 {
  matrix bC = e(b)
  matrix list bC, format(%21x)
}
* design truth under both id choices
qui gen byte insmp = 1
local TRUTHSMP insmp
local TRUTHI unit
local TRUTHW w
include "`GENDIR'/gen_truth.do"
local TRUTHI cohort
include "`GENDIR'/gen_truth.do"

* sgall: a one-level subgroup partition must reproduce the main row (invariance I13)
di as result _n "@@@@ SMOKE sgall"
local SHAPE canon_exact
qui include "`GENDIR'/gen_shapes.do"
capture noisily wooldid ync unit t ttre, att cluster(unit) subgroups(sgall)
di as result "@@@@ RESULT sgall rc=" _rc

* sgvary: a subgroup that changes within a unit must be refused or ruled on
di as result _n "@@@@ SMOKE sgvary"
capture noisily wooldid ync unit t ttre, att cluster(unit) subgroups(sgvary)
di as result "@@@@ RESULT sgvary rc=" _rc

* ---- does espre() move att?  (pre horizons leave the reference pool) ----
di as result _n "@@@@ SMOKE espre_moves_att"
local SHAPE canon
qui include "`GENDIR'/gen_shapes.do"
qui wooldid ync unit t ttre, att cluster(unit)
matrix b0 = e(b)
qui wooldid ync unit t ttre, att cluster(unit) espre(3) espost(3)
matrix b3 = e(b)
di as result "@@@@ ESPRE att espre0 = " %21x b0[1,1] "  espre3 = " %21x b3[1,1]
di as result "@@@@ ESPRE reldiff = " %12.3e (b3[1,1]-b0[1,1])/b0[1,1]

* ---- fast collapse vs nofastcollapse on a repeated cross-section ----
di as result _n "@@@@ SMOKE fcol_rcs"
local SHAPE rcs
qui include "`GENDIR'/gen_shapes.do"
qui wooldid ync unit t ttre [pw=w], att att_it cluster(unit) espre(2) espost(2)
matrix f1 = e(b)
local fc1 = e(wdid_fcol)
qui wooldid ync unit t ttre [pw=w], att att_it cluster(unit) espre(2) espost(2) nofastcollapse
matrix f0 = e(b)
di as result "@@@@ FCOL fcon=`fc1' collapsed = " %21x f1[1,1] "  full = " %21x f0[1,1]
di as result "@@@@ FCOL reldiff = " %12.3e (f1[1,1]-f0[1,1])/f0[1,1]

* ---- same with controls that vary WITHIN the cell (not in the collapse key) ----
di as result _n "@@@@ SMOKE fcol_rcs_ctrl"
qui wooldid y unit t ttre [pw=w], att cluster(unit) controls(x1 x2) espre(2) espost(2)
matrix g1 = e(b)
qui wooldid y unit t ttre [pw=w], att cluster(unit) controls(x1 x2) espre(2) espost(2) nofastcollapse
matrix g0 = e(b)
di as result "@@@@ FCOLC collapsed = " %21x g1[1,1] "  full = " %21x g0[1,1]
di as result "@@@@ FCOLC reldiff = " %12.3e (g1[1,1]-g0[1,1])/g0[1,1]

* ---- aweight vs pweight point estimates ----
di as result _n "@@@@ SMOKE awpw"
local SHAPE pw_survey
qui include "`GENDIR'/gen_shapes.do"
qui wooldid ync unit t ttre [aw=w], att cluster(unit)
matrix a1 = e(b)
qui wooldid ync unit t ttre [pw=w], att cluster(unit)
matrix p1 = e(b)
di as result "@@@@ AWPW aw = " %21x a1[1,1] "  pw = " %21x p1[1,1]
di as result "@@@@ AWPW reldiff = " %12.3e (a1[1,1]-p1[1,1])/p1[1,1]

di as result _n "EXTRADONE"
