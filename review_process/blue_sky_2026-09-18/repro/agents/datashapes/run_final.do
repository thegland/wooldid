clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/datashapes/gen"

* guard: gaisd addition must leave canon_exact bit-identical
local SHAPE canon_exact
local GOVR
qui include "`GENDIR'/gen_shapes.do"
qui wooldid ync unit t ttre, att att_i cluster(unit)
matrix bg = e(b)
di as result "@@@@ GUARD att = " %21x bg[1,1] " (want +1.3777777777781X+000)"

* i = COHORT exactness harness: gaisd 0 removes the unit effects
local SHAPE hetero_unbal
local GOVR "gysd 0 gaisd 0"
qui include "`GENDIR'/gen_shapes.do"
capture noisily wooldid ync cohort t ttre [aw=w], att att_i att_it att_itime cluster(cohort)
di as result "@@@@ RC cohortexact " _rc
if _rc == 0 {
  matrix bC = e(b)
  qui gen byte insmp = e(sample)
  local TRUTHSMP insmp
  local TRUTHI cohort
  local TRUTHW w
  include "`GENDIR'/gen_truth.do"
  di as result "@@@@ COHCMP att   est=" %14.9f bC[1,1] " truth=" %14.9f `TRUE_att'   " rel=" %10.2e (bC[1,1]-`TRUE_att')/`TRUE_att'
  di as result "@@@@ COHCMP att_i est=" %14.9f bC[1,3] " truth=" %14.9f `TRUE_att_i' " rel=" %10.2e (bC[1,3]-`TRUE_att_i')/`TRUE_att_i'
  drop insmp
}

* att_i vs att_itime separation needs weight mass to vary ACROSS cells within i
local SHAPE rcs
local GOVR "gysd 0 gtauhet 0.4"
qui include "`GENDIR'/gen_shapes.do"
capture noisily wooldid ync unit t ttre [aw=w], att att_i att_it att_itime cluster(unit)
di as result "@@@@ RC rcs4 " _rc
if _rc == 0 {
  matrix b4 = e(b)
  di as result "@@@@ FOUR att=" %12.8f b4[1,1] " att_it=" %12.8f b4[1,2] " att_i=" %12.8f b4[1,3] " att_itime=" %12.8f b4[1,4]
  qui gen byte insmp = e(sample)
  local TRUTHSMP insmp
  local TRUTHI unit
  local TRUTHW w
  include "`GENDIR'/gen_truth.do"
  di as result "@@@@ FOURCMP rel att=" %10.2e (b4[1,1]-`TRUE_att')/`TRUE_att' " att_it=" %10.2e (b4[1,2]-`TRUE_att_it')/`TRUE_att_it' " att_i=" %10.2e (b4[1,3]-`TRUE_att_i')/`TRUE_att_i' " att_itime=" %10.2e (b4[1,4]-`TRUE_att_itime')/`TRUE_att_itime'
}
di as result _n "FINALDONE"
