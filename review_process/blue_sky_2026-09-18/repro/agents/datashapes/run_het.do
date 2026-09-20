clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/datashapes/gen"

* regression guard: canon must be UNCHANGED by the gtauhet addition
local SHAPE canon_exact
qui include "`GENDIR'/gen_shapes.do"
qui count
di as result "@@@@ GUARD canon_exact rows=" r(N)
qui summarize truetau, meanonly
di as result "@@@@ GUARD truetau mean = " %21x r(mean)
capture noisily wooldid ync unit t ttre, att att_i cluster(unit)
matrix bg = e(b)
di as result "@@@@ GUARD att = " %21x bg[1,1] " att_i = " %21x bg[1,3]

foreach S in hetero_unit hetero_unbal {
  di as result _n "@@@@ SMOKE `S'"
  local SHAPE `S'
  local gysd 0
  qui include "`GENDIR'/gen_shapes.do"
  qui count
  di as result "@@@@ BUILT `S' rows=" r(N)
  capture noisily wooldid ync unit t ttre, att att_i att_it att_itime cluster(unit)
  di as result "@@@@ RESULT `S' i=unit rc=" _rc
  if _rc == 0 {
    matrix bU = e(b)
    matrix list bU
    qui gen byte insmp = e(sample)
    local TRUTHSMP insmp
    local TRUTHI unit
    local TRUTHW w
    include "`GENDIR'/gen_truth.do"
    drop insmp
  }
  capture noisily wooldid ync cohort t ttre, att att_i att_it att_itime cluster(cohort)
  di as result "@@@@ RESULT `S' i=cohort rc=" _rc
  if _rc == 0 {
    matrix bC = e(b)
    matrix list bC
  }
}
di as result _n "HETDONE"
