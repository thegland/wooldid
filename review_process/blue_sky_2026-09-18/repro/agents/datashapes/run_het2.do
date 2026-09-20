clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/datashapes/gen"

* guard: canon_exact unchanged after the GOVR edit
local SHAPE canon_exact
local GOVR
qui include "`GENDIR'/gen_shapes.do"
qui wooldid ync unit t ttre, att att_i cluster(unit)
matrix bg = e(b)
di as result "@@@@ GUARD att = " %21x bg[1,1] " (want +1.3777777777781X+000)"

* exact recovery on an UNBALANCED, WEIGHTED, WITHIN-COHORT-HETEROGENEOUS design
foreach S in hetero_unbal unbal admin {
  di as result _n "@@@@ EXACTTEST `S'"
  local SHAPE `S'
  local GOVR "gysd 0"
  if "`S'" == "admin" local GOVR "gysd 0 gcounts 40 3 2 4 2 gnevern 60"
  qui include "`GENDIR'/gen_shapes.do"
  qui count
  di as result "@@@@ BUILT `S' rows=" r(N)
  capture noisily wooldid ync unit t ttre [aw=w], att att_i att_it att_itime cluster(unit)
  di as result "@@@@ RC `S' i=unit " _rc
  if _rc == 0 {
    matrix bU = e(b)
    qui gen byte insmp = e(sample)
    local TRUTHSMP insmp
    local TRUTHI unit
    local TRUTHW w
    include "`GENDIR'/gen_truth.do"
    di as result "@@@@ CMP `S' att     est=" %14.9f bU[1,1] " truth=" %14.9f `TRUE_att'     " rel=" %10.2e (bU[1,1]-`TRUE_att')/`TRUE_att'
    di as result "@@@@ CMP `S' att_it  est=" %14.9f bU[1,2] " truth=" %14.9f `TRUE_att_it'  " rel=" %10.2e (bU[1,2]-`TRUE_att_it')/`TRUE_att_it'
    di as result "@@@@ CMP `S' att_i   est=" %14.9f bU[1,3] " truth=" %14.9f `TRUE_att_i'   " rel=" %10.2e (bU[1,3]-`TRUE_att_i')/`TRUE_att_i'
    di as result "@@@@ CMP `S' att_itm est=" %14.9f bU[1,4] " truth=" %14.9f `TRUE_att_itime' " rel=" %10.2e (bU[1,4]-`TRUE_att_itime')/`TRUE_att_itime'
    drop insmp
  }
  capture noisily wooldid ync cohort t ttre [aw=w], att att_i att_it att_itime cluster(cohort)
  di as result "@@@@ RC `S' i=cohort " _rc
  if _rc == 0 {
    matrix bC = e(b)
    di as result "@@@@ COH `S' att=" %14.9f bC[1,1] " att_i=" %14.9f bC[1,3]
  }
}
di as result _n "HET2DONE"
