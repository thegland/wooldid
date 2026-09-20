clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/datashapes/gen"

local ALL canon canon_exact unbal rotate admin tinytreated rcs rcs_big ///
          biennial quinq holeyear longT shortT fewG5 fewG10 manyG ///
          pw_survey w_spread3 wzero miss_y miss_x sg_gaps sg_tiny sg_big ///
          dose_spike dose_skew dose_lattice dose_cellconst ///
          clustonset firstperiod afterend bothends pretrend antic ///
          poisson poissep nevergrp datesm datesq datesd ///
          ed_ttrevary ed_ttremix ed_dup ed_singleton ed_zerovar ed_outlier ///
          ed_badw ed_collin ed_ctrlconst ed_ctrlistreat ed_never1 ed_nopre ///
          ed_ifflag ed_float ed_tscale ed_onecell ed_tsingle

foreach S of local ALL {
  di as result _n "@@@@ SMOKE `S'"
  local SHAPE `S'
  qui include "`GENDIR'/gen_shapes.do"

  local YV ync
  local OPT att att_i cluster(unit) espre(2) espost(2)
  local WT
  if "`S'" == "poisson" | "`S'" == "poissep" {
    local YV ycount
    local OPT att poisson cluster(unit)
  }
  if substr("`S'", 1, 5) == "dose_" local OPT att contreat(dose) cluster(unit)
  if "`S'" == "rcs" | "`S'" == "rcs_big" | "`S'" == "pw_survey" | "`S'" == "w_spread3" ///
     | "`S'" == "wzero" | "`S'" == "ed_badw" local WT [pw=w]

  capture noisily wooldid `YV' unit t ttre `WT', `OPT'
  local rc = _rc
  local bb = .
  local nn = .
  if `rc' == 0 {
    tempname bmat
    cap matrix `bmat' = e(b)
    if _rc == 0 local bb = `bmat'[1, 1]
    local nn = e(N)
  }
  di as result "@@@@ RESULT `S' rc=`rc' N=`nn' b1=" %12.6f `bb'
}

di as result _n "======== extra legs ========"

* i = cohort (the clean-coverage axis) and finer-than-i clustering
local SHAPE nevergrp
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE nevergrp_icohort"
capture noisily wooldid ync cohort t ttre, att att_i cluster(cohort) espre(2) espost(2)
di as result "@@@@ RESULT nevergrp_icohort rc=" _rc
di as result "@@@@ SMOKE nevergrp_finercluster"
capture noisily wooldid ync cohort t ttre, att att_i cluster(unit) espre(2) espost(2)
di as result "@@@@ RESULT nevergrp_finercluster rc=" _rc

* repeated cross-sections with a finer-than-i household cluster
local SHAPE rcs
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE rcs_hhcluster"
capture noisily wooldid ync unit t ttre [pw=w], att cluster(hh) espre(2) espost(2)
di as result "@@@@ RESULT rcs_hhcluster rc=" _rc

* subgroup legs
local SHAPE sg_big
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE sg_big_sub"
capture noisily wooldid ync unit t ttre, att subgroups(sgbig) cluster(unit) espre(1) espost(1)
di as result "@@@@ RESULT sg_big_sub rc=" _rc
di as result "@@@@ WDIDSGMATS: `e(wdid_cfxmats)'"

local SHAPE sg_gaps
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE sg_gaps_sub"
capture noisily wooldid ync unit t ttre, att subgroups(sggap) cluster(unit) espre(1) espost(1)
di as result "@@@@ RESULT sg_gaps_sub rc=" _rc

* the if-flag rule
local SHAPE ed_ifflag
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE ifflag_run"
capture noisily wooldid ync unit t ttre if flag, att cluster(unit)
di as result "@@@@ RESULT ifflag_run rc=" _rc

* controls that are pathological
local SHAPE ed_collin
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE collin_controls"
capture noisily wooldid ync unit t ttre, att controls(x1 xcollin) cluster(unit)
di as result "@@@@ RESULT collin_controls rc=" _rc

local SHAPE ed_ctrlistreat
qui include "`GENDIR'/gen_shapes.do"
di as result "@@@@ SMOKE ctrl_is_treat"
capture noisily wooldid ync unit t ttre, att controls(xpost) cluster(unit)
di as result "@@@@ RESULT ctrl_is_treat rc=" _rc

di as result _n "SMOKEDONE"
