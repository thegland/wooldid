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
  di as result _n "======== BUILD `S' ========"
  local SHAPE `S'
  capture noisily include "`GENDIR'/gen_shapes.do"
  if _rc != 0 {
    di as error "BUILDFAIL `S' rc=" _rc
    continue
  }
  qui count
  local nobs = r(N)
  tempvar tu
  qui bysort unit: gen byte `tu' = (_n == 1)
  qui count if `tu' == 1
  local nun = r(N)
  qui count if `tu' == 1 & ttre < .
  local ntre = r(N)
  qui count if `tu' == 1 & ttre >= .
  local nnev = r(N)
  qui levelsof t, local(tl)
  local nt : word count `tl'
  qui levelsof ttre, local(ol)
  qui summarize w, meanonly
  local wmin = r(min)
  local wmax = r(max)
  qui summarize y, meanonly
  local ymn = r(mean)
  qui count if y >= .
  local nymiss = r(N)
  di as text "BUILDOK `S' | rows=`nobs' units=`nun' tre=`ntre' nev=`nnev' T=`nt'" ///
     " onsets=`ol' wrange=`wmin'/`wmax' ymean=" %8.4f `ymn' " ymiss=`nymiss'"
  drop `tu'
}

di as result _n "======== MUTATORS on canon ========"
local MUTS dupall dupcohort halvew relabel reorder shift_t scale_t scale_w ///
           tofloat extmiss ifflag negw zerow dropmid
foreach M of local MUTS {
  local SHAPE canon
  qui include "`GENDIR'/gen_shapes.do"
  local MUT `M'
  capture noisily include "`GENDIR'/gen_mutate.do"
  if _rc != 0 di as error "MUTFAIL `M' rc=" _rc
  else di as text "MUTOK `M' rows=" _N
}
* cellmean needs a repeated cross-section
local SHAPE rcs
qui include "`GENDIR'/gen_shapes.do"
local MUT cellmean
capture noisily include "`GENDIR'/gen_mutate.do"
if _rc != 0 di as error "MUTFAIL cellmean rc=" _rc
else di as text "MUTOK cellmean rows=" _N

di as result _n "======== TRUTH on canon_exact ========"
local SHAPE canon_exact
qui include "`GENDIR'/gen_shapes.do"
qui gen byte insmp = 1
local TRUTHSMP insmp
local TRUTHI unit
local TRUTHW w
capture noisily include "`GENDIR'/gen_truth.do"
if _rc != 0 di as error "TRUTHFAIL rc=" _rc

di as result _n "ALLDONE"
