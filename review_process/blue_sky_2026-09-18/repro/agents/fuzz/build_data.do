clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local OUT    "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

file open mf using "`OUT'/manifest.csv", write replace
file write mf "name,nobs,nunits,ncohort,nnever,ndosepos,hascount,desc" _n

local SHAPES canon unbal rotate admin rcs biennial fewG10 pw_survey miss_y dose_spike dose_lattice dose_cellconst poisson nevergrp longT hetero_unit sg_gaps onsetgap clustonset
foreach S of local SHAPES {
  local SHAPE  "`S'"
  local GOVR   "gN 60 gT 8"
  if "`S'" == "longT"  local GOVR "gN 40 gT 14"
  if "`S'" == "rcs"    local GOVR "gN 30 gT 7"
  if "`S'" == "admin"  local GOVR "gN 60 gT 8"
  capture noisily include "`GENDIR'/gen_shapes.do"
  if _rc {
    dis "GENFAIL `S' rc=" _rc
    continue
  }
  qui count
  local nobs = r(N)
  qui levelsof unit, local(uu)
  local nunits : word count `uu'
  qui levelsof cohort if !missing(ttre), local(cc)
  local ncoh : word count `cc'
  qui count if missing(ttre)
  local nnev = r(N)
  local ndp = 0
  capture qui count if dose > 0 & !missing(dose)
  if !_rc local ndp = r(N)
  local hc = 0
  capture confirm variable ycount
  if !_rc local hc = 1
  qui save "`OUT'/ds_`S'.dta", replace
  file write mf "`S',`nobs',`nunits',`ncoh',`nnev',`ndp',`hc',`SHAPEDESC'" _n
  dis "BUILT `S' nobs=`nobs' nunits=`nunits' ncoh=`ncoh' nnev=`nnev' ndose=`ndp' cnt=`hc'"
}
file close mf
type "`OUT'/manifest.csv"
