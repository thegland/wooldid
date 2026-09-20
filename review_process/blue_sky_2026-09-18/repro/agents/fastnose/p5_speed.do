* p5 -- is fastnose 2-4x faster, as sthlp:288 claims?  Measure on realistic specs.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fastnose"
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gN 400
local gT 14
local gonsets "2003 2005 2007 2009 2011"
local gysd 0.30
local gdose unitconst
local gpois 1
include "`GEN'/gen_core.do"
dis "NOBS=" _N
qui save "`WD'/p5panel.dta", replace

local spA "y cohort t ttre , fe(unit) att att_it att_i att_itime esprelength(4) espostlength(5)"
local spB "y cohort t ttre , fe(unit) att att_it att_i att_itime esprelength(4) espostlength(5) subgroups(sg2)"
local spC "ycount cohort t ttre , poisson fe(unit) esprelength(4) espostlength(5)"
local spD "y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) att att_it esprelength(4) espostlength(5)"

timer clear
local k 0
foreach arm in A B C D {
  local ++k
  local tA = 2*`k' - 1
  local tF = 2*`k'
  use "`WD'/p5panel.dta", clear
  timer on `tA'
  cap noi qui wooldid `sp`arm''
  timer off `tA'
  dis "SPEED_`arm'_AN_RC=" _rc
  use "`WD'/p5panel.dta", clear
  timer on `tF'
  cap noi qui wooldid `sp`arm'' fastnose
  timer off `tF'
  dis "SPEED_`arm'_FN_RC=" _rc
}
timer list
qui timer list
local k 0
foreach arm in A B C D {
  local ++k
  local tA = 2*`k' - 1
  local tF = 2*`k'
  local ta = r(t`tA')
  local tf = r(t`tF')
  local rat = `ta'/`tf'
  dis "SPEEDUP_`arm': analytic=`ta' fastnose=`tf' ratio=`rat'"
}
dis "DONE_P5"
