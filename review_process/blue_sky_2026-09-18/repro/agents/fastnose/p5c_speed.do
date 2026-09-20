* p5c -- the 2-4x claim in the LARGE-K regime the help names (i = unit, so the number of
* cohort-by-period parameters margins must differentiate is large).  Prints the ratio
* after EVERY pair, so a deadline kill still leaves usable numbers.
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

local gN 140
local gT 13
local gonsets "2003 2005 2007 2009"
local gysd 0.30
local gpois 1
include "`GEN'/gen_core.do"
dis "NOBS=" _N
qui save "`WD'/p5cpanel.dta", replace

local spE "y unit t ttre , esprelength(4) espostlength(4)"
local spA "y unit t ttre , att att_it att_i att_itime esprelength(4) espostlength(4)"
local spS "y unit t ttre , att att_it esprelength(4) espostlength(4) subgroups(sg2)"
local spC "ycount unit t ttre , poisson esprelength(4) espostlength(4)"

timer clear
local k 0
foreach arm in E A S C {
  local ++k
  local tA = 2*`k' - 1
  local tF = 2*`k'
  use "`WD'/p5cpanel.dta", clear
  timer on `tA'
  cap noi qui wooldid `sp`arm''
  local rcA = _rc
  timer off `tA'
  use "`WD'/p5cpanel.dta", clear
  timer on `tF'
  cap noi qui wooldid `sp`arm'' fastnose
  local rcF = _rc
  timer off `tF'
  qui timer list
  local ta = r(t`tA')
  local tf = r(t`tF')
  dis "SPEEDPAIR_`arm': rcA=`rcA' rcF=`rcF' analytic=`ta' fastnose=`tf' ratio=" %6.3f `ta'/`tf'
}
dis "DONE_P5C"
