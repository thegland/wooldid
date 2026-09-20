* p5b -- the 2-4x claim, measured in the regime the help names: "specification search on
* LARGE panels, where the delta-method variance calculation inside margins dominates".
* i = unit (so K, the number of cohort-by-period parameters, is large).  Each arm is run
* twice, alternating AN/FN, and the MINIMUM of the two is used, because the machine is
* running 7-9 other Stata jobs.
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

local gN 200
local gT 16
local gonsets "2004 2006 2008 2010"
local gysd 0.30
local gdose unitconst
local gpois 1
include "`GEN'/gen_core.do"
dis "NOBS=" _N
qui save "`WD'/p5bpanel.dta", replace

local spA "y unit t ttre , att att_it att_i att_itime esprelength(5) espostlength(6)"
local spB "y unit t ttre , att att_it esprelength(5) espostlength(6) subgroups(sg2)"
local spC "ycount unit t ttre , poisson esprelength(5) espostlength(6)"
local spD "y unit t ttre , contreat(dose) contreatcontrols(untreatedZero) att att_it esprelength(5) espostlength(6)"
local spE "y unit t ttre , esprelength(5) espostlength(6)"

timer clear
local k 0
foreach arm in A B C E {
  local ++k
  forval rep = 1/2 {
    local tA = (`k'-1)*4 + `rep'
    local tF = (`k'-1)*4 + 2 + `rep'
    use "`WD'/p5bpanel.dta", clear
    timer on `tA'
    cap noi qui wooldid `sp`arm''
    timer off `tA'
    local rcA = _rc
    use "`WD'/p5bpanel.dta", clear
    timer on `tF'
    cap noi qui wooldid `sp`arm'' fastnose
    timer off `tF'
    dis "SPEED_`arm'_rep`rep'_rcA=`rcA' rcF=" _rc
  }
}
timer list
qui timer list
local k 0
foreach arm in A B C E {
  local ++k
  local a1 = r(t`=(`k'-1)*4+1')
  local a2 = r(t`=(`k'-1)*4+2')
  local f1 = r(t`=(`k'-1)*4+3')
  local f2 = r(t`=(`k'-1)*4+4')
  local am = min(`a1',`a2')
  local fm = min(`f1',`f2')
  dis "SPEEDUP_`arm': AN(`a1',`a2')->`am'  FN(`f1',`f2')->`fm'  ratio=" %6.3f `am'/`fm'
}
dis "DONE_P5B"
