*==========================================================================
* probeF.do -- metamorphic invariance of the permutation p-value
*              CONVENTIONS under transformations of the OUTCOME.
*   y -> c*y (c>0) must leave every p and t-ratio bit-identical.
*   y -> -y      must leave the equal-tailed, absolute and step-down
*                p-values bit-identical, and must swap permside(left)
*                with permside(right).
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

mkpanel , nunit(24) percl(2) nper(9)
save panelF2.dta, replace
qui su state
dis "NSTATES=" r(max)
qui levelsof ttre, missing
dis "ONSETS=" r(levels)

local CORE att att_i cluster(state) permute(shuffle) permreps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0) replace

*--------- permtstat, three outcome scales
foreach s in 1 1000 0.000001 {
  use panelF2.dta, clear
  qui replace y = y * `s'
  local tg = cond(`s'==1, "F1ts1", cond(`s'==1000, "F2tsK", "F3tsm"))
  cap noi wooldid y unit t ttre , `CORE' permtstat
  dis "RC_`tg'=" _rc
  cap dmp `tg'
}

*--------- permabsp, y and -y
use panelF2.dta, clear
cap noi wooldid y unit t ttre , `CORE' permabsp
dis "RC_F4abs=" _rc
cap dmp F4abs
use panelF2.dta, clear
qui replace y = -y
cap noi wooldid y unit t ttre , `CORE' permabsp
dis "RC_F5absneg=" _rc
cap dmp F5absneg

*--------- permside(left) on y vs permside(right) on -y
use panelF2.dta, clear
cap noi wooldid y unit t ttre , `CORE' permside(left)
dis "RC_F6left=" _rc
cap dmp F6left
use panelF2.dta, clear
qui replace y = -y
cap noi wooldid y unit t ttre , `CORE' permside(right)
dis "RC_F7rightneg=" _rc
cap dmp F7rightneg

*--------- permstepdown, y and -y
use panelF2.dta, clear
cap noi wooldid y unit t ttre , `CORE' permstepdown
dis "RC_F8sd=" _rc
cap dmp F8sd
use panelF2.dta, clear
qui replace y = -y
cap noi wooldid y unit t ttre , `CORE' permstepdown
dis "RC_F9sdneg=" _rc
cap dmp F9sdneg

*--------- simultaneousbands, y and -y
use panelF2.dta, clear
cap noi wooldid y unit t ttre , `CORE' simultaneousbands
dis "RC_F10sb=" _rc
cap dmp F10sb
use panelF2.dta, clear
qui replace y = -y
cap noi wooldid y unit t ttre , `CORE' simultaneousbands
dis "RC_F11sbneg=" _rc
cap dmp F11sbneg

*--------- default convention, y and -y and 1000*y (control for probeB's panel)
use panelF2.dta, clear
cap noi wooldid y unit t ttre , `CORE'
dis "RC_F12base=" _rc
cap dmp F12base
use panelF2.dta, clear
qui replace y = -y
cap noi wooldid y unit t ttre , `CORE'
dis "RC_F13neg=" _rc
cap dmp F13neg

dis "PROBEF_DONE"
