*==========================================================================
* probeG.do -- subgroups() under permute(): relabeling the SUBGROUP VALUES
*              in an order-REVERSING way must permute the per-subgroup
*              tables and nothing else.  A surface that keys on the sorted
*              INDEX rather than the value would cross the two subgroups.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

*--- 40 counties in 10 states; a county-level subgroup constant within the unit
mkpanel
qui gen byte sg = 1 + mod(state, 2)
save panelG.dta, replace

local CG att att_i cluster(state) subgroups(sg) permute(shuffle) permreps(30) rseed(20260917) dots(0) replace

*--- G0: sg in {1,2}
cap noi wooldid y unit t ttre , `CG'
dis "RC_G0=" _rc
cap dmp G0

*--- G1: sg relabelled 1<->2 (order reversing)
use panelG.dta, clear
qui replace sg = 3 - sg
cap noi wooldid y unit t ttre , `CG'
dis "RC_G1=" _rc
cap dmp G1

*--- G2: sg relabelled order-PRESERVING to {5,9}
use panelG.dta, clear
qui replace sg = cond(sg==1, 5, 9)
cap noi wooldid y unit t ttre , att att_i cluster(state) subgroups(sg) permute(shuffle) permreps(30) rseed(20260917) dots(0) replace
dis "RC_G2=" _rc
cap dmp G2

dis "PROBEG_DONE"
