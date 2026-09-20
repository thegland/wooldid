*==========================================================================
* probeB.do -- metamorphic invariance of permute() at fixed rseed:
*              DATA transformations (row duplication, if-vs-drop,
*              affine remap of (t,ttre), affine remap of y, weight scale).
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

local OPT att att_i cluster(state) permute(shuffle) permreps(60) rseed(20260917) esprelength(2) espostlength(3) dots(0) replace

mkpanel
save panelB.dta, replace

*---------------------------------------------------------------- B0 base (no weights)
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B0)
dis "RC_B0=" _rc
dmp B0
asgdmp B0

*---------------------------------------------------------------- B1 every row duplicated
use panelB.dta, clear
qui expand 2
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B1dup)
dis "RC_B1dup=" _rc
cap dmp B1dup
cap asgdmp B1dup

*---------------------------------------------------------------- B2 if-restricted
use panelB.dta, clear
cap noi wooldid y unit t ttre if state <= 9 , `OPT' permsave(ps_B2if)
dis "RC_B2if=" _rc
cap dmp B2if
cap asgdmp B2if

*---------------------------------------------------------------- B3 physically dropped
use panelB.dta, clear
qui drop if state > 9
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B3drop)
dis "RC_B3drop=" _rc
cap dmp B3drop
cap asgdmp B3drop

*---------------------------------------------------------------- B4 (t,ttre) shifted by +50
use panelB.dta, clear
qui replace t    = t + 50
qui replace ttre = ttre + 50
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B4tsh)
dis "RC_B4tsh=" _rc
cap dmp B4tsh
cap asgdmp B4tsh

*---------------------------------------------------------------- B5 y + 1000
use panelB.dta, clear
qui replace y = y + 1000
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B5yadd)
dis "RC_B5yadd=" _rc
cap dmp B5yadd
cap asgdmp B5yadd

*---------------------------------------------------------------- B6 y * 1000
use panelB.dta, clear
qui replace y = y * 1000
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B6ymul)
dis "RC_B6ymul=" _rc
cap dmp B6ymul
cap asgdmp B6ymul

*---------------------------------------------------------------- B7 y * (-1)
use panelB.dta, clear
qui replace y = -y
cap noi wooldid y unit t ttre , `OPT' permsave(ps_B7yneg)
dis "RC_B7yneg=" _rc
cap dmp B7yneg
cap asgdmp B7yneg

*---------------------------------------------------------------- B8 aweights, and B9 the same weights x 1000
use panelB.dta, clear
cap noi wooldid y unit t ttre [aw=w] , `OPT' permsave(ps_B8w)
dis "RC_B8w=" _rc
cap dmp B8w
cap asgdmp B8w

use panelB.dta, clear
qui replace w = w*1000
cap noi wooldid y unit t ttre [aw=w] , `OPT' permsave(ps_B9w1k)
dis "RC_B9w1k=" _rc
cap dmp B9w1k
cap asgdmp B9w1k

dis "PROBEB_DONE"
