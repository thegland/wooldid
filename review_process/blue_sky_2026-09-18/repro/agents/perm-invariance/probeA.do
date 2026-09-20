*==========================================================================
* probeA.do -- metamorphic invariance of permute() at fixed rseed:
*              ROW ORDER and LABEL transformations.
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

*---------------------------------------------------------------- base
mkpanel
save panelA.dta, replace
cap noi wooldid y unit t ttre , `OPT' permsave(ps_base)
dis "RC_base=" _rc
dmp base
asgdmp base

*---------------------------------------------------------------- rerun (same session, same data)
cap noi wooldid y unit t ttre , `OPT' permsave(ps_rerun)
dis "RC_rerun=" _rc
dmp rerun
asgdmp rerun

*---------------------------------------------------------------- row shuffle (deterministic, no RNG)
use panelA.dta, clear
qui gen double _rr = mod(_n*2654435761, 1000003)
sort _rr
qui drop _rr
cap noi wooldid y unit t ttre , `OPT' permsave(ps_rowsh)
dis "RC_rowsh=" _rc
dmp rowsh
asgdmp rowsh

*---------------------------------------------------------------- reverse row order
use panelA.dta, clear
qui gen long _n0 = -_n
sort _n0
qui drop _n0
cap noi wooldid y unit t ttre , `OPT' permsave(ps_revrow)
dis "RC_revrow=" _rc
dmp revrow
asgdmp revrow

*---------------------------------------------------------------- option list reordered
use panelA.dta, clear
cap noi wooldid y unit t ttre , replace dots(0) espostlength(3) rseed(20260917) esprelength(2) permreps(60) permute(shuffle) cluster(state) att_i att permsave(ps_optord)
dis "RC_optord=" _rc
dmp optord
asgdmp optord

*---------------------------------------------------------------- i relabelled, order preserving
use panelA.dta, clear
qui replace unit = unit*7 + 3
cap noi wooldid y unit t ttre , `OPT' permsave(ps_irelU)
dis "RC_irelU=" _rc
dmp irelU
asgdmp irelU

*---------------------------------------------------------------- i relabelled, order REVERSING
use panelA.dta, clear
qui replace unit = 1000 - unit
cap noi wooldid y unit t ttre , `OPT' permsave(ps_irelD)
dis "RC_irelD=" _rc
dmp irelD
asgdmp irelD

*---------------------------------------------------------------- cluster relabelled, order preserving
use panelA.dta, clear
qui replace state = state*13 + 5
cap noi wooldid y unit t ttre , `OPT' permsave(ps_crelU)
dis "RC_crelU=" _rc
dmp crelU
asgdmp crelU

*---------------------------------------------------------------- cluster relabelled, order REVERSING
use panelA.dta, clear
qui replace state = 100 - state
cap noi wooldid y unit t ttre , `OPT' permsave(ps_crelD)
dis "RC_crelD=" _rc
dmp crelD
asgdmp crelD

dis "PROBEA_DONE"
