*==========================================================================
* probeE2.do -- (b) permstrata under permexact on an INTERLEAVED design
*                   (ascending cluster order crosses stratum order, the
*                   NFPERM-02 geometry): p must be invariant to relabeling
*                   the strata and the clusters.
*               (c) SHARPENED missing-outcome probe: the ONLY cluster
*                   carrying a given onset level has no usable outcome.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

prog def mkstr
  clear
  set seed 2718
  set sortseed 2718
  qui set obs 16
  qui gen long unit  = _n
  qui gen long state = ceil(unit/2)
  qui gen byte strat = cond(mod(state,2)==1, 1, 2)
  qui gen double ai  = rnormal()
  qui gen double ttre = .
  qui replace ttre = 2004 if state == 1 | state == 3
  qui replace ttre = 2006 if state == 2 | state == 4
  qui expand 10
  qui bysort unit: gen int t = 1999 + _n
  qui gen double eps  = rnormal()*0.5
  qui gen byte   post = (t >= ttre) & !missing(ttre)
  qui gen double y    = ai + 0.2*(t-2000) + 0.8*post + eps
  qui drop eps
end

mkstr
save panelF.dta, replace
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permexact permstrata(strat) rseed(20260917) dots(0) replace permsave(ps_F0)
dis "RC_F0=" _rc
dis "REPS_F0=" e(perm_reps)
cap dmp F0
cap asgdmp F0

use panelF.dta, clear
qui replace strat = cond(strat==1, 9, 8)
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permexact permstrata(strat) rseed(20260917) dots(0) replace permsave(ps_F1strrev)
dis "RC_F1strrev=" _rc
dis "REPS_F1=" e(perm_reps)
cap dmp F1strrev
cap asgdmp F1strrev

use panelF.dta, clear
qui replace state = 100 - state
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permexact permstrata(strat) rseed(20260917) dots(0) replace permsave(ps_F2crelD)
dis "RC_F2crelD=" _rc
cap dmp F2crelD
cap asgdmp F2crelD

*==========================================================================
mkpanel
save panelJ.dta, replace
local CJ att att_i cluster(state) permute(shuffle) permreps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0) replace

use panelJ.dta, clear
qui replace y = . if state == 7
cap noi wooldid y unit t ttre , `CJ' permsave(ps_J0)
dis "RC_J0=" _rc
cap dmp J0
cap asgdmp J0

use panelJ.dta, clear
qui drop if state == 7
cap noi wooldid y unit t ttre , `CJ' permsave(ps_J1)
dis "RC_J1=" _rc
cap dmp J1
cap asgdmp J1

use panelJ.dta, clear
cap noi wooldid y unit t ttre if state != 7 , `CJ' permsave(ps_J2)
dis "RC_J2=" _rc
cap dmp J2
cap asgdmp J2

dis "PROBEE2_DONE"
