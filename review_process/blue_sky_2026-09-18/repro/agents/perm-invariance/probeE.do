*==========================================================================
* probeE.do -- (a) permonsetcohorts (shipped 2026-09-16, one day old):
*                  label/order invariance, and the documented "no-op when
*                  cluster() == i" claim.
*              (b) permstrata under permexact on an INTERLEAVED design
*                  (ascending cluster order crosses stratum order, the
*                  NFPERM-02 geometry): p must be invariant to relabeling
*                  the strata.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

*--- 40 units in 4 onset cohorts (i == the onset grouping); cluster = unit
prog def mkonc
  clear
  set seed 5150
  set sortseed 5150
  qui set obs 40
  qui gen long unit   = _n
  qui gen long cohort = ceil(unit/10)
  qui gen double ai   = rnormal()
  qui gen double ttre = .
  qui replace ttre = 2004 if cohort == 1
  qui replace ttre = 2006 if cohort == 2
  qui replace ttre = 2008 if cohort == 3
  qui expand 10
  qui bysort unit: gen int t = 1999 + _n
  qui gen double eps  = rnormal()*0.5
  qui gen byte   post = (t >= ttre) & !missing(ttre)
  qui gen double y    = ai + 0.2*(t-2000) + 0.8*post + eps
  qui drop eps
  qui order unit cohort t ttre y
end

mkonc
save panelE.dta, replace

local OE att att_i permute(shuffle) permreps(30) rseed(20260917) permonsetcohorts dots(0) replace

*--- E0 baseline: i = cohort, cluster = unit (FINER than i)
cap noi wooldid y cohort t ttre , `OE' cluster(unit) permsave(ps_E0)
dis "RC_E0=" _rc
dis "IREBUILT_E0=" e(perm_irebuilt)
cap dmp E0
cap asgdmp E0

*--- E1 i relabelled, order REVERSING
use panelE.dta, clear
qui replace cohort = 20 - cohort
cap noi wooldid y cohort t ttre , `OE' cluster(unit) permsave(ps_E1irelD)
dis "RC_E1irelD=" _rc
dis "IREBUILT_E1=" e(perm_irebuilt)
cap dmp E1irelD
cap asgdmp E1irelD

*--- E2 cluster relabelled, order preserving
use panelE.dta, clear
qui replace unit = unit*3 + 2
cap noi wooldid y cohort t ttre , `OE' cluster(unit) permsave(ps_E2crelU)
dis "RC_E2crelU=" _rc
cap dmp E2crelU
cap asgdmp E2crelU

*--- E3 row shuffle
use panelE.dta, clear
qui gen double _rr = mod(_n*2654435761, 1000003)
sort _rr
qui drop _rr
cap noi wooldid y cohort t ttre , `OE' cluster(unit) permsave(ps_E3rowsh)
dis "RC_E3rowsh=" _rc
cap dmp E3rowsh
cap asgdmp E3rowsh

*--- E4 / E5 the documented no-op: cluster() == i, with and without the option
use panelE.dta, clear
cap noi wooldid y cohort t ttre , att att_i permute(shuffle) permreps(30) rseed(20260917) permonsetcohorts cluster(cohort) dots(0) replace permsave(ps_E4onc)
dis "RC_E4onc=" _rc
dis "IREBUILT_E4=" e(perm_irebuilt)
cap dmp E4onc
cap asgdmp E4onc

use panelE.dta, clear
cap noi wooldid y cohort t ttre , att att_i permute(shuffle) permreps(30) rseed(20260917) cluster(cohort) dots(0) replace permsave(ps_E5noonc)
dis "RC_E5noonc=" _rc
cap dmp E5noonc
cap asgdmp E5noonc

*==========================================================================
* (b) INTERLEAVED strata under permexact
*==========================================================================
prog def mkstr
  clear
  set seed 2718
  set sortseed 2718
  qui set obs 16
  qui gen long unit  = _n
  qui gen long state = ceil(unit/2)
  qui gen byte  str1 = cond(mod(state,2)==1, 1, 2)
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
  qui order unit state str1 t ttre y
end

mkstr
save panelF.dta, replace
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permexact permstrata(str1) rseed(20260917) dots(0) replace permsave(ps_F0)
dis "RC_F0=" _rc
dis "REPS_F0=" e(perm_reps)
cap dmp F0
cap asgdmp F0

*--- stratum labels REVERSED (1<->2 swapped via 9 and 8)
use panelF.dta, clear
qui replace str1 = cond(str1==1, 9, 8)
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permexact permstrata(str1) rseed(20260917) dots(0) replace permsave(ps_F1strrev)
dis "RC_F1strrev=" _rc
dis "REPS_F1=" e(perm_reps)
cap dmp F1strrev
cap asgdmp F1strrev

*--- cluster ids reversed, strata unchanged (exhaustive enumeration => same p)
use panelF.dta, clear
qui replace state = 100 - state
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permexact permstrata(str1) rseed(20260917) dots(0) replace permsave(ps_F2crelD)
dis "RC_F2crelD=" _rc
cap dmp F2crelD
cap asgdmp F2crelD

*==========================================================================
* (c) SHARPENED missing-outcome probe: the ONLY cluster carrying a given
*     onset level has no usable outcome data.  The estimation sample then
*     contains no such cohort, but the permutation pool still contains the
*     onset and hands it to data-carrying clusters in replications.
*==========================================================================
mkpanel
save panelJ.dta, replace
local CJ att att_i cluster(state) permute(shuffle) permreps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0) replace

*--- J0: state 7 (the sole 2008 adopter) has y missing on every row, rows kept
use panelJ.dta, clear
qui replace y = . if state == 7
cap noi wooldid y unit t ttre , `CJ' permsave(ps_J0)
dis "RC_J0=" _rc
cap dmp J0
cap asgdmp J0

*--- J1: the same rows physically dropped
use panelJ.dta, clear
qui drop if state == 7
cap noi wooldid y unit t ttre , `CJ' permsave(ps_J1)
dis "RC_J1=" _rc
cap dmp J1
cap asgdmp J1

*--- J2: the same cluster excluded by the user's if instead
use panelJ.dta, clear
cap noi wooldid y unit t ttre if state != 7 , `CJ' permsave(ps_J2)
dis "RC_J2=" _rc
cap dmp J2
cap asgdmp J2

dis "PROBEE_DONE"
