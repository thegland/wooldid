* p10_stepdown.do -- SELF-CONTAINED.  Under `permexact' the step-down kernel
* (wooldid_permstep, ado:16290-16320) computes its p_raw / p_singlestep / p_stepdown as
* the PLAIN exact share  nex/n  whenever the enumeration is complete:
*     out[2, h] = (exk | memb ? nex/n : (1 + nex)/(n + 1))
* but the runtime footer and the help both describe that row unconditionally as the
* ADD-ONE share:
*   sthlp:514  "The stored kernel's p_raw row (and the footer's first column) is the
*               two-sided |t| share, add-one, in the step-down's own studentized metric"
* This probe prints the reported values at full precision so the convention is decidable:
* with M = 30 enumerated assignments and one exceedance, the plain share is 1/30 = .0333
* and the add-one share is 2/31 = .0645.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

clear
set obs 6
gen long id = _n
gen double ons = .
replace ons = 4 in 1
replace ons = 7 in 2
expand 10
bysort id: gen int t = _n
expand 4
gen double y = 0.6*id + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
compress

cap noi wooldid y id t ons, cluster(id) esprelength(2) espostlength(3) ///
    permute(shuffle) permexact permstepdown
dis "RC=" _rc
dis "M=" e(perm_npossible) "  completed=" e(perm_completed)
dis "1/M    = " 1/e(perm_npossible)
dis "2/(M+1)= " 2/(e(perm_npossible)+1)

foreach mm in perm_stepdown perm_stepdown_pre perm_stepdown_i perm_stepdown_i_pre {
  cap confirm matrix e(`mm')
  if _rc == 0 {
    dis "--- e(`mm') ---"
    mat li e(`mm'), format(%20.15f)
  }
  else dis "--- e(`mm') ABSENT ---"
}

* the same run under SAMPLED draws of the same size, where the add-one form does apply
cap noi wooldid y id t ons, cluster(id) esprelength(2) espostlength(3) ///
    permute(shuffle) permreps(30) permstepdown rseed(12345)
dis "RC2=" _rc
dis "completed=" e(perm_completed)
cap mat li e(perm_stepdown), format(%20.15f)
dis "ALLDONE"
