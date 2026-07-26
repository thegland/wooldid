* Targeted regression test for the E1b df fix: the two binary whole-sample
* pre-trend joint tests previously hardcoded df(`df_r') where every other margins
* call threads `marginsdof', so under unconditionalse they silently used the
* conditional (delta) VCE. After the fix they use df(`df_r') + vce(unconditional).
*
* This is a TOLERANCE-based suite, not a battery config: the unconditionalse path
* is not bit-reproducible run-to-run (same-build noise ~2.5e-14 on this binary
* config; see review_process/agent_brief.md), so an rtol-0 baseline would flap.
* Expected constants measured 2026-07-26 on the fixed build (commit introducing
* this file); parent-build values for the same config were joint_horiz =
* +1.8a51bceb69fcaX-003 (~0.192) and joint_horiz_i = +1.c33eddf4b1bf8X-002
* (~0.441) in %21x - the fix moves both materially, as intended.
*
* Run from test/:  "<stata>" -e do test_dffix.do   -> expect DFFIX_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* T1: the bite config - unconditional VCE now threads into the binary joint tests
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) unconditionalse jointtests esprelength(3) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "T1 FAIL: rc = " _rc
  local ++fails
}
else {
  if abs(e(joint_horiz)  - .2748732074342924) > 1e-6 {
    dis as error "T1 FAIL: joint_horiz = " %18.0g e(joint_horiz) " (expected ~.27487321)"
    local ++fails
  }
  if abs(e(joint_horiz_i) - .508946717132605) > 1e-6 {
    dis as error "T1 FAIL: joint_horiz_i = " %18.0g e(joint_horiz_i) " (expected ~.50894672)"
    local ++fails
  }
}

* T2: inertness - without unconditionalse the fix must change nothing
* (deterministic path; both scalars coincide on this balanced panel because the
* att and att_it weights are identical there)
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att att_it cluster(i) jointtests esprelength(3) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "T2 FAIL: rc = " _rc
  local ++fails
}
else {
  if abs(e(joint_horiz)  - .5189716315504179) > 1e-10 {
    dis as error "T2 FAIL: joint_horiz = " %18.0g e(joint_horiz) " (expected .51897163)"
    local ++fails
  }
  if abs(e(joint_horiz_i) - .5189716315504179) > 1e-10 {
    dis as error "T2 FAIL: joint_horiz_i = " %18.0g e(joint_horiz_i) " (expected .51897163)"
    local ++fails
  }
}

* T3: the contreat twins already threaded marginsdof - still working, still present
use "`DATA'/dta_cont", clear
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) unconditionalse jointtests esprelength(2) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "T3 FAIL: rc = " _rc
  local ++fails
}
else if "`e(joint_Choriz)'" == "" {
  dis as error "T3 FAIL: e(joint_Choriz) absent"
  local ++fails
}

if `fails' == 0 dis "DFFIX_ALL_PASS"
else dis as error "DFFIX_FAILS=`fails'"
exit, clear
