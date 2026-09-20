* p6 -- the documented fastnose contract, item by item (sthlp:286-297).
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

local gN 80
local gT 10
local gonsets "2004 2006 2008"
local gysd 0.30
include "`GEN'/gen_core.do"
qui save "`WD'/p6panel.dta", replace

dis "===== C1: bare fastnose, the e() contract ====="
use "`WD'/p6panel.dta", clear
ereturn clear
cap noi wooldid y unit t ttre , esprelength(3) espostlength(4) esfixedbaseperiod
dis "C1_RC=" _rc
dis "C1_properties=" e(properties)
dis "C1_vcvstate=" e(wdid_vcvstate)
cap confirm matrix e(V)
dis "C1_hasV_rc=" _rc

use "`WD'/p6panel.dta", clear
ereturn clear
cap noi wooldid y unit t ttre , esprelength(3) espostlength(4) esfixedbaseperiod fastnose
dis "C2_RC=" _rc
dis "C2_properties=" e(properties)
dis "C2_vcvstate=" e(wdid_vcvstate)
cap confirm matrix e(V)
dis "C2_hasV_rc=" _rc
cap confirm matrix e(b_es_att)
dis "C2_has_b_es_att_rc=" _rc
dis "C2_numpre_att=" e(numpre_att)
cap confirm matrix e(V_es_att)
dis "C2_has_V_es_att_rc=" _rc
dis "C2_df_r=" e(df_r)
dis "C2_MAINRESULTS:"
mat list e(wdidmainresults), format(%12.0g)

dis "--- consumers that should REFUSE ---"
cap noi test main_att = 0
dis "C2_test_rc=" _rc
cap noi lincom _b[main_att]
dis "C2_lincom_rc=" _rc
cap noi dis _se[main_att]
dis "C2_se_rc=" _rc
cap noi estat vce
dis "C2_estatvce_rc=" _rc
cap noi dis _b[main_att]
dis "C2_b_rc=" _rc

dis "--- est store / restore ---"
cap noi est store fnmod
dis "C2_eststore_rc=" _rc
cap noi est restore fnmod
dis "C2_estrestore_rc=" _rc
dis "C2_after_restore_properties=" e(properties)

dis "--- replay ---"
cap noi wooldid
dis "C2_replay_rc=" _rc

dis "===== C3: simultaneousbands + fastnose (accepted, caution, no bands) ====="
use "`WD'/p6panel.dta", clear
ereturn clear
cap noi wooldid y unit t ttre , esprelength(3) espostlength(4) esfixedbaseperiod fastnose simultaneousbands
dis "C3_RC=" _rc
dis "C3_esband_families=" e(esband_families)

dis "===== C4: makeplots + fastnose ====="
use "`WD'/p6panel.dta", clear
ereturn clear
cap noi wooldid y unit t ttre , esprelength(3) espostlength(4) fastnose makeplots
dis "C4_RC=" _rc

dis "===== C5: refusals ====="
foreach opt in unconditionalse impvalues jointtests adversarial "pretrendpower(0.5)" "bootstrap(standard) reps(5)" {
  use "`WD'/p6panel.dta", clear
  ereturn clear
  cap noi wooldid y unit t ttre , esprelength(3) espostlength(4) fastnose `opt'
  dis "C5_[`opt']_RC=" _rc
}

dis "===== C6: fastnose + subgroups with NO es window (sgmarginsdof stays nose) ====="
use "`WD'/p6panel.dta", clear
ereturn clear
cap noi wooldid y unit t ttre , subgroups(sg2)
dis "C6_AN_RC=" _rc
cap mat list e(pwsgs), format(%12.0g)
use "`WD'/p6panel.dta", clear
ereturn clear
cap noi wooldid y unit t ttre , subgroups(sg2) fastnose
dis "C6_FN_RC=" _rc
cap mat list e(pwsgs), format(%12.0g)

dis "DONE_P6"
