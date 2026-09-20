* p7 -- (a) permute() + permanalyticse + fastnose is ACCEPTED by wooldid_validate
*         (ado:5374-5412 refuses unconditionalse/impvalues/pretrendpower, and adversarial
*          /jointtests only outside the permute fast-mode default; permanalyticse is not
*          screened, and permtstat/permtci ARE, ado:5882).  With fastnose in force no
*          analytic SE exists, yet ado:8215 picks the permanalyticse header and ado:11679
*          prints "e(V) and e(V_es_*) hold the analytic covariance".
*      (b) does e(gteffects) carry real se/p/ci under fastnose?  Its cells are read off
*          the STAGE-1 e(V) (ado:1765-1773), which fastnose never touches.
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

local gN 60
local gT 9
local gonsets "2003 2005 2007"
local gysd 0.30
include "`GEN'/gen_core.do"
qui save "`WD'/p7panel.dta", replace

dis "===== G1: gteffects, analytic ====="
use "`WD'/p7panel.dta", clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) gteffects esprelength(2) espostlength(3)
dis "G1_RC=" _rc
cap mat list e(gteffects), format(%10.0g)

dis "===== G2: gteffects, fastnose ====="
use "`WD'/p7panel.dta", clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) gteffects esprelength(2) espostlength(3) fastnose
dis "G2_RC=" _rc
cap mat list e(gteffects), format(%10.0g)

dis "===== P1: permute fast-mode default ====="
use "`WD'/p7panel.dta", clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(20) esprelength(2) espostlength(2)
dis "P1_RC=" _rc
dis "P1_vcvstate=" e(wdid_vcvstate)
dis "P1_perm_analyticse=" e(perm_analyticse)

dis "===== P2: permute + permanalyticse (no fastnose) ====="
use "`WD'/p7panel.dta", clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(20) permanalyticse esprelength(2) espostlength(2)
dis "P2_RC=" _rc
dis "P2_vcvstate=" e(wdid_vcvstate)
cap confirm matrix e(V)
dis "P2_hasV_rc=" _rc
cap noi test main_att = 0
dis "P2_test_rc=" _rc

dis "===== P3: permute + permanalyticse + fastnose  (the suspect) ====="
use "`WD'/p7panel.dta", clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(20) permanalyticse fastnose esprelength(2) espostlength(2)
dis "P3_RC=" _rc
dis "P3_vcvstate=" e(wdid_vcvstate)
dis "P3_perm_analyticse=" e(perm_analyticse)
dis "P3_vcespec=" e(wdid_vcespec)
cap confirm matrix e(V)
dis "P3_hasV_rc=" _rc
cap noi test main_att = 0
dis "P3_test_rc=" _rc
cap mat list e(wdidmainresults), format(%12.0g)
cap mat list e(permkernel), format(%12.0g)

dis "===== P4: permute + permtstat + fastnose (should be REFUSED) ====="
use "`WD'/p7panel.dta", clear
ereturn clear
cap noi wooldid y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(20) permtstat fastnose esprelength(2) espostlength(2)
dis "P4_RC=" _rc
dis "DONE_P7"
