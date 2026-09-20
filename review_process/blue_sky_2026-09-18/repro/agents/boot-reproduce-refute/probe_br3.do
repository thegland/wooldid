* =================================================================
* probe_br3.do -- SKEPTIC repro of BR-3.
* Claim under test: WITHOUT cluster(), `bootstrap(standard) att` (count leg)
* and `bootstrap(standard) att att_i` (physical re-identification leg) draw
* DIFFERENT resamples at a common rseed(), so the reported att SE moves when
* a second estimand is added. WITH cluster() they agree.
* Evidence used: e(bootdist), the replication coefficient draws.
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_br3.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150

set obs 40
gen long unit = _n
gen int ttre = .
replace ttre = 2004 if unit <= 13
replace ttre = 2006 if unit > 13 & unit <= 26
gen double ufe = rnormal()
expand 8
bysort unit: gen int year = 2000 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.5, 0))
gen double y  = 1 + ufe + 0.1*(year-2000) + te + rnormal()*0.7
save "sk_b3.dta", replace

local RS 20260917
local RP 120

* ---- UNCLUSTERED: att alone (count leg) ----
use "sk_b3.dta", clear
cap noi wooldid y unit year ttre, att bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
dis "B3_U1_RC=" _rc
matrix U1 = e(bootdist)
dis "B3_U1_SE=" %21x sqrt(el(e(V),1,1))
dis "B3_U1_COLS=" colsof(U1) " ROWS=" rowsof(U1)
local u1nm : colnames U1
dis "B3_U1_NAMES=`u1nm'"

* ---- UNCLUSTERED: att + att_i (physical re-id leg) ----
use "sk_b3.dta", clear
cap noi wooldid y unit year ttre, att att_i bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
dis "B3_U2_RC=" _rc
matrix U2 = e(bootdist)
dis "B3_U2_SE=" %21x sqrt(el(e(V),1,1))
dis "B3_U2_COLS=" colsof(U2) " ROWS=" rowsof(U2)
local u2nm : colnames U2
dis "B3_U2_NAMES=`u2nm'"

* compare the att column replication by replication
matrix U1a = U1[1...,1]
matrix U2a = U2[1...,1]
dis "B3_UNCL_MRELDIF=" mreldif(U1a, U2a)
forvalues r = 1/5 {
  dis "B3_UNCL_REP`r'  count=" %12.8f el(U1a,`r',1) "   reid=" %12.8f el(U2a,`r',1)
}

* ---- CLUSTERED: att alone ----
use "sk_b3.dta", clear
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
dis "B3_C1_RC=" _rc
matrix C1 = e(bootdist)
dis "B3_C1_SE=" %21x sqrt(el(e(V),1,1))

* ---- CLUSTERED: att + att_i ----
use "sk_b3.dta", clear
cap noi wooldid y unit year ttre, att att_i cluster(unit) bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
dis "B3_C2_RC=" _rc
matrix C2 = e(bootdist)
dis "B3_C2_SE=" %21x sqrt(el(e(V),1,1))

matrix C1a = C1[1...,1]
matrix C2a = C2[1...,1]
dis "B3_CL_MRELDIF=" mreldif(C1a, C2a)
forvalues r = 1/5 {
  dis "B3_CL_REP`r'  count=" %12.8f el(C1a,`r',1) "   reid=" %12.8f el(C2a,`r',1)
}

dis "PROBE_BR3_DONE"
