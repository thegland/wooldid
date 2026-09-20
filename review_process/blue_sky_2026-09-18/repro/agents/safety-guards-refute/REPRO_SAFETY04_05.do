*==============================================================================
* SELF-CONTAINED REPRO -- SAFETY-04 (weight-sign screen on the RAW-if domain)
*                          SAFETY-05 (esfixedbaseperiod refusal under safety(off))
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*-------------------- SAFETY-04 --------------------
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2003 + mod(unit,3) if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t>=ttre & ttre<.)
gen double wgt = 0.5 + runiform()
gen double ymiss = y
replace ymiss = . if inlist(_n, 7, 8, 9, 101, 102)
replace wgt = -0.5 if missing(ymiss)
qui count if wgt < 0 & !missing(ymiss)
dis "NEG_WEIGHTS_ON_USABLE_ROWS=" r(N)   // 0: every negative weight is on a row markout drops

dis "--- D1 default safety: EXPECT r(411) even though no usable row has w<0"
cap noi wooldid ymiss unit t ttre [aw=wgt], att cluster(unit)
dis "D1_RC=" _rc
dis "--- D2 safety(off): runs"
cap noi wooldid ymiss unit t ttre [aw=wgt], att cluster(unit) safety(off)
dis "D2_RC=" _rc
if _rc==0 {
  mat d2=e(b)
  dis "D2_N=" e(N)
  dis "D2_ATT_HEX=" %21x d2[1,1]
}
dis "--- D3 default safety + the one-token workaround if !missing(ymiss)"
cap noi wooldid ymiss unit t ttre [aw=wgt] if !missing(ymiss), att cluster(unit)
dis "D3_RC=" _rc
if _rc==0 {
  mat d3=e(b)
  dis "D3_N=" e(N)
  dis "D3_ATT_HEX=" %21x d3[1,1]
  dis "D2_EQ_D3=" cond(d3[1,1]==d2[1,1],"YES","NO")
}

*-------------------- SAFETY-05 --------------------
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2000 if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t>=ttre & ttre<.)
save "_repro_s05.dta", replace

dis "--- E1 default safety: the named r(459)"
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(2) espostlength(2)
dis "E1_RC=" _rc
dis "--- E2 safety(off), same call: OBSERVED r(198) 'invalid syntax'"
use "_repro_s05.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(2) espostlength(2) safety(off)
dis "E2_RC=" _rc
dis "--- E3 default safety, no ES window options: the named r(459)"
use "_repro_s05.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod
dis "E3_RC=" _rc
dis "--- E4 safety(off), no ES window options: OBSERVED r(461) 'no observations in subpop()'"
use "_repro_s05.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod safety(off)
dis "E4_RC=" _rc
dis "REPRO_SAFETY04_05_DONE"
