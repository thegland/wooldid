* R4 -- MY OWN probes on the safety surface.
* BLOCK E: sthlp's refuse/caution/skip policy says safety(off) "has never suppressed a
*          disclosure about what was estimated - if the run changed your specification, it
*          says so whatever safety() is set to".  Under safety(off) the no-control-period
*          CAUTION is gated off but the period DROP still happens.  Is the change disclosed?
* BLOCK F: engine split for the si8 preflight AT THE DEFAULT CEILING + does a mid-run
*          r(907) under safety(off) leak wooldid tempvars into the user's dataset?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
dis "MAXVAR=" c(maxvar)

*=====================================================================
* BLOCK E -- unbalanced panel: the never-treated units exit after 2005, so
* 2006 and 2007 contain ONLY treated units and get dropped by the estimability rule.
*=====================================================================
clear
set obs 60
gen long unit = _n
gen byte tg = unit <= 20
expand 8
bysort unit: gen int t = 2000 + _n - 1
* never-treated units leave the panel after 2005 (attrition -- ordinary admin data)
drop if tg==0 & t > 2005
gen double ttre = .
replace ttre = 2003 + mod(unit,2) if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
qui count
dis "E_ROWS_IN=" r(N)
qui count if t > 2005
dis "E_ROWS_AFTER2005=" r(N)
save "r4e.dta", replace

dis "BLOCKE1_START default safety"
cap noi wooldid y unit t ttre, att cluster(unit)
dis "E1_RC=" _rc
if _rc==0 {
  mat e1 = e(b)
  dis "E1_N=" e(N)
  dis "E1_NDROP=" e(N_dropped)
  dis "E1_ATT_HEX=" %21x e1[1,1]
}
dis "BLOCKE1_END"

dis "BLOCKE2_START safety(off)"
use "r4e.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) safety(off)
dis "E2_RC=" _rc
if _rc==0 {
  mat e2 = e(b)
  dis "E2_N=" e(N)
  dis "E2_NDROP=" e(N_dropped)
  dis "E2_ATT_HEX=" %21x e2[1,1]
}
dis "BLOCKE2_END"

*=====================================================================
* BLOCK F -- 5200 units at the DEFAULT maxvar 5000: engine split + leak check
*=====================================================================
clear
set obs 5200
gen long unit = _n
gen byte tg = unit <= 20
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 + mod(unit,3) if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
qui ds
local varsbefore `r(varlist)'
dis "F_VARS_BEFORE=`varsbefore'"
qui count
dis "F_ROWS=" r(N)

dis "BLOCKF1_START default safety (expect the si8 refusal)"
cap noi wooldid y unit t ttre, att cluster(unit)
dis "F1_RC=" _rc
dis "BLOCKF1_END"

dis "BLOCKF2_START safety(off) reghdfe path"
cap noi wooldid y unit t ttre, att cluster(unit) safety(off)
dis "F2_RC=" _rc
if _rc==0 dis "F2_N=" e(N)
dis "BLOCKF2_END"

dis "BLOCKF3_START safety(off) customvce -> reg path (demand should be REAL here)"
cap noi wooldid y unit t ttre, att customvce(cluster unit) safety(off)
dis "F3_RC=" _rc
dis "BLOCKF3_END"

qui ds
local varsafter `r(varlist)'
dis "F_VARS_AFTER=`varsafter'"
dis "F_LEAKED=" cond("`varsbefore'"=="`varsafter'","NO","YES")
cap est dir
dis "F_ESTDIR_RC=" _rc
dis "R4_DONE"
