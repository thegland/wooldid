*=====================================================================
* s2_falsepos_default.do  -- SELF-CONTAINED.  ~3 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s2_falsepos_default.do
*
* Same claim as s1, at Stata's UNTOUCHED DEFAULT c(maxvar) = 5000, and with the
* two controls s1 did not reach.
*
* SHOWS: wooldid's variable-count preflight (ado:1371) is NOT a lower bound on
* the default estimation path.  It adds si8ilev + si8tlev, which reghdfe ABSORBS
* (ado:1476: absorb(i.`i' i.`t' `fe')) and which are therefore never design-matrix
* columns.  A 5200-unit administrative panel with 20 treated units is REFUSED with
* "requires at least 5542 variables" although the identical command estimates in
* about two seconds at that same maxvar.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid
dis as text "UNTOUCHED default c(maxvar) = " c(maxvar)

capture program drop mkpanel
program define mkpanel
  args nunits ntreat
  clear
  qui set obs `nunits'
  qui gen long unit = _n
  qui gen int ttre = .
  qui replace ttre = 2012 + mod(_n-1,4) if _n <= `ntreat'
  qui gen double ufe = rnormal()
  qui expand 12
  sort unit
  qui by unit: gen int year = 2007 + _n
  qui gen double yfe = 0.2*sin(year/2)
  qui gen byte post = (year >= ttre) & ttre < .
  qui gen double y = ufe + yfe + 1.5*post + rnormal()*0.5
end

*======== PART A: default maxvar 5000, 5200 units ====================
mkpanel 5200 20
dis as result "PANEL: 5200 units, rows=" _N ", dataset vars=" c(k)
tempfile BIG
qui save `BIG'

dis as text "=== A1: 5200 units, maxvar " c(maxvar) " (DEFAULT), DEFAULT safety ==="
cap noi wooldid y unit year ttre , att
dis as result "A1_RC=" _rc

dis as text "=== A2: same data, same DEFAULT maxvar, safety(off) ==="
timer clear 1
timer on 1
cap noi wooldid y unit year ttre , att safety(off)
local rcA2 = _rc
timer off 1
dis as result "A2_RC=" `rcA2'
if `rcA2' == 0 {
  mat A = e(wdidmainresults)
  dis as result "A2 att b = " A[1,1] "  se = " A[1,2] "  N = " e(N)
  scalar bA = A[1,1]
}
timer list 1

*======== PART B: raise maxvar, DEFAULT safety -> guard silent, same number ====
clear
clear matrix
estimates clear
set maxvar 8000
use "`BIG'", clear
dis as text "=== B1: same data/spec, DEFAULT safety, maxvar " c(maxvar) " ==="
cap noi wooldid y unit year ttre , att
local rcB = _rc
dis as result "B1_RC=" `rcB'
if `rcB' == 0 {
  mat B = e(wdidmainresults)
  dis as result "B1 att b = " B[1,1]
  cap dis as result "B1 == A2 bit-for-bit? " (B[1,1]==bA)
  cap dis as result "B1 att in hex " %21x B[1,1]
}

*======== PART C: how many columns did the fit ACTUALLY use? ==========
* Rebuild wooldid's own design by hand on the same data and count the columns
* reghdfe actually estimates.  This is the number si8demand claims to bound.
clear
clear matrix
estimates clear
set maxvar 8000
use "`BIG'", clear
set emptycells drop
gen long _i = cond(ttre<., unit, 0)
gen int  rel = year - ttre
egen long _t = group(rel) if ttre<.
replace _t = 0 if missing(_t)
gen byte tf = (year>=ttre) & ttre<.
cap noi reghdfe y i._i#i._t#i.tf , absorb(i.unit i.year) vce(robust)
dis as result "C_RC=" _rc
cap dis as result "ACTUAL columns in the wooldid-shaped design (colsof e(b)) = " colsof(e(b))

*======== PART D: control -- 500 units, same 20 treated, maxvar 2048 ===
clear
clear matrix
estimates clear
set maxvar 2048
mkpanel 500 20
dis as text "=== D1: 500 units, maxvar " c(maxvar) ", DEFAULT safety ==="
cap noi wooldid y unit year ttre , att
dis as result "D1_RC=" _rc

dis "S2_DONE"
