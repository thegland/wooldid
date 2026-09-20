*=====================================================================
* s1_falsepos.do  -- SELF-CONTAINED.  ~3 minutes.
*
*   cd <this directory>
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s1_falsepos.do
*
* CLAIM UNDER TEST (mine, not the prober's):
*   wooldid.ado:1371   si8demand = si8cells*(1+contreat) + si8ilev + si8tlev + 200
*   is documented (ado:1342-1348; review_process/open_defects.md "RETRACTED the same
*   day") as a LOWER BOUND whose whole justification is "zero false positives ... an
*   early refusal can only fire when the late r(907) was certain anyway."
*
*   On the DEFAULT estimation path the command is (ado:1476)
*      reghdfe y ... i._i#i._t#i.targetflag ... , absorb(i.`i' i.`t' `fe')
*   so the i and t levels are ABSORBED and are NOT design-matrix columns.
*   si8ilev + si8tlev is therefore not part of any bound -- it is additive slack
*   that grows one-for-one with the number of NEVER-TREATED units.  A panel with
*   many units and few treated cohorts trips the guard on slack alone and is
*   refused although it estimates fine at that very maxvar.
*
* DESIGN: 2200 hospitals, 12 years, 20 treated (onsets 2012-2015).  Ordinary
*   administrative panel.  i = unit.  maxvar 2048 == Stata/BE's FIXED ceiling,
*   which a BE user cannot raise at all.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid
dis as text "c(flavor)=" c(flavor) "  c(stata_version)=" c(stata_version) "  default c(maxvar)=" c(maxvar)

*---------------- data ----------------------------------------------
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

mkpanel 2200 20
dis as result "BIG PANEL: units=2200  rows=" _N "  vars in dataset=" c(k)
tempfile BIG
qui save `BIG'

mkpanel 500 20
tempfile SMALL
qui save `SMALL'

*---------------- ARM 1: 2200 units, maxvar 2048, DEFAULT safety -----
clear all
set maxvar 2048
use "`BIG'", clear
dis as text "=== ARM 1: 2200 units, maxvar " c(maxvar) ", DEFAULT safety ==="
cap noi wooldid y unit year ttre , att
dis as result "ARM1_RC=" _rc

*---------------- ARM 2: IDENTICAL data, IDENTICAL maxvar, safety(off)
dis as text "=== ARM 2: same data, same maxvar " c(maxvar) ", safety(off) ==="
timer clear 1
timer on 1
cap noi wooldid y unit year ttre , att safety(off)
local rc2 = _rc
timer off 1
dis as result "ARM2_RC=" `rc2'
if `rc2' == 0 {
  mat A = e(wdidmainresults)
  dis as result "ARM2 att b = " A[1,1] "   se = " A[1,2]
  dis as result "ARM2 e(N) = " e(N)
  scalar b2 = A[1,1]
}
timer list 1

*---------------- ARM 3: same data, default safety, maxvar 3000 ------
clear
set maxvar 3000
use "`BIG'", clear
dis as text "=== ARM 3: same spec, DEFAULT safety, maxvar " c(maxvar) " ==="
cap noi wooldid y unit year ttre , att
local rc3 = _rc
dis as result "ARM3_RC=" `rc3'
if `rc3' == 0 {
  mat B = e(wdidmainresults)
  dis as result "ARM3 att b = " B[1,1]
  cap dis as result "ARM3 matches ARM2 exactly? " (B[1,1]==b2)
}

*---------------- ARM 4: control -- 500 units, same 20 treated -------
clear
set maxvar 2048
use "`SMALL'", clear
dis as text "=== ARM 4: 500 units, maxvar " c(maxvar) ", DEFAULT safety ==="
cap noi wooldid y unit year ttre , att
dis as result "ARM4_RC=" _rc

dis "S1_DONE"
