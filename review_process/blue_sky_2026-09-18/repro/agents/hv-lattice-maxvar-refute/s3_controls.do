*=====================================================================
* s3_controls.do  -- SELF-CONTAINED.  ~3 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s3_controls.do
*
* The controls for the false-positive result in s1/s2.  Same 5200-unit panel,
* same seed, so the data is identical to s2's.
*   B: maxvar 8000, DEFAULT safety -> guard silent, run completes, and the att
*      must equal s2's safety(off) att BIT FOR BIT (1.3745078...).  If it does,
*      safety(off) changed nothing about the estimate and the refusal in s2 was
*      purely a variable-count decision.
*   C: count the columns the wooldid-shaped design ACTUALLY uses.
*   D: 500-unit control -- only si8ilev changes, guard must go silent at 2048.
*=====================================================================
clear all
set more off
set processors 2
set maxvar 8000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

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

mkpanel 5200 20
tempfile BIG
qui save `BIG'

*---- B: DEFAULT safety at maxvar 8000 --------------------------------
dis as text "=== B1: 5200 units, maxvar " c(maxvar) ", DEFAULT safety ==="
cap noi wooldid y unit year ttre , att
dis as result "B1_RC=" _rc
if _rc == 0 {
  mat B = e(wdidmainresults)
  dis as result "B1 att b = " B[1,1]
  dis as result "B1 att hex = " %21x B[1,1]
  dis as result "s2 A2 att hex was 1.3745078 -> compare the hex printed above"
}

*---- C: how many columns does the design ACTUALLY have? ---------------
use "`BIG'", clear
set emptycells drop
qui gen long _i = cond(ttre<., unit, 0)
qui gen int  rel = year - ttre
qui egen long _t = group(rel) if ttre<.
qui replace _t = 0 if missing(_t)
qui gen byte tf = (year>=ttre) & ttre<.
dis as text "=== C: wooldid-shaped design, hand-built ==="
cap noi reghdfe y i._i#i._t#i.tf , absorb(i.unit i.year) vce(robust)
dis as result "C_RC=" _rc
cap dis as result "ACTUAL columns (colsof e(b)) = " colsof(e(b))
cap dis as result "ACTUAL rank (e(rank)) = " e(rank)
dis as result "si8demand claimed 'at least 5542 variables' for this specification."

*---- D: 500-unit control ---------------------------------------------
* (maxvar stays 8000 here; the point is the PREDICTED number, printed below.)
mkpanel 500 20
qui gen int rel = year - ttre
tempvar tg
qui egen byte `tg' = tag(unit rel) if (year>=ttre) & ttre<.
qui count if `tg'==1
local cells = r(N)
dis as result "D: 500-unit panel -> si8cells would be about " `cells' ///
              ", si8ilev=500, si8tlev=12, so si8demand ~ " `cells'+500+12+200
dis as result "   (2200-unit panel gave 2542; 5200-unit gave 5542; the ONLY term that moved is si8ilev)"
dis as text "=== D1: 500 units, DEFAULT safety, maxvar 8000 ==="
cap noi wooldid y unit year ttre , att
dis as result "D1_RC=" _rc

dis "S3_DONE"
