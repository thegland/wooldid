*==========================================================================
* p7_overflow.do -- the 33-character e(b) name that ado:12283 warns about:
*   pretreatment_SG<8 digits>_att_itime.  Clean error, silent truncation, or
*   a COLLISION between two 8-digit subgroup values?  SELF-CONTAINED.
*   (p5 RUN S did NOT reach this: with no esprelength/espostlength no
*    pretreatment_* row is built, so the longest name was only 25 chars.)
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop dumpn
prog def dumpn
  args tag
  local nn : colnames e(b)
  local k = 0
  foreach n of local nn {
    local k = `k' + 1
    di "NAME|`tag'|`k'|`n'|len=" strlen("`n'")
  }
  di "NCOL|`tag'|`k'"
end

clear
set obs 60
gen long unit = _n
gen double u_fe = rnormal()
gen int ttre = .
replace ttre = 2003 if _n <= 20
replace ttre = 2006 if inrange(_n,21,36)
gen double coheff = 0
replace coheff = 1.0 if ttre < .
expand 10
bysort unit: gen int year = 1999 + _n
gen double post = (year >= ttre) & ttre < .
gen double y = u_fe + 0.2*(year-2000) + post*coheff + 0.20*rnormal()
save "p7_data.dta", replace

* two 8-digit subgroup values differing ONLY in the last digit:
*   pretreatment_SG10000001_att_itime  = 33 characters (Stata's cap is 32)
di "=== RUN W: 8-digit subgroup values, att_itime + PRETREATMENT rows ==="
gen long sg8 = cond(mod(unit,2)==0, 10000001, 10000002)
cap noi wooldid y unit year ttre, att att_itime subgroups(sg8) esprelength(1) espostlength(1) robust
di "W_RC=" _rc
if _rc == 0 dumpn "W"

di "=== RUN X: 7-digit values, same request (32-char name, at the cap) ==="
use "p7_data.dta", clear
gen long sg7 = cond(mod(unit,2)==0, 1000001, 1000002)
cap noi wooldid y unit year ttre, att att_itime subgroups(sg7) esprelength(1) espostlength(1) robust
di "X_RC=" _rc
if _rc == 0 dumpn "X"

di "PROBE7_DONE"
