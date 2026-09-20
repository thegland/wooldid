*==========================================================================
* p5_semantics_and_overflow.do
*  (1) is the single shared ES_att_it block actually DIFFERENT from ES_att?
*      (i.e. is the name a user gets when asking for att_i alone carrying a
*      number they would not have got from ES_att?)
*  (2) 8-digit subgroup values: clean error, or a silently TRUNCATED /
*      COLLIDING e(b) name?
*  SELF-CONTAINED.
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
end
cap prog drop dumpb
prog def dumpb
  args tag
  local nn : colnames e(b)
  foreach n of local nn {
    di "BVAL|`tag'|`n'|" %21x _b[`n']
  }
end

* ---- DGP with VERY unequal cohort sizes so att (cell-weighted) and the
*      equal-weighted family must differ ----
clear
set obs 60
gen long unit = _n
gen double u_fe = rnormal()
gen int ttre = .
replace ttre = 2003 if _n <= 30          // big cohort, 30 units
replace ttre = 2005 if inrange(_n,31,33) // tiny cohort, 3 units
replace ttre = 2007 if inrange(_n,34,36) // tiny cohort, 3 units
* _n 37..60 never treated
gen double coheff = .
replace coheff = 0.5 if ttre == 2003
replace coheff = 3.0 if ttre == 2005
replace coheff = 3.0 if ttre == 2007
replace coheff = 0   if ttre >= .
expand 10
bysort unit: gen int year = 1999 + _n
gen double post = (year >= ttre) & ttre < .
gen double y = u_fe + 0.2*(year-2000) + post*coheff + 0.20*rnormal()
save "p5_data.dta", replace

di "=== RUN R1: att + att_it + att_i + att_itime together ==="
cap noi wooldid y unit year ttre, att att_it att_i att_itime esprelength(1) espostlength(1) robust
di "R1_RC=" _rc
if _rc == 0 {
  dumpn "R1"
  dumpb "R1"
}

di "=== RUN R2: att_i ALONE (what a user who wants att_i actually gets) ==="
use "p5_data.dta", clear
cap noi wooldid y unit year ttre, att_i esprelength(1) espostlength(1) robust
di "R2_RC=" _rc
if _rc == 0 {
  dumpn "R2"
  dumpb "R2"
}

di "=== RUN R3: att ALONE (for the ES_att comparison) ==="
use "p5_data.dta", clear
cap noi wooldid y unit year ttre, att esprelength(1) espostlength(1) robust
di "R3_RC=" _rc
if _rc == 0 dumpb "R3"

* ================= 8-digit subgroup values (name overflow) ==============
di "=== RUN S: 8-digit subgroup values, att_itime requested ==="
use "p5_data.dta", clear
gen long sg8 = cond(mod(unit,2)==0, 10000001, 10000002)
cap noi wooldid y unit year ttre, att att_itime subgroups(sg8) robust
di "S_RC=" _rc
if _rc == 0 {
  dumpn "S"
  dumpb "S"
}

di "=== RUN T: 7-digit subgroup values, att_itime requested ==="
use "p5_data.dta", clear
gen long sg7 = cond(mod(unit,2)==0, 1000001, 1000002)
cap noi wooldid y unit year ttre, att att_itime subgroups(sg7) robust
di "T_RC=" _rc
if _rc == 0 {
  dumpn "T"
}

di "PROBE5_DONE"
