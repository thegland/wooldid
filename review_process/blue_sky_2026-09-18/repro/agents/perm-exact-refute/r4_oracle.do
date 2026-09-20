* r4_oracle.do -- INDEPENDENT BRUTE-FORCE ORACLE for the permexact p-values.
* Same 8-cluster stratified design as r3_enum.do (permexact is deterministic, so this
* rerun reproduces it exactly). Writes the replication draws AND every reported
* permutation quantity to CSV so an outside oracle can recompute them from scratch.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 8
gen long unit = _n
gen byte strat = mod(unit, 2)
gen double onset = .
replace onset = 4 if unit == 1
replace onset = 4 if unit == 3
replace onset = 6 if unit == 5
replace onset = 5 if unit == 2
gen double ue = rnormal()*0.5
expand 9
bysort unit: gen int year = _n
gen double te = rnormal()*0.3
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.15*year + ue + te + 0.8*post
capture erase r4_asg.dta
capture noisily wooldid y unit year onset, att espostlength(2) cluster(unit) ///
  permute(shuffle) permexact permstrata(strat) permstepdown permsave(r4_asg) replace dots(0) level(95)
dis "R4_RC=" _rc
dis "R4_perm_reps=" e(perm_reps)

foreach m in perm_p perm_p_lo perm_p_hi perm_ci perm_sd perm_n perm_comparator_b perm_stepdown b {
  cap confirm matrix e(`m')
  if _rc == 0 {
    preserve
    clear
    matrix MM = e(`m')
    svmat double MM, names(col)
    export delimited using "r4_`m'.csv", replace nolabel
    restore
    dis "R4_WROTE `m'"
    matrix list e(`m'), format(%21x)
  }
  else dis "R4_MISSING `m'"
}
use r4_asg.dta, clear
export delimited using "r4_asg.csv", replace nolabel
dis "DONE_R4"
