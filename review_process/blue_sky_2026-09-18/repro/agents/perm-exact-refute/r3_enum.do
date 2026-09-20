* r3_enum.do -- INDEPENDENT AUDIT of the permexact enumeration itself.
* wooldid computes the possible-assignment count from log-factorials (ado:13655-13665)
* and the enumeration from a separate recursive routine (wooldid_permx_enum, ado:13736);
* its own comment says "nothing reconciles them" and that a disagreement yields either
* all-control replications or a non-exhaustive enumeration still flagged exact.
* This probe writes the assignments with permsave() and audits them against a
* brute-force enumeration in Python.
* Design: 8 clusters, strata INTERLEAVED with cluster order (stratum = mod(unit,2)),
*   stratum 1 = units 1,3,5,7 with onsets {4,4,6,.}  -> 4!/(2!1!1!) = 12
*   stratum 0 = units 2,4,6,8 with onsets {5,.,.,.}  -> 4!/(1!3!)   =  4
*   M = 48
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
tempfile d
save `d', replace

capture erase r3_asg.dta
dis "===== stratified permexact ====="
capture noisily wooldid y unit year onset, att espostlength(2) cluster(unit) ///
  permute(shuffle) permexact permstrata(strat) permsave(r3_asg) replace dots(0)
dis "R3_RC=" _rc
dis "R3_perm_exact=" e(perm_exact)
dis "R3_perm_reps=" e(perm_reps)
cap confirm matrix e(perm_n)
if _rc==0 matrix list e(perm_n)
cap confirm matrix e(perm_p)
if _rc==0 matrix list e(perm_p), format(%20.17f)

use r3_asg.dta, clear
desc, fullnames
export delimited using "r3_asg.csv", replace nolabel
dis "R3_NASG_ROWS=" _N
dis "DONE_R3"
