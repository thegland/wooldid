* R8_onsetcohorts.do -- SELF-CONTAINED.  My own probe of permonsetcohorts (the T1
* machinery that known_bydesign section 7 flags as one day old and unassessed).
* i = the onset grouping (a cohort id), cluster() = the finer real unit.
* Checks:
*  (1) without permonsetcohorts, cluster finer than i is refused r(459)  [documented]
*  (2) with it, the run completes and e(perm_irebuilt)=1
*  (3) the OBSERVED point estimates are unchanged by adding permute()+permonsetcohorts
*  (4) at cluster()==i the option is a documented no-op: e(perm_irebuilt)=0 AND the
*      p-values must equal the plain permute() run's on the same rseed
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606060
set sortseed 606060

set obs 60
gen int unit = _n
gen double a_i = rnormal(0, 1.0)
gen byte cohort = mod(_n, 4) + 1          /* 1 = never treated, 2/3/4 = onset groups */
gen int onset = .
replace onset = 2003 if cohort == 2
replace onset = 2005 if cohort == 3
replace onset = 2007 if cohort == 4
expand 10
bysort unit: gen int year = 2000 + _n
gen double e2 = rnormal(0, 1.0)
bysort unit (year): replace e2 = .4*e2[_n-1] + e2 if _n > 1
gen int rel = year - onset
gen double eff = 0
replace eff = 0.45 + 0.15*rel if onset < . & rel >= 0
gen double y = 4 + a_i + .08*(year-2000) + eff + e2
drop rel eff

*===== (0) the plain analytic run on i = cohort, cluster = unit
ereturn clear
cap noi wooldid y cohort year onset, att cluster(unit)
dis "P0_RC=" _rc
matrix A0 = e(wdidmainresults)
matrix list A0, format(%21x)

*===== (1) permute() with cluster finer than i, WITHOUT permonsetcohorts
ereturn clear
cap noi wooldid y cohort year onset, att cluster(unit) permute(shuffle) permreps(40) rseed(808)
dis "P1_RC=" _rc

*===== (2)(3) permute() + permonsetcohorts
ereturn clear
cap noi wooldid y cohort year onset, att cluster(unit) permute(shuffle) permreps(40) rseed(808) permonsetcohorts
dis "P2_RC=" _rc
dis "P2_irebuilt=" e(perm_irebuilt)
dis "P2_npossible=" e(perm_npossible_str)
matrix A2 = e(wdidmainresults)
matrix list A2, format(%21x)
dis "P2_estimate_matches_P0=" (A2[1,1] == A0[1,1])

*===== (4) at cluster()==i the option must be a no-op
ereturn clear
cap noi wooldid y cohort year onset, att cluster(cohort) permute(shuffle) permreps(40) rseed(909)
dis "P3_RC=" _rc
cap confirm matrix e(wdidmainresults)
if _rc == 0 {
  matrix A3 = e(wdidmainresults)
  matrix list A3, format(%21x)
}
ereturn clear
cap noi wooldid y cohort year onset, att cluster(cohort) permute(shuffle) permreps(40) rseed(909) permonsetcohorts
dis "P4_RC=" _rc
dis "P4_irebuilt=" e(perm_irebuilt)
cap confirm matrix e(wdidmainresults)
if _rc == 0 {
  matrix A4 = e(wdidmainresults)
  matrix list A4, format(%21x)
}
cap matrix DNOOP = A4 - A3
cap matrix list DNOOP, format(%21x)
dis "R8 done"
