* R9_perm_ifpop.do -- SELF-CONTAINED.  Does the permute() randomization design use the
* SAME population as the estimation sample when the user supplies an if-restriction?
* Oracle: with G ever/never-treated clusters split k per onset group, the number of
* distinct shuffle assignments is the multinomial coefficient G!/(prod k_g!).  wooldid
* reports it in e(perm_npossible_str), so the count is a direct read-out of which
* clusters the design believes it has.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242

set obs 60
gen int firm = _n
gen double a_i = rnormal(0, 1.1)
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
* region: firms 1-40 in region 1, 41-60 in region 2, EXACTLY 10 per onset group in region 1
gen byte region = cond(_n <= 40, 1, 2)
* a bare 0/1/2 flag for the BS-37 "==1 only" rule
gen byte flag = cond(_n <= 40, 1, 2)
expand 10
bysort firm: gen int year = 1997 + _n
gen double e2 = rnormal(0, 1.0)
bysort firm (year): replace e2 = .4*e2[_n-1] + e2 if _n > 1
gen int rel = year - onset
gen double eff = 0
replace eff = 0.5 + 0.2*rel if onset < . & rel >= 0
gen double y = 4 + a_i + .1*(year-1997) + eff + e2
drop rel eff

* how many firms per onset group inside region 1?
preserve
bysort firm: keep if _n == 1
tab grp region
restore

* ORACLES computed in mata
mata:
  lg = lnfactorial(60) - 4*lnfactorial(15)
  printf("ORACLE_full = %g\n", exp(lg))
  lg2 = lnfactorial(40) - 4*lnfactorial(10)
  printf("ORACLE_region1 = %g\n", exp(lg2))
end

*===== (1) unrestricted
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(40) rseed(11)
dis "U_RC=" _rc
dis "U_nclusters=" e(nclusters)
dis "U_N=" e(N)
dis "U_npossible=" e(perm_npossible_str)

*===== (2) if region == 1
ereturn clear
cap noi wooldid y firm year onset if region == 1, att cluster(firm) permute(shuffle) permreps(40) rseed(11)
dis "R_RC=" _rc
dis "R_nclusters=" e(nclusters)
dis "R_N=" e(N)
dis "R_npossible=" e(perm_npossible_str)

*===== (3) bare 0/1/2 flag -- BS-37 says only rows where the expression is EXACTLY 1 are kept
ereturn clear
cap noi wooldid y firm year onset if flag, att cluster(firm) permute(shuffle) permreps(40) rseed(11)
dis "F_RC=" _rc
dis "F_nclusters=" e(nclusters)
dis "F_N=" e(N)
dis "F_npossible=" e(perm_npossible_str)
dis "F_Ndropped=" e(N_dropped)
dis "R9 done"
