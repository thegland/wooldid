* R6_perm_oracle.do -- SELF-CONTAINED.  My OWN probes of the permute() surface.
* (1) Does level() actually reach the permutation shift interval, and is the interval
*     the documented "estimate minus the equal-tailed quantiles of the permutation
*     distribution"?  Verified by hand against e(permdist_b).
* (2) nofastcollapse under permute() claims "no answer changes".  Same rseed, both ways.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150

set obs 44
gen int u = _n
gen double a_i = rnormal(0, 1.2)
gen byte g = mod(_n, 4)
gen int onset = .
replace onset = 2003 if g == 1
replace onset = 2005 if g == 2
replace onset = 2007 if g == 3
expand 10
bysort u: gen int year = 2000 + _n
gen double e2 = rnormal(0, 1.0)
bysort u (year): replace e2 = .45*e2[_n-1] + e2 if _n > 1
gen int rel = year - onset
gen double eff = 0
replace eff = 0.5 + 0.2*rel if onset < . & rel >= 0
gen double y = 5 + a_i + .1*(year-2000) + eff + e2
drop rel eff

*=========== (1a) default level 95
ereturn clear
cap noi wooldid y u year onset, att cluster(u) permute(shuffle) permreps(120) rseed(31337)
dis "L95_RC=" _rc
matrix W95 = e(wdidmainresults)
dis "L95_colnames=" "`: colnames W95'"
dis "L95_permlevel=" e(perm_level)
matrix D95 = e(permdist_b)
dis "L95_draws_rows=" rowsof(D95) " cols=" colsof(D95)
dis "L95_draw_colnames=" "`: colnames D95'"
matrix list W95, format(%21x)

*=========== (1b) level 90 on the same seed
ereturn clear
cap noi wooldid y u year onset, att cluster(u) permute(shuffle) permreps(120) rseed(31337) level(90)
dis "L90_RC=" _rc
matrix W90 = e(wdidmainresults)
dis "L90_colnames=" "`: colnames W90'"
dis "L90_permlevel=" e(perm_level)
matrix D90 = e(permdist_b)
matrix list W90, format(%21x)

*=========== hand oracle: shift interval from the draws
* rebuild lb/ub as  b - q_{hi}(draws)  and  b - q_{lo}(draws)
foreach L in 95 90 {
  local alpha = (100-`L')/200
  preserve
  clear
  qui svmat double D`L', names(col)
  local nm : colfullnames D`L'
  local first : word 1 of `nm'
  qui ds
  local v1 : word 1 of `r(varlist)'
  qui count if !missing(`v1')
  local n = r(N)
  _pctile `v1', p(`=100*`alpha'' `=100*(1-`alpha')')
  local qlo = r(r1)
  local qhi = r(r2)
  dis "ORACLE`L' n=`n' qlo=`qlo' qhi=`qhi'"
  restore
  matrix WW = W`L'
  dis "ORACLE`L' b=" WW[1,1] " reported_lb=" WW[1,5] " reported_ub=" WW[1,6]
  dis "ORACLE`L' hand_lb=" WW[1,1]-`qhi' " hand_ub=" WW[1,1]-`qlo'
}

*=========== (2) nofastcollapse invariance under permute()
ereturn clear
cap noi wooldid y u year onset, att cluster(u) permute(shuffle) permreps(60) rseed(2024)
dis "FC0_RC=" _rc
matrix F0 = e(wdidmainresults)
matrix list F0, format(%21x)
ereturn clear
cap noi wooldid y u year onset, att cluster(u) permute(shuffle) permreps(60) rseed(2024) nofastcollapse
dis "FC1_RC=" _rc
matrix F1 = e(wdidmainresults)
matrix list F1, format(%21x)
matrix DIFF = F1 - F0
matrix list DIFF, format(%21x)
dis "---- R6 done ----"
