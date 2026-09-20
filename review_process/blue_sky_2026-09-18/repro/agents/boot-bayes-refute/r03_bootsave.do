* r03: bb-02.  sthlp:557 says the bootsave() _rep=0 row under bootalpha(size) holds
* 1/N with "N the number of observations in the estimation sample".  Code (ado:437)
* divides by wdb_r_nrows = the MAP-FRAME row count (user-if x nonmissing cluster),
* set at ado:8949.  Build a panel with item nonresponse on y and see which it is.
* ALSO: check the bootalpha(equal) row (1/G) and check that the replicate multipliers
* really do sum to 1 over the MAP-FRAME rows rather than over e(sample) rows.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

set seed 515151
set sortseed 515151
clear
set obs 30
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2006 if st <=  8
replace onset = 2009 if st >  8 & st <= 16
replace onset = 2012 if st > 16 & st <= 22
expand 16
bysort st: gen int year = 2000 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.04*(year-2000) + 0.8*post + rnormal()*0.6
* ordinary item nonresponse on the outcome, 15%
gen double u = runiform()
replace y = . if u < 0.15
qui count
local nrows = r(N)
qui count if !missing(y)
local nyok = r(N)
dis "ROWS_IN_DATA=`nrows'   ROWS_WITH_Y=`nyok'"
tempfile dd
save `dd', replace

cap erase r03size.dta
cap noi qui wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(20) rseed(2468) dots(0) bootsave(r03size) replace
dis "SIZE_RC=" _rc
local eN = e(N)
local G  = e(N_clust)
dis "e(N)=`eN'   e(N_clust)=`G'"
global ENN = `eN'
global NR  = `nrows'
global GG  = `G'

use r03size.dta, clear
qui ds w_*
local wv `r(varlist)'
local nw : word count `wv'
dis "w_ columns=`nw'"
qui su `: word 1 of `wv'' if _rep==0, meanonly
local obsrow = r(mean)
dis "OBSROW=" %21x `obsrow' "   =" %14.10f `obsrow'
dis "1/nrows=" %21x 1/$NR "   =" %14.10f 1/$NR "   MATCH_NROWS=" (reldif(`obsrow', 1/$NR) < 1e-13)
dis "1/e(N) =" %21x 1/$ENN "   =" %14.10f 1/$ENN "   MATCH_eN   =" (reldif(`obsrow', 1/$ENN) < 1e-13)
* are all observed-row w_ columns identical?
local allsame 1
foreach v of local wv {
  qui su `v' if _rep==0, meanonly
  if reldif(r(mean), `obsrow') > 1e-13 local allsame 0
}
dis "OBSROW_FLAT_ACROSS_CLUSTERS=`allsame'"
* do the replicate multipliers sum to 1 over the map-frame rows?
* per-cluster map-frame row count is 16 for every state here (balanced), so
* sum_c n_c * w_c = 16 * sum_c w_c
qui su `: word 1 of `wv'' if _rep==1, meanonly
local w1 = r(mean)
qui gen double __tot = 0 if _rep==1
foreach v of local wv {
  qui replace __tot = __tot + `v' if _rep==1
}
qui su __tot if _rep==1, meanonly
dis "SUM_w_over_clusters(rep1)=" %18.16f r(mean) "   times 16 = " %18.16f 16*r(mean)
* the ESTIMATED rows per cluster differ (missing y), so recompute the estimated-row sum
use `dd', clear
qui gen byte ok = !missing(y)
qui gegen double nyc = total(ok), by(st)
qui gegen byte tg = tag(st)
qui keep if tg
qui sort st
mkmat st nyc, matrix(NYC)
use r03size.dta, clear
local tot2 = 0
forvalues j = 1/`nw' {
  local v : word `j' of `wv'
  qui su `v' if _rep==1, meanonly
  local tot2 = `tot2' + NYC[`j',2]*r(mean)
}
dis "SUM n_estimated_c * w_c (rep1) = " %18.16f `tot2'

*----- bootalpha(equal) leg: observed row should be 1/G
use `dd', clear
cap erase r03equal.dta
cap noi qui wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(20) rseed(2468) dots(0) bootsave(r03equal) replace
dis "EQUAL_RC=" _rc
use r03equal.dta, clear
qui ds w_*
local wv2 `r(varlist)'
qui su `: word 1 of `wv2'' if _rep==0, meanonly
dis "EQUAL OBSROW=" %14.10f r(mean) "   1/G=" %14.10f 1/$GG "   MATCH=" (reldif(r(mean),1/$GG)<1e-13)
dis "R03_DONE"
