* p2_size_degen.do
*  A) degenerate: ONE live cluster padded to 6 by dead ones -> what SE comes back?
*  B) bootalpha(size) with PARTIAL within-cluster missingness -> SE shrunk by sqrt(f)?
*  C) bootsave(): do the w_ columns include dead clusters, and does the _rep=0 row
*     hold 1/N of the ESTIMATION SAMPLE as the help promises?
*  D) complete-data triggers: a singleton unit, and rows outside the ES window.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootsample"

*================================================== A: degenerate 1 live cluster
* 6 states. State 6 alone is live and holds BOTH treated and never-treated units,
* so the run is estimable with a single cluster. States 1-5 report no outcome.
clear
set obs 6
gen int state = _n
gen double sfx = rnormal()*0.4
expand 4
bysort state: gen int unit = (state-1)*4 + _n
gen double ufx = rnormal()*0.3
gen int onset = .
replace onset = 2004 if mod(unit,2)==1
expand 10
bysort unit: gen int year = 1999 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .10*(year-1999) + sfx + ufx + .5*post + rnormal()*0.6
tempfile a0
save `a0', replace
replace y = . if state <= 5
qui count if !missing(y)
dis "A_LIVEROWS=" r(N)
dis "=== A1: 6 clusters in the frame, 1 in the estimation sample ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(100) bootdf(clusters) rseed(31) dots(0)
dis "A1_RC=" _rc
dis "A1_N=" e(N)
dis "A1_BOOTDF=" e(boot_df)
cap mat va = e(V)
cap mat ba = e(b)
cap dis "A1_B=" ba[1,1]
cap dis "A1_SE=" sqrt(va[1,1])
cap dis "A1_SE_HEX=" %21x sqrt(va[1,1])
cap mat li e(boot_ci_percentile)

dis "=== A2: same data, dead rows dropped (oracle) ==="
use `a0', clear
keep if state == 6
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(100) bootdf(clusters) rseed(31) dots(0)
dis "A2_RC=" _rc

*================================================== B: bootalpha(size), f = 1/2
* 14 states x 2 units x 12 years. y missing for a random half of the cells
* (an outcome from a survey with partial coverage). Every cluster keeps ~half
* its rows, so NO cluster is dead: the only thing that moves is the Dirichlet
* concentration n_c, which the engine reads off the UNSCREENED frame.
clear
set seed 8080
set obs 14
gen int state = _n
gen double sfx = rnormal()*0.4
gen int onset = .
replace onset = 2003 + mod(_n-1,4) if _n <= 8
expand 2
bysort state: gen int unit = (state-1)*2 + _n
gen double ufx = rnormal()*0.3
expand 12
bysort unit: gen int year = 1997 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .08*(year-1997) + sfx + ufx + .5*post + rnormal()*0.6
gen double u = runiform()
tempfile b0
save `b0', replace
replace y = . if u < 0.5
qui count
local btot = r(N)
qui count if !missing(y)
dis "B_ROWS_FULL=" `btot'
dis "B_ROWS_LIVE=" r(N)
dis "B_F=" r(N)/`btot'

dis "=== B1: FULL DATA, bootalpha(size) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) bootalpha(size) reps(400) bootdf(normal) rseed(555) dots(0) bootsave(b1save) replace
dis "B1_RC=" _rc
dis "B1_N=" e(N)
mat vb1 = e(V)
mat bb1 = e(b)
dis "B1_B=" bb1[1,1]
dis "B1_SE=" sqrt(vb1[1,1])

dis "=== B2: DEAD ROWS DROPPED, bootalpha(size) (oracle) ==="
use `b0', clear
drop if u < 0.5
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) bootalpha(size) reps(400) bootdf(normal) rseed(555) dots(0) bootsave(b2save) replace
dis "B2_RC=" _rc
dis "B2_N=" e(N)
mat vb2 = e(V)
mat bb2 = e(b)
dis "B2_B=" bb2[1,1]
dis "B2_SE=" sqrt(vb2[1,1])
dis "B_SE_RATIO_FULL_OVER_ORACLE=" sqrt(vb1[1,1])/sqrt(vb2[1,1])

dis "=== B3: same comparison under the DEFAULT bootalpha(equal) ==="
use `b0', clear
replace y = . if u < 0.5
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(400) bootdf(normal) rseed(556) dots(0)
mat ve1 = e(V)
dis "B3A_SE=" sqrt(ve1[1,1])
use `b0', clear
drop if u < 0.5
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(400) bootdf(normal) rseed(556) dots(0)
mat ve2 = e(V)
dis "B3B_SE=" sqrt(ve2[1,1])
dis "B3_SE_RATIO=" sqrt(ve1[1,1])/sqrt(ve2[1,1])

*================================================== C: bootsave identities
use b1save.dta, clear
qui ds w_*
dis "C_FULL_NWCOLS=" `: word count `r(varlist)''
dis "C_FULL_WCOLNAMES=" "`r(varlist)'"
qui su _rep
list _rep in 1/1
* the help: "_rep = 0 row holds ... 1/N, N the number of observations in the
* estimation sample", and "it is the rows, sum over c of n_c times w_c, that sum to 1"
qui sum _rep
foreach v of varlist w_* {
  qui sum `v' if _rep == 0
  local c_obs`v' = r(mean)
}
qui sum w_1 if _rep == 0
dis "C_FULL_REP0_W=" r(mean)
dis "C_FULL_REP0_INV=" 1/r(mean)
use b2save.dta, clear
qui sum w_1 if _rep == 0
dis "C_ORACLE_REP0_W=" r(mean)
dis "C_ORACLE_REP0_INV=" 1/r(mean)

*================================================== D: complete-data triggers
* D1 singleton unit: every value present, but the unit is observed in ONE period.
clear
set seed 606
set obs 10
gen int state = _n
gen double sfx = rnormal()*0.4
gen int onset = .
replace onset = 2003 + mod(_n-1,3) if _n <= 6
expand 2
bysort state: gen int unit = (state-1)*2 + _n
gen double ufx = rnormal()*0.3
expand 8
bysort unit: gen int year = 2000 + _n
gen byte post = (onset<. & year>=onset)
gen double y = 1 + .08*(year-2000) + sfx + ufx + .5*post + rnormal()*0.6
* state 11: one unit, one period, complete data
expand 2 if _n == 1
replace state = 11 if _n == _N
replace unit  = 99 if _n == _N
replace onset = .  if _n == _N
replace year  = 2004 if _n == _N
replace y     = 2.0 if _n == _N
tempfile d0
save `d0', replace
dis "=== D1: singleton cluster, COMPLETE data ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(60) bootdf(clusters) rseed(9) dots(0)
dis "D1_RC=" _rc
dis "D1_N=" e(N)
dis "D1_BOOTDF=" e(boot_df)
dis "=== D2: same panel without the singleton row ==="
use `d0', clear
drop if state == 11
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(60) bootdf(clusters) rseed(9) dots(0)
dis "D2_RC=" _rc
dis "D2_N=" e(N)
dis "D2_BOOTDF=" e(boot_df)

dis "P2_ALLDONE"
