* ============================================================================
* repro_nocaution.do -- SELF-CONTAINED, one draw, ~1 minute.
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_nocaution.do
* The analytic path says NOTHING when a treated cohort-by-period cell holds only three
* clusters and only six treated clusters exist in total, while the bootstrap engine - which
* sthlp:532 says is "not a remedy for few treated clusters" - cautions about exactly that.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 606061
set sortseed 606061
clear
set obs 40
gen long state = _n
gen double s_g = rnormal()*0.8
gen double ttre = cond(state<=2, 5, cond(state<=4, 7, .))
expand 4
bysort state: gen long unit = (state-1)*4 + _n
gen double a_u = rnormal()*0.5
expand 12
bysort unit: gen int t = _n
tempvar fst vv
bysort state t: gen byte `fst' = (_n==1)
gen double `vv' = rnormal()*0.4 if `fst'
bysort state t (`fst'): replace `vv' = `vv'[_N]
gen double v_st = `vv'
gen double e = rnormal()*0.4
gen double rel = t - ttre
gen double tau = 0
replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
gen double y = s_g + v_st + a_u + 0.2*t + 0.5*cos(t) + tau + e
gen double coh = cond(ttre<., ttre, 99)
save "_nocau.dta", replace

di "### geometry: 40 clusters(state); i = coh; each treated (coh,period) cell holds 2 states"
di "### only 4 of the 40 states are ever treated.  TRUE att = 1.1535714285714"
tab ttre, m
egen tag = tag(state coh)
tab coh if tag & coh<99
drop tag

di ""
di "=========== 1. THE DEFAULT ANALYTIC PATH -- verbatim, nothing suppressed ==========="
use "_nocau.dta", clear
cap noi wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
di "RC1=" _rc
di ""
di ">>> Read the output above: no caution about clusters per cell, none about treated"
di ">>> clusters, and the joint pre-trend p-value is printed as an ordinary test."
di ""
di "=========== 2. THE SAME RUN UNDER bootstrap(standard) ==========="
di ">>> sthlp:532 says of this engine: 'Neither mode is a remedy for few treated clusters.'"
use "_nocau.dta", clear
cap noi wooldid y coh t ttre, att cluster(state) bootstrap(standard) reps(29) rseed(11)
di "RC2=" _rc
di ""
di "=========== 3. cluster AT i, for contrast: the one caution that does exist ==========="
use "_nocau.dta", clear
cap noi wooldid y state t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
di "RC3=" _rc
di "REPRO_NOCAUTION_DONE"
