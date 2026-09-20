* mkdata_R.do -- my own ragged administrative-shaped panel for the FP1 refutation
* 7 states (clusters), i = onset cohort (2 treated cohorts + never), cluster(state)
* strictly finer than i.  Ragged observation windows and uneven caseloads.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771131
set sortseed 771131

clear
set obs 7
gen int state = _n
* cohort membership: 3 early adopters, 2 late adopters, 2 never
gen byte coh = 1
replace coh = 2 if inlist(state,4,5)
replace coh = 3 if inlist(state,6,7)
gen int onset = .
replace onset = 2005 if coh==1
replace onset = 2009 if coh==2
* ragged windows
gen int t0 = 1998
replace t0 = 2000 if inlist(state,2,6)
replace t0 = 1999 if inlist(state,3,7)
replace t0 = 2001 if state==4
gen int t1 = 2015
replace t1 = 2014 if state==2
replace t1 = 2013 if state==4
replace t1 = 2012 if state==6
* state level and trend heterogeneity (administrative panels have both)
gen double a_s = .
replace a_s = 1.20 in 1
replace a_s = 0.40 in 2
replace a_s = 2.10 in 3
replace a_s = -0.60 in 4
replace a_s = 1.70 in 5
replace a_s = 0.10 in 6
replace a_s = -1.30 in 7
gen double g_s = .
replace g_s =  0.030 in 1
replace g_s = -0.015 in 2
replace g_s =  0.055 in 3
replace g_s = -0.040 in 4
replace g_s =  0.010 in 5
replace g_s =  0.025 in 6
replace g_s = -0.005 in 7

expand t1-t0+1
bys state: gen int t = t0 + _n - 1
drop t0 t1

* a genuine COHORT-level covariate: the cohort's statutory generosity index.
* constant within coh by construction, which is exactly the class FP1 names.
gen double cohbase = .
replace cohbase = 0.80 if coh==1
replace cohbase = 1.45 if coh==2
replace cohbase = 2.10 if coh==3
* an "adoption wave" grouping (early adopters vs the rest)
gen byte wavegrp = (coh==1)

* outcome: log caseload per capita, cohort-level heterogeneity in LEVEL and TREND
gen double post = (onset<. & t>=onset)
set seed 771131
gen double eps = rnormal()*0.45
gen double y = a_s + g_s*(t-2006) + 0.35*cohbase + 0.012*cohbase*(t-2006) ///
               + 0.10*rnormal() + eps + 0.60*post
* uneven caseloads (an analytic weight a real admin panel would carry)
gen double cases = exp(3 + 0.8*a_s + 0.2*rnormal())

label var y "log caseload per capita"
compress
save "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute/panelR.dta", replace
tab coh state
dis "ROWS=" _N
