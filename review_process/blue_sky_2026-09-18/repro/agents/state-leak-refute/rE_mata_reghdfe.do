* rE_mata_reghdfe.do -- SKEPTIC arm for state-leak-03: REGHDFE ONLY, own process. No wooldid.
clear all
set more off
set processors 2
set seed 20260917
set sortseed 20260917
clear
set obs 40
gen int st = _n
gen int onset = .
replace onset = 2009 if st <= 10
replace onset = 2011 if st > 10 & st <= 20
replace onset = 2013 if st > 20 & st <= 30
gen double sfe = rnormal()
expand 12
bys st: gen int yr = 2004 + _n
gen byte post = (yr >= onset) & onset < .
gen double y = sfe + 0.03*(yr-2004) + 0.40*post + rnormal()*0.5

mata: mata clear
dis "=====PRE MATA DIR====="
mata: mata describe
dis "=====PRE SCALAR DIR====="
scalar dir
mata: k = 12345
mata: st_numscalar("KBEFORE", k)
dis "K_BEFORE=" KBEFORE

qui reghdfe y post, absorb(st yr) cluster(st)

dis "=====POST MATA DIR (reghdfe only)====="
mata: mata describe
dis "=====POST SCALAR DIR====="
scalar dir
cap mata: st_numscalar("KAFTER", k)
dis "KAFTER_RC=" _rc
dis "K_AFTER=" KAFTER
dis "REGHDFE_CLOBBERED_K=" (KAFTER != 12345)
dis "ALLDONE"
