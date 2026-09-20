* rF_mata_gtools.do -- SKEPTIC arm for state-leak-03: GTOOLS ONLY, own process. No wooldid.
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
gen double sfe = rnormal()
expand 12
bys st: gen int yr = 2004 + _n
gen double y = sfe + rnormal()*0.5

mata: mata clear
dis "=====PRE MATA DIR====="
mata: mata describe
dis "=====PRE SCALAR DIR====="
scalar dir
mata: k = 12345
mata: st_numscalar("KBEFORE", k)
dis "K_BEFORE=" KBEFORE

qui gegen double ybar = mean(y), by(st)
preserve
qui gcollapse (mean) y, by(st yr) fast
restore

dis "=====POST MATA DIR (gtools only)====="
mata: mata describe
dis "=====POST SCALAR DIR====="
scalar dir
cap mata: st_numscalar("KAFTER", k)
dis "KAFTER_RC=" _rc
dis "K_AFTER=" KAFTER
dis "GTOOLS_CLOBBERED_K=" (KAFTER != 12345)
dis "ALLDONE"
