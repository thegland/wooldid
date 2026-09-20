* rD_mata_wooldid.do -- SKEPTIC arm for state-leak-03: WOOLDID ONLY, own process.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

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

* start from a provably empty mata workspace, and plant a user variable named k
mata: mata clear
dis "=====PRE MATA DIR====="
mata: mata describe
dis "=====PRE SCALAR DIR====="
scalar dir
mata: k = 12345
mata: kk_check = 1
dis "PLANTED"
mata: st_numscalar("KBEFORE", k)
dis "K_BEFORE=" KBEFORE

qui wooldid y st yr onset, cluster(st) att

dis "=====POST MATA DIR (wooldid only)====="
mata: mata describe
dis "=====POST SCALAR DIR====="
scalar dir
cap mata: st_numscalar("KAFTER", k)
dis "KAFTER_RC=" _rc
dis "K_AFTER=" KAFTER
dis "WOOLDID_CLOBBERED_K=" (KAFTER != 12345)
dis "ALLDONE"
