clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
clear
set obs 48
gen int unit = _n
gen byte grp  = mod(_n,4)
gen int onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i = rnormal()
expand 10
bys unit: gen int year = 1999 + _n
tab year onset, missing
tab grp, missing
gen byte post = !missing(onset) & year >= onset
gen double tau = 0.8*post*(year - onset + 1)
gen double x1 = rnormal()
gen double y = a_i + 0.10*(year-2000) + tau + 0.4*x1 + 0.30*rnormal()
dis "--- try 1: i=unit ---"
cap noisily wooldid y unit year onset, att
dis "D_RC1=" _rc
dis "--- try 2: i=grp ---"
cap noisily wooldid y grp year onset, att
dis "D_RC2=" _rc
dis "--- try 3: i=unit verbose ---"
cap noisily wooldid y unit year onset, att verbose
dis "D_RC3=" _rc
dis "DIAG DONE"
