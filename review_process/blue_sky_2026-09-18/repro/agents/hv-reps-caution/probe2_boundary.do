*==============================================================================
* probe2_boundary.do -- (a) is the caution CORRECT at its .05 boundary when it
* is reachable?  (b) does the defect fire at a realistic reps() just under the
* documented 40?  (c) is the permute() counterpart really unstrandable?
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

clear
set obs 100
gen long unit = _n
gen double _u = runiform()
sort _u unit
gen long _rk = _n
gen double ttre = .
qui replace ttre = 2003 if _rk <= 20
qui replace ttre = 2005 if _rk > 20 & _rk <= 40
qui replace ttre = 2007 if _rk > 40 & _rk <= 60
drop _u _rk
sort unit
gen double a_i = rnormal()
expand 10
bysort unit: gen int t = 1999 + _n
gen byte post = (ttre < . & t >= ttre)
gen double y0 = a_i + 0.10*(t - 2000) + rnormal()*0.5
compress

dis "WDLEG_B1_R39_CLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(39) rseed(7) dots(0) cluster(unit)
dis "WDRCB1=" _rc

dis "WDLEG_B2_R40_CLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(40) rseed(7) dots(0) cluster(unit)
dis "WDRCB2=" _rc

dis "WDLEG_B3_R39_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(39) rseed(7) dots(0)
dis "WDRCB3=" _rc

dis "WDLEG_B4_PERMUTE_NOCLUS"
cap noi wooldid y0 unit t ttre, permute(shuffle) permreps(20) rseed(7) dots(0)
dis "WDRCB4=" _rc

dis "WDLEG_B5_PERMUTE_CLUS_R20"
cap noi wooldid y0 unit t ttre, permute(shuffle) permreps(20) rseed(7) dots(0) cluster(unit)
dis "WDRCB5=" _rc

dis "WDLEG_B_ALLDONE"
