*==============================================================================
* probe3_waldjoint.do -- does a JOINT p actually display under low-reps
* bootdisplay(wald), and is its 1/(reps+1) floor disclosed anywhere when there
* is no cluster()?  18 is the largest reps with 1/(reps+1) = .0526 > .05.
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

dis "WDLEG_C1_WALDJT_R18_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(18) rseed(7) dots(0) bootdisplay(wald) jointtests esprelength(2) espostlength(2)
dis "WDRCC1=" _rc

dis "WDLEG_C2_WALDJT_R18_CLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(18) rseed(7) dots(0) bootdisplay(wald) jointtests esprelength(2) espostlength(2) cluster(unit)
dis "WDRCC2=" _rc

dis "WDLEG_C3_WALDJT_R10_NOCLUS"
cap noi wooldid y0 unit t ttre, bootstrap(standard) reps(10) rseed(7) dots(0) bootdisplay(wald) jointtests esprelength(2) espostlength(2)
dis "WDRCC3=" _rc

dis "WDLEG_C_ALLDONE"
