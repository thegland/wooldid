clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

which wooldid
dis "BUILD_OK"

*--- realistic-ish unbalanced state-year panel -------------------------------
clear
set obs 30
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2006 if st <= 7
replace onset = 2009 if st > 7 & st <= 14
replace onset = 2012 if st > 14 & st <= 20
* st 21-30 never treated
expand 16
bysort st: gen int year = 1999 + _n
* unbalanced: unequal cluster row counts (entry years staggered)
gen double u = runiform()
bysort st: gen double u1 = u[1]
drop if year < 2000 + floor(u1*4)
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.04*(year-2000) + 1.3*post + rnormal()*0.6
gen double pop = exp(rnormal()*0.5)*100
tab st if year==2010
qui bysort st: gen long nc = _N if _n==1
qui su nc
dis "CLUSTER_ROWCOUNT min=" r(min) " max=" r(max) " N=" _N
drop nc u u1
save p01_panel, replace

*--- baseline bayes run with bootsave --------------------------------------
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(50) rseed(777) bootsave(p01_eq) replace dots(0)
dis "RUN_EQ_RC=" _rc
mat li e(b)
dis "nclus=" e(boot_nclus) " alpha=" e(boot_alpha)
ereturn list

*--- can I call the mata draw function directly? ----------------------------
cap noi mata: mata describe wooldid_bootdw()
dis "MATADESC_RC=" _rc
