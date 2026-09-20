* self-contained DGP used by every probe here: staggered panel, 2 treated cohorts + never-treated
clear all
set seed 20260917
set sortseed 20260917
set obs 60
gen long uid = _n
gen byte grp = mod(_n,3)          /* 0,1,2 */
gen ttre = .
replace ttre = 5 if grp==1
replace ttre = 8 if grp==2
gen double ufe = rnormal()
expand 12
bysort uid: gen int t = _n
gen double y = ufe + 0.10*t + rnormal()
replace y = y + 0.8*(t - ttre + 1) if ttre < . & t >= ttre
gen double w = 1 + runiform()
gen double x1 = rnormal()
label var y "outcome"
compress
save "dat.dta", replace
