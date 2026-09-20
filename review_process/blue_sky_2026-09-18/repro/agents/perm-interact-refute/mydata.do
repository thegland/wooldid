* mydata.do -- my own DGP, written from scratch for the refutation wave.
* A realistic administrative-style staggered panel: 60 firms x 12 years (1998-2009),
* three treatment cohorts plus a never-treated group, unit and year fixed effects,
* persistent unit-level noise, heterogeneous dynamic effects growing with exposure.
clear all
set more off
set seed 771103
set sortseed 771103

set obs 60
gen int firm = _n
gen double fe_i = rnormal(0, 1.4)
* cohort assignment: 15 never-treated, 15 each in 2002, 2004, 2006
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
* a firm-level covariate that is plausible as a control
gen double size = 3 + rnormal(0, .8)
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
* AR(1)-ish persistence within firm
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
gen double y = 2.0 + fe_i + fe_t + 0.4*size + effect + eps
label var y "outcome"
compress
save "mypanel.dta", replace
summ y
