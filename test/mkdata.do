* Build deterministic test datasets for the wooldid regression battery
clear all
set more off
local OUT "data"
cap mkdir "`OUT'"

* ---- binary-treatment staggered panel (medium)
clear
set seed 20260719
set sortseed 20260719
set obs 80
gen i = _n
gen ttre = 2005 if i <= 25
replace ttre = 2008 if i > 25 & i <= 50
gen sg = 1 + mod(i,2)
gen ife = rnormal()
expand 12
bysort i: gen t = 1999 + _n
gen w = 0.25 + rchi2(3)
gen eff = cond(ttre==2005, 1.0, 0.6)
gen y = 2 + ife + 0.05*(t-2000) + eff*(t >= ttre & ttre != .) + rnormal(0, 0.5)
gen ycount = round(exp(0.2*y))
save "`OUT'/dta_binary", replace

* ---- tiny degenerate panel
clear
set seed 111
set sortseed 111
set obs 8
gen i = _n
gen ttre = 4 if i <= 3
gen sg = 1 + mod(i,2)
gen ife = rnormal()
expand 6
bysort i: gen t = _n
gen w = 1 + runiform()
gen y = ife + 0.1*t + 1*(t >= ttre & ttre != .) + rnormal(0,0.3)
save "`OUT'/dta_tiny", replace

* ---- simple-example binary data (repo example, verbatim construction)
clear
set seed 1234
set sortseed 1234
set obs 100
gen i = mod(_n,10) + 1
bysort i: gen t = _n
gen ttre = 3 if i <= 3
replace ttre = 7 if i >= 4 & i <= 6
gen y = rnormal(0,.11)
replace y = y + ttre if t >= ttre
save "`OUT'/dta_simple1", replace

* ---- continuous-treatment data (repo example construction, seeded)
clear
set seed 5678
set sortseed 5678
set obs 10
gen i = _n
expand 10
bysort i: gen t = _n
gen ttre = 5 if i <= 5
gen treated = t >= ttre
expand 10
gen intensity = rnormal(10,.5)
gen y = rnormal(0,.01) + intensity*2 + 3*i + t*-2 + treated*intensity
gen intensity_uz = intensity
replace intensity_uz = 0 if ttre == .
gen w = 0.5 + runiform()
save "`OUT'/dta_cont", replace

di "MKDATA_DONE"
exit, clear
