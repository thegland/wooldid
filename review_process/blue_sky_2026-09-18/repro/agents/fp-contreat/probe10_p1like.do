clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 121212
set sortseed 121212
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat"

* ---- a P1-shaped fixture: 12 states, 5-12 counties each, 20 years ---------
clear
set obs 12
gen int st = _n
gen int nc = 5 + mod(_n*7,8)          /* 5..12 counties */
gen double ttre = .
qui replace ttre =  8 if inlist(st,1,2,3)
qui replace ttre = 12 if inlist(st,4,5,6)
qui replace ttre = 15 if inlist(st,7,8)
* states 9-12 never treated
expand nc
bysort st: gen int cwithin = _n
gen long unit = _n
expand 20
bysort unit: gen int year = _n
* unit-level lognormal dose with 10% time noise
bysort unit (year): gen double ubase = exp(rnormal(0,0.6)) if _n==1
bysort unit (year): replace ubase = ubase[1]
gen double dose_it = ubase*(1 + 0.10*rnormal())
bysort unit (year): gen double ufe = rnormal(0,1) if _n==1
bysort unit (year): replace ufe = ufe[1]
gen double y_dose = 3 + ufe + 0.2*year + 0.5*dose_it + rnormal(0,1)
qui replace y_dose = y_dose + 0.4*dose_it + 0.3 if ttre<. & year>=ttre
qui compress
dis "N=" _N
save "`D'/p1like.dta", replace
tab st if year==1

foreach spec in "lattice(all 3) cluster(unit)" "lattice(all 3) robust" "cluster(unit)" "robust" "lattice(all 2) robust" "lattice(all 5) robust" {
  use "`D'/p1like.dta", clear
  dis "########## SPEC: `spec'"
  cap noisily wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) `spec' esprelength(2) espostlength(2)
  dis "RC=" _rc
  dis "N=" e(N) " df_r=" e(df_r) " rss=" e(rss) " rmse=" e(rmse) " r2=" e(r2)
  cap mat li e(wdidmainresults)
}
dis "DONE"
