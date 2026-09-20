* =================================================================
* probe_br4c.do -- BR-4 realism: do the COMMON routes out of e(sample)
* also inflate the bootstrap's resampling universe?
*   route S: singleton units (one period only; reghdfe drops singletons)
*   route C: units whose controls() variable is missing on every row
*   route Z: the same 9 estimation units with nothing extra (the truth)
* All three must give the same e(N)/e(N_clust). If e(boot_df) and the two
* few-cluster cautions differ, the bootstrap universe is not e(sample).
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_br4c.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606060
set sortseed 606060

* --- 9 core units, 10 periods, plus 6 extra units used differently per route ---
set obs 15
gen long unit = _n
gen byte extra = unit > 9
gen int ttre = .
replace ttre = 2004 if inlist(unit,1,2)
replace ttre = 2006 if inlist(unit,3,4)
replace ttre = 2004 if inlist(unit,10,11,12,13)
gen double ufe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.6, 0))
gen double y  = 1 + ufe + 0.1*(year-1999) + te + rnormal()*0.6
gen double x  = rnormal()
save "sk_b4c_raw.dta", replace

* route Z: extras absent
use "sk_b4c_raw.dta", clear
drop if extra
save "sk_b4c_Z.dta", replace

* route S: extras kept but observed in ONE period each (singletons)
use "sk_b4c_raw.dta", clear
bysort unit (year): drop if extra & _n > 1
save "sk_b4c_S.dta", replace

* route C: extras kept in full, but the control x is missing on every extra row
use "sk_b4c_raw.dta", clear
replace x = . if extra
save "sk_b4c_C.dta", replace

* route Zc: route Z with the same controls(x) call, for a like-for-like comparison
use "sk_b4c_raw.dta", clear
drop if extra
save "sk_b4c_Zc.dta", replace

foreach arm in Z S {
  use "sk_b4c_`arm'.dta", clear
  dis "########## ROUTE `arm' (no controls) ##########"
  cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(30) rseed(4242) dots(0)
  dis "C_`arm'_RC=" _rc
  dis "C_`arm'_N=" e(N) "  NCLUST=" e(N_clust) "  BOOTDF=" e(boot_df) "  DFR=" e(df_r)
  dis "C_`arm'_ATT=" %21x el(e(b),1,1)
  dis "C_`arm'_FAILED=" e(boot_failed) "  COMPLETED=" e(boot_completed)
}

foreach arm in Zc C {
  use "sk_b4c_`arm'.dta", clear
  dis "########## ROUTE `arm' (controls(x)) ##########"
  cap noi wooldid y unit year ttre, att controls(x) cluster(unit) bootstrap(standard) reps(30) rseed(4242) dots(0)
  dis "C_`arm'_RC=" _rc
  dis "C_`arm'_N=" e(N) "  NCLUST=" e(N_clust) "  BOOTDF=" e(boot_df) "  DFR=" e(df_r)
  dis "C_`arm'_ATT=" %21x el(e(b),1,1)
  dis "C_`arm'_FAILED=" e(boot_failed) "  COMPLETED=" e(boot_completed)
}

dis "PROBE_BR4C_DONE"
