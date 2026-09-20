* p05 — bootstrap() with cluster() STRICTLY FINER than i, and CROSS-CUTTING i.
* (the two-regime bootstrap unit rule, ado:8987-8998 / 9323-9349, shipped 2026-09-16)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

* 40 counties nested in 10 states; i = state (the cohort), cluster = county (finer)
clear
set obs 40
gen int county = _n
gen int st = ceil(county/4)
gen int onset = .
replace onset = 5 if st <= 3
replace onset = 7 if st > 3 & st <= 6
replace onset = 9 if st > 6 & st <= 8
gen byte xcut = mod(county,7) + 1
gen double cfe = rnormal()
expand 12
bysort county: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + cfe + 0.05*t + te + rnormal()*0.8
save "p05_data.dta", replace

local BASE esfixedbaseperiod esprelength(3) espostlength(3)

foreach CL in "st" "county" "xcut" "" {
  local tag = cond("`CL'"=="", "none", "`CL'")
  dis "##### CLUSTER=`tag' #####"
  cap noi wooldid y st t onset, `BASE' `=cond("`CL'"=="","","cluster(`CL')")'
  dis "AN_`tag'_RC=" _rc
  if _rc==0 {
    matrix VA = e(V)
    dis "AN_`tag'_SE=" %12.6g sqrt(VA[1,1]) " DFR=" e(df_r)
  }
  cap noi wooldid y st t onset, `BASE' `=cond("`CL'"=="","","cluster(`CL')")' bootstrap(standard) reps(200) rseed(31) dots(0)
  dis "BO_`tag'_RC=" _rc
  if _rc==0 {
    matrix VB = e(V)
    matrix VA2 = e(V_analytic)
    dis "BO_`tag'_SEBOOT=" %12.6g sqrt(VB[1,1]) " SEANA=" %12.6g sqrt(VA2[1,1]) " RATIO=" %8.4f sqrt(VB[1,1]/VA2[1,1])
    dis "BO_`tag'_COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed) " NV=" e(boot_N_V) " DFR=" e(df_r) " BOOTDF=" e(boot_df)
    dis "BO_`tag'_VCV=`e(wdid_vcvstate)'"
  }
  * and the cell-weighted estimand, which takes the PHYSICAL resampling leg
  cap noi wooldid y st t onset, `BASE' att_it `=cond("`CL'"=="","","cluster(`CL')")' bootstrap(standard) reps(200) rseed(31) dots(0)
  dis "BOIT_`tag'_RC=" _rc
  if _rc==0 {
    matrix VB2 = e(V)
    matrix VA3 = e(V_analytic)
    dis "BOIT_`tag'_SEBOOT=" %12.6g sqrt(VB2[1,1]) " SEANA=" %12.6g sqrt(VA3[1,1]) " RATIO=" %8.4f sqrt(VB2[1,1]/VA3[1,1])
    dis "BOIT_`tag'_COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed)
  }
  * boottstat takes the physical leg under either regime
  cap noi wooldid y st t onset, `BASE' `=cond("`CL'"=="","","cluster(`CL')")' bootstrap(standard) reps(60) rseed(31) dots(0) boottstat
  dis "BOTS_`tag'_RC=" _rc
  if _rc==0 dis "BOTS_`tag'_COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed)
}
dis "P05_DONE"
