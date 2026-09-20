*=============================================================================
* p08 — the very bottom of the blessed reps() range, and pretrendpower() on a
* pre-period covariance block with FEWER replications than pre-periods.
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 40
gen int state = _n
gen int onset = .
replace onset = 9  if state <= 10
replace onset = 11 if state > 10 & state <= 20
replace onset = 13 if state > 20 & state <= 30
gen double sfe = rnormal()
expand 20
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + sfe + 0.05*t + te + rnormal()*0.8
save "p08_data.dta", replace

local BASE cluster(state) esfixedbaseperiod esprelength(6) espostlength(6)

* ---- A. reps 1,2,3: what does a "timing" run report? ----
foreach R in 1 2 3 5 {
  dis "##### A reps(`R') #####"
  cap noi wooldid y state t onset, `BASE' bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "A`R'_RC=" _rc
  if _rc == 0 {
    dis "A`R'_COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed) " NV=" e(boot_N_V) " VPOST=" e(boot_v_posted)
    matrix M = e(wdidmainresults)
    matrix list M, format(%12.6g)
    matrix BV = e(V)
    mata: V=st_matrix("BV"); V=(V+V')/2; st_numscalar("rk",rank(V))
    dis "A`R'_VRANK=" rk " of " colsof(BV)
    dis "A`R'_SEATT=" %12.6g sqrt(BV[1,1])
    cap noi matrix list e(boot_ci_percentile), format(%12.6g)
  }
  dis "##### END A reps(`R') #####"
}

* ---- B. pretrendpower with 5 pre-period coefficients and reps 3 / 4 / 6 / 400 ----
local PBASE cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) pretrendpower(50 80)
dis "##### B analytic #####"
cap noi wooldid y state t onset, `PBASE'
dis "BAN_RC=" _rc " NPRE=" e(roth_npre_att) " KEFF=" e(roth_keff_att)
dis "BAN_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " PREP=" e(roth_prep_att) " W=" e(roth_W_att)
cap noi matrix list e(roth_power_att), format(%12.6g)
foreach R in 3 4 6 400 {
  dis "##### B reps(`R') #####"
  cap noi wooldid y state t onset, `PBASE' bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "B`R'_RC=" _rc " COMPLETED=" e(boot_completed)
  dis "B`R'_NPRE=" e(roth_npre_att) " KEFF=" e(roth_keff_att) " ROTHOK=" e(roth_ok)
  dis "B`R'_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
  dis "B`R'_PREP=" e(roth_prep_att) " W=" e(roth_W_att)
  cap noi matrix list e(roth_power_att), format(%12.6g)
  dis "##### END B reps(`R') #####"
}
dis "P08_DONE"
