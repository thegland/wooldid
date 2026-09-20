*=============================================================================
* p10 — SELF-CONTAINED REPRO.
* CLAIM: under bootstrap(), pretrendpower()'s Roth diagnostics are computed on the
* bootstrap event-study covariance with NO check that the covariance has enough
* replications behind it.  The only guard is invsym's rank test on the k x k pre block,
* which stops firing as soon as completed reps > k.  At reps = k+1 the run prints a
* pre-trend slope, a level() confidence interval, a joint pre-trend p and two
* minimum-detectable-trend columns with NO caution of any kind, and those numbers are
* arbitrary: they move by an order of magnitude with the replication count and with the
* seed, and the slope's own point estimate moves (Roth's slope is GLS in the covariance).
* Run: stata-mp -e do p10_ptpboundary.do
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

local PB cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) pretrendpower(80)

dis "##### ANALYTIC #####"
cap noi wooldid y state t onset, `PB'
dis "ANA_NPRE=" e(roth_npre_att) " KEFF=" e(roth_keff_att) " SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att)
dis "ANA_CI=" e(roth_slope_ll_att) " " e(roth_slope_ul_att) " PREP=" e(roth_prep_att)
matrix PA = e(roth_power_att)
matrix list PA, format(%12.6g)

* boundary sweep: k = 6 pre coefficients, so KEFF stops being clipped once completed > 6
foreach R in 7 8 10 15 25 50 100 400 {
  dis "##### reps(`R') seed 4242 #####"
  cap noi wooldid y state t onset, `PB' bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "S`R'_RC=" _rc " COMPLETED=" e(boot_completed)
  dis "S`R'_NPRE=" e(roth_npre_att) " KEFF=" e(roth_keff_att)
  dis "S`R'_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
  dis "S`R'_PREP=" e(roth_prep_att) " W=" e(roth_W_att) " NISSIZE=" e(roth_nissize_att)
  cap noi matrix list e(roth_power_att), format(%12.6g)
  dis "##### END reps(`R') #####"
}

* seed sensitivity at the boundary: same reps, five seeds
foreach S in 4242 777 99 1234 55555 {
  dis "##### reps(10) seed `S' #####"
  cap noi wooldid y state t onset, `PB' bootstrap(standard) reps(10) rseed(`S') dots(0)
  dis "Z`S'_KEFF=" e(roth_keff_att) " SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att) " PREP=" e(roth_prep_att)
  cap noi matrix list e(roth_power_att), format(%12.6g)
  dis "##### END reps(10) seed `S' #####"
}
dis "P10_DONE"
