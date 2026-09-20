* p14 — does BI-1 survive bootdf(normal) (chi2 reference), bootstrap(bayes), boottstat,
*       and an unclustered bootstrap?  Same data, 6 pre coefficients, reps(8).
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
local V1 "bootdf(normal)"
local V2 "bootdf(3)"
local V3 "boottstat"
local V4 "bootdisplay(bc)"
forval k = 1/4 {
  dis "##### variant `k': [`V`k''] #####"
  cap noi wooldid y state t onset, `PB' bootstrap(standard) reps(8) rseed(4242) dots(0) `V`k''
  dis "V`k'_RC=" _rc " COMPLETED=" e(boot_completed) " DFR=" e(df_r) " BOOTDF=" e(boot_df)
  dis "V`k'_KEFF=" e(roth_keff_att) " NPRE=" e(roth_npre_att) " W=" e(roth_W_att) " PREP=" e(roth_prep_att)
  dis "V`k'_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
  cap noi matrix list e(roth_power_att), format(%12.6g)
  dis "##### END variant `k' #####"
}
dis "##### variant 5: bayes #####"
cap noi wooldid y state t onset, `PB' bootstrap(bayes) reps(8) rseed(4242) dots(0)
dis "V5_RC=" _rc " COMPLETED=" e(boot_completed) " DFR=" e(df_r)
dis "V5_KEFF=" e(roth_keff_att) " W=" e(roth_W_att) " PREP=" e(roth_prep_att)
dis "V5_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
cap noi matrix list e(roth_power_att), format(%12.6g)
dis "##### variant 6: no cluster() #####"
cap noi wooldid y state t onset, esfixedbaseperiod esprelength(6) espostlength(4) pretrendpower(80) bootstrap(standard) reps(8) rseed(4242) dots(0)
dis "V6_RC=" _rc " COMPLETED=" e(boot_completed) " DFR=" e(df_r) " BOOTDF=" e(boot_df)
dis "V6_KEFF=" e(roth_keff_att) " W=" e(roth_W_att) " PREP=" e(roth_prep_att)
dis "P14_DONE"
