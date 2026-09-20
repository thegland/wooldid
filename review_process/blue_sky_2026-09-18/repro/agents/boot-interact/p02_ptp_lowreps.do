* p02 — pretrendpower() under bootstrap() at LOW reps: is the Roth block silent about
* a covariance estimated from a handful of draws?  Same data, same design, only reps moves.
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
replace onset = 5 if state <= 10
replace onset = 7 if state > 10 & state <= 20
replace onset = 9 if state > 20 & state <= 30
gen double sfe = rnormal()
expand 12
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + sfe + 0.05*t + te + rnormal()*0.8
save "p02_data.dta", replace

local OPTS cluster(state) esfixedbaseperiod esprelength(4) espostlength(4) pretrendpower(50 80)

dis "===== ANALYTIC ====="
cap noi wooldid y state t onset, `OPTS'
dis "AN_RC=" _rc
dis "AN_SLOPE=" %21x e(roth_slope_att) " SE=" %21x e(roth_slope_se_att)
dis "AN_SLOPE_D=" e(roth_slope_att) " SE_D=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
dis "AN_PREP=" e(roth_prep_att) " W=" e(roth_W_att) " KEFF=" e(roth_keff_att) " NPRE=" e(roth_npre_att)
matrix AN_P = e(roth_power_att)
matrix list AN_P

foreach S in 4242 777 99 1234 55555 {
  dis "===== BOOT reps(8) seed(`S') ====="
  cap noi wooldid y state t onset, `OPTS' bootstrap(standard) reps(8) rseed(`S') dots(0)
  dis "B8_`S'_RC=" _rc
  dis "B8_`S'_COMPLETED=" e(boot_completed) " DFR=" e(df_r) " BOOTDF=" e(boot_df)
  dis "B8_`S'_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
  dis "B8_`S'_PREP=" e(roth_prep_att) " W=" e(roth_W_att) " KEFF=" e(roth_keff_att) " NPRE=" e(roth_npre_att)
  cap noi matrix list e(roth_power_att)
}

foreach S in 4242 777 {
  dis "===== BOOT reps(300) seed(`S') ====="
  cap noi wooldid y state t onset, `OPTS' bootstrap(standard) reps(300) rseed(`S') dots(0)
  dis "B300_`S'_RC=" _rc
  dis "B300_`S'_COMPLETED=" e(boot_completed) " DFR=" e(df_r) " BOOTDF=" e(boot_df)
  dis "B300_`S'_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
  dis "B300_`S'_PREP=" e(roth_prep_att) " W=" e(roth_W_att) " KEFF=" e(roth_keff_att)
  cap noi matrix list e(roth_power_att)
}
dis "P02_DONE"
