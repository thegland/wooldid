* s8_boot.do -- SELF-CONTAINED. Unseeded probe: pretrendpower() under bootstrap().
* sthlp:789 says that under bootstrap() the diagnostics run on the BOOTSTRAP covariance.
* Question: what reference distribution does the joint pre-trend test then use, and does
* it know how many replications the covariance was built from?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 6101993
set sortseed 6101993
which wooldid
clear
set obs 16
gen long st = _n
gen byte coh = 1 + floor((st-1)/4)
gen int ttre = .
replace ttre = 6  if coh == 1
replace ttre = 8  if coh == 2
replace ttre = 10 if coh == 3
gen double alpha_s = rnormal()*0.8
expand 12
bysort st: gen int t = _n
sort st t
gen double e = rnormal()
forvalues k = 2/12 {
  qui replace e = 0.5*e[_n-1] + sqrt(.75)*rnormal() if t == `k'
}
gen double y = 2 + alpha_s + 0.12*t - 0.02*t*t + 0.8*(ttre<. & t>=ttre) + e

dis "=========== analytic ==========="
cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "AN_RC=" _rc "  df_r=" e(df_r) "  W=" e(roth_W_att) "  prep=" e(roth_prep_att)
cap mat list e(roth_power_att), format(%12.6g)

foreach R in 39 199 {
  dis "=========== bootstrap(standard) reps(`R') default bootdf ==========="
  cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80) bootstrap(standard) reps(`R') rseed(4242)
  dis "B`R'_RC=" _rc
  dis "B`R' df_r=" e(df_r) "  boot_df=" e(boot_df) "  boot_completed=" e(boot_completed) "  boot_N_V=" e(boot_N_V)
  dis "B`R' roth_W=" e(roth_W_att) "  prep=" e(roth_prep_att) "  keff=" e(roth_keff_att) "  slope_se=" e(roth_slope_se_att)
  cap mat list e(roth_power_att), format(%12.6g)
  dis "=========== bootstrap(standard) reps(`R') bootdf(normal) ==========="
  cap noi wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80) bootstrap(standard) reps(`R') rseed(4242) bootdf(normal)
  dis "BN`R'_RC=" _rc
  dis "BN`R' df_r=" e(df_r) "  boot_df=" e(boot_df)
  dis "BN`R' roth_W=" e(roth_W_att) "  prep=" e(roth_prep_att) "  slope_se=" e(roth_slope_se_att)
  cap mat list e(roth_power_att), format(%12.6g)
}
dis "DONE"
