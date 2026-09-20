clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---- simple staggered DGP, iid errors, no pre-trend ----
local N 80
local T 10
set obs `N'
gen long unit = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 5 if mod(unit,5)==1
replace onset = 7 if mod(unit,5)==2
replace onset = 9 if mod(unit,5)==3
expand `T'
bysort unit: gen int t = _n
gen double f_t = 0.1*t
gen byte D = (onset<. & t>=onset)
gen double y = a_i + f_t + 1.0*D + rnormal()

timer clear 1
timer on 1
wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(50 80)
timer off 1
timer list 1

dis "COLNAMES b_es_att: " 
mat list e(b_es_att)
mat list e(V_es_att), format(%21x)
dis "roth_slope_att = " %21x e(roth_slope_att)
dis "roth_slope_se_att = " %21x e(roth_slope_se_att)
dis "roth_W_att = " %21x e(roth_W_att)
dis "roth_prep_att = " %21x e(roth_prep_att)
dis "roth_hbar_att = " %21x e(roth_hbar_att)
dis "roth_npre_att = " e(roth_npre_att)
dis "roth_keff_att = " e(roth_keff_att)
dis "roth_nissize_att = " %21x e(roth_nissize_att)
dis "df_r = " e(df_r)
dis "N = " e(N)
mat list e(roth_power_att), format(%21x)
dis "estimands: " e(wdid_estimands)
dis "vcvstate: " e(wdid_vcvstate)
ereturn list
