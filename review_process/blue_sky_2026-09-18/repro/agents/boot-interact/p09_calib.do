*=============================================================================
* p09 — (A) calibration: does Stata's own bootstrap prefix behave the same way?
*       (B) the realistic no-caution zone: reps(50) with ~90 reported coefficients,
*           where 2/reps = .04 < .05 so NO low-reps caution fires at all.
*=============================================================================
clear all
set more off
set processors 2
set maxvar 20000
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
gen byte grp = mod(state,3)
expand 20
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + sfe + 0.05*t + te + rnormal()*0.8
save "p09_data.dta", replace

* ---------------- (A) calibration: Stata's own bootstrap prefix ----------------
dis "##### A: stata bootstrap prefix #####"
cap noi bootstrap, reps(8) cluster(state) seed(4242) nodots: regress y i.t i.state
dis "SBP_RC=" _rc
if _rc == 0 {
  matrix SV = e(V)
  mata: V=st_matrix("SV"); V=(V+V')/2; st_numscalar("rk",rank(V))
  dis "SBP_RANK=" rk " of " colsof(SV) " REPS=" e(N_reps)
  cap noi test (2.t=0) (3.t=0) (4.t=0) (5.t=0) (6.t=0) (7.t=0)
  dis "SBP_TEST6_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
}
dis "##### END A #####"

* ---------------- (B) large k, reps(50): does ANY caution fire? ----------------
local BIG cluster(state) esfixedbaseperiod esprelength(8) espostlength(8) subgroups(grp)
dis "##### B analytic #####"
cap noi wooldid y state t onset, `BIG'
dis "BAN_RC=" _rc " K=" colsof(e(b))
matrix BB = e(b)
local bn : colnames BB
local esall
foreach nm of local bn {
  if strpos("`nm'","ES_")>0 local esall `esall' `nm'
}
local tcall
foreach nm of local esall {
  local tcall `tcall' (`nm' = 0)
}
dis "BAN_NES=" wordcount("`esall'")
cap noi test `tcall'
dis "BAN_TESTALL_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)

dis "##### B reps(50) - CAUTION WATCH START #####"
cap noi wooldid y state t onset, `BIG' bootstrap(standard) reps(50) rseed(4242) dots(0)
dis "##### B reps(50) - CAUTION WATCH END #####"
dis "B50_RC=" _rc " K=" colsof(e(b)) " COMPLETED=" e(boot_completed) " NV=" e(boot_N_V)
matrix BV = e(V)
mata: V=st_matrix("BV"); V=(V+V')/2; st_numscalar("rk",rank(V)); ev=symeigenvalues(V)
mata: st_numscalar("emin",min(ev))
dis "B50_VRANK=" rk " of " colsof(BV) " EMIN=" %12.4e emin
cap noi test `tcall'
dis "B50_TESTALL_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop) " NCON=" wordcount("`esall'")

* ---------------- (C) REALISTIC reps: 15 pre-period coefficients, reps(60) ----------------
* 2/60 = .033 < .05 so NO low-reps caution fires at all.
local CB cluster(state) esfixedbaseperiod esprelength(16) espostlength(4) pretrendpower(80) simultaneousbands
dis "##### C analytic #####"
cap noi wooldid y state t onset, `CB'
dis "CAN_RC=" _rc " NPRE=" e(roth_npre_att) " KEFF=" e(roth_keff_att) " DFR=" e(df_r)
dis "CAN_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " PREP=" e(roth_prep_att) " W=" e(roth_W_att)
cap noi matrix list e(roth_power_att), format(%12.6g)
dis "CAN_BANDCRIT=" e(esband_crit_att)
foreach R in 60 300 {
  dis "##### C reps(`R') - CAUTION WATCH START #####"
  cap noi wooldid y state t onset, `CB' bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "##### C reps(`R') - CAUTION WATCH END #####"
  local n = e(boot_completed)
  local k = e(roth_keff_att)
  local W = e(roth_W_att)
  local ph = .
  if `n' > `k' local ph = Ftail(`k', `n'-`k', (`n'-`k')/((`n'-1)*`k')*`W')
  dis "C`R'_RC=" _rc " COMPLETED=" `n' " KEFF=" `k' " NPRE=" e(roth_npre_att) " DFR=" e(df_r) " W=" `W'
  dis "C`R'_PREP_WOOLDID=" e(roth_prep_att) "  PREP_HOTELLING=" `ph' "  INFL=" (`n'-1)/(`n'-`k')
  dis "C`R'_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
  cap noi matrix list e(roth_power_att), format(%12.6g)
  dis "C`R'_BANDCRIT=" e(esband_crit_att) " N=" e(esband_n_att)
}
dis "P09_DONE"
