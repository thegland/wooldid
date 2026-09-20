*=============================================================================
* p09b — SELF-CONTAINED REPRO for the ZERO-CAUTION zone.
* (a) k = 72 reported coefficients, reps(50): rank(e(V)) = 49, NO caution at all,
*     a hand `test` over the 64 event-study coefficients returns p ~ 1e-33.
* (b) 11 pre-period coefficients, reps(60): pretrendpower() prints a joint pre-trend p
*     of .4984 where the Hotelling reference gives .6431, again with NO caution.
* Run: stata-mp -e do p09b_zerocaution.do
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

local BIG cluster(state) esfixedbaseperiod esprelength(8) espostlength(8) subgroups(grp)
dis "##### (a) analytic #####"
cap noi wooldid y state t onset, `BIG'
dis "AAN_RC=" _rc " K=" colsof(e(b))
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
dis "AAN_NES=" wordcount("`esall'")
cap noi test `tcall'
dis "AAN_TESTALL p=" r(p) " df=" r(df) " drop=" r(drop)
dis "##### (a) reps(50) CAUTION WATCH START #####"
cap noi wooldid y state t onset, `BIG' bootstrap(standard) reps(50) rseed(4242) dots(0)
dis "##### (a) reps(50) CAUTION WATCH END #####"
dis "A50_RC=" _rc " K=" colsof(e(b)) " COMPLETED=" e(boot_completed)
matrix BV = e(V)
mata: V=st_matrix("BV"); V=(V+V')/2; st_numscalar("rk",rank(V)); ev=symeigenvalues(V)
mata: st_numscalar("emin",min(ev))
dis "A50_VRANK=" rk " of " colsof(BV) " EMIN=" %12.4e emin
cap noi test `tcall'
dis "A50_TESTALL p=" r(p) " df=" r(df) " drop=" r(drop) " NCON=" wordcount("`esall'")

local CB cluster(state) esfixedbaseperiod esprelength(16) espostlength(4) pretrendpower(80) simultaneousbands
dis "##### (b) analytic #####"
cap noi wooldid y state t onset, `CB'
dis "BAN_RC=" _rc " NPRE=" e(roth_npre_att) " KEFF=" e(roth_keff_att) " DFR=" e(df_r)
dis "BAN_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " PREP=" e(roth_prep_att) " W=" e(roth_W_att) " BAND=" e(esband_crit_att)
dis "##### (b) reps(60) CAUTION WATCH START #####"
cap noi wooldid y state t onset, `CB' bootstrap(standard) reps(60) rseed(4242) dots(0)
dis "##### (b) reps(60) CAUTION WATCH END #####"
local n = e(boot_completed)
local k = e(roth_keff_att)
local W = e(roth_W_att)
local ph = .
if `n' > `k' local ph = Ftail(`k', `n'-`k', (`n'-`k')/((`n'-1)*`k')*`W')
dis "B60_RC=" _rc " COMPLETED=" `n' " KEFF=" `k' " NPRE=" e(roth_npre_att) " DFR=" e(df_r) " W=" `W'
dis "B60_PREP_WOOLDID=" e(roth_prep_att) "  PREP_HOTELLING=" `ph' "  INFL=" (`n'-1)/(`n'-`k')
dis "B60_SLOPE=" e(roth_slope_att) " SE=" e(roth_slope_se_att) " LL=" e(roth_slope_ll_att) " UL=" e(roth_slope_ul_att)
dis "B60_BAND=" e(esband_crit_att) " N=" e(esband_n_att)
dis "P09B_DONE"
