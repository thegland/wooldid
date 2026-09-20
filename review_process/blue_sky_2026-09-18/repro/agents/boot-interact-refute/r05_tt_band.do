* r05_tt_band.do -- SELF-CONTAINED. BOOTX-4 (timetrends x cluster(i) x bootstrap)
* and BOOTX-5 (simultaneousbands critical value vs draw count), plus the BOOTX-3
* hand-test / e(V)-rank leg.  My own data.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

* ---------- design B: i = state, cluster(state) (cluster ON the cohort id)
clear
set obs 40
gen state = _n
gen onset = .
replace onset = 6 + mod(_n,5) if mod(_n,2)==0
expand 16
bysort state: gen t = _n
bysort state (t): gen double ueff = rnormal(0,.4) if _n==1
bysort state (t): replace ueff = ueff[_n-1] if _n>1
bysort state (t): gen double utr  = rnormal(0,.05) if _n==1
bysort state (t): replace utr = utr[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .1*rel, 0)
gen double y = 2 + .05*t + ueff + utr*t + te + e
drop rel te
tempfile B
save `B'

dis "======== A. timetrends x cluster(i): ANALYTIC"
cap noi wooldid y state t onset, cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) timetrends
dis "TT_A_RC=" _rc
dis "TT_A_vcvstate=" "`e(wdid_vcvstate)'"
cap confirm matrix e(V)
dis "TT_A_confirmV_rc=" _rc
cap confirm matrix e(V_es_att)
dis "TT_A_confirmVes_rc=" _rc
tempname M
cap matrix `M' = e(wdidmainresults)
dis "TT_A_mainres_rc=" _rc
cap matrix list `M'

dis "======== A2. timetrends x cluster(i): BOOTSTRAP reps(60)"
cap noi wooldid y state t onset, cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) timetrends bootstrap(standard) reps(60) rseed(4242)
dis "TT_B_RC=" _rc
dis "TT_B_vcvstate=" "`e(wdid_vcvstate)'"
dis "TT_B_v_posted=" e(boot_v_posted) " completed=" e(boot_completed)
cap confirm matrix e(V)
dis "TT_B_confirmV_rc=" _rc
cap confirm matrix e(V_boot)
dis "TT_B_confirmVboot_rc=" _rc
cap confirm matrix e(V_analytic)
dis "TT_B_confirmVanalytic_rc=" _rc
cap confirm matrix e(V_es_att)
dis "TT_B_confirmVes_rc=" _rc
cap matrix `M' = e(wdidmainresults)
cap matrix list `M'
cap noi test ES_att_post0 = 0
dis "TT_B_handtest_rc=" _rc
cap noi lincom _b[ES_att_post0]
dis "TT_B_lincom_rc=" _rc

dis "======== B. simultaneousbands: analytic vs bootstrap draw counts"
use `B', clear
local OPTB "cluster(state) esfixedbaseperiod esprelength(5) espostlength(5) simultaneousbands"
cap noi wooldid y state t onset, `OPTB'
dis "SB_A_RC=" _rc
dis "SBRES|ANALYTIC|crit=" %12.8g e(esband_crit_att) "|n=" e(esband_n_att) "|engine=" "`e(esband_engine)'"
foreach R in 8 20 40 100 400 {
  qui cap noi wooldid y state t onset, `OPTB' bootstrap(standard) reps(`R') rseed(5150)
  dis "SB_`R'_RC=" _rc
  dis "SBRES|BOOT`R'|crit=" %12.8g e(esband_crit_att) "|n=" e(esband_n_att) "|engine=" "`e(esband_engine)'" "|completed=" e(boot_completed)
}

dis "======== C. BOOTX-3: rank of e(V) and a hand k-dim test"
use `B', clear
local OPTC "cluster(state) esfixedbaseperiod esprelength(5) espostlength(5) jointtests"
cap noi wooldid y state t onset, `OPTC'
dis "C_A_RC=" _rc
tempname VV
matrix `VV' = e(V)
mata: st_local("rk", strofreal(rank(st_matrix("`VV'"))))
dis "CRES|ANALYTIC|colsofb=" colsof(e(b)) "|rankV=`rk'|dfr=" e(df_r)
cap noi test ES_att_pre6 ES_att_pre5 ES_att_pre4 ES_att_pre3 ES_att_pre2
dis "C_A_test_rc=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
foreach R in 8 20 60 400 {
  qui cap noi wooldid y state t onset, `OPTC' bootstrap(standard) reps(`R') rseed(5150)
  dis "C_`R'_RC=" _rc
  matrix `VV' = e(V)
  mata: st_local("rk", strofreal(rank(st_matrix("`VV'"))))
  mata: st_local("mine", strofreal(min(symeigenvalues(makesymmetric(st_matrix("`VV'")))), "%12.4e"))
  local nb = e(boot_completed)
  local nv = e(boot_N_V)
  cap noi test ES_att_pre6 ES_att_pre5 ES_att_pre4 ES_att_pre3 ES_att_pre2
  dis "CRES|BOOT`R'|colsofb=" colsof(e(b)) "|rankV=`rk'|mineig=`mine'|completed=`nb'|NV=`nv'|test_rc=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
}
dis "DONE"
