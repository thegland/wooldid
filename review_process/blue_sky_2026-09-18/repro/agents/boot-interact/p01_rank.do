* p01 — reps ladder: rank of posted e(V), and what `test` does with it.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------- realistic staggered panel: 40 states x 5 units x 12 periods ----------
clear
set obs 200
gen long unit = _n
gen int state = ceil(unit/5)
gen int onset = .
replace onset = 5 if state <= 10
replace onset = 7 if state > 10 & state <= 20
replace onset = 9 if state > 20 & state <= 30
gen double ufe = rnormal()
sort state unit
by state: gen double sfe = rnormal() if _n==1
by state: replace sfe = sfe[1]
expand 12
bysort unit: gen int t = _n
gen double tfe = 0.05*t
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + ufe + sfe + tfe + te + rnormal()*0.8

compress
save "p01_data.dta", replace

* ---------- analytic baseline ----------
cap noi wooldid y unit t onset, cluster(state) esfixedbaseperiod esprelength(4) espostlength(6)
dis "AN_RC=" _rc
dis "AN_K=" colsof(e(b))
dis "AN_VCVSTATE=" "`e(wdid_vcvstate)'"
dis "AN_DFR=" e(df_r)
matrix AN_b = e(b)
local anames : colnames AN_b
dis "AN_NAMES: `anames'"

* pick three event-study columns for the test
local c1
local c2
local c3
local j 0
foreach nm of local anames {
  local ++j
  if strpos("`nm'","post1")>0 & "`c1'"=="" local c1 `nm'
  if strpos("`nm'","post2")>0 & "`c2'"=="" local c2 `nm'
  if strpos("`nm'","post3")>0 & "`c3'"=="" local c3 `nm'
}
dis "TESTCOLS: `c1' | `c2' | `c3'"

cap noi test `c1' = 0
dis "AN_TEST1_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
cap noi test (`c1'=0) (`c2'=0) (`c3'=0)
dis "AN_TEST3_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop) " dropped=" "`r(dropped)'"

* ---------- reps ladder ----------
foreach R in 8 19 21 40 200 {
  dis "===== REPS=`R' ====="
  cap noi wooldid y unit t onset, cluster(state) esfixedbaseperiod esprelength(4) espostlength(6) bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "B`R'_RC=" _rc
  dis "B`R'_K=" colsof(e(b))
  dis "B`R'_VCVSTATE=" "`e(wdid_vcvstate)'"
  dis "B`R'_COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed)
  dis "B`R'_DFR=" e(df_r) " DFRAN=" e(df_r_analytic) " BOOTDF=" e(boot_df)
  dis "B`R'_VPOSTED=" e(boot_v_posted) " NV=" e(boot_N_V)
  matrix BV = e(V)
  mata: V = st_matrix("BV")
  mata: st_numscalar("rk", rank(V))
  mata: ev = symeigenvalues((V+V')/2)
  mata: st_numscalar("emin", min(ev)); st_numscalar("emax", max(ev))
  dis "B`R'_RANK=" rk " of " colsof(BV) " EMIN=" %12.4e emin " EMAX=" %12.4e emax
  * provenance
  cap confirm matrix e(V_analytic)
  dis "B`R'_HAS_V_ANALYTIC_RC=" _rc
  cap confirm matrix e(V_es_analytic_att)
  dis "B`R'_HAS_VESANA_RC=" _rc
  cap noi dis "B`R'_MRELDIF_VES=" mreldif(e(V_es_att), e(V_es_analytic_att))
  dis "B`R'_VESFAILED=" "`e(boot_ves_failed)'"
  * the tests
  cap noi test `c1' = 0
  dis "B`R'_TEST1_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
  cap noi test (`c1'=0) (`c2'=0) (`c3'=0)
  dis "B`R'_TEST3_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop) " dropped=" "`r(dropped)'"
  cap noi lincom `c1' - `c2'
  dis "B`R'_LINCOM_RC=" _rc " est=" r(estimate) " se=" r(se) " p=" r(p)
}

dis "P01_DONE"
