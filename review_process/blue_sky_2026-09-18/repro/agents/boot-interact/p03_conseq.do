* p03 — CONSEQUENCE SURFACE of a rank-deficient bootstrap e(V)/e(V_es_*) at low reps.
* Long-ish event study so that #coefficients > reps.  Everything a user would plausibly
* do next: test, lincom, honestdid, simultaneousbands, esttab-style sqrt(diag).
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
save "p03_data.dta", replace

local OPTS cluster(state) esfixedbaseperiod esprelength(6) espostlength(6)

dis "===== ANALYTIC ====="
cap noi wooldid y state t onset, `OPTS' simultaneousbands
dis "AN_RC=" _rc
dis "AN_K=" colsof(e(b))
matrix ANB = e(b)
local anames : colnames ANB
dis "AN_NAMES: `anames'"
local esall
foreach nm of local anames {
  if strpos("`nm'","ES_att_")>0 local esall `esall' `nm'
}
dis "AN_NES=" wordcount("`esall'") " :: `esall'"
cap noi test `=subinstr("`esall'"," "," = 0) (",.)' = 0
dis "AN_TESTALL_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
dis "AN_BANDCRIT=" e(esband_crit_att) " ENGINE=`e(esband_engine)'"
matrix ANVES = e(V_es_att)
mata: V=st_matrix("ANVES"); V=(V+V')/2; st_numscalar("rk",rank(V))
dis "AN_VES_RANK=" rk " of " colsof(ANVES)

foreach R in 8 12 2000 {
  dis "===== REPS=`R' ====="
  cap noi wooldid y state t onset, `OPTS' simultaneousbands bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "C`R'_RC=" _rc
  dis "C`R'_K=" colsof(e(b))
  dis "C`R'_COMPLETED=" e(boot_completed) " NV=" e(boot_N_V) " DFR=" e(df_r)
  matrix BV = e(V)
  mata: V=st_matrix("BV"); V=(V+V')/2; st_numscalar("rk",rank(V)); ev=symeigenvalues(V); st_numscalar("emin",min(ev))
  dis "C`R'_V_RANK=" rk " of " colsof(BV) " EMIN=" %12.4e emin
  matrix BVES = e(V_es_att)
  mata: V=st_matrix("BVES"); V=(V+V')/2; st_numscalar("rk2",rank(V)); ev=symeigenvalues(V); st_numscalar("emin2",min(ev))
  dis "C`R'_VES_RANK=" rk2 " of " colsof(BVES) " EMIN=" %12.4e emin2
  * (a) a hand joint test over the WHOLE event study
  cap noi test `=subinstr("`esall'"," "," = 0) (",.)' = 0
  dis "C`R'_TESTALL_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop) " dropped=`r(dropped)'"
  * (b) the natural hand pre-trend test
  local espre
  foreach nm of local esall {
    if strpos("`nm'","pre")>0 local espre `espre' `nm'
  }
  cap noi test `=subinstr("`espre'"," "," = 0) (",.)' = 0
  dis "C`R'_TESTPRE_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop) " n=" wordcount("`espre'")
  * (c) simultaneous band critical value
  dis "C`R'_BANDCRIT=" e(esband_crit_att) " N=" e(esband_n_att) " ENGINE=`e(esband_engine)'"
  * (d) honestdid on the bootstrap event-study covariance
  matrix HB = e(b_es_att)
  matrix HV = e(V_es_att)
  local hn : colnames HB
  local npre 0
  local j 0
  local preidx
  local postidx
  foreach nm of local hn {
    local ++j
    if substr("`nm'",1,3)=="pre" local preidx `preidx' `j'
    if substr("`nm'",1,4)=="post" local postidx `postidx' `j'
  }
  local p1 : word 1 of `postidx'
  dis "C`R'_HONEST_PRE=`preidx' POST1=`p1'"
  cap noi honestdid, b(HB) vcov(HV) pre(`preidx') post(`p1') mvec(0 0.5 1) parallel(0)
  dis "C`R'_HONESTDID_RC=" _rc
  cap noi matrix list r(table)
}
dis "P03_DONE"
