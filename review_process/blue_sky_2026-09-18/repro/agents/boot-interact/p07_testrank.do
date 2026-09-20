*=============================================================================
* p07 — SELF-CONTAINED REPRO.
* CLAIM: after bootstrap() with reps() < the number of reported coefficients, the
* posted e(V) is SINGULAR (rank = completed replications - 1) and nothing in the run
* says so.  A hand `test` over more coefficients than that rank returns a p-value with
* full degrees of freedom, no dropped constraint and no warning.  wooldid's OWN joint
* test on the same coefficients IS guarded (eigen screen) and is withheld.
* Run: stata-mp -e do p07_testrank.do   (needs the pinned build on the adopath)
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

local BASE cluster(state) esfixedbaseperiod esprelength(6) espostlength(6) jointtests

dis "===== ANALYTIC ====="
qui wooldid y state t onset, `BASE'
matrix ANB = e(b)
local anames : colnames ANB
local esall
local espre
foreach nm of local anames {
  if strpos("`nm'","ES_att_")>0 local esall `esall' `nm'
  if strpos("`nm'","ES_att_pre")>0 local espre `espre' `nm'
}
local tcall
foreach nm of local esall {
  local tcall `tcall' (`nm' = 0)
}
local tcpre
foreach nm of local espre {
  local tcpre `tcpre' (`nm' = 0)
}
dis "NES=" wordcount("`esall'") " NPRE=" wordcount("`espre'")
cap noi test `tcall'
dis "AN_TESTALL_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
cap noi test `tcpre'
dis "AN_TESTPRE_RC=" _rc " p=" r(p) " df=" r(df) " drop=" r(drop)
dis "AN_JOINT_HORIZ=" e(joint_horiz) " JOINT_CHORIZ=" e(joint_Choriz)

foreach R in 4 8 12 30 400 {
  dis "===== REPS=`R' ====="
  cap noi wooldid y state t onset, `BASE' bootstrap(standard) reps(`R') rseed(4242) dots(0)
  dis "R`R'_RC=" _rc
  dis "R`R'_COMPLETED=" e(boot_completed) " NV=" e(boot_N_V) " K=" colsof(e(b)) " DFR=" e(df_r)
  matrix BV = e(V)
  mata: V=st_matrix("BV"); V=(V+V')/2; st_numscalar("rk",rank(V)); ev=symeigenvalues(V)
  mata: st_numscalar("emin",min(ev)); st_numscalar("emax",max(ev))
  dis "R`R'_VRANK=" rk " of " colsof(BV) " EMIN=" %12.4e emin " EMAX=" %12.4e emax
  cap noi test `tcall'
  dis "R`R'_TESTALL_RC=" _rc " p=" %21x r(p) " pd=" r(p) " df=" r(df) " drop=" r(drop) " dropped=[`r(dropped)']"
  cap noi test `tcpre'
  dis "R`R'_TESTPRE_RC=" _rc " p=" %21x r(p) " pd=" r(p) " df=" r(df) " drop=" r(drop)
  * wooldid's own guarded joint pre-trend test on the SAME coefficients
  dis "R`R'_JOINT_HORIZ=" e(joint_horiz) " JOINT_CHORIZ=" e(joint_Choriz) " JOINT_PRE=" e(joint_prehoriz)
  dis "R`R'_JTABSENT=[`e(boot_jt_absent)'] JTFAM=[`e(boot_jt_family)']"
  cap noi matrix list e(boot_jointWobs), format(%12.6g)
}
dis "P07_DONE"
