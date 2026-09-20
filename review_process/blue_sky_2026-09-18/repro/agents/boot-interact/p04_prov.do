* p04 — COVARIANCE PROVENANCE MATRIX after a bootstrap run.
* Which posted objects carry bootstrap inference and which are still analytic?
* Same data, same seed, run twice (analytic / bootstrap); compare every SE-bearing surface.
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
gen byte grp = mod(state,3)
expand 12
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely + 0.10*grp, 0)
gen double y = 1 + sfe + 0.05*t + te + rnormal()*0.8
save "p04_data.dta", replace

local OPTS cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) subgroups(grp) jointtests adversarial gteffects

dis "===== ANALYTIC ====="
cap noi wooldid y state t onset, `OPTS'
dis "AN_RC=" _rc
dis "AN_VCVSTATE=`e(wdid_vcvstate)'"
foreach M in wdidmainresults pwsgs gteffects escomp adversarial_line_att {
  cap confirm matrix e(`M')
  if _rc == 0 {
    matrix AN_`M' = e(`M')
    dis "AN_HAS_`M'=1 dim=" rowsof(AN_`M') "x" colsof(AN_`M')
    matrix list AN_`M', format(%12.6g)
  }
  else dis "AN_HAS_`M'=0"
}
dis "AN_ADVP_ATT=" e(adversarial_pval_att)
cap noi matrix list e(wdidsgresults_0), format(%12.6g)
dis "AN_ELIST_START"
ereturn list
dis "AN_ELIST_END"

dis "===== BOOTSTRAP ====="
cap noi wooldid y state t onset, `OPTS' bootstrap(standard) reps(200) rseed(4242) dots(0)
dis "BO_RC=" _rc
dis "BO_VCVSTATE=`e(wdid_vcvstate)'"
dis "BO_COMPLETED=" e(boot_completed)
foreach M in wdidmainresults pwsgs gteffects escomp adversarial_line_att {
  cap confirm matrix e(`M')
  if _rc == 0 {
    matrix BO_`M' = e(`M')
    dis "BO_HAS_`M'=1 dim=" rowsof(BO_`M') "x" colsof(BO_`M')
    matrix list BO_`M', format(%12.6g)
    cap noi dis "MRELDIF_`M'=" mreldif(AN_`M', BO_`M')
  }
  else dis "BO_HAS_`M'=0"
}
dis "BO_ADVP_ATT=" e(adversarial_pval_att) " ADVABSENT=[`e(boot_adv_absent)']"
cap noi matrix list e(wdidsgresults_0), format(%12.6g)
dis "BO_ELIST_START"
ereturn list
dis "BO_ELIST_END"
dis "P04_DONE"
