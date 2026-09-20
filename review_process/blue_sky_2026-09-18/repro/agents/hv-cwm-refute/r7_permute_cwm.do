* r7_permute_cwm.do -- closes the gap r6 left: does permute() carry the multiplier into
* its replications, or is the permutation p computed for the UNWEIGHTED estimand while the
* point estimate is the weighted one?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
which wooldid

clear
set obs 24
gen int county = _n
gen int slot   = mod(county-1,4)+1
gen int onset = .
replace onset = 2011 if slot == 1
replace onset = 2013 if slot == 2
gen double pop = 1
replace pop = 50 if slot == 1
replace pop = 20 if slot == 2
replace pop = 4  if slot == 3
gen double cfe = rnormal(0,1)
gen double teff = 1.0 + 0.05*pop + rnormal(0,0.3)
expand 10
bysort county: gen int year = 2007 + _n
expand 3
bysort county year: gen int estab = _n
gen double sw = exp(rnormal(1.0,0.5))
gen byte post = (onset<. & year>=onset)
gen double y = 3 + cfe + 0.11*(year-2007) + teff*post + rnormal(0,0.8)
save "panel_r7.dta", replace

capture program drop showse
program define showse
  syntax , TAG(string)
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    dis "`tag'_NOMAIN"
    exit
  }
  matrix WM = e(wdidmainresults)
  dis "==== MAINRESULTS `tag' ===="
  matrix list WM, format(%20.14g)
  dis "==== END `tag' ===="
end

use "panel_r7.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county) permute(shuffle) permreps(80) rseed(4242) dots(0)
dis "Q0_RC=" _rc
if _rc == 0 showse , tag(Q0_perm_nocwm)

use "panel_r7.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county) permute(shuffle) permreps(80) rseed(4242) dots(0) customweightmultiplier(pop)
dis "Q1_RC=" _rc
if _rc == 0 showse , tag(Q1_perm_cwm)

dis "R7_DONE"
