* r3_note_honesty.do -- MY OWN angle.
* The only runtime signal a customweightmultiplier() user gets is
*   "Note: all Estimates Shown Apply <var> Weights Post-Estimation"  (ado:11369/11529).
* Is that note TRUE of every number it sits above?  Probe the display paths where
* wooldid documents that it aggregates differently: impvalues, subgroups(), and the
* ame family under contreat().
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 480915
set sortseed 480915
which wooldid

clear
set obs 20
gen int county = _n
gen int onset = .
replace onset = 2011 if inrange(county,1,4)
replace onset = 2013 if inrange(county,5,8)
replace onset = 2015 if inrange(county,9,11)
gen byte region = mod(county,2)
gen double pop = 1
replace pop = 40 if inlist(county,1,6)
replace pop = 18 if inlist(county,2,9)
replace pop = 6  if inlist(county,5,12,15)
replace pop = 3  if inlist(county,3,10,17)
gen double cfe = rnormal(0,1)
gen double teff = 1.2 + 0.055*pop + rnormal(0,0.3)
expand 10
bysort county: gen int year = 2007 + _n
expand 4
bysort county year: gen int estab = _n
gen double sw = exp(rnormal(1.1,0.55))
gen double yfe = 0.12*(year-2007)
gen byte post = (onset<. & year>=onset)
* continuous dose, exactly zero outside treated cohort-periods (untreatedZero)
gen double dose = 0
replace dose = 0.5 + 0.35*(year-onset) + runiform() if post
gen double y = 3 + cfe + yfe + teff*post + 0.40*dose + rnormal(0,0.8)
compress
save "panel_r3.dta", replace

capture program drop showb
program define showb
  syntax , TAG(string)
  cap confirm matrix e(wdidmainresults)
  if _rc != 0 {
    dis "`tag'_NOMAIN"
    exit
  }
  matrix WM = e(wdidmainresults)
  local rn : rownames WM
  local k = 1
  foreach r of local rn {
    dis "`tag'|" "`r'" "|b=" %21x el(WM,`k',1) " (" %18.10f el(WM,`k',1) ")"
    local ++k
  }
end

*=================== impvalues x cwm ===================
use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county) impvalues
dis "I0_RC=" _rc
showb , tag(I0_nocwm)

use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county) impvalues customweightmultiplier(pop)
dis "I1_RC=" _rc
showb , tag(I1_cwm)

* the analytic (non-impvalues) counterparts, for contrast
use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county)
dis "I2_RC=" _rc
showb , tag(I2_nocwm_noimp)

use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county) customweightmultiplier(pop)
dis "I3_RC=" _rc
showb , tag(I3_cwm_noimp)

*=================== subgroups x cwm ===================
use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it subgroups(region)
dis "SG0_RC=" _rc
showb , tag(SG0_nocwm)

use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it subgroups(region) customweightmultiplier(pop)
dis "SG1_RC=" _rc
showb , tag(SG1_cwm)

*=================== contreat: does cwm move BOTH att and ame? ===================
use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it contreat(dose) contreatcontrols(untreatedZero)
dis "CT0_RC=" _rc
showb , tag(CT0_nocwm)

use "panel_r3.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it contreat(dose) contreatcontrols(untreatedZero) customweightmultiplier(pop)
dis "CT1_RC=" _rc
showb , tag(CT1_cwm)

dis "R3_DONE"
