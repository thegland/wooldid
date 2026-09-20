* r4_impvalues_cwm.do -- does customweightmultiplier() reach the impvalues estimates?
* r3 found bit-identical estimates with and without cwm under impvalues, but there the
* multiplier was constant WITHIN the cluster, where an equal-weight cross-cluster
* aggregation makes it cancel arithmetically.  Here the multiplier varies WITHIN the
* cluster (counties of very different population inside one state), so a genuine
* application of the multiplier MUST move the numbers.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31337
set sortseed 31337
which wooldid

* 6 states x 4 counties = 24 counties, 10 years, 3 establishments per county-year.
* In each state: one county treated 2011, one treated 2013, two never treated.
clear
set obs 24
gen int county = _n
gen int state  = ceil(county/4)
gen int slot   = mod(county-1,4)+1
gen int onset = .
replace onset = 2011 if slot == 1
replace onset = 2013 if slot == 2
* population varies sharply WITHIN state
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
compress
save "panel_r4.dta", replace

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

*--- J0: impvalues, cluster(state), NO cwm
use "panel_r4.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(state) impvalues
dis "J0_RC=" _rc
showb , tag(J0_imp_nocwm)

*--- J1: impvalues, cluster(state), cwm(pop) varying WITHIN state
use "panel_r4.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(state) impvalues customweightmultiplier(pop)
dis "J1_RC=" _rc
showb , tag(J1_imp_cwmpop)

*--- J2: impvalues, cluster(state), cwm varying WITHIN the (i,t) CELL
use "panel_r4.dta", clear
gen double epop = pop*(1 + 0.8*(estab-2))
cap noi wooldid y county year onset [aw=sw], att att_it cluster(state) impvalues customweightmultiplier(epop)
dis "J2_RC=" _rc
showb , tag(J2_imp_cwmepop)

*--- J3/J4: the SAME two runs without impvalues, to show cwm is not inert generally
use "panel_r4.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(state)
dis "J3_RC=" _rc
showb , tag(J3_noimp_nocwm)

use "panel_r4.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(state) customweightmultiplier(pop)
dis "J4_RC=" _rc
showb , tag(J4_noimp_cwmpop)

*--- J5/J6: event-study rows under impvalues, with and without cwm
use "panel_r4.dta", clear
cap noi wooldid y county year onset [aw=sw], att cluster(state) impvalues esprelength(2) espostlength(3)
dis "J5_RC=" _rc
cap matrix list e(b_es_att), format(%21x)

use "panel_r4.dta", clear
cap noi wooldid y county year onset [aw=sw], att cluster(state) impvalues esprelength(2) espostlength(3) customweightmultiplier(pop)
dis "J6_RC=" _rc
cap matrix list e(b_es_att), format(%21x)

dis "R4_DONE"
