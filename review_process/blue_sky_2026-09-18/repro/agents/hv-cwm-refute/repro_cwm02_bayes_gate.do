* repro_cwm02_bayes_gate.do -- SELF-CONTAINED check of hv-cwm-02.
* Claim under test: bootstrap(bayes)'s Dirichlet draw travels as the SAMPLE weight, so it
* cancels out of the cell-weighted estimand weights, while customweightmultiplier() travels
* as a numerator-only multiplier and does not.  The two are therefore not inconsistent.
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_cwm02_bayes_gate.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909090
set sortseed 909090
which wooldid

clear
set obs 20
gen int county = _n
gen int onset = .
replace onset = 2011 if inrange(county,1,4)
replace onset = 2013 if inrange(county,5,8)
replace onset = 2015 if inrange(county,9,11)
gen double pop = 1
replace pop = 40 if inlist(county,1,6)
replace pop = 18 if inlist(county,2,9)
replace pop = 6  if inlist(county,5,12,15)
gen double cfe = rnormal(0,1)
gen double teff = 1.2 + 0.055*pop + rnormal(0,0.3)
expand 10
bysort county: gen int year = 2007 + _n
expand 4
bysort county year: gen int estab = _n
gen double sw = exp(rnormal(1.1,0.55))
gen byte post = (onset<. & year>=onset)
gen double y = 3 + cfe + 0.12*(year-2007) + teff*post + rnormal(0,0.8)
* a multiplier that varies WITHIN the (county,year) cell -- finer than the cell
gen double epop = pop*(1 + 0.8*(estab-2.5))
replace epop = 0.05 if epop <= 0
save "panel_cwm02.dta", replace

*--- K1: the multiplier folded into the SAMPLE WEIGHT (the route the Dirichlet draw takes).
*        w_att_it must then be identically 1: sum_cell(w*m)/sum_cell(w*m) = 1.
use "panel_cwm02.dta", clear
gen double wfold = sw*epop
cap noi wooldid y county year onset [aw=wfold], att att_it gteffects
dis "K1_RC=" _rc
matrix G1 = e(gteffects)
preserve
clear
svmat double G1, names(col)
qui sum w_att_it
dis "K1_w_att_it  min=" %20.17f r(min) "  max=" %20.17f r(max) "   (sample-weight route: CANCELS)"
restore

*--- K2: the SAME multiplier as customweightmultiplier() -- the numerator-only route.
use "panel_cwm02.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it gteffects customweightmultiplier(epop)
dis "K2_RC=" _rc
matrix G2 = e(gteffects)
preserve
clear
svmat double G2, names(col)
qui sum w_att_it
dis "K2_w_att_it  min=" %20.17f r(min) "  max=" %20.17f r(max) "   (cwm route: does NOT cancel)"
restore

*--- K3: the gate itself -- bootstrap(bayes) + a cell-weighted estimand + a cwm must refuse.
use "panel_cwm02.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it cluster(county) bootstrap(bayes) reps(20) customweightmultiplier(pop)
dis "K3_RC=" _rc

*--- K4: bootstrap(bayes) + att only + a cwm must proceed.
use "panel_cwm02.dta", clear
cap noi wooldid y county year onset [aw=sw], att cluster(county) bootstrap(bayes) reps(20) customweightmultiplier(pop)
dis "K4_RC=" _rc

dis "CWM02_DONE"
