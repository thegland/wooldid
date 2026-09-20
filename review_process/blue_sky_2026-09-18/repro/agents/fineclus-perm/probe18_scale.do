* probe18: the headline use case at realistic scale, plus the opt-in multiplicity surfaces
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"
set seed 5150
set sortseed 5150

* 40 states, 5 adoption cohorts (4 treated + never), 12 years, ragged caseloads
clear
set obs 40
gen int state = _n
gen double sfx = rnormal()*2
gen int nper = 1 + mod(state,4)
expand 12
bysort state: gen int t = _n
expand nper
bysort state t: gen int unit = _n
gen double ufx = rnormal()*0.5
gen double eps = rnormal()
gen int coh = 1 + mod(state,5)
gen double onset = cond(coh==1,4,cond(coh==2,6,cond(coh==3,8,cond(coh==4,10,.))))
gen byte D = (onset<. & t>=onset)
gen double y = 2 + sfx + ufx + 0.2*t + 0.8*D + eps
save "dataD.dta", replace
tab coh

* analytic clustered run (cluster at state, finer than i)
use "dataD.dta", clear
cap noi wooldid y coh t onset, att esprelength(3) espostlength(3) cluster(state)
dis "D_ANALYTIC rc=" _rc " b=" _b[main_att] " se=" _se[main_att]

timer clear 1
timer on 1
use "dataD.dta", clear
cap noi wooldid y coh t onset, att esprelength(3) espostlength(3) cluster(state) permute(shuffle) permreps(199) rseed(2026) permonsetcohorts dots(50)
timer off 1
qui timer list 1
dis "D_PERM rc=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att) " b=" _b[main_att] " secs=" r(t1)
cap mat PN = e(perm_n)
cap mat list PN
cap mat PP = e(perm_p)
cap mat list PP

* multiplicity surfaces under the rebuild
use "dataD.dta", clear
cap noi wooldid y coh t onset, att esprelength(3) espostlength(3) cluster(state) permute(shuffle) permreps(99) rseed(2026) permonsetcohorts permstepdown simultaneousbands dots(50)
dis "D_SD rc=" _rc " done=" e(perm_completed) " fail=" e(perm_failed)
cap mat SD = e(perm_p_stepdown)
cap mat list SD
ereturn list
