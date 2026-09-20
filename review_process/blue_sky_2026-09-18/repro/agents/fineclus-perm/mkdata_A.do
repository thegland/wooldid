* mkdata_A.do -- design A: 5 states, i = onset-cohort, cluster = state
* cohort 1 (onset 4): states 1,2 ; cohort 2 (onset 6): state 3 ; cohort 3 (NT): states 4,5
clear all
set more off
set seed 424242
set sortseed 424242
set obs 5
gen int state = _n
gen double sfx = rnormal()*2
expand 8
bysort state: gen int t = _n
gen double tfx = 0.3*t
gen double eps = rnormal()
gen int coh = cond(inlist(state,1,2),1,cond(state==3,2,3))
gen double onset = cond(coh==1,4,cond(coh==2,6,.))
gen byte D = (onset<. & t>=onset)
gen double y = 1 + sfx + tfx + 1.5*D + eps
keep state t coh onset y D
order state coh onset t y D
sort state t
save "dataA.dta", replace
