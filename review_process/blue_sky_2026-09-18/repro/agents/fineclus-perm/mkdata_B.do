* design B: RAGGED panel so that finer-than-i FE and controls actually bite.
* 6 states (clusters), i = coh (onset grouping), cluster = state
* coh1 (onset 4) = s1,s2 ; coh2 (onset 6) = s3,s4 ; coh3 (never) = s5,s6
* M = 6!/(2!2!2!) = 90 assignments
clear all
set more off
set seed 909090
set sortseed 909090
set obs 6
gen int state = _n
gen double sfx = rnormal()*2
gen int nper = cond(state==1,2,cond(state==2,3,cond(state==3,1,cond(state==4,2,cond(state==5,3,1)))))
gen int tlo  = cond(state==2,2,cond(state==5,3,1))
gen int thi  = cond(state==3,7,8)
expand 8
bysort state: gen int t = _n
drop if t < tlo | t > thi
expand nper
bysort state t: gen int unit = _n
gen double ufx = rnormal()
gen double xcov = rnormal() + 0.4*t
gen double eps = rnormal()
gen int coh = cond(inlist(state,1,2),1,cond(inlist(state,3,4),2,3))
gen double onset = cond(coh==1,4,cond(coh==2,6,.))
gen byte D = (onset<. & t>=onset)
gen double y = 1 + sfx + ufx + 0.3*t + 0.7*xcov + 1.5*D + eps
gen int cohcopy = coh
gen double cohctrl = 10*coh + 0.5*coh*coh
gen double statectrl = sfx
keep state unit t coh cohcopy onset y D xcov cohctrl statectrl
sort state t unit
save "dataB.dta", replace
tab state coh
