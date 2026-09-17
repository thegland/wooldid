set linesize 255
clear all
do "./w_timer.ado"
clear
set obs 500
gen long unit = _n
gen coh = ceil(unit/63)
gen st  = ceil(unit/10)
gen ttre = cond(coh==8, ., 4 + coh)
expand 20
bysort unit: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
tempfile p
save `p'
qui count
local N = r(N)
use `p', clear
cap qui wooldid y coh t ttre, att cluster(st)

* --- ANALYTIC: one worker entry
timer clear
use `p', clear
timer on 1
cap qui wooldid y coh t ttre, att cluster(st)
timer off 1
qui timer list
di "DIRECT analytic N=`N' total=" r(t1) " F13=" r(t91) " si8=" r(t92)

* --- BOOTSTRAP reps(50): 51 worker entries
timer clear
use `p', clear
timer on 1
cap noi qui wooldid y coh t ttre, att_it cluster(st) bootstrap(standard) reps(50) rseed(4242)
timer off 1
qui timer list
di "DIRECT boot50 N=`N' total=" r(t1) " F13=" r(t91) " si8=" r(t92)

* --- BOOTSTRAP reps(150)
timer clear
use `p', clear
timer on 1
cap noi qui wooldid y coh t ttre, att_it cluster(st) bootstrap(standard) reps(150) rseed(4242)
timer off 1
qui timer list
di "DIRECT boot150 N=`N' total=" r(t1) " F13=" r(t91) " si8=" r(t92)
