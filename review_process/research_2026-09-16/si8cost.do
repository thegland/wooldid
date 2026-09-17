* usage: stata-mp -b do si8cost.do <nunits> <ncoh> <nt>
args nunits ncoh nt
set linesize 255
set maxvar 20000
clear all
do "./w_si8timer.ado"
clear
set obs `nunits'
gen long unit = _n
gen coh = ceil(unit/(`nunits'/`ncoh'))
gen st  = ceil(unit/10)
gen ttre = cond(coh==`ncoh', ., 4 + mod(coh,6))
expand `nt'
bysort unit: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
tempfile p
save `p'
qui count
local N = r(N)
use `p', clear
cap qui wooldid y coh t ttre, att cluster(st)

timer clear
use `p', clear
timer on 1
cap noi qui wooldid y coh t ttre, att cluster(st)
timer off 1
qui timer list
di "SI8COST N=`N' units=`nunits' ncoh=`ncoh' nt=`nt' rc=" _rc " total=" r(t1) " block=" r(t84) " tag_i=" r(t81) " tag_iXt=" r(t82) " tag_t=" r(t83)
