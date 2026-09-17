* usage: stata-mp -b do time_arm.do <A|B> <nunits> <reps>
args arm nunits reps
set linesize 255
clear all
do "./w`arm'.ado"
clear
set obs `nunits'
gen long unit = _n
gen coh = ceil(unit/(`nunits'/8))
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

* warm: neither arm should pay first-call costs
use `p', clear
cap qui wooldid y coh t ttre, att cluster(st)

* ARM 1 - the headline the author asked for: guard share of a TYPICAL analytic run.
* safety(off) disables more than the two blocks under study, so this is an UPPER BOUND
* on guard share, not an exact attribution.
timer clear
use `p', clear
timer on 1
cap qui wooldid y coh t ttre, att cluster(st)
timer off 1
use `p', clear
timer on 2
cap qui wooldid y coh t ttre, att cluster(st) safety(off)
timer off 2

* ARM 2 - the repeated cost: bootstrap, where the blocks run once per replication
use `p', clear
timer on 3
cap noi qui wooldid y coh t ttre, att_it cluster(st) bootstrap(standard) reps(`reps') rseed(4242)
timer off 3
local rc = 0
local rc = _rc
qui timer list
di "TIMING arm=`arm' N=`N' reps=`reps' rc=`rc' analytic_safetyon=" r(t1) " analytic_safetyoff=" r(t2) " boot=" r(t3)
