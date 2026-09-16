clear all
do "./w_si8.ado"
clear
set obs 200
gen long unit = _n
gen st = ceil(unit/10)
gen ttre = cond(mod(unit,4)==0, ., 5 + mod(unit,6))
expand 20
bysort unit: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
di "=== the fixture that raised Stata's own r(907) ==="
cap noi wooldid y unit t ttre, att cluster(st)
di "rc=" _rc
