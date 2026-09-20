* genA.do -- Design A: unit-level, unbalanced, 3 cohorts, 8 periods
clear
set obs 30
gen int i = _n
gen double ui = rnormal()*2
gen int ttre = .
replace ttre = 4 if i <= 8
replace ttre = 6 if i > 8 & i <= 16
expand 8
bysort i: gen int t = _n
gen double ut = 0.3*t + 0.05*t*t
gen double relyr = t - ttre
gen double tau = 0
replace tau = 1 + 0.4*relyr if relyr >= 0 & relyr < .
gen double y = ui + ut + tau + rnormal()*0.5
* unbalance: drop scattered (i,t) rows, but never a unit's whole record
drop if inlist(i,3,11,19) & t == 5
drop if inlist(i,5,14,22) & t == 2
drop if inlist(i,7,16,27) & t == 7
drop if i == 25 & inlist(t,3,4)
* sample weights, varying within unit across t
gen double w = 0.5 + runiform()*2
* custom weight multiplier, varying by (i,t)
gen double cwm = 0.25 + runiform()*3
keep i t ttre y w cwm
order i t ttre y w cwm
sort i t
