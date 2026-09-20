clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 1
clear
qui set obs 10
gen long ind = _n
qui expand 4
bysort ind: gen byte uw = _n
gen long unit = _n
gen int ttre = .
qui replace ttre = 4 if uw == 1
qui replace ttre = 6 if uw == 2
gen double tau = 0.15*ind + 0.30*uw if !missing(ttre)
qui replace tau = 0 if missing(tau)
qui expand 8
bysort unit: gen int t = _n
gen double y = tau*(t>=ttre & ttre<.) + 0.04*t + rnormal()
wooldid y unit t ttre, att att_i att_it
mat b = e(b)
dis "COLNAMES: " "`: colnames b'"
mat list b
dis "MAINRES:"
mat list e(wdidmainresults)
