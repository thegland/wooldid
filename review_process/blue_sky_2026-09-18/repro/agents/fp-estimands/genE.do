* genE.do -- Design E: strictly-positive outcome with WIDELY DISPERSED unit levels
* (the shape of firm revenue / county spending / hospital caseload), for probing
* semielasticity, whose reported quantity divides the effect by a prediction.
clear
set obs 24
gen int i = _n
gen double ui = 20 + 200*runiform()
gen int ttre = .
replace ttre = 4 if i <= 8
replace ttre = 6 if i > 8 & i <= 16
expand 8
bysort i: gen int t = _n
gen double relyr = t - ttre
gen double tau = 0
replace tau = 5 if relyr >= 0 & relyr < .
gen double y = ui + 2*t + tau + rnormal()*1.5
gen double w = 0.5 + runiform()*2
keep i t ttre y w
order i t ttre y w
sort i t
