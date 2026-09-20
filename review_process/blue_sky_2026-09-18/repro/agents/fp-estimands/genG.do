* genG.do -- Design G: Poisson counts with dispersed unit levels, for the
* poisson / poisexpresults engine-split contrast.
clear
set obs 24
gen int i = _n
gen double ai = 1.0 + 2.0*runiform()
gen int ttre = .
replace ttre = 4 if i <= 8
replace ttre = 6 if i > 8 & i <= 16
expand 8
bysort i: gen int t = _n
gen double mu = exp(ai + 0.05*t + 0.3*((t >= ttre) & ttre < .))
gen int y = rpoisson(mu)
keep i t ttre y
order i t ttre y
sort i t
