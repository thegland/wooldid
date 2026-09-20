* genI.do -- multi-row-per-cell admin panel WITH a covariate that varies inside the
* (i,t) cell, plus a cluster that crosses i, for the cell-collapse invariance probe.
clear
set obs 16
gen int i = _n
gen double ui = rnormal()*2
gen int ttre = .
replace ttre = 3 if i <= 5
replace ttre = 5 if i > 5 & i <= 10
expand 7
bysort i: gen int t = _n
gen int nrec = 2 + int(runiform()*4)
expand nrec
gen double x  = rnormal()                  // varies WITHIN the (i,t) cell
gen double relyr = t - ttre
gen double tau = 0
replace tau = 0.8 + 0.4*relyr if relyr >= 0 & relyr < .
gen double y = ui + 0.3*t + 0.5*x + tau + rnormal()*0.6
gen double w = 0.3 + runiform()*3
gen byte sg = mod(i,2)
gen int cl = 1 + mod(_n,9)                 // a cluster that CROSSES i
keep i t ttre y w x sg cl
order i t ttre y w x sg cl
sort i t
