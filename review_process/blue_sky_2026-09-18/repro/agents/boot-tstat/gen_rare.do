* gen_rare.do -- one early cohort (state 1) supplies the long post horizons.
clear
set obs 20
gen int state = _n
gen int onset = .
replace onset = 2002 in 1
replace onset = 2005 if inrange(state,2,6)
replace onset = 2006 if inrange(state,7,11)
replace onset = 2007 if inrange(state,12,15)
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand 10
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1.0 + 0.25*(t - onset)
qui replace tau = 0 if !post
gen double eps = rnormal()*0.7
gen double y = a_i + u_i + tshock + tau*post + eps
gen int ttre = onset
