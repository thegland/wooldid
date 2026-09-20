*=======================================================================
* gen_skep.do -- MY OWN state-policy panel (skeptic's, not the prober's).
*   36 states (cluster), 4 counties each (rows), 2006-2019 (14 periods).
*   5 single-state early adopters, three multi-state cohorts, 13 never
*   treated.  3 regions cross the timing.  AR(1) county noise, a state
*   random effect, dynamic heterogeneous effects, NO pre-trend.
* Out: state county t ttre reg y
*=======================================================================
clear
set seed 71351
set sortseed 71351
local NS 36
local NC 4
local T0 2006
local T1 2019
local TT = `T1' - `T0' + 1

qui set obs `NS'
gen int state = _n
gen double ttre = .
qui replace ttre = 2009 in 1
qui replace ttre = 2010 in 2
qui replace ttre = 2011 in 3
qui replace ttre = 2012 in 4
qui replace ttre = 2013 in 5
qui replace ttre = 2016 in 6/11
qui replace ttre = 2017 in 12/17
qui replace ttre = 2018 in 18/23
* 24..36 never treated
gen byte reg = mod(state - 1, 3) + 1
gen double a_state = rnormal()*0.8

expand `NC'
bysort state: gen int cnum = _n
gen long county = state*100 + cnum
gen double a_unit = a_state + rnormal()*0.6

expand `TT'
bysort county: gen int t = `T0' + _n - 1
sort county t

* common time shocks
tempname fm
mata: st_matrix("`fm'", rnormal(`TT', 1, 0, .4))
gen double f_t = .
forvalues j = 1/`TT' {
    qui replace f_t = `fm'[`j',1] if t == `T0' + `j' - 1
}
* AR(1) county noise
gen double e = .
by county: replace e = rnormal(0,.5) if _n == 1
by county: replace e = 0.45*e[_n-1] + rnormal(0,.45) if _n > 1

gen int h = t - ttre
gen byte post = (h >= 0) & !missing(h)
gen double truetau = cond(post == 1, 0.35 + 0.07*h + 0.10*(reg - 2), 0)
gen double y = a_unit + f_t + truetau + e

keep state county t ttre reg y truetau post h
sort county t
qui compress
qui count
dis as text "gen_skep: rows=" as result r(N) as text " states=`NS' counties=`=`NS'*`NC'' T=`T0'-`T1'"
