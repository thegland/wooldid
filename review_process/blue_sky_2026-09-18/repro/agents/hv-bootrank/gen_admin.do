*=======================================================================
* gen_admin.do -- a REALISTIC US-state-style staggered adoption panel.
*   30 states (the cluster), 6 counties each (the unit i), 2005-2019.
*   Adoption is genuinely staggered: six single-state early adopters,
*   three large late cohorts, ten never-treated states.  Four census
*   regions cross the adoption timing.  No pre-trends: EVERY pre-period
*   event-study coefficient has a true value of exactly 0.
* Out: state county t ttre reg y truetau post h
*=======================================================================
clear
set seed 4471
set sortseed 4471

local NS 30
local NC 6
local T0 2005
local T1 2019

qui set obs `NS'
gen int state = _n

* onsets: 6 single-state early adopters, then three big cohorts, then never
gen double ttre = .
qui replace ttre = 2008 in 1
qui replace ttre = 2009 in 2
qui replace ttre = 2010 in 3
qui replace ttre = 2011 in 4
qui replace ttre = 2012 in 5
qui replace ttre = 2013 in 6
qui replace ttre = 2016 in 7/11
qui replace ttre = 2017 in 12/16
qui replace ttre = 2018 in 17/20
* 21..30 never treated

* four regions, crossing adoption timing (state 1,5,9,... in region 1 etc.)
gen byte reg = mod(state - 1, 4) + 1

gen double a_state = rnormal()
expand `NC'
bysort state: gen int cnum = _n
gen long county = state*100 + cnum
gen double a_unit = a_state + rnormal()*0.7

expand `=`T1'-`T0'+1'
bysort county: gen int t = `T0' + _n - 1

tempname fm
mata: st_matrix("`fm'", rnormal(`=`T1'-`T0'+1', 1, 0, .5))
gen double f_t = .
forvalues j = 1/`=`T1'-`T0'+1' {
    qui replace f_t = `fm'[`j',1] if t == `T0' + `j' - 1
}

gen int h = t - ttre
gen byte post = (h >= 0) & !missing(h)
* a dynamic, heterogeneous effect - the realistic case
gen double truetau = cond(post == 1, 0.4 + 0.10*h + 0.15*reg, 0)
gen double y = a_unit + f_t + truetau + rnormal(0, 0.45)

keep state county t ttre reg y truetau post h
order state county t ttre reg y truetau post h
sort county t
qui compress
qui count
dis as text "gen_admin: rows=" as result r(N) as text " counties=`=`NS'*`NC'' states=`NS' T=`T0'-`T1'"
