clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---- Shape B: state-year panel, sharp null (NO treatment effect at all) ----
* G=10 states, T=8, onsets staggered, 4 never-treated states.
local G  10
local T  8
clear
set obs `G'
gen long st = _n
* onset assignment: states 1-2 -> 2004, 3-4 -> 2005, 5-6 -> 2006, 7-10 never
gen double ttre = .
replace ttre = 2004 if inrange(st,1,2)
replace ttre = 2005 if inrange(st,3,4)
replace ttre = 2006 if inrange(st,5,6)
gen double ai = rnormal()
expand `T'
bys st: gen int t = 2000 + _n - 1
bys t: gen double ft = 0
* period shock: same across states within t
tempvar junk
gen double e = rnormal()
* add a common time effect
sort t st
by t: replace ft = rnormal() if _n==1
by t: replace ft = ft[1]
gen double y = ai + ft + e
compress
list in 1/12

timer clear
timer on 1
qui wooldid y st t ttre, cluster(st)
timer off 1
dis "PLAIN_b=" _b[main_att] " se=" _se[main_att]

timer on 2
cap noi wooldid y st t ttre, cluster(st) permute(shuffle) permreps(79) rseed(123)
dis "PERM_RC=" _rc
timer off 2
dis "perm_p_att=" e(perm_p_att)

timer on 3
cap noi wooldid y st t ttre, cluster(st) bootstrap(standard) reps(99) rseed(123)
dis "BOOT_RC=" _rc
timer off 3
mat lb = e(boot_p_percentile)
mat lc = e(boot_p_bc)
mat li = e(boot_ci_percentile)
mat list lb
mat list lc
mat list li

timer on 4
cap noi wooldid y st t ttre, cluster(st) bootstrap(standard) reps(99) rseed(123) boottstat
dis "BOOTT_RC=" _rc
timer off 4

timer on 5
cap noi wooldid y st t ttre, cluster(st) bootstrap(bayes) reps(99) rseed(123)
dis "BAYES_RC=" _rc
timer off 5

timer list
