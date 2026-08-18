* mkdata_admin.do — realistic administrative-style panel for the comprehensive suite.
* Usage: stata -e do mkdata_admin.do <outdir>
* Builds dta_s_admin.dta: ~2,000 counties (FIPS-style ids = state*1000+county) x
* years 2005-2019; staggered adoption with UNEQUAL cohorts (2009 ~40%, 2012 ~25%,
* 2015 ~10%, never ~25%); lognormal population weights spanning ~3e2..1e7; ~4% of
* counties enter late and ~4% exit early (unbalanced); ~2% scattered missing y
* (REGULAR missing only — the extended-missing divergence class is deliberately
* exercised by dta_s_extmiss, not here); near-collinear controls x1/x2 with ~1%
* missing x1; a rare-count outcome for poisson; region (4) subgroups; state (48)
* cluster/fe variable. Deterministic: fixed seed AND sortseed.
local outdir `1'
if "`outdir'" == "" local outdir .

clear all
set seed 20260818
set sortseed 20260818

set obs 2000
gen state = 1 + mod(_n - 1, 48)
gen county3 = 1 + floor((_n - 1) / 48)
gen double fips = state * 1000 + county3
gen region = 1 + mod(state, 4)
gen double u = runiform()
gen ttre = cond(u < .40, 2009, cond(u < .65, 2012, cond(u < .75, 2015, .)))
gen double pop = exp(rnormal() * 1.5 + 9.2)
replace pop = 300 if pop < 300
replace pop = 1e7 if pop > 1e7
gen entry = 2005
gen exit = 2019
replace entry = 2005 + ceil(runiform() * 3) if runiform() < .04
* late adopters slightly more likely to enter late (entry correlated w/ timing)
replace entry = 2005 + ceil(runiform() * 3) if ttre == 2015 & runiform() < .06
replace exit = 2019 - ceil(runiform() * 3) if runiform() < .04
gen double ceff = rnormal() * .3

expand 15
bysort fips: gen year = 2004 + _n
drop if year < entry | year > exit

gen double treatnow = (ttre < . & year >= ttre)
gen double relyr = year - ttre if ttre < .
gen double effect = 0
replace effect = .5 + .1 * (ttre == 2012) - .2 * (ttre == 2015) + .05 * min(relyr, 4) if treatnow == 1
gen double y = 10 + .3 * region + ceff + .08 * (year - 2004) + effect + rnormal() * .8
gen double x1 = .4 * ceff + rnormal()
gen double x2 = .98 * x1 + rnormal() * .05
replace x1 = . if runiform() < .01
replace y = . if runiform() < .02
gen double ycount = rpoisson(exp(-1.5 + .2 * treatnow + .1 * region))
drop u entry exit ceff treatnow relyr effect

* ---- deterministic self-test (hard asserts; abort = no sentinel printed) ----
assert inrange(_N, 27000, 30000)
qui count if ttre == 2009
assert inrange(r(N) / _N, .35, .45)
qui count if ttre == .
assert inrange(r(N) / _N, .20, .30)
qui count if missing(y)
assert inrange(r(N) / _N, .01, .03)
qui count if y == .a | y == .b
assert r(N) == 0
assert inrange(fips, 1001, 48999)
qui sum pop
assert r(min) >= 300 & r(max) <= 1e7
qui sum pop, detail
assert r(max) / r(min) > 1000
bysort fips (year): gen __nyr = _N
qui sum __nyr
assert r(min) >= 8 & r(max) == 15
qui count if __nyr < 15
assert r(N) > 0
drop __nyr
bysort fips (ttre): assert ttre[1] == ttre[_N]

save "`outdir'/dta_s_admin.dta", replace
di as result "MKDATA_ADMIN_SELFTEST_PASS"
exit, clear
