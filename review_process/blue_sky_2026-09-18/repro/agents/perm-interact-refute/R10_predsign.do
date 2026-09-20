* R10_predsign.do -- SELF-CONTAINED.
* PART 1 (decides SK-01): in the semielasticity SUCCESS band with min(y) < 0, are the
*   model's fitted values ever non-positive?  If they are all positive, the reported
*   semi-elasticity is a legitimate quantity and SK-01 must be withdrawn.
*   Oracle model = wooldid's own specification: TWFE absorbing i and t, with one dummy per
*   treated (cohort, relative-period) cell and everything untreated in the base cell.
* PART 2: the premise of the rc-0/no-e(b) mechanism at ado:4051-4085 -- does `ereturn post'
*   really refuse a coefficient vector containing a missing value?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

set obs 60
gen int firm = _n
gen double fe_i = rnormal(0, 1.4)
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
gen double ycross = fe_i + fe_t + effect + eps
drop effect
forvalues c = 0/6 {
  gen double yshift`c' = ycross + `c'
}

* wooldid's design for an att-only run: one cell per (onset, relative period >= 0)
gen byte treated = (onset < . & year >= onset)
gen long cellid = 0
qui egen long tmpcell = group(onset rel) if treated == 1
qui replace cellid = tmpcell if treated == 1
drop tmpcell

foreach c in 0 1 2 3 4 5 {
  qui summ yshift`c'
  local mn = r(min)
  qui reghdfe yshift`c' ib0.cellid, absorb(firm year) resid
  qui predict double fit`c', xbd
  qui count if fit`c' <= 0
  local np = r(N)
  qui count if fit`c' <= 0 & treated == 1
  local npt = r(N)
  qui summ fit`c'
  dis "PRED c=`c' miny=`mn' fitmin=" r(min) " fitmax=" r(max) " n_fit_nonpos=`np' n_fit_nonpos_treated=`npt' of " _N
  * how close to zero does the fitted value get, relative to its spread?
  qui summ fit`c', detail
  dis "PRED c=`c' fit_p1=" r(p1) " fit_p50=" r(p50)
}

* cross-check that this specification reproduces wooldid's att on yshift2
qui reghdfe yshift2 ib0.cellid, absorb(firm year)
qui margins , dydx(cellid) nose noestimcheck
dis "XCHECK: wooldid att on yshift2 was 1.02-ish (fitted from the ladder hyperbola)"
qui margins , eydx(cellid) subpop(if treated==1) nose noestimcheck
matrix MM = r(b)
dis "XCHECK eydx cols=" colsof(MM)

*========================= PART 2
clear
set obs 10
gen double z = rnormal()
qui reg z
matrix BB = (1.5, .)
matrix colnames BB = a b
matrix rownames BB = y
matrix VV = I(2)
matrix colnames VV = a b
matrix rownames VV = a b
cap ereturn post BB VV
dis "POST_b_with_missing_and_V_rc=" _rc
cap ereturn post BB
dis "POST_b_with_missing_only_rc=" _rc
cap ereturn post
dis "POST_bare_rc=" _rc
cap confirm matrix e(b)
dis "POST_bare_has_b=" (_rc==0)
dis "R10 done"
