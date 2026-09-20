* RF03_semielasticity.do
* SELF-CONTAINED repro for the skeptic's PARTIAL verdict on PX-03.
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do RF03_semielasticity.do
*
* SHOWS TWO THINGS:
*  (A) the SUCCESS band. With min(y) = -2.42 wooldid reports a semi-elasticity with se, t,
*      p and a 95% interval and NO caution -- while (section C) 14 of the 720 FITTED values
*      are non-positive, three of them inside the treated-and-post subpopulation the att is
*      averaged over. eydx divides by the prediction, so those three rows contribute ratios
*      of the wrong sign to the reported "% change in y".
*      -> contradicts PX-03's stated trigger ("non-positive outcome fails"), and is the
*         separate finding SK-01.
*  (B) the FAILURE band. Push the same outcome two units lower and margins dies with a
*      bare, undiagnosed message; wooldid never names the condition in its own words and
*      leaves e(cmd) as reghdfe with tempvar coefficient names.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
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
drop rel effect
forvalues c = 0/6 {
  gen double yshift`c' = ycross + `c'
}

dis "===== (A) SUCCESS BAND: a reported semi-elasticity on a sign-crossing outcome ====="
forvalues c = 2/4 {
  qui summ yshift`c'
  dis "--- shift `c': min(y)=" r(min) " max(y)=" r(max) " share non-positive:"
  qui count if yshift`c' <= 0
  dis "    n_nonpositive=" r(N) " of " _N
  ereturn clear
  cap noi wooldid yshift`c' firm year onset, att cluster(firm) semielasticity
  dis "SHIFT`c'_RC=" _rc
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    matrix list WR
  }
}

dis "===== (B) FAILURE BAND: the same outcome two units lower ====="
forvalues c = 0/1 {
  qui summ yshift`c'
  dis "--- shift `c': min(y)=" r(min)
  ereturn clear
  cap noi wooldid yshift`c' firm year onset, att cluster(firm) semielasticity
  dis "SHIFT`c'_RC=" _rc
  dis "SHIFT`c'_ecmd=" "`e(cmd)'"
  ereturn clear
  cap noi wooldid yshift`c' firm year onset, att cluster(firm) semielasticity fastnose
  dis "SHIFT`c'_FN_RC=" _rc
  dis "SHIFT`c'_FN_ecmd=" "`e(cmd)'"
  ereturn clear
  cap noi wooldid yshift`c' firm year onset, att cluster(firm) semielasticity permute(shuffle) permreps(30) rseed(99)
  dis "SHIFT`c'_PERM_RC=" _rc
}
dis "===== (C) the fitted values behind (A): wooldid's own att-only specification ====="
* one dummy per treated (cohort, relative period) cell, everything untreated in the base
* cell, absorbing i and t -- the model wooldid fits for an att-only run.
cap drop rel
gen int rel = year - onset
gen byte treated = (onset < . & year >= onset)
gen long cellid = 0
qui egen long tmpcell = group(onset rel) if treated == 1
qui replace cellid = tmpcell if treated == 1
drop tmpcell
forvalues c = 0/5 {
  qui reghdfe yshift`c' ib0.cellid, absorb(firm year) resid
  cap drop fit`c'
  qui predict double fit`c', xbd
  qui summ yshift`c'
  local mn = r(min)
  qui count if fit`c' <= 0
  local np = r(N)
  qui count if fit`c' <= 0 & treated == 1
  local npt = r(N)
  qui summ fit`c'
  dis "PRED shift=`c' miny=`mn' fitmin=" r(min) " n_fit_nonpos=`np' n_fit_nonpos_in_att_subpop=`npt' of " _N
}

dis "---- RF03 done ----"
