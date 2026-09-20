* R3_semi_oracle.do -- SELF-CONTAINED. What does wooldid's semielasticity actually divide by,
* and is "the outcome is non-positive" the right trigger?
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
gen double size = 3 + rnormal(0, .8)
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
gen double ycross = fe_i + fe_t + effect + eps
gen double y12  = ycross + 12
gen double y120 = ycross + 120
gen double shock = 0
replace shock = -30 if mod(_n, 197) == 0
gen double ymix = y12 + shock
gen byte treatpost = (onset < . & year >= onset)

*---------------------------------------------------------------- oracle inputs
summ y12 if treatpost == 1
local m12 = r(mean)
summ y120 if treatpost == 1
local m120 = r(mean)
dis "ORACLE_mean_y12_treatpost=" `m12'
dis "ORACLE_mean_y120_treatpost=" `m120'

cap prog drop shw
prog def shw
  args tag
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    dis "`tag'_att=" WR[1,1]
  }
  else dis "`tag'_att=NOMATRIX"
  cap confirm matrix e(b)
  dis "`tag'_hasb=" (_rc==0)
  dis "`tag'_ecmd=" "`e(cmd)'"
end

*=========== 1. dydx (plain) on y12 and y120 -- level-invariant by construction
ereturn clear
cap noi wooldid y12 firm year onset, att cluster(firm)
dis "D12_RC=" _rc
shw D12
ereturn clear
cap noi wooldid y120 firm year onset, att cluster(firm)
dis "D120_RC=" _rc
shw D120

*=========== 2. semielasticity on y12 and y120
ereturn clear
cap noi wooldid y12 firm year onset, att cluster(firm) semielasticity
dis "S12_RC=" _rc
shw S12
ereturn clear
cap noi wooldid y120 firm year onset, att cluster(firm) semielasticity
dis "S120_RC=" _rc
shw S120

*=========== 3. fitted-value signs: is the PREDICTION ever non-positive?
* refit the analogous absorbing model by hand and look at both counterfactual predictions
foreach v in ycross y12 ymix {
  qui reghdfe `v' i.treatpost, absorb(firm year) 
  qui predict double fit_`v', xbd
  qui predict double xbonly_`v', xb
  qui count if fit_`v' <= 0
  dis "FIT_`v'_nonpos_xbd=" r(N)
  qui count if xbonly_`v' <= 0
  dis "FIT_`v'_nonpos_xb=" r(N)
  qui summ fit_`v'
  dis "FIT_`v'_xbd_min=" r(min) " max=" r(max)
  qui summ xbonly_`v'
  dis "FIT_`v'_xb_min=" r(min) " max=" r(max)
}

*=========== 4. does margins on a reghdfe fit divide by xb or xbd? direct measurement
qui reghdfe y12 i.treatpost, absorb(firm year)
qui margins , dydx(treatpost) subpop(if treatpost==1) noestimcheck nose
matrix MD = r(b)
dis "MARG_dydx=" MD[1,2]
qui margins , eydx(treatpost) subpop(if treatpost==1) noestimcheck nose
matrix ME = r(b)
dis "MARG_eydx=" ME[1,2]
qui summ fit_y12 if treatpost==1
dis "MARG_meanfit_xbd=" r(mean)
qui summ xbonly_y12 if treatpost==1
dis "MARG_meanfit_xb=" r(mean)
dis "MARG_dydx_over_meanxbd=" MD[1,2]/`m12'
dis "---- R3 done ----"
