*======================================================================
* d04 -- the permutation p and the permutation INTERVAL disagree.
* On MC draw 8 of my G=10 sharp-null series, wooldid reports
*     e(perm_p_att)  = .025      (the observed statistic beats ALL 79 draws)
*     e(perm_lb_att) = -1.2956 , e(perm_ub_att) = 1.1604
* i.e. a 95% interval that contains 0 with enormous margin, while the p
* rejects at 2.5%.  sthlp:512 calls the default equal-tailed p "the exact
* dual of the equal-tailed shift interval".
* This file rebuilds that draw and prints every ingredient.
* Self-contained.  stata-mp -e do d04_perm_dual.do
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

* --- exactly the generator mc_main.do used, at its draw-8 seed ---
set seed 12010008
set sortseed 12010008
clear
set obs 10
gen long st = _n
gen double u = runiform()
sort u
gen long rk = _n
gen double ttre = .
replace ttre = 2004 if rk <= 2
replace ttre = 2005 if rk >  2 & rk <= 4
replace ttre = 2006 if rk >  4 & rk <= 6
gen double ai = rnormal()
gen double esd = 1
drop u rk
expand 8
bysort st: gen int t = 1999 + _n
sort t st
by t: gen double ft = 0.5*rnormal() if _n==1
by t: replace ft = ft[1]
sort st t
by st: gen double e = rnormal()*esd if _n==1
by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(1-0.25)*esd if _n>1
gen double y = ai + ft + e
keep st t ttre y
order st t ttre y
list in 1/8
tab ttre, missing

ereturn clear
cap noi wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(700066)
dis "RC=" _rc
dis "p_att   = " e(perm_p_att)    "   lb=" e(perm_lb_att)    " ub=" e(perm_ub_att)
dis "p_att_i = " e(perm_p_att_i)  "   lb=" e(perm_lb_att_i)  " ub=" e(perm_ub_att_i)
dis "bandlo_att=" e(perm_bandlo_att) " bandhi_att=" e(perm_bandhi_att) " sd_att=" e(perm_sd_att)
dis "completed=" e(perm_completed) " failed=" e(perm_failed) " blank=" e(perm_blank)

mat CB = e(perm_comparator_b)
mat OB = e(b)
dis "comparator_b:"
mat list CB
dis "e(b):"
mat list OB

mat D = e(permdist_b)
dis "permdist_b is " rowsof(D) " x " colsof(D)
local cn : colnames D
dis "cols: `cn'"

* --- rebuild the p and the interval BY HAND from the posted draws ---
clear
qui svmat double D, names(dd)
local K = colsof(D)
qui gen long _r = _n
forvalues k = 1/`K' {
  local nm : word `k' of `cn'
  qui count if dd`k' < .
  local n = r(N)
  qui sum dd`k', detail
  local sd = r(sd)
  qui count if dd`k' <= CB[1,`k'] & dd`k' < .
  local m0 = r(N)
  qui count if dd`k' >= CB[1,`k'] & dd`k' < .
  local m1 = r(N)
  local phand = min(1, 2*min((1+`m0')/(`n'+1), (1+`m1')/(`n'+1)))
  * 2.5 / 97.5 order statistics of the draws
  qui _pctile dd`k' if dd`k' < ., percentiles(2.5 97.5)
  local qlo = r(r1)
  local qhi = r(r2)
  dis "col `k' (`nm'): n=`n' comparator=" CB[1,`k'] " sd=`sd'"
  dis "   m0=`m0' m1=`m1'   p_by_hand=`phand'"
  dis "   draw q2.5=`qlo'  q97.5=`qhi'   min=" %9.4f r(min)
  qui sum dd`k'
  dis "   draw min=" %9.4f r(min) "  max=" %9.4f r(max) "  mean=" %9.4f r(mean)
  dis "   shift interval by hand = [" %9.4f (CB[1,`k'] - `qhi') "," %9.4f (CB[1,`k'] - `qlo') "]"
}
