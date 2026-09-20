*==========================================================================
* r2.do — MECHANISM, outside wooldid entirely.
*
* wooldid's collapse hands margins one row per (i,t,targetflag) cell carrying the
* cell's summed weight and the FIRST row's covariate values.  NFC-01 says the
* resulting VARIANCE moves by ~1e-6 relative when a continuous control varies
* inside the cell.  This file reproduces that geometry with plain reghdfe +
* margins, and adds the one thing the prober never computed: the EXACT
* delta-method answer.
*
*   The estimand is an average of coefficients:  AME = c'b  with c = the treated
*   rows' weight shares over the (cohort x relative-period) cells.  Because the
*   model is linear and the control does not interact with the treatment, the
*   exact Jacobian is c and the exact variance is c'Vc — no numerical
*   differentiation needed.  Anything margins reports away from c'Vc is its own
*   numerical-Jacobian error, present with or without wooldid.
*
* ARMS
*   ORACLE     c'Vc, exact
*   FULL       margins on the full micro data          (= wooldid's nofastcollapse)
*   COLL1      margins on one row per cell, FIRST row's age   (= wooldid's default)
*   COLLMEAN   margins on one row per cell, cell-MEAN age     (the candidate remedy)
*   NOAGE      the same three with age dropped from the model (the floor)
*==========================================================================
clear all
set more off
set processors 2
set seed 31337
set sortseed 31337

*------------------------------------------------------------------ data (own)
local NS 30
local NI 25
local T  8
set obs `NS'
gen int state = _n
gen double u_s = rnormal()
gen int ttre = .
replace ttre = 2014 if state <= 10
replace ttre = 2016 if state > 10 & state <= 20
expand `NI'
bysort state: gen int pid = _n
gen double u_i = rnormal()
expand `T'
bysort state pid: gen int year = 2010 + _n
gen double age = 25 + 20*runiform()
gen byte tf = (ttre < . & year >= ttre)
gen double y = 2 + u_s + u_i + 0.05*(year-2010) + 0.4*tf + 0.01*age + 0.5*rnormal()

* (cohort, relative period) cell id for treated rows; 0 for every untreated row
gen int rel = year - ttre
gen int cellid = 0
replace cellid = 100*(ttre-2000) + (rel+10) if tf==1
gen double wv = 1
compress
save "_r2.dta", replace

*------------------------------------------------------------------ helper
* exact delta-method variance of the weight-share average of the cell coefficients
prog def oracle, rclass
  syntax , [NOAGE]
  tempname bb VV cvec av vv2
  matrix `bb' = e(b)
  matrix `VV' = e(V)
  matrix `cvec' = J(1, colsof(`bb'), 0)
  qui sum wv if tf==1, meanonly
  local wtot = r(sum)
  levelsof cellid if tf==1, local(cids)
  foreach k of local cids {
    qui sum wv if tf==1 & cellid==`k', meanonly
    local sh = r(sum)/`wtot'
    local cn = colnumb(`bb', "`k'.cellid#c.tf")
    if `cn' >= . {
      dis as error "ORACLE: coefficient `k'.cellid#c.tf not found"
      exit 111
    }
    matrix `cvec'[1,`cn'] = `sh'
  }
  matrix `av' = `cvec' * `bb''
  matrix `vv2' = `cvec' * `VV' * `cvec''
  return scalar b = `av'[1,1]
  return scalar V = `vv2'[1,1]
end

*------------------------------------------------------------------ run
foreach spec in AGE NOAGE {
  local ctrl c.age
  if "`spec'"=="NOAGE" local ctrl

  use "_r2.dta", clear
  qui reghdfe y c.tf#ib0.cellid `ctrl', absorb(state year)
  estimates store MOD
  oracle
  scalar oB = r(b)
  scalar oV = r(V)

  * FULL: margins over every row, as wooldid does under nofastcollapse
  qui margins [aw=wv], dydx(tf) subpop(if tf==1) noestimcheck
  matrix rB = r(b)
  matrix rV = r(V)
  scalar fB = rB[1,1]
  scalar fV = rV[1,1]

  * COLL1: one row per (state, year, tf) with summed weight and the FIRST row's age
  use "_r2.dta", clear
  sort state year tf pid
  by state year tf: gen long _ord = _n
  by state year tf: egen double _wsum = total(wv)
  keep if _ord==1
  replace wv = _wsum
  qui estimates restore MOD
  qui margins [aw=wv], dydx(tf) subpop(if tf==1) noestimcheck noesample
  matrix rB = r(b)
  matrix rV = r(V)
  scalar c1B = rB[1,1]
  scalar c1V = rV[1,1]

  * COLLMEAN: identical, except the surviving row carries the cell MEAN of age
  use "_r2.dta", clear
  sort state year tf pid
  by state year tf: gen long _ord = _n
  by state year tf: egen double _wsum = total(wv)
  by state year tf: egen double _agem = mean(age)
  replace age = _agem
  keep if _ord==1
  replace wv = _wsum
  qui estimates restore MOD
  qui margins [aw=wv], dydx(tf) subpop(if tf==1) noestimcheck noesample
  matrix rB = r(b)
  matrix rV = r(V)
  scalar cmB = rB[1,1]
  scalar cmV = rV[1,1]

  dis _n as text "{hline 88}"
  dis as text "SPEC = `spec'    (oracle = exact c'Vc from the same reghdfe fit)"
  dis as text "{hline 88}"
  dis as result "  oracle   V = " %20.18g oV "   se = " %20.18g sqrt(oV)
  dis as result "  FULL     V = " %20.18g fV "   rel dev vs oracle = " %9.2e abs(fV-oV)/abs(oV)
  dis as result "  COLL1    V = " %20.18g c1V "   rel dev vs oracle = " %9.2e abs(c1V-oV)/abs(oV)
  dis as result "  COLLMEAN V = " %20.18g cmV "   rel dev vs oracle = " %9.2e abs(cmV-oV)/abs(oV)
  dis as result "  FULL vs COLL1     rel = " %9.2e abs(fV-c1V)/abs(fV)
  dis as result "  FULL vs COLLMEAN  rel = " %9.2e abs(fV-cmV)/abs(fV)
  dis as result "  b:  oracle " %20.18g oB "  FULL " %20.18g fB "  COLL1 " %20.18g c1B
  dis as result "  b rel dev: FULL " %9.2e abs(fB-oB)/abs(oB) "   COLL1 " %9.2e abs(c1B-oB)/abs(oB)
}

*---- does the deviation track the SIZE of the control's contribution to xb?
* same data, age replaced by k*age for k = 1, 10, 100 (so beta*x is unchanged in
* the fit but |xb| is not, since beta rescales exactly) -- and separately a control
* whose COEFFICIENT is large.
dis _n as text "{hline 88}"
dis as text "SCALE SWEEP: control contribution to the linear prediction"
dis as text "{hline 88}"
foreach mult in 0.01 1 100 10000 {
  use "_r2.dta", clear
  gen double xs = age*`mult'
  * make the DGP contribution of the control large or small by rescaling y's own
  * dependence: refit y on the scaled control so beta*x ~ 0.35*`mult' in size
  replace y = y + 0.01*xs
  qui reghdfe y c.tf#ib0.cellid c.xs, absorb(state year)
  estimates store M2
  qui predict double _xb, xb
  qui sum _xb, meanonly
  local xbbar = r(mean)
  oracle
  scalar oV = r(V)
  qui margins [aw=wv], dydx(tf) subpop(if tf==1) noestimcheck
  matrix rV = r(V)
  scalar fV = rV[1,1]
  use "_r2.dta", clear
  gen double xs = age*`mult'
  replace y = y + 0.01*xs
  sort state year tf pid
  by state year tf: gen long _ord = _n
  by state year tf: egen double _wsum = total(wv)
  keep if _ord==1
  replace wv = _wsum
  qui estimates restore M2
  qui margins [aw=wv], dydx(tf) subpop(if tf==1) noestimcheck noesample
  matrix rV = r(V)
  scalar cV = rV[1,1]
  dis as result "  mult=" %9.0g `mult' "  mean xb=" %12.4g `xbbar' ///
     "  FULLvsORACLE=" %9.2e abs(fV-oV)/oV "  COLLvsORACLE=" %9.2e abs(cV-oV)/oV ///
     "  FULLvsCOLL=" %9.2e abs(fV-cV)/abs(fV)
}

dis "R2_DONE"
