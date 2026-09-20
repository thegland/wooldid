*======================================================================
* r08_boottstat -- SELF-CONTAINED.  MCBP-BOOTTSTAT claims boottstat at G=10 is
* correctly sized (.025) but buys it with 95% intervals whose median width is
* 2.7x, p90 7x and max 22x what a calibrated interval needs.  Independent test
* on my own data: 40 draws, cluster(st)==i, sharp null, reps(99).
* The comparator is computed IN THIS FILE from the same 40 draws: the MC sd of
* the att estimates times 2*invttail(9,.025).
* Run:  stata-mp -e do r08_boottstat.do           (~20 min)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpan
prog def mkpan
  syntax , G(integer) T(integer) SEEDX(integer) [KPER(integer 2) RHO(real 0.4) NTRTC(integer 3) NSUB(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `g'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  gen double ttre = .
  forvalues c = 1/`ntrtc' {
    local lo = (`c'-1)*`kper' + 1
    local hi = `c'*`kper'
    local on = 2003 + 2*`c'
    replace ttre = `on' if rk >= `lo' & rk <= `hi'
  }
  gen double ai = rnormal()*1.2
  drop u rk
  expand `t'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.6*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal() if _n==1
  by st: replace e = `rho'*e[_n-1] + rnormal()*sqrt(1-`rho'^2) if _n>1
  if `nsub' > 1 {
    expand `nsub'
    bysort st t: gen int sub = _n
    gen double v = rnormal()*0.8
  }
  else gen double v = 0
  gen double y = ai + ft + e + v
  keep st t ttre y
  order st t ttre y
end

* pull one main-table row out of e(wdidmainresults): b se stat p lo hi
cap prog drop getrow
prog def getrow, rclass
  syntax , ROWN(string)
  return scalar b  = .
  return scalar se = .
  return scalar p  = .
  return scalar lo = .
  return scalar hi = .
  cap mat TT = e(wdidmainresults)
  cap confirm matrix TT
  if _rc exit
  local rn : rownames TT
  local re : roweq TT
  local nr : word count `rn'
  forvalues z = 1/`nr' {
    local nz : word `z' of `rn'
    local ez : word `z' of `re'
    if "`nz'"=="`rown'" & "`ez'"=="main" {
      return scalar b  = TT[`z',1]
      return scalar se = TT[`z',2]
      return scalar p  = TT[`z',4]
      return scalar lo = TT[`z',5]
      return scalar hi = TT[`z',6]
      continue, break
    }
  }
end

tempname P
postfile `P' int draw str6 arm double b double se double p double lo double hi double ncomp int rc using r08_out.dta, replace
forvalues d = 1/40 {
  local sx = 770000 + `d'
  local rs = 660000 + 23*`d'
  mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(1)
  tempfile dd
  qui save "`dd'", replace
  foreach a in TSTAT PLAIN {
    qui use "`dd'", clear
    ereturn clear
    if "`a'"=="TSTAT" cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(standard) boottstat reps(99) rseed(`rs')
    if "`a'"=="PLAIN" cap noi qui wooldid y st t ttre, cluster(st) att
    local rc = _rc
    local nc = .
    if `rc'==0 cap local nc = e(boot_completed)
    qui getrow , rown(att)
    post `P' (`d') ("`a'") (r(b)) (r(se)) (r(p)) (r(lo)) (r(hi)) (`nc') (`rc')
  }
  if mod(`d',10)==0 dis "... r08 draw `d'"
}
postclose `P'
use r08_out.dta, clear
gen byte rej = (p <= .05+1e-9) if p<.
gen byte cov = (lo <= 0 & hi >= 0) if lo<.
gen double wid = hi - lo
qui sum b if arm=="PLAIN"
local truesd = r(sd)
local need = 2*invttail(9,.025)*`truesd'
dis "R08 MC sd of att over these draws = " %8.5f `truesd'
dis "R08 width a calibrated 95% t(9) interval needs = " %8.5f `need'
foreach a in TSTAT PLAIN {
  qui sum rej if arm=="`a'"
  local rj = r(mean)
  local n = r(N)
  qui sum cov if arm=="`a'"
  local cv = r(mean)
  qui sum wid if arm=="`a'", detail
  dis "R08 `a': n=`n' rej@5%=" %6.3f `rj' " cov=" %6.3f `cv' ///
      " width p50=" %8.4f r(p50) " mean=" %8.4f r(mean) " p90=" %8.4f r(p90) " max=" %9.4f r(max)
  dis "R08 `a': p50/need=" %6.2f `=r(p50)/`need'' "  max/need=" %6.2f `=r(max)/`need''
}
dis "--- p ECDF TSTAT ---"
foreach x in .05 .10 .20 .30 {
  qui count if arm=="TSTAT" & p <= `x'+1e-9 & p<.
  local c1 = r(N)
  qui count if arm=="TSTAT" & p<.
  dis "  TSTAT P(p<=`x') = " %6.3f `=`c1'/r(N)'
}
qui sum ncomp if arm=="TSTAT"
dis "R08 completions min=" r(min) " mean=" r(mean)
count if rc != 0
dis "R08 NONZERO_RC=" r(N)
dis "R08_DONE"
