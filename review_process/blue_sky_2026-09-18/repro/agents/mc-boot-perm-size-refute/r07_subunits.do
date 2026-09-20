*======================================================================
* r07_subunits -- SELF-CONTAINED.  NEW (not seeded).  sthlp:251 documents the
* analytic att SE at cluster() at-or-above i as "roughly 0.66 to 0.75 of the
* reference".  Does that hold when the (i,t) cell holds SEVERAL observations,
* i.e. microdata with i = state and cluster(state) — the commonest applied shape?
* i = state throughout; cluster(st) == i throughout; sharp null; 250 draws.
* nsub = 1, 3, 10 observations per state-year.  Also esfixedbaseperiod beside it.
* Run:  stata-mp -e do r07_subunits.do            (~12 min)
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
postfile `P' int draw int nsub str3 spec double b double se double p double nn int rc using r07_out.dta, replace
forvalues d = 1/250 {
  local sx = 660000 + `d'
  foreach ns in 1 3 10 {
    mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(`ns')
    tempfile dd
    qui save "`dd'", replace
    foreach s in D0 D2 {
      qui use "`dd'", clear
      ereturn clear
      if "`s'"=="D0" cap noi qui wooldid y st t ttre, cluster(st) att
      if "`s'"=="D2" cap noi qui wooldid y st t ttre, cluster(st) att esfixedbaseperiod
      local rc = _rc
      local nn = .
      if `rc'==0 cap local nn = e(N)
      qui getrow , rown(att)
      post `P' (`d') (`ns') ("`s'") (r(b)) (r(se)) (r(p)) (`nn') (`rc')
    }
  }
  if mod(`d',50)==0 dis "... r07 draw `d'"
}
postclose `P'
use r07_out.dta, clear
gen byte rej = (p <= .05+1e-9) if p<.
dis "=== R07 analytic att SE calibration vs observations per (i,t) cell, cluster(st)==i ==="
foreach ns in 1 3 10 {
  foreach s in D0 D2 {
    qui sum b if nsub==`ns' & spec=="`s'"
    local sd = r(sd)
    local n = r(N)
    qui sum se if nsub==`ns' & spec=="`s'"
    local ms = r(mean)
    qui sum rej if nsub==`ns' & spec=="`s'"
    local rj = r(mean)
    qui sum nn if nsub==`ns' & spec=="`s'"
    dis "R07 nsub=`ns' `s': n=`n' N=" %6.0f r(mean) " sd(b)=" %7.4f `sd' " mse=" %7.4f `ms' " se/sd=" %6.3f `=`ms'/`sd'' " rej@5%=" %6.3f `rj'
  }
}
count if rc != 0
dis "R07 NONZERO_RC=" r(N)
dis "R07_DONE"
