*======================================================================
* d05 -- WHICH LEVER fixes the analytic att standard error at cluster()==i:
*        the event-study WINDOW, esfixedbaseperiod, or both?
* Same sharp-null data, 25 draws, four plain (analytic) specifications.
* Self-contained.  stata-mp -e do d05_es_lever.do     (~1 minute)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
cap prog drop mkdgp
prog def mkdgp
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `ngrp'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  local nc = round(`ngrp'*0.2)
  gen double ttre = .
  replace ttre = 2004 if rk <= `nc'
  replace ttre = 2005 if rk >    `nc' & rk <= 2*`nc'
  replace ttre = 2006 if rk >  2*`nc' & rk <= 3*`nc'
  gen double ai = rnormal()
  drop u rk
  expand `nper'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75) if _n>1
  gen double y = ai + ft + e
  keep st t ttre y
end
tempname P
postfile `P' int draw str6 spec double b double se double p double lo double hi int rc ///
  using d05_out.dta, replace
forvalues d = 1/25 {
  local sx = 2000000 + 10*100000 + 10*1000 + `d'
  mkdgp , ngrp(10) nper(8) seedx(`sx')
  tempfile dat
  qui save "`dat'", replace
  foreach s in S0 S1 S2 S3 {
    qui use "`dat'", clear
    ereturn clear
    cap mat drop TT
    if "`s'"=="S0" cap noi qui wooldid y st t ttre, cluster(st) att
    if "`s'"=="S1" cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3)
    if "`s'"=="S2" cap noi qui wooldid y st t ttre, cluster(st) att esfixedbaseperiod
    if "`s'"=="S3" cap noi qui wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3) esfixedbaseperiod
    local rc = _rc
    local b = .
    local se = .
    local p = .
    local lo = .
    local hi = .
    if `rc' == 0 {
      cap mat TT = e(wdidmainresults)
      cap confirm matrix TT
      if !_rc {
        local rn : rownames TT
        local re : roweq TT
        local nr : word count `rn'
        local j = 0
        forvalues z = 1/`nr' {
          local nz : word `z' of `rn'
          local ez : word `z' of `re'
          if `j'==0 & "`nz'"=="att" & "`ez'"=="main" local j = `z'
        }
        if `j' > 0 {
          local b = TT[`j',1]
          local se = TT[`j',2]
          local p = TT[`j',4]
          local lo = TT[`j',5]
          local hi = TT[`j',6]
        }
      }
    }
    post `P' (`d') ("`s'") (`b') (`se') (`p') (`lo') (`hi') (`rc')
  }
}
postclose `P'
use d05_out.dta, clear
gen byte rej = (p <= .05 + 1e-9) if p < .
dis "S0 = att only | S1 = +window | S2 = +esfixedbaseperiod only | S3 = +both"
foreach s in S0 S1 S2 S3 {
  qui sum b if spec=="`s'"
  local sd = r(sd)
  local nb = r(N)
  qui sum se if spec=="`s'"
  local mse = r(mean)
  qui sum rej if spec=="`s'"
  dis "`s': n=`nb'  sd(b)=" %7.4f `sd' "  mean se=" %7.4f `mse' "  se/sd=" %6.3f `=`mse'/`sd'' "  rej@5%=" %5.3f r(mean)
}
count if rc != 0
dis "NONZERO_RC=" r(N)
tab spec rc
