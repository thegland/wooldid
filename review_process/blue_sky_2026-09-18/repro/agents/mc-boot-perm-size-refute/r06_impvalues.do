*======================================================================
* r06_impvalues -- SELF-CONTAINED.  A path NOBODY in this wave measured:
* impvalues (Ibragimov-Muller), the textbook FEW-CLUSTER remedy, at G=10 and
* G=6 with cluster()==i under a TRUE SHARP NULL.  Compared against the default
* analytic att on the identical draws.  300 draws, analytic only (fast).
* Run:  stata-mp -e do r06_impvalues.do          (~12 min)
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
postfile `P' int draw str6 arm double b double se double p int rc using r06_out.dta, replace
forvalues d = 1/300 {
  local sx = 910000 + `d'
  mkpan , g(10) t(10) seedx(`sx') kper(2) ntrtc(3) nsub(1)
  tempfile d10
  qui save "`d10'", replace
  mkpan , g(6) t(10) seedx(`=`sx'+400000') kper(1) ntrtc(3) nsub(1)
  tempfile d6
  qui save "`d6'", replace
  foreach a in G10AN G10IM G6AN G6IM {
    local f = cond(substr("`a'",1,3)=="G10", "`d10'", "`d6'")
    qui use "`f'", clear
    ereturn clear
    if "`a'"=="G10AN" cap noi qui wooldid y st t ttre, cluster(st) att
    if "`a'"=="G10IM" cap noi qui wooldid y st t ttre, cluster(st) att impvalues
    if "`a'"=="G6AN"  cap noi qui wooldid y st t ttre, cluster(st) att
    if "`a'"=="G6IM"  cap noi qui wooldid y st t ttre, cluster(st) att impvalues
    local rc = _rc
    qui getrow , rown(att)
    post `P' (`d') ("`a'") (r(b)) (r(se)) (r(p)) (`rc')
  }
}
postclose `P'
use r06_out.dta, clear
gen byte rej = (p <= .05+1e-9) if p<.
dis "=== R06 impvalues (Ibragimov-Muller) size under a TRUE SHARP NULL, cluster()==i ==="
foreach a in G10AN G10IM G6AN G6IM {
  qui sum b if arm=="`a'"
  local sd = r(sd)
  local n = r(N)
  qui sum se if arm=="`a'"
  local ms = r(mean)
  qui sum rej if arm=="`a'"
  dis "R06 `a': n=`n' sd(b)=" %7.4f `sd' " mse=" %7.4f `ms' " se/sd=" %6.3f `=`ms'/`sd'' " rej@5%=" %6.3f r(mean) " (nrej=" r(N) ")"
}
foreach a in G10IM G6IM {
  dis "--- p ECDF `a' ---"
  foreach x in .05 .10 .20 .50 {
    qui count if arm=="`a'" & p <= `x'+1e-9 & p<.
    local c1 = r(N)
    qui count if arm=="`a'" & p<.
    dis "  `a' P(p<=`x') = " %6.3f `=`c1'/r(N)'
  }
}
count if rc != 0
dis "R06 NONZERO_RC=" r(N)
tab arm rc
dis "R06_DONE"
