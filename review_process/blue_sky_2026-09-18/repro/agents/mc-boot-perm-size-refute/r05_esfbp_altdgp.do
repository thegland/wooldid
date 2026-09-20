*======================================================================
* r05_esfbp_altdgp -- SELF-CONTAINED.  MCBP-ESFBP said esfixedbaseperiod moves
* the analytic att se/sd from ~0.71 to 0.93-1.01 and rejection from .24 to .08
* at G=10 with cluster()==i.  That was measured on 25 draws.  Here the SAME
* question on a panel whose ONSET STRUCTURE matches theirs (T=8, three
* CONSECUTIVE onset years, 2 states each, AR(1) rho .5) at 300 draws, plus a
* second panel shape (T=12, 4 cohorts) to see whether any "fix" generalizes.
* Run:  stata-mp -e do r05_esfbp_altdgp.do          (~12 min)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkcons
prog def mkcons
  syntax , G(integer) T(integer) SEEDX(integer) [KPER(integer 2) RHO(real 0.5) NTRTC(integer 3) FIRST(integer 2004)]
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
    local on = `first' + `c' - 1
    replace ttre = `on' if rk >= `lo' & rk <= `hi'
  }
  gen double ai = rnormal()
  drop u rk
  expand `t'
  bysort st: gen int t = (`first' - (`t' - `ntrtc')) + _n - 1
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal() if _n==1
  by st: replace e = `rho'*e[_n-1] + rnormal()*sqrt(1-`rho'^2) if _n>1
  gen double y = ai + ft + e
  keep st t ttre y
  order st t ttre y
end

cap prog drop getatt
prog def getatt, rclass
  return scalar b  = .
  return scalar se = .
  return scalar p  = .
  cap mat TT = e(wdidmainresults)
  cap confirm matrix TT
  if _rc exit
  local rn : rownames TT
  local re : roweq TT
  local nr : word count `rn'
  forvalues z = 1/`nr' {
    local nz : word `z' of `rn'
    local ez : word `z' of `re'
    if "`nz'"=="att" & "`ez'"=="main" {
      return scalar b  = TT[`z',1]
      return scalar se = TT[`z',2]
      return scalar p  = TT[`z',4]
      continue, break
    }
  }
end

tempname P
postfile `P' int draw str4 shape str4 spec double b double se double p int rc using r05_out.dta, replace
forvalues d = 1/300 {
  local sx = 810000 + `d'
  * shape C8 : T=8, onsets 2004 2005 2006, 2 states each, 4 never-treated  (their shape)
  mkcons , g(10) t(8) seedx(`sx') kper(2) ntrtc(3) first(2004) rho(0.5)
  tempfile d8
  qui save "`d8'", replace
  * shape C12: T=12, onsets 2004..2007, 2 states each, 12 states, 4 never-treated
  mkcons , g(12) t(12) seedx(`=`sx'+500000') kper(2) ntrtc(4) first(2004) rho(0.5)
  tempfile d12
  qui save "`d12'", replace
  foreach sh in C8 C12 {
    local f = cond("`sh'"=="C8", "`d8'", "`d12'")
    foreach s in P0 P2 {
      qui use "`f'", clear
      ereturn clear
      if "`s'"=="P0" cap noi qui wooldid y st t ttre, cluster(st) att
      if "`s'"=="P2" cap noi qui wooldid y st t ttre, cluster(st) att esfixedbaseperiod
      local rc = _rc
      qui getatt
      post `P' (`d') ("`sh'") ("`s'") (r(b)) (r(se)) (r(p)) (`rc')
    }
  }
}
postclose `P'
use r05_out.dta, clear
gen byte rej = (p <= .05+1e-9) if p<.
dis "=== R05 esfixedbaseperiod as an SE 'lever', 300 draws, cluster(st)==i, sharp null ==="
foreach sh in C8 C12 {
  foreach s in P0 P2 {
    qui sum b if shape=="`sh'" & spec=="`s'"
    local sd = r(sd)
    local n = r(N)
    qui sum se if shape=="`sh'" & spec=="`s'"
    local ms = r(mean)
    qui sum rej if shape=="`sh'" & spec=="`s'"
    dis "R05 `sh' `s': n=`n' sd(b)=" %7.4f `sd' " mse=" %7.4f `ms' " se/sd=" %6.3f `=`ms'/`sd'' " rej@5%=" %6.3f r(mean)
  }
}
count if rc != 0
dis "R05 NONZERO_RC=" r(N)
dis "R05_DONE"
