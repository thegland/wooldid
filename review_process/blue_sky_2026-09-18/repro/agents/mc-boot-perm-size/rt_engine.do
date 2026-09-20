*=======================================================================
* repro_engine_size.do -- SELF-CONTAINED.  Size of every wooldid inference
* engine under a TRUE SHARP NULL at G = 10 clusters, i = cluster.
*
*   stata-mp -e do repro_engine_size.do          (~50 min at NDRAW 60)
*
* WHAT IT SHOWS.  y is generated with NO treatment effect at any horizon, so
* every estimand's truth is 0 and a 5% test must reject 5% of the time.
* Measured here (see the summary at the end):
*   permute(shuffle)      ~ 5-8%   (correct; randomization inference is EXACT
*                                   under this design, so this is a theorem)
*   analytic cluster se   ~ 20-25% (the documented at-or-above-i understatement)
*   bootstrap(standard)   ~ 15-20%
*   bootstrap(bayes)      ~ 30-40% (WORSE than the analytic se it is offered
*                                   as a remedy for; sthlp:255 names
*                                   "bootstrap() or permute()" as the remedies)
* The point estimate is BITWISE identical across the engines, so the difference
* is entirely in the reported dispersion.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAW 2

* ---------------- sharp-null DGP: 10 states x 8 years, no effect -------------
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
  gen double y = ai + ft + e         /* <-- no treatment term anywhere */
  keep st t ttre y
end

tempname P
postfile `P' int draw str16 arm str6 est double b double se double p double lo double hi int rc ///
  using rt_engine_out.dta, replace

forvalues d = 1/`NDRAW' {
  local sx = 2000000 + 10*100000 + 10*1000 + `d'
  local rs = 700000 + `d'*7 + 10
  mkdgp , ngrp(10) nper(8) seedx(`sx')
  tempfile dat
  qui save "`dat'", replace
  foreach a in analytic permute boot_standard boot_bayes {
    qui use "`dat'", clear
    ereturn clear
    cap mat drop TT
    * bootstrap(bayes) is refused with att_i under default safety, so that arm is att only
    if "`a'" == "analytic"      cap noi qui wooldid y st t ttre, cluster(st) att att_i
    if "`a'" == "permute"       cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs')
    if "`a'" == "boot_standard" cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs')
    if "`a'" == "boot_bayes"    cap noi qui wooldid y st t ttre, cluster(st) att       bootstrap(bayes)    reps(79) rseed(`rs')
    local rc = _rc
    cap mat TT = e(wdidmainresults)
    foreach ee in att att_i {
      local b = .
      local se = .
      local p = .
      local lo = .
      local hi = .
      if `rc' == 0 {
        cap confirm matrix TT
        if !_rc {
          local rn : rownames TT
          local re : roweq TT
          local nr : word count `rn'
          local j = 0
          forvalues z = 1/`nr' {
            local nz : word `z' of `rn'
            local ez : word `z' of `re'
            if `j'==0 & "`nz'"=="`ee'" & "`ez'"=="main" local j = `z'
          }
          if `j' > 0 {
            local b  = TT[`j',1]
            local se = TT[`j',2]
            local p  = TT[`j',4]
            local lo = TT[`j',5]
            local hi = TT[`j',6]
          }
        }
      }
      post `P' (`d') ("`a'") ("`ee'") (`b') (`se') (`p') (`lo') (`hi') (`rc')
    }
  }
  if mod(`d',10)==0 dis "... draw `d' done"
}
postclose `P'

use rt_engine_out.dta, clear
gen byte rej  = (p  <= 0.05 + 1e-9) if p  < .
gen byte excl = (lo >  0 | hi < 0)  if lo < . & hi < .
dis " "
dis "=== SIZE OF EACH ENGINE UNDER A TRUE SHARP NULL, G=10, nominal 5% ==="
table arm est, stat(mean rej) stat(mean excl) stat(mean se) stat(count p)
dis " "
dis "=== reported se vs the TRUE Monte-Carlo sd of the estimator ==="
qui levelsof arm, local(A)
foreach ee in att att_i {
  qui sum b if arm=="analytic" & est=="`ee'"
  local truesd = r(sd)
  dis "`ee': true MC sd of the estimate across draws = " %7.4f `truesd'
  foreach a of local A {
    qui sum se if arm=="`a'" & est=="`ee'"
    if r(N) > 0 dis "   `a': mean reported se " %7.4f r(mean) "   ratio to truth " %6.3f `=r(mean)/`truesd''
  }
}
dis " "
dis "=== point estimates identical across engines? ==="
keep if est=="att"
qui reshape wide b se p lo hi rc, i(draw) j(arm) string
gen double d1 = abs(banalytic - bboot_standard)
gen double d2 = abs(banalytic - bboot_bayes)
gen double d3 = abs(banalytic - bpermute)
sum d1 d2 d3
