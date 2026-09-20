*======================================================================
* SELF-CONTAINED REPRO -- size of a wooldid inference engine under a
* TRUE SHARP NULL (no treatment effect at any horizon).
*
*   Run with:  stata-mp -e do <thisfile>.do
*   Nothing outside this file is needed except the pinned wooldid build.
*
* EDIT THESE THREE LINES to select the arm:
local NDRAW  100                  /* Monte Carlo draws                */
local GG     10                   /* number of states = clusters = i  */
local ARM    "perm"               /* perm | permt | boot | boottstat | bayes | analytic */
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

* ---- DGP: staggered state-year panel, NO treatment effect anywhere ----
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
  gen double y = ai + ft + e          /* <- no treatment term at all */
  keep st t ttre y
end

tempname P
postfile `P' int draw double p double lo double hi int rc using repro_out.dta, replace
forvalues d = 1/`NDRAW' {
  local sx = 2000000 + `GG'*1000 + `d'
  local rs = 700000 + `d'*7 + `GG'
  mkdgp , ngrp(`GG') nper(8) seedx(`sx')
  ereturn clear
  if "`ARM'" == "perm"      cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs')
  if "`ARM'" == "permt"     cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs') permtstat
  if "`ARM'" == "boot"      cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs')
  if "`ARM'" == "boottstat" cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs') boottstat
  if "`ARM'" == "bayes"     cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(bayes) reps(79) rseed(`rs')
  if "`ARM'" == "analytic"  cap noi qui wooldid y st t ttre, cluster(st) att att_i
  local rc = _rc
  local p = .
  local lo = .
  local hi = .
  if `rc' == 0 {
    if inlist("`ARM'","perm","permt") {
      local p  = e(perm_p_att)
      local lo = e(perm_lb_att)
      local hi = e(perm_ub_att)
    }
    else {
      cap mat drop TT
      if "`ARM'" == "analytic"  cap mat TT = e(wdidmainresults)
      else if "`ARM'" == "boottstat" cap mat TT = e(boot_tpctl_mainresults)
      else cap mat TT = e(boot_pctl_mainresults)
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
          local p  = TT[`j',4]
          local lo = TT[`j',5]
          local hi = TT[`j',6]
        }
      }
    }
  }
  post `P' (`d') (`p') (`lo') (`hi') (`rc')
}
postclose `P'
use repro_out.dta, clear
gen byte rej  = (p <= 0.05 + 1e-9) if p < .
gen byte excl = (lo > 0 | hi < 0)  if lo < . & hi < .
dis "=== ARM=`ARM'  G=`GG'  draws=`NDRAW'  TRUE EFFECT = 0 EVERYWHERE ==="
sum rej excl
count if rc != 0
dis "NONZERO_RC=" r(N)
