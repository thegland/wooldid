* v02: do the e(perm_p_<est>) / e(perm_lb_) / e(perm_ub_) scalars I harvest equal
* the numbers the PRINTED table carries?  (If not, the MC measures the wrong thing.)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
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
forvalues d = 1/3 {
  local sx = 2000000 + 10*1000 + `d'
  mkdgp , ngrp(10) nper(8) seedx(`sx')
  ereturn clear
  cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(555)
  mat M = e(wdidmainresults)
  local rn : rownames M
  local re : roweq M
  local nr : word count `rn'
  forvalues z = 1/`nr' {
    local nz : word `z' of `rn'
    local ez : word `z' of `re'
    if "`ez'" == "main" {
      dis "d`d' `nz'  TABLE p=" %21x M[`z',4] " lb=" %21x M[`z',5] " ub=" %21x M[`z',6]
      dis "d`d' `nz'  SCALR p=" %21x e(perm_p_`nz') " lb=" %21x e(perm_lb_`nz') " ub=" %21x e(perm_ub_`nz')
    }
  }
}

* --- permtstat: does e(perm_p_att) still equal the printed table's p? ---
mkdgp , ngrp(10) nper(8) seedx(2000010001)
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(555) permtstat
dis "TS_rc=" _rc "  perm_stat=" e(perm_stat)
mat M = e(wdidmainresults)
local rn : rownames M
local re : roweq M
local nr : word count `rn'
forvalues z = 1/`nr' {
  local nz : word `z' of `rn'
  local ez : word `z' of `re'
  if "`ez'" == "main" {
    dis "TS `nz'  TABLE p=" %21x M[`z',4] "  SCALR p=" %21x e(perm_p_`nz')
    dis "TS `nz'  TABLE lb=" %21x M[`z',5] "  SCALR lb=" %21x e(perm_lb_`nz')
  }
}

* --- PO design: is the permonsetcohorts comparator the observed fit? ---
clear
set seed 3000001
set sortseed 3000001
set obs 12
gen long st = _n
gen double u = runiform()
sort u
gen long rk = _n
gen double ttre = .
gen long coh = 4
replace ttre = 2004 if rk <= 3
replace coh  = 1    if rk <= 3
replace ttre = 2005 if rk > 3 & rk <= 6
replace coh  = 2    if rk > 3 & rk <= 6
replace ttre = 2006 if rk > 6 & rk <= 9
replace coh  = 3    if rk > 6 & rk <= 9
gen double ai = rnormal()
drop u rk
expand 8
bysort st: gen int t = 1999 + _n
sort t st
by t: gen double ft = 0.5*rnormal() if _n==1
by t: replace ft = ft[1]
sort st t
by st: gen double e = rnormal()
by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75) if _n>1
gen double y = ai + ft + e
keep st coh t ttre y
save v02_po.dta, replace
ereturn clear
cap noi wooldid y coh t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(555) permonsetcohorts
dis "PO_rc=" _rc "  irebuilt=" e(perm_irebuilt) "  npossible=" e(perm_npossible) "  completed=" e(perm_completed) " failed=" e(perm_failed)
mat CB = e(perm_comparator_b)
mat BB = e(b)
dis "PO comparator vs e(b):"
mat list CB
mat list BB
dis "PO p_att=" e(perm_p_att) " p_att_i=" e(perm_p_att_i)

* --- EX design: enumeration size and exact flag ---
clear
set seed 3000001
set sortseed 3000001
set obs 8
gen long st = _n
gen double u = runiform()
sort u
gen long rk = _n
gen double ttre = .
replace ttre = 2004 if rk <= 2
replace ttre = 2006 if rk == 3
gen double ai = rnormal()
drop u rk
expand 8
bysort st: gen int t = 1999 + _n
sort t st
by t: gen double ft = 0.5*rnormal() if _n==1
by t: replace ft = ft[1]
sort st t
by st: gen double e = rnormal()
by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75) if _n>1
gen double y = ai + ft + e
keep st t ttre y
ereturn clear
cap noi wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permexact rseed(555)
dis "EX_rc=" _rc "  exact=" e(perm_exact) " npossible=" e(perm_npossible) " completed=" e(perm_completed) " failed=" e(perm_failed)
dis "EX p_att=" e(perm_p_att) " p_att_i=" e(perm_p_att_i)
