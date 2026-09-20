*=== d01_specdiag.do : why does A0 report att=-0.295 and A5 att=+0.462 on the
*    SAME draw?  Candidates: (a) esprelength/espostlength change the estimate,
*    (b) an earlier wooldid call mutates the data for the next one.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

*---------------- DGP: SHARP NULL state-year panel (no treatment effect) -----
* ngrp   states (= i = cluster);  nper periods 2000..2000+nper-1
* onsets 2004/2005/2006 to 20% of states each; 40% never treated
* y = a_i + f_t + e_it ,  e_it AR(1) rho=.5 within state (Bertrand-Duflo-
* Mullainathan serial correlation).  hetv scales the error sd of TREATED
* states only (hetv=1 -> homogeneous).
cap prog drop mkdgp
prog def mkdgp
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1) HETV(real 1) NPERCELL(integer 1)]
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
  gen double esd = cond(ttre < ., `hetv', 1)
  drop u rk
  expand `nper'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()*esd if _n==1
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(1-0.25)*esd if _n>1
  gen double y = ai + ft + e
  if `npercell' > 1 {
    expand `npercell'
    bysort st t: gen int sub = _n
    gen long unit = st*1000 + sub
    replace y = y + rnormal()*0.5
  }
  keep st t ttre y
  order st t ttre y
end

local sx = 2000000 + 100*100000 + 10*1000 + 1
mkdgp , ngrp(10) nper(8) seedx(`sx')
save d01_data.dta, replace
count
dis "NOBS=" r(N)
codebook ttre, compact

* 1. plain, no ES
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i
dis "S1_rc=" _rc
mat A = e(wdidmainresults)
mat list A

* 2. plain, WITH ES window
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) jointtests
dis "S2_rc=" _rc
mat B = e(wdidmainresults)
mat list B

* 3. plain, no ES AGAIN on the same in-memory data (data mutation test)
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i
dis "S3_rc=" _rc
mat C = e(wdidmainresults)
mat list C

* 4. reload pristine data and repeat 1
use d01_data.dta, clear
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i
dis "S4_rc=" _rc
mat D = e(wdidmainresults)
mat list D

* 5. ES window WIDE enough to cover every realized horizon
use d01_data.dta, clear
ereturn clear
cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(9) espostlength(9)
dis "S5_rc=" _rc
mat E = e(wdidmainresults)
mat list E

dis "ATT_noES  = " %21x A[1,1]
dis "ATT_ES33  = " %21x B[1,1]
dis "ATT_noES2 = " %21x C[1,1]
dis "ATT_fresh = " %21x D[1,1]
dis "ATT_ES99  = " %21x E[1,1]

* 6. NOISY run with ES window 3/3 to see any caution
use d01_data.dta, clear
ereturn clear
cap noi wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) jointtests
dis "S6_rc=" _rc
