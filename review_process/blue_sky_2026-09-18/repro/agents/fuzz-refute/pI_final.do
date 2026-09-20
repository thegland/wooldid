*=====================================================================
* pI: closing probes.
*  (1) permute() + semielasticity on a zero-straddling outcome
*  (2) what the shift ladder does to the FITTED values (cell means)
*  (3) a harder sort-order test for FUZZ-08 (ill-conditioned absorb,
*      weights over 3 orders of magnitude, cluster VCE)
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set seed 606
set sortseed 606
qui set obs 60
gen long unit = _n
gen double a_i = rnormal(0,0.5)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.4 & u <= 0.6
qui replace ttre = 2005 if u > 0.6 & u <= 0.8
qui replace ttre = 2007 if u > 0.8
qui expand 10
bysort unit: gen int t = 1999 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double g = a_i + 0.3*rnormal() + 0.6*D + 0.2*x1
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i
tempfile gd
qui save "`gd'", replace

* (1) permute leg, legal permreps
cap noi wooldid g i t ttre , att cluster(i) semielasticity permute(shuffle) permreps(24) rseed(3) dots(0)
dis "I_PERM_SEMI_RC=" _rc

* (2) fitted-value (saturated cell mean) minimum along the shift ladder
foreach sh in 0 1 2 {
  use "`gd'", clear
  qui replace g = g + `sh'
  qui egen double cm = mean(g), by(icoh t D)
  qui sum cm
  local cmin = r(min)
  qui sum g
  local ymin = r(min)
  cap noi qui wooldid g i t ttre , att cluster(i) semielasticity
  dis "I_SHIFT_`sh'  ymin=" %9.4f `ymin' "  min_cellmean=" %9.4f `cmin' "  rc=" _rc
}

* (3) harder sort test
clear
set seed 1234567
set sortseed 1234567
qui set obs 80
gen long unit = _n
gen double a_i = rnormal(0,3)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.30 & u <= 0.53
qui replace ttre = 2005 if u > 0.53 & u <= 0.76
qui replace ttre = 2007 if u > 0.76
qui expand 12
bysort unit: gen int t = 1997 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double x2 = x1 + 0.001*rnormal()
gen double fwt = exp(rnormal(0,1.5))
gen double y = a_i + 0.7*(t-2000) + 2.5*D + 0.4*x1 - 0.2*x2 + rnormal()
gen long i = unit
gen byte region = mod(unit,10)+1
tempfile hard
qui save "`hard'", replace

cap prog drop shot2
prog def shot2
  syntax , tag(string)
  cap noi qui wooldid y i t ttre [aw=fwt], att cluster(region) controls(x1 x2) esprelength(3) level(90) safety(off)
  local rc = _rc
  if `rc' != 0 {
    dis "SHOT2 `tag' rc=`rc'"
    exit
  }
  local b = e(b)[1,1]
  mat V = e(V)
  local s = sqrt(V[1,1])
  mat M = e(wdidmainresults)
  dis "SHOT2 `tag' b=" %21x `b' " se=" %21x `s'
  dis "SHOT2 `tag' t=" %21x M[1,3] " p=" %21x M[1,4] " lb=" %21x M[1,5]
end
use "`hard'", clear
sort i t
shot2 , tag(H_FWD1)
use "`hard'", clear
sort i t
shot2 , tag(H_FWD2)
use "`hard'", clear
gsort -i -t
shot2 , tag(H_REV1)
use "`hard'", clear
gsort -i -t
shot2 , tag(H_REV2)
use "`hard'", clear
set seed 99
gen double sh = runiform()
sort sh
drop sh
shot2 , tag(H_SHUF)
* regtol sweep on the reversed sort
foreach rt in 9 12 13 {
  use "`hard'", clear
  sort i t
  cap noi qui wooldid y i t ttre [aw=fwt], att cluster(region) controls(x1 x2) esprelength(3) level(90) safety(off) regtol(`rt')
  local bf = e(b)[1,1]
  mat Vf = e(V)
  local sf = sqrt(Vf[1,1])
  use "`hard'", clear
  gsort -i -t
  cap noi qui wooldid y i t ttre [aw=fwt], att cluster(region) controls(x1 x2) esprelength(3) level(90) safety(off) regtol(`rt')
  local br = e(b)[1,1]
  mat Vr = e(V)
  local sr = sqrt(Vr[1,1])
  dis "I_REGTOL `rt' relb=" %12.4e (abs(`bf'-`br')/abs(`bf')) "  relse=" %12.4e (abs(`sf'-`sr')/abs(`sf'))
}
dis "=================== DONE pI ==================="
