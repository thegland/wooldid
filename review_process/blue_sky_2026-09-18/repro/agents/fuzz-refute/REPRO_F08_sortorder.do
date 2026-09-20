clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ============ FUZZ-08: does re-sorting move the analytic SEs? ============
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
drop u a_i
tempfile hard
qui save "`hard'", replace
cap prog drop shot2
prog def shot2
  syntax , tag(string)
  cap noi qui wooldid y i t ttre [aw=fwt], att cluster(region) controls(x1 x2) esprelength(3) level(90) safety(off)
  local b = e(b)[1,1]
  mat V = e(V)
  local s = sqrt(V[1,1])
  mat M = e(wdidmainresults)
  dis "`tag' b=" %21x `b' "  se=" %21x `s' "  t=" %21x M[1,3]
end
* PRISTINE TWICE first -- rule out plain nondeterminism
use "`hard'", clear
sort i t
shot2 , tag(FWD1)
use "`hard'", clear
sort i t
shot2 , tag(FWD2)
use "`hard'", clear
gsort -i -t
shot2 , tag(REV1)
use "`hard'", clear
gsort -i -t
shot2 , tag(REV2)
use "`hard'", clear
set seed 99
gen double sh = runiform()
sort sh
drop sh
shot2 , tag(SHUF)
* magnitudes, and whether regtol() moves them
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
  dis "regtol(`rt') rel|db|=" %12.4e (abs(`bf'-`br')/abs(`bf')) "  rel|dse|=" %12.4e (abs(`sf'-`sr')/abs(`sf'))
}
* attribution: reghdfe alone, forward vs reverse
use "`hard'", clear
sort i t
qui reghdfe y x1 [aw=fwt], absorb(i t) vce(robust) tol(1e-9)
scalar bf2 = _b[x1]
scalar sf2 = _se[x1]
use "`hard'", clear
gsort -i -t
qui reghdfe y x1 [aw=fwt], absorb(i t) vce(robust) tol(1e-9)
dis "reghdfe alone: b identical=" (bf2==_b[x1]) "  se identical=" (sf2==_se[x1])
