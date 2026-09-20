*=====================================================================
* pG: FUZZ-08 -- do analytic SEs move when the data are re-sorted?
* Step 1 PRISTINE TWICE (same sort, same call, twice).
* Step 2 reverse sort through wooldid.
* Step 3 reghdfe directly, forward vs reverse sort.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz-refute/gen_mypanel.do"
local mN 40
local mseed 90210
include "`GEN'"
gen double fwt = 1 + mod(_n,5)/2
tempfile base
qui save "`base'", replace

cap prog drop shot
prog def shot
  syntax , tag(string)
  cap noi qui wooldid y i t ttre [aw=fwt], att robust controls(x1) esprelength(3) level(90) safety(off)
  local rc = _rc
  if `rc' != 0 {
    dis "SHOT `tag' rc=`rc'"
    exit
  }
  local b = e(b)[1,1]
  mat V = e(V)
  local s = sqrt(V[1,1])
  mat M = e(wdidmainresults)
  dis "SHOT `tag' b=" %21x `b' " se=" %21x `s' " t=" %21x M[1,3] " lb=" %21x M[1,5]
end

* --- pristine twice, identical sort ---
use "`base'", clear
sort i t
shot , tag(P1)
use "`base'", clear
sort i t
shot , tag(P2)
use "`base'", clear
sort i t
shot , tag(P3)

* --- reverse sort ---
use "`base'", clear
gsort -i -t
shot , tag(R1)
use "`base'", clear
gsort -i -t
shot , tag(R2)

* --- random shuffle of rows ---
use "`base'", clear
set seed 4321
gen double sh = runiform()
sort sh
drop sh
shot , tag(S1)

* --- id relabelling (monotone, order preserving) ---
use "`base'", clear
qui replace i = i*7 + 3
sort i t
shot , tag(IDX)

* --- y rescaled by 1000 (should scale b and se exactly) ---
use "`base'", clear
qui replace y = y*1000
sort i t
cap noi qui wooldid y i t ttre [aw=fwt], att robust controls(x1) esprelength(3) level(90) safety(off)
local b = e(b)[1,1]
mat V = e(V)
local s = sqrt(V[1,1])
mat M = e(wdidmainresults)
dis "SHOT YSC b/1000=" %21x (`b'/1000) " se/1000=" %21x (`s'/1000) " t=" %21x M[1,3]

* --- reghdfe directly: forward vs reverse ---
use "`base'", clear
sort i t
qui reghdfe y x1 [aw=fwt], absorb(i t) vce(robust) tol(1e-9)
mat Vf = e(V)
scalar bf = _b[x1]
scalar sf = _se[x1]
use "`base'", clear
gsort -i -t
qui reghdfe y x1 [aw=fwt], absorb(i t) vce(robust) tol(1e-9)
mat Vr = e(V)
scalar br = _b[x1]
scalar sr = _se[x1]
dis "RH_FWD b=" %21x bf " se=" %21x sf
dis "RH_REV b=" %21x br " se=" %21x sr
dis "RH_SAME_b=" (bf==br) "  RH_SAME_se=" (sf==sr)
dis "RH_RELDIFF_se=" %12.4e (abs(sf-sr)/abs(sf))
dis "=================== DONE pG ==================="
