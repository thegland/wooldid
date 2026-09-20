*=========================================================================
* p4_mech.do -- mechanism probes.
*  (A) mcci-03: is it esprelength() alone that moves the att, or espostlength()
*      too?  Is e(escomp_refN) returned, and does it shrink with the request?
*  (B) mcci-01 boundary: impvalues with i = unit and cluster(unit) -- the
*      geometry where each cluster owns its own coefficients.  If the IM se
*      is sane there, the defect is specifically "cluster finer than i".
*  (C) mcci-01 boundary: cluster COARSER than i (region nesting cohorts).
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p4_mech.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkdat
prog def mkdat
  args nunit sd
  clear
  set seed `sd'
  set sortseed `sd'
  set obs `nunit'
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen byte reg = mod(_n-1,5)+1
  gen double ttre = .
  replace ttre = 2004 if cg==1
  replace ttre = 2006 if cg==2
  replace ttre = 2008 if cg==3
  gen double ui   = rnormal(0,0.70)
  gen double hmul = 1 + rnormal(0,0.25)
  expand 12
  bys unit: gen int t = 2000 + _n - 1
  gen double e = .
  bys unit (t): replace e = rnormal(0,0.55) if _n==1
  bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
  gen double truetau = 0
  replace truetau = 1.00*hmul if cg==1 & t>=ttre
  replace truetau = 1.15*hmul if cg==2 & t>=ttre
  replace truetau = 0.85*hmul if cg==3 & t>=ttre
  gen double y = ui + 0.25*(t-2000) + truetau + e
end

mkdat 60 4411
tempfile dA
qui save "`dA'"

di as txt _n "########## A : which ES option moves the att ? ##########"
* spec list: none / post only / pre only / both / both+summarystats
local k = 0
foreach spec in none post2 pre2 both2 sumstat {
  use "`dA'", clear
  local opt ""
  if "`spec'"=="post2"   local opt "espostlength(2)"
  if "`spec'"=="pre2"    local opt "esprelength(2) espostlength(-999)"
  if "`spec'"=="both2"   local opt "esprelength(2) espostlength(2)"
  if "`spec'"=="sumstat" local opt "esprelength(2) espostlength(2) summarystats"
  di as txt _n "=== A spec=`spec' : `opt'"
  cap noi wooldid y cg t ttre, att cluster(unit) `opt'
  di "A_RC_`spec'=" _rc
  if _rc==0 {
    mat T = e(wdidmainresults)
    di "A_`spec'_b="  %12.9f T[1,1]
    di "A_`spec'_se=" %12.9f T[1,2]
    di "A_`spec'_N="  e(N)
    local rr = .
    cap local rr = e(escomp_refN)
    di "A_`spec'_refN=`rr'"
    local rw = .
    cap local rw = e(escomp_refW)
    di "A_`spec'_refW=`rw'"
  }
}

di as txt _n "########## B : impvalues with i = unit, cluster(unit) ##########"
use "`dA'", clear
di as txt _n "=== B1 i=unit cluster(unit) impvalues"
cap noi wooldid y unit t ttre, att cluster(unit) impvalues level(95)
di "B1_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B1_b="   %12.9f T[1,1]
  di "B1_se="  %12.9f T[1,2]
  di "B1_dfr=" e(df_r)
}
use "`dA'", clear
di as txt _n "=== B2 i=unit cluster(unit) analytic"
cap noi wooldid y unit t ttre, att cluster(unit) level(95)
di "B2_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B2_b="   %12.9f T[1,1]
  di "B2_se="  %12.9f T[1,2]
  di "B2_dfr=" e(df_r)
}

di as txt _n "########## C : cluster COARSER than i (reg nests nothing, crosses cg) ##########"
* reg has 5 levels and CROSSES cg (each reg holds units of every cohort)
use "`dA'", clear
di as txt _n "=== C1 i=cg cluster(reg) impvalues"
cap noi wooldid y cg t ttre, att cluster(reg) impvalues level(95)
di "C1_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "C1_b="   %12.9f T[1,1]
  di "C1_se="  %12.9f T[1,2]
  di "C1_dfr=" e(df_r)
}
use "`dA'", clear
di as txt _n "=== C2 i=cg cluster(reg) analytic"
cap noi wooldid y cg t ttre, att cluster(reg) level(95)
di "C2_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "C2_b="   %12.9f T[1,1]
  di "C2_se="  %12.9f T[1,2]
}

di as txt _n "ALLDONE_P4"
