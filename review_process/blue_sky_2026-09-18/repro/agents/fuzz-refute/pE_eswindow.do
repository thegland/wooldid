*=====================================================================
* pE: FUZZ-05 -- does esprelength() reprice the headline ATT?
* Discriminating test: ZERO-NOISE designs where the true ATT is known
* exactly and every pre-treatment effect is identically zero.  If the
* reported att is invariant to esprelength() there, the 5-11% move the
* prober measured is sampling variation across specifications, not a
* repricing defect.  Plus the level()-restore poison persistence check.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*---------------------------------------------------------------------
* 0. level() restore poison: does the LEAKED WOOLDID_LEVEL0 keep firing?
*---------------------------------------------------------------------
cap prog drop mkexact
prog def mkexact
  syntax , n(integer) [gaps(real 0) unitfe(real 1) ramp(real 0) noise(real 0) sd(integer 11)]
  clear
  set seed `sd'
  set sortseed `sd'
  qui set obs `n'
  gen long unit = _n
  gen double a_i = `unitfe'*rnormal()
  gen double u = runiform()
  gen double ttre = .
  qui replace ttre = 2003 if u > 0.35 & u <= 0.5667
  qui replace ttre = 2005 if u > 0.5667 & u <= 0.7833
  qui replace ttre = 2007 if u > 0.7833
  qui expand 10
  bysort unit: gen int t = 1999 + _n
  gen double ft = 0.4*(t-2000) + 0.1*(t-2000)^2
  gen byte D = (ttre < . & t >= ttre)
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 2 + `ramp'*rel if D==1
  gen double y = a_i + ft + tau + `noise'*rnormal()
  if `gaps' > 0 {
    gen double g = runiform()
    qui drop if g < `gaps' & t > 2000 & t < 2008
    drop g
  }
  gen long i = unit
  egen long icoh = group(ttre)
  qui replace icoh = 99 if ttre == .
  drop u a_i ft D
  * design truth: simple mean of tau over treated (i,t) cells, att weighting
  qui sum tau if rel >= 0 & ttre < .
  c_local truth = r(mean)
  c_local ncells = r(N)
end

cap prog drop runwin
prog def runwin
  syntax , tag(string) ivar(string) [pre(integer 0) post(integer 0) cl(string)]
  local o
  if `pre'  > 0 local o `o' esprelength(`pre')
  if `post' > 0 local o `o' espostlength(`post')
  local v robust
  if "`cl'" != "" local v cluster(`cl')
  cap noi qui wooldid y `ivar' t ttre , att `v' `o' safety(off)
  local rc = _rc
  if `rc' == 0 {
    local b  = e(b)[1,1]
    mat V = e(V)
    local s  = sqrt(V[1,1])
    local rn = e(escomp_refN)
    dis "ESW `tag' pre=`pre' post=`post' rc=0 att=" %21x `b' " att10=" %15.10f `b' " se=" %12.8f `s' " N=" e(N) " refN=" `rn'
  }
  else dis "ESW `tag' pre=`pre' post=`post' rc=`rc'"
end

*---- level poison persistence ----
mkexact , n(40) noise(0.3)
tempfile lvl
qui save "`lvl'", replace
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 99.9
cap noi qui wooldid y i t ttre , att robust level(95) safety(off)
dis "E0_FIRST_RC=" _rc "  LEVEL0=[${WOOLDID_LEVEL0}]  c(level)=" c(level)
cap noi qui wooldid y i t ttre , att robust safety(off)
dis "E0_SECOND_RC=" _rc "  (no level() asked for at all)  e(N)=" e(N) "  LEVEL0=[${WOOLDID_LEVEL0}]"
cap noi qui wooldid y i t ttre , att robust safety(off)
dis "E0_THIRD_RC=" _rc "  LEVEL0=[${WOOLDID_LEVEL0}]"
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 95

*---------------------------------------------------------------------
* E1 balanced, ZERO noise, homogeneous effect, i = unit
*---------------------------------------------------------------------
mkexact , n(60) noise(0) ramp(0) gaps(0)
dis "E1_TRUTH=" %15.10f `truth' " ncells=" `ncells'
foreach p in 0 1 2 3 {
  runwin , tag(E1) ivar(i) pre(`p') post(0)
}
runwin , tag(E1) ivar(i) pre(0) post(3)
runwin , tag(E1) ivar(i) pre(2) post(3)

*---------------------------------------------------------------------
* E2 balanced, ZERO noise, DYNAMIC ramp effect, i = unit
*---------------------------------------------------------------------
mkexact , n(60) noise(0) ramp(0.5) gaps(0)
dis "E2_TRUTH=" %15.10f `truth' " ncells=" `ncells'
foreach p in 0 1 2 3 {
  runwin , tag(E2) ivar(i) pre(`p') post(0)
}
runwin , tag(E2) ivar(i) pre(0) post(3)

*---------------------------------------------------------------------
* E3 UNBALANCED (25% interior gaps), ZERO noise, ramp, i = unit
*---------------------------------------------------------------------
mkexact , n(80) noise(0) ramp(0.5) gaps(0.25)
dis "E3_TRUTH=" %15.10f `truth' " ncells=" `ncells'
foreach p in 0 1 2 3 {
  runwin , tag(E3) ivar(i) pre(`p') post(0)
}

*---------------------------------------------------------------------
* E4 UNBALANCED, ZERO noise, ramp, i = COHORT, unit FE present (the
*    aggregation the author calls the user's own omitted-variable problem)
*---------------------------------------------------------------------
mkexact , n(80) noise(0) ramp(0.5) gaps(0.25) unitfe(1)
dis "E4_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(E4) ivar(icoh) pre(`p') post(0)
}

*---------------------------------------------------------------------
* E5 UNBALANCED, ZERO noise, ramp, i = COHORT, NO unit FE
*---------------------------------------------------------------------
mkexact , n(80) noise(0) ramp(0.5) gaps(0.25) unitfe(0)
dis "E5_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(E5) ivar(icoh) pre(`p') post(0)
}

*---------------------------------------------------------------------
* E6 NOISY realistic: how big is the move relative to its own SE?
*---------------------------------------------------------------------
mkexact , n(150) noise(1.0) ramp(0.5) gaps(0.15) sd(777)
dis "E6_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(E6) ivar(i) pre(`p') post(0) cl(i)
}
dis "=================== DONE pE ==================="
