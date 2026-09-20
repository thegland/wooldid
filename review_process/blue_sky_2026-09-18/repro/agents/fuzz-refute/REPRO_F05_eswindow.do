clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ============ FUZZ-05: does esprelength() reprice the headline ATT? ============
* Discriminating test: ZERO noise, exact parallel trends, effects known exactly.
* Every pre-treatment effect is identically zero, so the reported att must not
* depend on how many pre periods get their own dummy.
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
    gen double gg = runiform()
    qui drop if gg < `gaps' & t > 2000 & t < 2008
    drop gg
  }
  gen long i = unit
  egen long icoh = group(ttre)
  qui replace icoh = 99 if ttre == .
  drop u a_i ft D
  qui sum tau if rel >= 0 & ttre < .
  c_local truth = r(mean)
end
cap prog drop runwin
prog def runwin
  syntax , tag(string) ivar(string) [pre(integer 0) cl(string)]
  local o
  if `pre' > 0 local o esprelength(`pre')
  local v robust
  if "`cl'" != "" local v cluster(`cl')
  cap noi qui wooldid y `ivar' t ttre , att `v' `o' safety(off)
  local b = e(b)[1,1]
  mat V = e(V)
  local s = sqrt(V[1,1])
  dis "`tag' esprelength(`pre') att=" %21x `b' " = " %15.10f `b' "  se=" %12.8f `s' "  N=" e(N)
end
* A: zero noise, homogeneous effect, i = unit  -> truth 2
mkexact , n(60) noise(0) ramp(0) gaps(0)
dis "A_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(A) ivar(i) pre(`p')
}
* B: zero noise, DYNAMIC effect, unbalanced panel, i = unit
mkexact , n(80) noise(0) ramp(0.5) gaps(0.25)
dis "B_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(B) ivar(i) pre(`p')
}
* C: zero noise, unbalanced, i = COHORT, unit FE present (coarse-i aggregation)
mkexact , n(80) noise(0) ramp(0.5) gaps(0.25) unitfe(1)
dis "C_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(C) ivar(icoh) pre(`p')
}
* D: same but with NO unit heterogeneity -> invariance returns
mkexact , n(80) noise(0) ramp(0.5) gaps(0.25) unitfe(0)
dis "D_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(D) ivar(icoh) pre(`p')
}
* E: NOISY -- compare the move with its own SE
mkexact , n(150) noise(1.0) ramp(0.5) gaps(0.15) sd(777)
dis "E_TRUTH=" %15.10f `truth'
foreach p in 0 1 2 3 {
  runwin , tag(E) ivar(i) pre(`p') cl(i)
}
