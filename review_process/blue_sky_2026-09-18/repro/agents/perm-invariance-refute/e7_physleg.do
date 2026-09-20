*==========================================================================
* e7_physleg.do -- SELF-CONTAINED.  Companion to e5_boot.do.
* e5 measured the REWEIGHTING legs (att under bootstrap(standard) = cluster
* frequency reweighting; att under bootstrap(bayes) = Dirichlet).  This file
* measures the PHYSICAL-RESAMPLE legs, where the scale-invariance argument does
* NOT apply: a cell-weighted estimand (att_i) under bootstrap(standard), and
* boottstat.  Does a phantom cluster move those?
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do e7_physleg.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

clear
set seed 31337
set sortseed 31337
qui set obs 27
qui gen long unit  = _n
qui gen long state = ceil(unit/3)
qui gen double ai  = rnormal()
qui gen double ttre = .
qui replace ttre = 2004 if state <= 2
qui replace ttre = 2006 if state == 3 | state == 4
qui expand 8
qui bysort unit: gen int t = 1999 + _n
qui gen double eps = rnormal()*0.4
qui gen byte post = (t >= ttre) & !missing(ttre)
qui gen double y  = ai + 0.15*(t-2000) + 0.9*post + eps
qui drop eps
tempfile base
save `base'

prog def brpt
  args tag
  dis "RES `tag' rc=" _rc " N=" e(N) " Nclust=" e(N_clust) " b=" %12.8f e(wdidmainresults)[1,1] " se=" %12.8f e(wdidmainresults)[1,2]
end

foreach sd in 101 202 303 {
  use `base', clear
  qui replace y = . if state == 9
  cap noi wooldid y unit t ttre , att_i cluster(state) bootstrap(standard) reps(400) rseed(`sd') dots(0)
  brpt A_atti_`sd'
  use `base', clear
  qui drop if state == 9
  cap noi wooldid y unit t ttre , att_i cluster(state) bootstrap(standard) reps(400) rseed(`sd') dots(0)
  brpt B_atti_`sd'
}
foreach sd in 101 202 {
  use `base', clear
  qui replace y = . if state == 9
  cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) boottstat reps(200) rseed(`sd') dots(0)
  brpt A_bts_`sd'
  use `base', clear
  qui drop if state == 9
  cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) boottstat reps(200) rseed(`sd') dots(0)
  brpt B_bts_`sd'
}
dis "DONE_E7"
