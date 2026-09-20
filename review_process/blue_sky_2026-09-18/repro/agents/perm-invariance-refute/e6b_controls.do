*==========================================================================
* e6b_controls.do -- SELF-CONTAINED.  Companion to e6_bootscreens.do.
* Does the SAME suppression of bootstrap()'s few-TREATED-cluster caution occur
* when the phantom treated cluster is killed by a missing CONTROL rather than a
* missing outcome?  (smp2, ado:6461, is the same predicate either way.)
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do e6b_controls.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

clear
set seed 771
set sortseed 771
qui set obs 36
qui gen long unit  = _n
qui gen long state = ceil(unit/3)
qui gen double ai  = rnormal()
qui gen double ttre = .
qui replace ttre = cond(mod(state,2)==1, 2004, 2006) if state <= 6
qui expand 8
qui bysort unit: gen int t = 1999 + _n
qui gen double x1  = rnormal()
qui gen double eps = rnormal()*0.4
qui gen byte post  = (t >= ttre) & !missing(ttre)
qui gen double y   = ai + 0.15*(t-2000) + 0.5*x1 + 0.9*post + eps
qui drop eps
tempfile b
save `b'

dis "=== D1: state 6 treated, x1 missing on every row of it (6 design-treated) ==="
use `b', clear
qui replace x1 = . if state == 6
cap noi wooldid y unit t ttre , att controls(x1) cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES D1 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)

dis "=== D2: same estimation sample, state 6's rows dropped (5 design-treated) ==="
use `b', clear
qui drop if state == 6
cap noi wooldid y unit t ttre , att controls(x1) cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES D2 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)

dis "DONE_E6B"
