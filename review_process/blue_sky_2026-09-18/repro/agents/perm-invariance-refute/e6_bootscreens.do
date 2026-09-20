*==========================================================================
* e6_bootscreens.do -- SELF-CONTAINED.  NEW probe (not seeded).
*
* CLAIM UNDER TEST: bootstrap()'s cluster-count screens (wooldid.ado:6460-6500)
* count clusters over  smp2 = ((if)==1) & !missing(cluster)  -- the same
* phantom-inclusive population permute() uses -- and NOT over the estimation
* sample.  Consequences:
*   (a) the "<10 clusters are available for the bootstrap" caution is SUPPRESSED
*       by a cluster that contributes no estimable row;
*   (b) the "<6 treated clusters are available for resampling" caution -- the
*       few-treated-cluster warning the help calls the important one -- is
*       likewise SUPPRESSED by a treated cluster with no estimable row;
*   (c) the refusal whose own message says "bootstrap() requires at least 2
*       clusters IN THE ESTIMATION SAMPLE" is passed by a sample holding one.
* sthlp:430-431 documents the phantom-inclusive population for permute()'s
* randomization design and its possible-assignment count.  It says nothing
* about bootstrap()'s screens, and the two-cluster message names the estimation
* sample explicitly.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do e6_bootscreens.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

prog def mk
  * args: nstates ntreated  -> nstates states x 3 counties x 8 years,
  * the first `ntreated' states treated (alternating 2004 / 2006)
  args ns nt
  clear
  set seed 771
  set sortseed 771
  qui set obs `=3*`ns''
  qui gen long unit  = _n
  qui gen long state = ceil(unit/3)
  qui gen double ai  = rnormal()
  qui gen double ttre = .
  qui replace ttre = cond(mod(state,2)==1, 2004, 2006) if state <= `nt'
  qui expand 8
  qui bysort unit: gen int t = 1999 + _n
  qui gen double eps = rnormal()*0.4
  qui gen byte post  = (t >= ttre) & !missing(ttre)
  qui gen double y   = ai + 0.15*(t-2000) + 0.9*post + eps
  qui drop eps
end

*==================== (a) total-cluster caution ============================
* 9 real clusters (3 treated) + 1 phantom never-treated cluster = 10 in design
dis "=== A1: 9 real + 1 phantom (design count 10) ==="
mk 10 3
qui replace y = . if state == 10
cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES A1 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)

dis "=== A2: the same 9 real clusters, phantom rows dropped ==="
mk 10 3
qui drop if state == 10
cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES A2 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)

*==================== (b) few-TREATED-cluster caution =======================
* 12 states: 1-5 treated with data, 6 treated PHANTOM, 7-12 never-treated.
* design treated-cluster count 6 (no caution); estimable treated clusters 5.
dis "=== B1: 5 real treated + 1 phantom treated (design treated count 6) ==="
mk 12 6
qui replace y = . if state == 6
cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES B1 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)

dis "=== B2: the same 5 real treated, phantom rows dropped ==="
mk 12 6
qui drop if state == 6
cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES B2 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)

*==================== (c) the two-cluster refusal ===========================
* one real cluster (half its units treated) + one phantom cluster.
dis "=== C1: 1 real cluster + 1 phantom (design count 2) ==="
clear
set seed 771
set sortseed 771
qui set obs 9
qui gen long unit = _n
qui gen long state = cond(unit <= 6, 1, 2)
qui gen double ai = rnormal()
qui gen double ttre = .
qui replace ttre = 2004 if inlist(unit,1,2,3)
qui expand 8
qui bysort unit: gen int t = 1999 + _n
qui gen double eps = rnormal()*0.4
qui gen byte post = (t >= ttre) & !missing(ttre)
qui gen double y  = ai + 0.15*(t-2000) + 0.9*post + eps
qui drop eps
tempfile c
save `c'
use `c', clear
qui replace y = . if state == 2
cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES C1 rc=" _rc " N=" e(N) " Nclust=" e(N_clust)
if _rc == 0 dis "RES C1 b=" %21x e(wdidmainresults)[1,1] " se=" %21x e(wdidmainresults)[1,2] " lo=" %21x e(wdidmainresults)[1,5] " hi=" %21x e(wdidmainresults)[1,6]

dis "=== C2: the same estimation sample, phantom rows dropped ==="
use `c', clear
qui drop if state == 2
cap noi wooldid y unit t ttre , att cluster(state) bootstrap(standard) reps(60) rseed(7) dots(0)
dis "RES C2 rc=" _rc

dis "DONE_E6"
