*==========================================================================
* e5_boot.do -- SELF-CONTAINED.  Tests PI-1's UNMEASURED extension that
* bootstrap() "shares the design population" and therefore shares the defect.
*
* Theory says the shared predicate has NO consequence under bootstrap(bayes):
* the normalized sub-vector of a Dirichlet over the real clusters is again
* Dirichlet over their own shapes, and aweights are scale-free, so the
* replication distribution is EXACTLY invariant to a phantom cluster -- only
* the RNG stream moves.  Under bootstrap(standard) a drawn phantom is a wasted
* draw, so some inflation is expected.  Three seeds per arm to separate the
* structural effect from Monte-Carlo noise.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do e5_boot.do
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
qui gen double eps  = rnormal()*0.4
qui gen byte   post = (t >= ttre) & !missing(ttre)
qui gen double y    = ai + 0.15*(t-2000) + 0.9*post + eps
qui drop eps
tempfile base
save `base'
* state 9 is the phantom (never-treated, y missing on every row)

prog def brpt
  args tag
  dis "RES `tag' rc=" _rc " N=" e(N) " Nclust=" e(N_clust) " b=" %12.8f e(wdidmainresults)[1,1] " se=" %12.8f e(wdidmainresults)[1,2] " lo=" %12.8f e(wdidmainresults)[1,5] " hi=" %12.8f e(wdidmainresults)[1,6]
end

foreach md in bayes standard {
  foreach sd in 101 202 303 {
    * A: phantom present, y missing on every one of its rows
    use `base', clear
    qui replace y = . if state == 9
    cap noi wooldid y unit t ttre , att cluster(state) bootstrap(`md') reps(400) rseed(`sd') dots(0) safety(off)
    brpt A_`md'_`sd'
    * B: phantom rows dropped (identical estimation sample)
    use `base', clear
    qui drop if state == 9
    cap noi wooldid y unit t ttre , att cluster(state) bootstrap(`md') reps(400) rseed(`sd') dots(0) safety(off)
    brpt B_`md'_`sd'
  }
}
dis "DONE_E5"
