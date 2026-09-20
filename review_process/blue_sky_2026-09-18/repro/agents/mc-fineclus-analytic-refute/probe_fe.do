* ============================================================================
* probe_fe.do -- follow-up to probe_unbal2.do.  Does fe(state) repair the
* composition bias that the coarse-i remedy geometry introduces, and what does
* it cost on the standard-error side?  ZERO NOISE: truth is exactly 1.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
cap program drop mkpanel
program define mkpanel
  clear
  qui set obs 40
  gen long state = _n
  gen double s_g = 2*mod(state,7) - 6
  gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
  qui expand 12
  bysort state: gen int t = _n
  gen double tau = cond(ttre<. & t>=ttre, 1, 0)
  gen double y = s_g + 0.2*t + tau
  gen double coh = cond(ttre<., ttre, 99)
  gen double g   = cond(ttre<., ttre, 0)
  set seed 4242
  qui gen double u = runiform()
end

di "===== ragged (8% MCAR), zero noise, truth = 1 ====="
qui mkpanel
qui drop if u<0.08
cap noi qui wooldid y state t ttre, att cluster(state)
di "A) i=state cluster(state)            rc=" _rc " att=" %18.12f _b[main_att] " se=" %12.8f _se[main_att]
qui mkpanel
qui drop if u<0.08
cap noi qui wooldid y coh t ttre, att cluster(state)
di "B) i=coh   cluster(state)            rc=" _rc " att=" %18.12f _b[main_att] " se=" %12.8f _se[main_att]
qui mkpanel
qui drop if u<0.08
cap noi qui wooldid y coh t ttre, att cluster(state) fe(state)
di "C) i=coh   cluster(state) fe(state)  rc=" _rc " att=" %18.12f _b[main_att] " se=" %12.8f _se[main_att]
qui mkpanel
qui drop if u<0.08
cap noi qui wooldid y coh t ttre, att cluster(state) controls(i.state)
di "D) i=coh   cluster(state) ctrl i.st  rc=" _rc " att=" %18.12f _b[main_att] " se=" %12.8f _se[main_att]

di "===== does wooldid say ANYTHING when i has sub-structure? (verbatim output) ====="
qui mkpanel
qui drop if u<0.08
cap noi wooldid y coh t ttre, att cluster(state)
di "PROBE_FE_DONE"
