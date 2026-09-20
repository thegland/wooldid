* ============================================================================
* probe_ar1.do -- SELF-CONTAINED. Why does MC-FC-1's K=10 rung (0.912 SE/SD) differ
* from mine (1.013) and from the project's own 2026-09-14 figure (0.972) at the SAME
* geometry?  The one structural difference in their generator is AR(1) errors within
* unit.  Serial correlation gives CR1 an extra finite-sample downward bias that has
* nothing to do with clusters-per-cell.  This file turns AR(1) on and off, K=10 fixed.
* 40 states x 4 units x 12 periods; states 1-10 adopt at t=5, 11-20 at t=7, 21-40 never.
* i = coh, cluster(state): cluster STRICTLY FINER than i, 10 states per treated cell.
* TRUE att = 1.1535714285714286 ; TRUE ES post h = 1 + .05h.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAWS = 120

cap program drop argen
program define argen
  args rho
  clear
  qui set obs 40
  gen long state = _n
  gen double s_g = rnormal()*0.8
  gen double ttre = cond(state<=10, 5, cond(state<=20, 7, .))
  qui expand 4
  bysort state: gen long unit = (state-1)*4 + _n
  gen double a_u = rnormal()*0.5
  qui expand 12
  bysort unit: gen int t = _n
  tempvar fst vv eps
  bysort state t: gen byte `fst' = (_n==1)
  gen double `vv' = rnormal()*0.4 if `fst'
  bysort state t (`fst'): replace `vv' = `vv'[_N]
  gen double v_st = `vv'
  gen double `eps' = rnormal()*0.4
  bysort unit (t): gen double e = `eps'[1]
  by unit (t): replace e = `rho'*e[_n-1] + sqrt(1-`rho'*`rho')*`eps'[_n] if _n > 1
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_st + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  keep y state unit t ttre coh
end

tempname M
postfile `M' int draw str4 arm str16 qty double b double se using "ar1_raw.dta", replace
forval d = 1/`NDRAWS' {
  foreach r in 0 5 {
    local rho = `r'/10
    set seed `=440000 + 29*`d''
    set sortseed `=440000 + 29*`d''
    qui argen `rho'
    cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod
    if _rc == 0 {
      foreach q in main_att ES_att_post0 ES_att_post1 ES_att_post2 ES_att_post3 {
        local bb = .
        local ss = .
        cap local bb = _b[`q']
        cap local ss = _se[`q']
        post `M' (`d') ("r`r'") ("`q'") (`bb') (`ss')
      }
    }
  }
  if mod(`d',30)==0 di "AR1_PROGRESS draw `d'"
}
postclose `M'

use "ar1_raw.dta", clear
gen double truth = .
replace truth = 1.1535714285714286 if qty=="main_att"
replace truth = 1.00 if qty=="ES_att_post0"
replace truth = 1.05 if qty=="ES_att_post1"
replace truth = 1.10 if qty=="ES_att_post2"
replace truth = 1.15 if qty=="ES_att_post3"
local tc = invttail(39,0.025)
di as txt ""
di as txt " K = 10 clusters per treated cohort-by-period cell, cluster() FINER than i."
di as txt " arm   qty                n      SD_b   mean_SE   SE/SD   cov95"
foreach a in r0 r5 {
  foreach q in main_att ES_att_post0 ES_att_post1 ES_att_post2 ES_att_post3 {
    qui count if arm=="`a'" & qty=="`q'" & se<.
    local n = r(N)
    qui su b if arm=="`a'" & qty=="`q'" & se<.
    local sd = r(sd)
    qui su se if arm=="`a'" & qty=="`q'" & se<.
    local ms = r(mean)
    qui count if arm=="`a'" & qty=="`q'" & se<. & abs(b-truth)<=`tc'*se
    di as txt " `a'   " %-16s "`q'" " " as res %4.0f `n' "  " %8.4f `sd' "  " %8.4f `ms' "  " %6.3f `ms'/`sd' "  " %6.3f r(N)/`n'
  }
}
di as txt ""
di as txt " r0 = iid errors (my mc1_ladder DGP).  r5 = AR(1) rho .5 within unit (theirs)."
di "PROBE_AR1_DONE"
