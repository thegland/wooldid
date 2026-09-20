* ============================================================================
* mc2_noau.do -- MECHANISM TEST, CLEAN CONTROL for mc2_mech.do: the unit random effect
* is removed (a_u = 0), so cluster(unit) is exactly the right level and arm F should be
* correctly sized. Then E/F isolates the (Kcell-1)/Kcell loss with no offsetting bias.
* Original header follows.
* mc2_mech.do -- MECHANISM TEST.  DGP has NO state-level dependence whatever:
* units are independent, y = a_u + 0.2t + .5cos t + tau + e.  `state' is only a
* grouping label.  Three arms on the SAME draws, all with the SAME point est:
*   E  i=coh,  cluster(state)  -> treated (coh,rel) cell holds  3 clusters  (OVER-clustered)
*   F  i=coh,  cluster(unit)   -> treated (coh,rel) cell holds 12 clusters  (correct level)
*   G  i=unit, cluster(unit)   -> cluster AT i                              (K=1)
* If the seeded mechanism is right, SE/SD ~ sqrt((Kcell-1)/Kcell):
*   E -> 0.816,  F -> 0.957,  G -> far below.
* That would mean OVER-clustering (normally conservative) is ANTI-conservative here.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAWS = 120

cap program drop ugen
program define ugen
  clear
  qui set obs 40
  gen long state = _n
  gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
  qui expand 4
  bysort state: gen long unit = (state-1)*4 + _n
  gen double a_u = 0   // NO unit random effect: both cluster arms are then clean
  qui expand 12
  bysort unit: gen int t = _n
  gen double e = rnormal()*0.6
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  keep y state unit t ttre coh
end

cap program drop grab2
program define grab2
  args M draw arm nclus
  local qs main_att ES_att_post0 ES_att_post1 ES_att_post2 ES_att_post3 ES_att_pre3 ES_att_pre2
  foreach q of local qs {
    local bb = .
    local ss = .
    cap local bb = _b[`q']
    cap local ss = _se[`q']
    post `M' (`draw') ("`arm'") ("`q'") (`bb') (`ss') (`nclus')
  }
end

tempname M
postfile `M' int draw str4 arm str16 qty double b double se double nclus using "mc2b_raw.dta", replace
tempfile dd
forval d = 1/`NDRAWS' {
  set seed `=520000 + 61*`d''
  set sortseed `=520000 + 61*`d''
  qui ugen
  qui save `"`dd'"', replace
  cap qui wooldid y coh t ttre, att cluster(state) esprelength(3) espostlength(4) esfixedbaseperiod
  if _rc == 0 grab2 `M' `d' "E" `=e(N_clust)'
  qui use `"`dd'"', clear
  cap qui wooldid y coh t ttre, att cluster(unit) esprelength(3) espostlength(4) esfixedbaseperiod
  if _rc == 0 grab2 `M' `d' "F" `=e(N_clust)'
  qui use `"`dd'"', clear
  cap qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(4) esfixedbaseperiod
  if _rc == 0 grab2 `M' `d' "G" `=e(N_clust)'
  qui use `"`dd'"', clear
  if mod(`d',20)==0 di "MC2B_PROGRESS draw `d'"
}
postclose `M'

use "mc2b_raw.dta", clear
gen double truth = .
replace truth = 1.1535714285714286 if qty=="main_att"
replace truth = 1.00 if qty=="ES_att_post0"
replace truth = 1.05 if qty=="ES_att_post1"
replace truth = 1.10 if qty=="ES_att_post2"
replace truth = 1.15 if qty=="ES_att_post3"
replace truth = 0 if inlist(qty,"ES_att_pre3","ES_att_pre2")
di as txt ""
di as txt " arm  qty              nclus    n     SD_b   mean_SE   SE/SD   cov95"
foreach a in E F G {
  qui su nclus if arm=="`a'"
  local nc = r(mean)
  local tc = invttail(`nc'-1, 0.025)
  foreach q in main_att ES_att_post0 ES_att_post1 ES_att_post2 ES_att_post3 ES_att_pre3 ES_att_pre2 {
    qui count if arm=="`a'" & qty=="`q'" & se<.
    local n = r(N)
    if `n' > 0 {
      qui su b if arm=="`a'" & qty=="`q'" & se<.
      local sd = r(sd)
      qui su se if arm=="`a'" & qty=="`q'" & se<.
      local ms = r(mean)
      qui count if arm=="`a'" & qty=="`q'" & se<. & abs(b-truth) <= `tc'*se
      local cv = r(N)/`n'
      di as txt " `a'   " %-16s "`q'" " " as res %5.0f `nc' " " %4.0f `n' " " %8.4f `sd' "  " %8.4f `ms' "  " %6.3f `ms'/`sd' "  " %6.3f `cv'
    }
  }
}
di as txt ""
di as txt "PREDICTIONS if 'clusters inside the saturated cell' is the driver:"
di as txt "  E (3 clusters in cell)  -> " %6.3f sqrt(2/3)
di as txt "  F (12 clusters in cell) -> " %6.3f sqrt(11/12)
di "MC2B_DONE"
