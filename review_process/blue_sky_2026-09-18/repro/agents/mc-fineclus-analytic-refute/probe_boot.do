* ============================================================================
* probe_boot.do -- sthlp:255 lists THREE remedies for the at-or-above-i
* understatement: "a finer cluster(), or bootstrap() or permute()".  The seeded
* finding MC-FC-1 tested only the first and concluded the remedy is incomplete
* at K = 3 clusters per treated cell.  Does remedy #2 close the gap there?
*
* DGP: 40 states, ONE row per state-year, T=12.  states 1-3 adopt at t=5,
* 4-6 at t=7, 7-40 never.  y = s_g + v_st + .2t + .5cos t + tau + e.
* i = coh (cluster STRICTLY FINER than i), cluster(state), K = 3 states in each
* treated cohort-by-period cell.  True att = 1.153571428571...
* Arms on the SAME draws: analytic clustered SE, and bootstrap(standard).
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAWS  = 60
local NBOOT   = 30
local REPS    = 99

cap program drop bgen
program define bgen
  clear
  qui set obs 40
  gen long state = _n
  gen double s_g = rnormal()*0.8
  gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
  qui expand 12
  bysort state: gen int t = _n
  gen double v_st = rnormal()*0.4
  gen double e = rnormal()*0.4
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_st + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  keep y state t ttre coh
end

tempname M
postfile `M' int draw double b double se_an double se_bt double lb double ub double nrep ///
   using "boot_raw.dta", replace
tempfile dd
timer clear
forval d = 1/`NDRAWS' {
  set seed `=880000 + 37*`d''
  set sortseed `=880000 + 37*`d''
  qui bgen
  qui save `"`dd'"', replace
  local bb = .
  local sa = .
  local sb = .
  local lb = .
  local ub = .
  local nr = .
  timer on 1
  cap qui wooldid y coh t ttre, att cluster(state)
  timer off 1
  if _rc == 0 {
    local bb = _b[main_att]
    local sa = _se[main_att]
  }
  if `d' <= `NBOOT' {
    qui use `"`dd'"', clear
    timer on 2
    cap qui wooldid y coh t ttre, att cluster(state) bootstrap(standard) reps(`REPS') rseed(`=3000+`d'')
    local brc = _rc
    timer off 2
    if `brc' == 0 {
      local sb = _se[main_att]
      tempname R
      cap matrix `R' = e(wdidmainresults)
      cap local lb = `R'[rownumb(`R',"main_att"), colnumb(`R',"lb_95ci")]
      cap local ub = `R'[rownumb(`R',"main_att"), colnumb(`R',"ub_95ci")]
      cap local nr = e(N_reps)
    }
    di "BOOTDRAW `d' rc=`brc' se_an=`sa' se_bt=`sb'"
    if `d' == 1 & `brc' == 0 {
      di "---- e() after first bootstrap draw ----"
      ereturn list
    }
  }
  post `M' (`d') (`bb') (`sa') (`sb') (`lb') (`ub') (`nr')
  qui use `"`dd'"', clear
}
postclose `M'
timer list

use "boot_raw.dta", clear
local TRUTH = 1.1535714285714286
local tc = invttail(39,0.025)
qui su b if se_an<.
di as txt "ALL `=r(N)' draws: empirical SD of b = " as res %8.5f r(sd)
local sdall = r(sd)
qui su b if se_bt<.
local sdb = r(sd)
di as txt "Bootstrapped subset (`=r(N)' draws): SD of b = " as res %8.5f `sdb'
qui su se_an if se_an<.
di as txt "  mean analytic SE (all draws)     = " as res %8.5f r(mean) as txt "   SE/SD = " as res %6.3f r(mean)/`sdall'
qui su se_an if se_bt<.
local msa = r(mean)
qui su se_bt if se_bt<.
local msb = r(mean)
di as txt "  on the bootstrapped subset: mean analytic SE = " as res %8.5f `msa' as txt "  -> SE/SD = " as res %6.3f `msa'/`sdb'
di as txt "                              mean bootstrap SE = " as res %8.5f `msb' as txt "  -> SE/SD = " as res %6.3f `msb'/`sdb'
qui count if se_an<. & abs(b-`TRUTH')<=`tc'*se_an
qui count if se_bt<.
local nb = r(N)
qui count if se_bt<. & abs(b-`TRUTH')<=`tc'*se_an
di as txt "  coverage (analytic, on boot subset)  = " as res %6.3f r(N)/`nb'
qui count if se_bt<. & abs(b-`TRUTH')<=1.959963985*se_bt
di as txt "  coverage (bootstrap SE, normal)      = " as res %6.3f r(N)/`nb'
qui count if lb<. & `TRUTH'>=lb & `TRUTH'<=ub
di as txt "  coverage (bootstrap printed 95% CI)  = " as res %6.3f r(N)/`nb'
qui su nrep
di as txt "  completed reps: mean " as res r(mean) as txt " min " as res r(min)
di "PROBE_BOOT_DONE"
