* ============================================================================
* mc1_ladder.do -- INDEPENDENT re-test of MC-FC-1 / MC-FC-2.
* My own generator (skgen), my own arms, my own truth. Long-format postfile.
* Arms, all on the SAME draw:
*   A  i=coh,   cluster(state), esfixedbaseperiod   [cluster STRICTLY FINER than i]
*   B  i=coh,   cluster(state), pooled reference    [phi lever: less saturation]
*   C  i=state, cluster(state), esfixedbaseperiod   [cluster AT i  = K=1 endpoint]
* K in {3,10} = number of states per treated wave (= clusters per treated cell).
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local NDRAWS = 120

cap program drop skgen
program define skgen
  args K sv sa se0 ss
  clear
  qui set obs 40
  gen long state = _n
  gen double s_g = rnormal()*`ss'
  gen double ttre = cond(state<=`K', 5, cond(state<=2*`K', 7, .))
  qui expand 4
  bysort state: gen long unit = (state-1)*4 + _n
  gen double a_u = rnormal()*`sa'
  qui expand 12
  bysort unit: gen int t = _n
  tempvar fst vv
  bysort state t: gen byte `fst' = (_n==1)
  gen double `vv' = rnormal()*`sv' if `fst'
  bysort state t (`fst'): replace `vv' = `vv'[_N]
  gen double v_st = `vv'
  gen double e = rnormal()*`se0'
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_st + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  keep y state unit t ttre coh
end

cap program drop grab
program define grab
  args M draw K arm
  local qs main_att pre_att ES_att_pre5 ES_att_pre4 ES_att_pre3 ES_att_pre2 ES_att_post0 ES_att_post1 ES_att_post2 ES_att_post3
  local jp = e(joint_horiz)
  if "`jp'" == "" local jp = .
  foreach q of local qs {
    local bb = .
    local ss = .
    cap local bb = _b[`q']
    cap local ss = _se[`q']
    post `M' (`draw') (`K') ("`arm'") ("`q'") (`bb') (`ss') (`jp')
  }
end

tempname M
postfile `M' int draw int K str4 arm str16 qty double b double se double jp ///
   using "mc1_raw.dta", replace

tempfile dd
foreach K in 3 10 {
  forval d = 1/`NDRAWS' {
    set seed `=310000 + 41*`d' + 7*`K''
    set sortseed `=310000 + 41*`d' + 7*`K''
    qui skgen `K' 0.4 0.5 0.4 0.8
    qui save `"`dd'"', replace

    cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
    if _rc == 0 grab `M' `d' `K' "A"
    qui use `"`dd'"', clear

    cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) jointtests
    if _rc == 0 grab `M' `d' `K' "B"
    qui use `"`dd'"', clear

    cap qui wooldid y state t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
    if _rc == 0 grab `M' `d' `K' "C"
    qui use `"`dd'"', clear

    if mod(`d',20)==0 di "MC1_PROGRESS K=`K' draw `d'"
  }
}
postclose `M'

* ---------------- analysis ----------------
use "mc1_raw.dta", clear
gen double truth = .
replace truth = 1.1535714285714286 if qty=="main_att"
replace truth = 0 if inlist(qty,"pre_att","ES_att_pre5","ES_att_pre4","ES_att_pre3","ES_att_pre2")
replace truth = 1.00 if qty=="ES_att_post0"
replace truth = 1.05 if qty=="ES_att_post1"
replace truth = 1.10 if qty=="ES_att_post2"
replace truth = 1.15 if qty=="ES_att_post3"
local tc = invttail(39,0.025)
gen byte cov = (abs(b-truth) <= `tc'*se) if se<.

di as txt ""
di as txt "  K  arm  qty                 n    mean_b     SD_b   mean_SE   SE/SD   cov95"
levelsof K, local(Ks)
foreach K of local Ks {
  foreach a in A B C {
    foreach q in main_att pre_att ES_att_pre5 ES_att_pre4 ES_att_pre3 ES_att_pre2 ES_att_post0 ES_att_post1 ES_att_post2 ES_att_post3 {
      qui count if K==`K' & arm=="`a'" & qty=="`q'" & se<.
      local n = r(N)
      if `n' > 0 {
        qui su b if K==`K' & arm=="`a'" & qty=="`q'" & se<.
        local mb = r(mean)
        local sd = r(sd)
        qui su se if K==`K' & arm=="`a'" & qty=="`q'" & se<.
        local ms = r(mean)
        qui su cov if K==`K' & arm=="`a'" & qty=="`q'"
        local cv = r(mean)
        di as txt "  " %2.0f `K' "  `a'   " %-16s "`q'" "  " as res %4.0f `n' ///
           "  " %8.4f `mb' " " %8.4f `sd' "  " %8.4f `ms' "  " %6.3f `ms'/`sd' "  " %6.3f `cv'
      }
    }
  }
}
di as txt ""
di as txt "JOINT PRE-TREND TEST SIZE (true pretrend is exactly zero); nominal 5%:"
di as txt "  K  arm   n   rej05   rej10"
foreach K of local Ks {
  foreach a in A B C {
    qui count if K==`K' & arm=="`a'" & qty=="main_att" & jp<.
    local n = r(N)
    if `n' > 0 {
      qui count if K==`K' & arm=="`a'" & qty=="main_att" & jp<0.05
      local r5 = r(N)
      qui count if K==`K' & arm=="`a'" & qty=="main_att" & jp<0.10
      local r10 = r(N)
      di as txt "  " %2.0f `K' "  `a'  " as res %4.0f `n' "  " %6.3f `r5'/`n' "  " %6.3f `r10'/`n'
    }
  }
}
di as txt ""
di as txt "sqrt((K-1)/K): K=3 -> " %6.4f sqrt(2/3) "   K=10 -> " %6.4f sqrt(9/10)
di "MC1_DONE"
