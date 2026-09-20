* ============================================================================
* probe_unbal2.do -- SELF-CONTAINED, ZERO NOISE, exact arithmetic.
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_unbal2.do
*
* wooldid.sthlp:255 ("What level to cluster at"): "The intended practice is to
* cluster BELOW the level of i - at a level whose units are nested inside
* cohorts, so that a treated cohort-by-period cell contains several clusters",
* and sthlp:255 again names "a finer cluster()" as THE remedy for the
* at-or-above-i standard-error understatement.
*
* In a state-policy panel where a whole wave of states adopts on the same date
* the ONLY way to reach that geometry is to hand wooldid the ONSET WAVE as `i'
* instead of the state. The help presents that as an INFERENCE choice.
*
* This file shows it is not: on an unbalanced panel it silently moves the POINT
* ESTIMATE. The DGP has NO noise at all and a homogeneous effect of EXACTLY 1,
* so any departure from 1.000000 is arithmetic, not sampling variation.
*   y = s_g + 0.2*t + 1*(t >= ttre)
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 991155
set sortseed 991155

cap program drop mkpanel
program define mkpanel
  clear
  qui set obs 40
  gen long state = _n
  * deterministic, permanent, LARGE state level differences
  gen double s_g = 2*mod(state,7) - 6
  gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
  qui expand 12
  bysort state: gen int t = _n
  gen double tau = cond(ttre<. & t>=ttre, 1, 0)
  gen double y = s_g + 0.2*t + tau          // NO error term at all
  gen double coh = cond(ttre<., ttre, 99)
  gen double g   = cond(ttre<., ttre, 0)
end

di "########## CASE 0: fully balanced (control) ##########"
foreach iv in state coh {
  qui mkpanel
  cap noi qui wooldid y `iv' t ttre, att robust
  di "BAL i=`iv'  rc=" _rc "  att = " %18.12f _b[main_att]
}

di "########## CASE 1: one treated state stops reporting after t=8 ##########"
di "   (state 4, a wave-2 adopter, attrits; every other cell is complete)"
foreach iv in state coh {
  qui mkpanel
  qui drop if state==4 & t>8
  cap noi qui wooldid y `iv' t ttre, att robust
  di "ATTRIT i=`iv'  rc=" _rc "  att = " %18.12f _b[main_att]
  scalar A_`iv' = _b[main_att]
}
di "   TRUE att = 1.000000000000 ; error(i=state) = " %18.12f (A_state-1) ///
   "  error(i=coh) = " %18.12f (A_coh-1)

di "########## CASE 2: one treated state enters the data late ##########"
foreach iv in state coh {
  qui mkpanel
  qui drop if state==1 & t<5
  cap noi qui wooldid y `iv' t ttre, att robust
  di "LATE i=`iv'  rc=" _rc "  att = " %18.12f _b[main_att]
  scalar B_`iv' = _b[main_att]
}
di "   TRUE att = 1.000000000000 ; error(i=state) = " %18.12f (B_state-1) ///
   "  error(i=coh) = " %18.12f (B_coh-1)

di "########## CASE 3: ordinary ragged administrative panel ##########"
di "   ~8% of state-years missing at random across ALL 40 states"
foreach iv in state coh {
  qui mkpanel
  set seed 4242
  qui gen double u = runiform()
  qui drop if u < 0.08
  cap noi qui wooldid y `iv' t ttre, att robust
  di "RAGGED i=`iv'  rc=" _rc "  att = " %18.12f _b[main_att]
  scalar C_`iv' = _b[main_att]
}
di "   TRUE att = 1.000000000000 ; error(i=state) = " %18.12f (C_state-1) ///
   "  error(i=coh) = " %18.12f (C_coh-1)

di "########## CASE 4: external oracles on the CASE-1 data ##########"
qui mkpanel
qui drop if state==4 & t>8
cap noi qui csdid y, ivar(state) time(t) gvar(g)
di "csdid rc=" _rc
cap noi qui estat simple
di "csdid simple ATT = " %18.12f _b[ATT]
qui mkpanel
qui drop if state==4 & t>8
qui gen byte D = (ttre<. & t>=ttre)
cap noi qui did_imputation y state t ttre
di "did_imputation rc=" _rc "  tau = " %18.12f _b[tau]
di "PROBE_UNBAL2_DONE"
