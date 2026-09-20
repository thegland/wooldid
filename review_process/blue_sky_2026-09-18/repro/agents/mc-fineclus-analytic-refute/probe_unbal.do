* ============================================================================
* probe_unbal.do -- SELF-CONTAINED.  Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_unbal.do
* QUESTION: wooldid.sthlp:255 tells the user to cluster BELOW the level of i.
* In a state-policy panel the only way to get there is to make `i' the ONSET
* COHORT instead of the state.  Is that a variance-only change?
* On a BALANCED panel the point estimate is unchanged (prior agents measured
* 6e-14).  This probe asks what happens on an UNBALANCED one.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 445501
set sortseed 445501

* --- 40 states, T=12.  states 1-3 adopt at t=5, 4-6 at t=7, 7-40 never. -------
* --- state level effects s_g are LARGE and permanent; tau = 1.0 exactly. ------
clear
set obs 40
gen long state = _n
gen double s_g = rnormal()*3
gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
expand 12
bysort state: gen int t = _n
gen double e = rnormal()*0.3
gen double tau = cond(ttre<. & t>=ttre, 1, 0)
gen double y = s_g + 0.2*t + tau + e
gen double coh = cond(ttre<., ttre, 99)
save "_unbal_full.dta", replace

di "=================== 1. BALANCED PANEL ==================="
foreach iv in state coh {
  use "_unbal_full.dta", clear
  cap noi qui wooldid y `iv' t ttre, att robust
  di "BAL_`iv'_RC=" _rc "  b_att=" %20.15f _b[main_att] "  se=" %12.8f _se[main_att]
  scalar b_bal_`iv' = _b[main_att]
}
di "BAL_DIFF = " %20.15f (b_bal_state - b_bal_coh)

di "=================== 2. UNBALANCED: state 1 enters the data at t=5 ========"
di "    (a wholly ordinary administrative-data shape: one unit reports late)"
foreach iv in state coh {
  use "_unbal_full.dta", clear
  qui drop if state==1 & t<5
  cap noi qui wooldid y `iv' t ttre, att robust
  di "UNB_`iv'_RC=" _rc "  b_att=" %20.15f _b[main_att] "  se=" %12.8f _se[main_att]
  scalar b_unb_`iv' = _b[main_att]
}
di "UNB_DIFF (state minus coh) = " %20.15f (b_unb_state - b_unb_coh)
di "TRUE att = 1.000000000000000"
di "BIAS i=state = " %20.15f (b_unb_state - 1)
di "BIAS i=coh   = " %20.15f (b_unb_coh   - 1)
use "_unbal_full.dta", clear
qui su s_g if state==1
di "s_g(state 1) = " %12.6f r(mean)
qui su s_g if state<=3
di "mean s_g(wave 1) = " %12.6f r(mean)

di "=================== 3. UNBALANCED: attrition instead of late entry ======="
foreach iv in state coh {
  use "_unbal_full.dta", clear
  qui drop if state==4 & t>8
  cap noi qui wooldid y `iv' t ttre, att robust
  di "ATR_`iv'_RC=" _rc "  b_att=" %20.15f _b[main_att]
  scalar b_atr_`iv' = _b[main_att]
}
di "ATR_DIFF = " %20.15f (b_atr_state - b_atr_coh)
di "ATR BIAS i=state = " %20.15f (b_atr_state - 1)
di "ATR BIAS i=coh   = " %20.15f (b_atr_coh   - 1)

di "=================== 4. BALANCED but 4 sub-units/state, ragged ============"
clear
set obs 40
gen long state = _n
gen double s_g = rnormal()*3
gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
expand 4
bysort state: gen long unit = (state-1)*4 + _n
gen double a_u = rnormal()*1.5
expand 12
bysort unit: gen int t = _n
gen double e = rnormal()*0.3
gen double tau = cond(ttre<. & t>=ttre, 1, 0)
gen double y = s_g + a_u + 0.2*t + tau + e
gen double coh = cond(ttre<., ttre, 99)
* ragged: within state 1, unit 1 leaves after t=6, unit 4 arrives at t=4
drop if unit==1 & t>6
drop if unit==4 & t<4
save "_unbal_sub.dta", replace
foreach iv in state coh {
  use "_unbal_sub.dta", clear
  cap noi qui wooldid y `iv' t ttre, att robust
  di "SUB_`iv'_RC=" _rc "  b_att=" %20.15f _b[main_att]
  scalar b_sub_`iv' = _b[main_att]
}
di "SUB_DIFF = " %20.15f (b_sub_state - b_sub_coh)
di "SUB BIAS i=state = " %20.15f (b_sub_state - 1)
di "SUB BIAS i=coh   = " %20.15f (b_sub_coh   - 1)
di "PROBE_UNBAL_DONE"
