* ============================================================================
* repro_joint_finer.do -- SELF-CONTAINED.  Run with:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_joint_finer.do
* Takes about 5 minutes.  Needs only the pinned wooldid build (csdid optional).
*
* CLAIM UNDER TEST.  wooldid.sthlp:743 (jointtests) says the joint pre-trend test is
* "severely oversized" when cluster() is at or above i, prints a caution there, and names
* the fixes: "use bootstrap-based clustered standard errors or randomization inference
* (see permute()), OR CLUSTER STRICTLY FINER THAN THE COHORT-BY-PERIOD CELLS."  The runtime
* caution at wooldid.ado:1443 ends with the same clause: "...as does clustering strictly
* finer than the cohort-by-period cells."
*
* THIS FILE clusters strictly finer than the cohort-by-period cells and measures the size
* of the test against a DGP whose pre-trend is EXACTLY ZERO.
*
* DGP.  40 states x 4 units x 12 periods.  With K states per adopting wave:
*   states 1..K adopt at t=5, K+1..2K adopt at t=7, the remaining 40-2K never adopt.
*   y = s_g + v_st + a_u + .2 t + .5 cos t + tau + e
*     s_g ~ N(0,.8^2) state RE ; v_st ~ N(0,.4^2) state-by-year shock ;
*     a_u ~ N(0,.5^2) unit RE ; e ~ N(0,.4^2) ; tau = 1 + .05*rel for rel >= 0.
*   tau is ZERO in every pre-treatment period, so the joint pre-trend null is TRUE.
* States are independent, so cluster(state) is the correct level, and there are 40 of them
* in every rung.
*
* ARMS, on the SAME draws:
*   FINER  i = coh   , cluster(state)  -> K states inside each treated cohort-by-period
*                                         cell: STRICTLY FINER than the cell.  No caution.
*   AT_I   i = state , cluster(state)  -> cluster AT i: the documented-bad regime, caution.
* ============================================================================
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAWS = 80

cap program drop skgen
program define skgen
  args K
  clear
  qui set obs 40
  gen long state = _n
  gen double s_g = rnormal()*0.8
  gen double ttre = cond(state<=`K', 5, cond(state<=2*`K', 7, .))
  qui expand 4
  bysort state: gen long unit = (state-1)*4 + _n
  gen double a_u = rnormal()*0.5
  qui expand 12
  bysort unit: gen int t = _n
  tempvar fst vv
  bysort state t: gen byte `fst' = (_n==1)
  gen double `vv' = rnormal()*0.4 if `fst'
  bysort state t (`fst'): replace `vv' = `vv'[_N]
  gen double v_st = `vv'
  gen double e = rnormal()*0.4
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_st + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  gen double g   = cond(ttre<., ttre, 0)
  keep y state unit t ttre coh g
end

tempname M
postfile `M' int draw int K str7 arm double jp using "jf_raw.dta", replace
tempfile dd
foreach K in 3 10 {
  forval d = 1/`NDRAWS' {
    set seed `=770000 + 53*`d' + 11*`K''
    set sortseed `=770000 + 53*`d' + 11*`K''
    qui skgen `K'
    qui save `"`dd'"', replace
    local p1 = .
    local p2 = .
    local p3 = .
    cap qui wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
    if _rc == 0 {
      local p1 = e(joint_horiz)
      if `d' == 1 {
        di "---- rank / conditioning of the pre block of e(V_es_att), K=`K', draw 1 ----"
        cap matrix VV = e(V_es_att)
        cap matrix list VV
        cap matrix PB = VV[1..4, 1..4]
        cap matrix list PB
        cap noi mata: PBm = st_matrix("PB"); ("rank", rank(PBm)); ("cond", cond(PBm)); ("det", det(PBm))
      }
    }
    qui use `"`dd'"', clear
    cap qui wooldid y state t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
    if _rc == 0 local p2 = e(joint_horiz)
    qui use `"`dd'"', clear
    local p4 = .
    cap qui wooldid y coh t ttre, att cluster(state) esprelength(2) espostlength(4) esfixedbaseperiod jointtests
    if _rc == 0 local p4 = e(joint_horiz)
    qui use `"`dd'"', clear
    cap qui csdid y, ivar(unit) time(t) gvar(g) cluster(state) long2
    if _rc == 0 {
      cap qui estat pretrend
      if _rc == 0 local p3 = r(p)
    }
    qui use `"`dd'"', clear
    post `M' (`d') (`K') ("FINER") (`p1')
    post `M' (`d') (`K') ("AT_I")  (`p2')
    post `M' (`d') (`K') ("CSDID") (`p3')
    post `M' (`d') (`K') ("FINER2") (`p4')
    if mod(`d',50)==0 di "JF_PROGRESS K=`K' draw `d'"
  }
}
postclose `M'

use "jf_raw.dta", clear
di as txt ""
di as txt "SIZE OF THE JOINT PRE-TREND TEST.  The true pre-trend is EXACTLY ZERO, so these"
di as txt "are rejection rates under the null.  Nominal 5% and 10%."
di as txt ""
di as txt "  FINER  = i=coh, cluster(state), esprelength(4)  (4 pre horizons; pre5 is"
di as txt "           identified off ONE cohort, so it is the ragged one)"
di as txt "  FINER2 = i=coh, cluster(state), esprelength(2)  (2 pre horizons, BOTH cohorts)"
di as txt ""
di as txt "  K   arm       n    reject@5%   reject@10%   caution printed by wooldid?"
foreach K in 3 10 {
  foreach a in FINER FINER2 AT_I CSDID {
    qui count if K==`K' & arm=="`a'" & jp<.
    local n = r(N)
    qui count if K==`K' & arm=="`a'" & jp<0.05
    local r5 = r(N)/`n'
    qui count if K==`K' & arm=="`a'" & jp<0.10
    local r10 = r(N)/`n'
    local cau = cond("`a'"=="AT_I", "yes", cond("`a'"=="CSDID", "n/a", "NO"))
    di as txt "  " %2.0f `K' "  " %-7s "`a'" "  " as res %4.0f `n' "     " %6.3f `r5' "      " %6.3f `r10' "        `cau'"
  }
}
di as txt ""
di as txt "sthlp:743 names 'cluster strictly finer than the cohort-by-period cells' as a fix"
di as txt "for the at-or-above-i oversizing, and ado:1443's caution repeats the clause."
di as txt "At K=3 that geometry rejects a TRUE null at well over 50%, and prints nothing."
di "REPRO_JF_DONE"
