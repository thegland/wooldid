* ================================================================
* repro_power.do -- SELF-CONTAINED demonstration that a pretrendpower()
*   minimum detectable trend does (or does not) have the stated power.
*
*   Design: 16-state policy panel, 10 years, state-level outcome,
*           AR(1) errors within state.  i = onset-cohort code,
*           fe(state), cluster(state).  The cluster is STRICTLY FINER
*           than i, so wooldid prints no caution on this geometry.
*
*   Step 1  fit once with pretrendpower(80); read off mdt_wald(80).
*   Step 2  simulate REPS panels that carry a linear pre-trend of
*           exactly that slope, and count how often the joint
*           pre-trend test rejects at 5%.
*   The rejection rate should be 0.80.
*
*   Set  local REPS  smaller for a quick look; 300 gives a Clopper-
*   Pearson half-width of about 0.045 at p = 0.8.
*   Run: stata-mp -e do repro_power.do  ; read repro_power.log
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local REPS 300

capture program drop mkstates
program define mkstates
  args delta
  clear
  set obs 16
  gen long st = _n
  gen byte coh = 4   // never-treated cohort code; wooldid requires positive i
  replace coh = 1 if st <= 4
  replace coh = 2 if st > 4 & st <= 8
  replace coh = 3 if st > 8 & st <= 12
  gen int ttre = .
  replace ttre = 5 if coh == 1
  replace ttre = 7 if coh == 2
  replace ttre = 9 if coh == 3
  expand 10
  bysort st: gen int t = _n
  sort st t
  gen double f_t = 0.1*t + 0.3*sin(t)
  gen byte D = (ttre<. & t>=ttre)
  gen double relyr = t - ttre
  * a linear pre-trend of slope `delta' in relative time, treated units only:
  * it makes E[b_es(pre r)] = delta*(r - reference) exactly
  gen double trendterm = cond(ttre<., relyr + 1, 0)
  gen double e = rnormal()
  forvalues k = 2/10 {
    qui replace e = 0.6*e[_n-1] + sqrt(0.64)*rnormal() if t == `k'
  }
  gen double y = f_t + 1.0*D + `delta'*trendterm + e
end

* ---------- Step 1: what does wooldid say is detectable at 80%? ----------
* ONE fit first, exactly what a user sees ...
set seed 20260917
set sortseed 20260917
mkstates 0
wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
mat RP = e(roth_power_att)
dis "STEP1a one fit: mdt_wald(80) = " RP[1,2] "   (df_r = " e(df_r) ")"
* ... but the reported MDT is proportional to a 15-df variance estimate, so it
* is itself noisy (CV ~ 23% on this design).  Target the MEDIAN over 25 fits so
* the demonstration is about the SYSTEMATIC shortfall, not one unlucky draw.
tempname mm
matrix `mm' = J(25, 1, .)
forvalues r = 1/25 {
  set seed `= 61000 + `r''
  mkstates 0
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
  if _rc == 0 {
    mat RR = e(roth_power_att)
    matrix `mm'[`r', 1] = RR[1,2]
  }
}
preserve
  clear
  qui svmat double `mm', names(m)
  qui sum m1, detail
  global DSTAR = r(p50)
  dis "STEP1b median mdt_wald(80) over 25 fits = " r(p50) "  (min " r(min) ", max " r(max) ")"
restore
dis "STEP1 target slope = $DSTAR"

* ---------- Step 2: simulate at exactly that slope ----------
local rej 0
local rejnis 0
local good 0
forvalues r = 1/`REPS' {
  set seed `= 81000 + `r''
  mkstates $DSTAR
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc == 0 {
    mat B = e(b_es_att)
    mat V = e(V_es_att)
    mat S = V[1..3,1..3]
    mat S = (S + S')/2
    mat bp = (B[1,1] \ B[1,2] \ B[1,3])
    mat Si = invsym(S)
    mat Q = bp'*Si*bp
    local keff = 3 - diag0cnt(Si)
    local p = Ftail(`keff', e(df_r), Q[1,1]/`keff')
    if `p' < .05 local ++rej
    local tc = invttail(e(df_r), .025)
    local hit 0
    forvalues j = 1/3 {
      if abs(B[1,`j'])/sqrt(S[`j',`j']) > `tc' local hit 1
    }
    local rejnis = `rejnis' + `hit'
    local ++good
  }
}
dis "STEP2 replications with a trend of exactly mdt_wald(80) = `good'"
dis "STEP2 joint pre-trend test rejects at 5%: " %8.4f `rej'/`good'  "   CLAIMED POWER = 0.80"
dis "STEP2 NIS pre-test rejects:               " %8.4f `rejnis'/`good'
dis "DONE_REPRO_POWER"
