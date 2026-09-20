* ================================================================
* repro_designE.do -- SELF-CONTAINED (no sibling files).
*   On a 16-state policy panel -- i = onset-cohort code with state
*   effects through fe() (the design sthlp:1334 recommends for
*   pretrendpower) and cluster(state) -- wooldid prints NO caution,
*   the joint pre-trend test is twice its nominal size, and the
*   reported minimum detectable trend does not have 80% power.
*   Step 1  null replications -> the size of the joint pre-trend test
*           and the median reported mdt_wald(80)
*   Step 2  replications carrying a trend of exactly that median ->
*           power under the nominal cutoff AND under a size-corrected
*           cutoff taken from step 1's own null distribution
*   Run: stata-mp -e do repro_designE.do ; read repro_designE.log
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local R 200

capture program drop mkpanel
program define mkpanel
  args delta
  clear
  set obs 16
  gen long st = _n
  gen byte coh = 4
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
  gen double trendterm = cond(ttre<., relyr + 1, 0)
  gen double e = rnormal()
  forvalues k = 2/10 {
    qui replace e = 0.6*e[_n-1] + sqrt(0.64)*rnormal() if t == `k'
  }
  gen double y = f_t + 1.0*D + `delta'*trendterm + e
end

* the statistic wooldid_ptp itself computes, recomputed from the posted
* e(b_es_att)/e(V_es_att)/e(df_r) (verified to 1-2 ulp against e(roth_*))
capture program drop wstat
program define wstat, rclass
  tempname B V S bp Si Q
  matrix `B' = e(b_es_att)
  matrix `V' = e(V_es_att)
  matrix `S' = `V'[1..3,1..3]
  matrix `S' = (`S' + `S'')/2
  matrix `bp' = (`B'[1,1] \ `B'[1,2] \ `B'[1,3])
  matrix `Si' = invsym(`S')
  matrix `Q' = `bp''*`Si'*`bp'
  local keff = 3 - diag0cnt(`Si')
  return scalar F = `Q'[1,1]/`keff'
  return scalar keff = `keff'
  matrix `Q' = (-3 \ -2 \ -1)'*`Si'*(-3 \ -2 \ -1)
  return scalar mdt80 = sqrt(npnF(`keff', e(df_r), invF(`keff', e(df_r), .95), .2)/`Q'[1,1])
end

file open fh using "repro_designE_null.csv", write text replace
file write fh "rep,F,mdt80" _n
file close fh
forvalues r = 1/`R' {
  set seed `= 51000 + `r''
  mkpanel 0
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc == 0 {
    wstat
    file open fh using "repro_designE_null.csv", write text append
    file write fh "`r',`r(F)',`r(mdt80)'" _n
    file close fh
  }
}
preserve
  qui import delimited using "repro_designE_null.csv", clear varnames(1)
  qui count if f > invF(3,15,.95)
  local nrej = r(N)
  qui count
  local n0 = r(N)
  qui sum mdt80, detail
  global DSTAR = r(p50)
  qui _pctile f, p(95)
  global CADJ = r(r1)
restore
dis "STEP1 null replications = `n0'"
dis "STEP1 size of the joint pre-trend test at 5% = " %8.4f `nrej'/`n0'   "   (nominal .05)"
dis "STEP1 nominal F(3,15) 5% critical value = " invF(3,15,.95)
dis "STEP1 empirical 95th percentile of the null statistic = $CADJ"
dis "STEP1 median reported mdt_wald(80) = $DSTAR"

local rej 0
local rejadj 0
local good 0
forvalues r = 1/`R' {
  set seed `= 52000 + `r''
  mkpanel $DSTAR
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc == 0 {
    wstat
    if r(F) > invF(3,15,.95) local ++rej
    if r(F) > $CADJ local ++rejadj
    local ++good
  }
}
dis "STEP2 replications at a trend of exactly the reported mdt_wald(80) = `good'"
dis "STEP2 power, NOMINAL cutoff       = " %8.4f `rej'/`good'    "   CLAIMED 0.80"
dis "STEP2 power, SIZE-CORRECTED cutoff = " %8.4f `rejadj'/`good'
dis "DONE_REPRO_DESIGNE"
