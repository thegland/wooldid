* ================================================================
* repro_designB.do -- SELF-CONTAINED (no sibling files).
*   At cluster() == i the pretrendpower() minimum detectable trend
*   LOOKS calibrated and is not: the apparent 80% power is bought
*   with a ~28% false-positive rate.
*   Step 1  null replications -> the size of the joint pre-trend test
*           and the median reported mdt_wald(80)
*   Step 2  replications carrying a trend of exactly that median ->
*           power under the nominal cutoff AND under a size-corrected
*           cutoff taken from step 1's own null distribution
*   Run: stata-mp -e do repro_designB.do ; read repro_designB.log
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local R 150

capture program drop mkpanel
program define mkpanel
  args delta
  clear
  set obs 80
  gen long unit = _n
  gen double a_i = rnormal()
  gen int onset = .
  replace onset = 5 if mod(unit,5)==1
  replace onset = 7 if mod(unit,5)==2
  replace onset = 9 if mod(unit,5)==3
  expand 10
  bysort unit: gen int t = _n
  sort unit t
  gen double f_t = 0.1*t + 0.3*sin(t)
  gen byte D = (onset<. & t>=onset)
  gen double relyr = t - onset
  gen double trendterm = cond(onset<., relyr + 1, 0)
  gen double e = rnormal()
  forvalues k = 2/10 {
    qui replace e = 0.6*e[_n-1] + sqrt(0.64)*rnormal() if t == `k'
  }
  gen double y = a_i + f_t + 1.0*D + `delta'*trendterm + e
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

file open fh using "repro_designB_null.csv", write text replace
file write fh "rep,F,mdt80" _n
file close fh
forvalues r = 1/`R' {
  set seed `= 91000 + `r''
  mkpanel 0
  cap qui wooldid y unit t onset, cluster(unit) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc == 0 {
    wstat
    file open fh using "repro_designB_null.csv", write text append
    file write fh "`r',`r(F)',`r(mdt80)'" _n
    file close fh
  }
}
preserve
  qui import delimited using "repro_designB_null.csv", clear varnames(1)
  qui count if f > invF(3,79,.95)
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
dis "STEP1 nominal F(3,79) 5% critical value = " invF(3,79,.95)
dis "STEP1 empirical 95th percentile of the null statistic = $CADJ"
dis "STEP1 median reported mdt_wald(80) = $DSTAR"

local rej 0
local rejadj 0
local good 0
forvalues r = 1/`R' {
  set seed `= 92000 + `r''
  mkpanel $DSTAR
  cap qui wooldid y unit t onset, cluster(unit) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc == 0 {
    wstat
    if r(F) > invF(3,79,.95) local ++rej
    if r(F) > $CADJ local ++rejadj
    local ++good
  }
}
dis "STEP2 replications at a trend of exactly the reported mdt_wald(80) = `good'"
dis "STEP2 power, NOMINAL cutoff       = " %8.4f `rej'/`good'    "   CLAIMED 0.80"
dis "STEP2 power, SIZE-CORRECTED cutoff = " %8.4f `rejadj'/`good'
dis "DONE_REPRO_DESIGNB"
