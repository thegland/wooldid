* ================================================================
* p04_simA.do -- DOES pretrendpower()'s minimum detectable trend
*                ACTUALLY HAVE THE STATED POWER?   DESIGN A.
* Balanced staggered panel, iid homoskedastic normal errors,
* default (OLS/"ols") VCE -- the geometry where the joint pre-trend
* F test should be EXACT, so this is the cleanest possible test of
* the mdt_wald column.
*   arm 0: delta = 0            -> size of the pre-test, and the
*                                  mdt the package reports
*   arm 1: delta = median mdt80 -> achieved power should be .80
*   arm 2: delta = median mdt50 -> achieved power should be .50
* Results appended line-by-line to CSVs so a kill leaves evidence.
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

capture program drop runarm
program define runarm
  syntax , delta(real) reps(int) out(string) [clus(string)]
  file open fh using "`out'", write text replace
  file write fh "rep,rc,dfr,slope,se,W,p,mdtw50,mdtw80,nisrej,b1,b2,b3,s11,s22,s33,s12,s13,s23" _n
  file close fh
  forvalues r = 1/`reps' {
    set seed `=20260000 + `r''
    qui replace y = a_i + f_t + 1.0*D + `delta'*trendterm + rnormal()
    cap qui wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) `clus'
    local rc = _rc
    local dfr = .
    local slope = .
    local se = .
    local W = .
    local p = .
    local m50 = .
    local m80 = .
    local nis = .
    local b1 = .
    local b2 = .
    local b3 = .
    local s11 = .
    local s22 = .
    local s33 = .
    local s12 = .
    local s13 = .
    local s23 = .
    if `rc' == 0 {
      local dfr = e(df_r)
      mat B = e(b_es_att)
      mat V = e(V_es_att)
      mat S = V[1..3,1..3]
      mat S = (S + S')/2
      mat bp = (B[1,1] \ B[1,2] \ B[1,3])
      mat mm = (-3 \ -2 \ -1)
      mat Si = invsym(S)
      mat Q1 = mm'*Si*mm
      mat Q2 = mm'*Si*bp
      mat Q3 = bp'*Si*bp
      local c2 = Q1[1,1]
      local slope = Q2[1,1]/`c2'
      local se = 1/sqrt(`c2')
      local W = Q3[1,1]
      local p = Ftail(3,`dfr',`W'/3)
      local m50 = sqrt(npnF(3,`dfr',invF(3,`dfr',.95),1-.50)/`c2')
      local m80 = sqrt(npnF(3,`dfr',invF(3,`dfr',.95),1-.80)/`c2')
      local tc = invttail(`dfr',.025)
      local nis = 0
      forvalues j = 1/3 {
        if abs(B[1,`j'])/sqrt(S[`j',`j']) > `tc' local nis = 1
      }
      local b1 = B[1,1]
      local b2 = B[1,2]
      local b3 = B[1,3]
      local s11 = S[1,1]
      local s22 = S[2,2]
      local s33 = S[3,3]
      local s12 = S[1,2]
      local s13 = S[1,3]
      local s23 = S[2,3]
    }
    file open fh using "`out'", write text append
    file write fh "`r',`rc',`dfr',`slope',`se',`W',`p',`m50',`m80',`nis',`b1',`b2',`b3',`s11',`s22',`s33',`s12',`s13',`s23'" _n
    file close fh
  }
end

* ---- fixed design ----
set obs 80
gen long unit = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 5 if mod(unit,5)==1
replace onset = 7 if mod(unit,5)==2
replace onset = 9 if mod(unit,5)==3
expand 10
bysort unit: gen int t = _n
gen double f_t = 0.1*t + 0.3*sin(t)
gen byte D = (onset<. & t>=onset)
gen double relyr = t - onset
gen double trendterm = cond(onset<., relyr + 1, 0)
gen double y = 0

local R 300

timer clear
timer on 1
runarm , delta(0) reps(`R') out("A0_null.csv")
timer off 1
dis "ARM0_DONE"
timer list 1

preserve
  qui import delimited using "A0_null.csv", clear varnames(1)
  qui sum mdtw80, detail
  global D80 = r(p50)
  qui sum mdtw50, detail
  global D50 = r(p50)
restore
dis "MEDIAN_MDT80 = $D80"
dis "MEDIAN_MDT50 = $D50"

timer on 2
runarm , delta($D80) reps(`R') out("A1_pow80.csv")
timer off 2
dis "ARM1_DONE"
timer list 2

timer on 3
runarm , delta($D50) reps(`R') out("A2_pow50.csv")
timer off 3
dis "ARM2_DONE"
timer list 3
dis "ALL_DONE"
