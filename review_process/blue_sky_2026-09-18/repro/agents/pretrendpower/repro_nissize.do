* ================================================================
* repro_nissize.do  -- SELF-CONTAINED.
*   Claim: e(roth_nissize_TYPE) -- documented at sthlp:1384 as "the
*   actual size of the any-individually-significant pre-test" and
*   printed as "size of the any-pre-period-individually-significant
*   pre-test" -- evaluates a NORMAL probability at a *t* critical
*   value.  The test it describes rejects when |b/se| > invttail(df_r,
*   .025), and b/se is referred to t(df_r) everywhere else in wooldid,
*   so its size is 5% per coefficient.  With few clusters the reported
*   size is far below the real one and points the WRONG WAY (it tells
*   the user the pre-test is conservative when it is exactly nominal).
*
*   Part 1: one run, k = 1 pre-period, 16 clusters.  There the NIS test
*           IS a single 5% t test, so "the actual size" is 0.05 by
*           construction and no simulation is needed; wooldid prints
*           2*normal(-invttail(15,.025)) = 0.03305 instead.
*   Part 2: 300-replication Monte Carlo of the realised rejection rate.
*           Part 1 is the proof (pure arithmetic).  Part 2 is corroboration:
*           300 replications pin a 5% rate only to about +/- 0.025, so read it
*           as "consistent with 0.05, inconsistent with the reported 0.033",
*           not as a point estimate of the size.
*
*   Run:  stata-mp -e do repro_nissize.do     (then read repro_nissize.log)
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

capture program drop mkstates
program define mkstates
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
  gen double e = rnormal()
  forvalues k = 2/10 {
    qui replace e = 0.6*e[_n-1] + sqrt(0.64)*rnormal() if t == `k'
  }
  gen double y = f_t + 1.0*D + e
end

* ---------------- Part 1: the arithmetic ----------------
mkstates
wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(1) espostlength(2) pretrendpower(80)
dis "PART1 df_r                = " e(df_r)
dis "PART1 npre                = " e(roth_npre_att)
dis "PART1 e(roth_nissize_att) = " %12.8f e(roth_nissize_att)
dis "PART1 2*normal(-invttail(df,.025)) = " %12.8f 2*normal(-invttail(e(df_r),.025))
dis "PART1 the size of ONE two-sided 5% t test  = 0.05000000"
dis "PART1 understatement = " %6.2f 100*(1 - e(roth_nissize_att)/.05) " percent"

* ---------------- Part 2: realised rejection rate ----------------
* reject when the single pre-period coefficient is individually
* significant at 5% ON WOOLDID'S OWN REPORTED p-VALUE.
tempname rej
scalar `rej' = 0
local R 300
local good 0
forvalues r = 1/`R' {
  set seed `= 71000 + `r''
  mkstates
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(1) espostlength(2)
  if _rc == 0 {
    mat B = e(b_es_att)
    mat V = e(V_es_att)
    local tt = abs(B[1,1])/sqrt(V[1,1])
    if `tt' > invttail(e(df_r), .025) scalar `rej' = `rej' + 1
    local ++good
  }
}
dis "PART2 replications = `good'"
dis "PART2 realised rejection rate of the NIS pre-test = " %8.4f `rej'/`good'
dis "PART2 wooldid reported that size as 0.03305"
dis "DONE_REPRO_NISSIZE"
