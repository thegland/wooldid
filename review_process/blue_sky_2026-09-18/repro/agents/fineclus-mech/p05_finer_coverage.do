clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"

*----------------------------------------------------------------------
* P05 -- SIZE of the NEW finer-than-i bootstrap arm (keep-i regime).
*
* i  = cohort group (12 values)   <- the estimand definition
* cluster() = firm (96 values)    <- STRICTLY FINER than i
* True effect = 0 everywhere, so the true att_i is 0.
* Dependence is at the FIRM level (a firm random effect constant over
* time), so cluster(firm) is the correct level and a correctly
* implemented cluster bootstrap must reject H0: att_i = 0 about 5% of
* the time at the 5% level.
*----------------------------------------------------------------------
prog def getatti
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  local j 0
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if "`this'" == "main_att_i" local j `q'
  }
  if `j' == 0 {
    c_local ATTI = .
    c_local ATTISE = .
  }
  else {
    c_local ATTI = `b'[1,`j']
    c_local ATTISE = sqrt(`V'[`j',`j'])
  }
end

local NMC   80
local REPS 79

tempname RES
mat `RES' = J(`NMC', 4, .)

timer clear 9
timer on 9
forvalues m = 1/`NMC' {
  clear
  set seed `=5000 + `m''
  set sortseed `=5000 + `m''
  set obs 96
  gen long firm = _n
  gen long cgrp = ceil(_n/8)
  gen int onset = .
  replace onset = 2004 if cgrp <= 3
  replace onset = 2006 if cgrp >= 4 & cgrp <= 6
  replace onset = 2008 if cgrp >= 7 & cgrp <= 9
  gen double fe_ = rnormal()*1.0
  expand 8
  bysort firm (cgrp): gen int year = 2001 + _n - 1
  gen double e = rnormal()*0.5
  gen byte D = (year >= onset) & !missing(onset)
  * TRUE EFFECT ZERO
  gen double y = fe_ + 0.05*(year-2001) + 0*D + e

  cap qui wooldid y cgrp year onset, att_i cluster(firm)
  if _rc == 0 {
    getatti
    mat `RES'[`m',1] = `ATTI'
    mat `RES'[`m',2] = `ATTISE'
  }
  cap qui wooldid y cgrp year onset, att_i cluster(firm) bootstrap(standard) reps(`REPS') rseed(`=900+`m'') dots(0)
  if _rc == 0 {
    getatti
    mat `RES'[`m',3] = `ATTI'
    mat `RES'[`m',4] = `ATTISE'
  }
  if mod(`m',10) == 0 {
    timer off 9
    timer list 9
    dis "P05 progress m=`m'"
    timer on 9
  }
}
timer off 9

clear
qui svmat double `RES'
rename `RES'1 b_an
rename `RES'2 se_an
rename `RES'3 b_bs
rename `RES'4 se_bs
gen byte rej_an = abs(b_an/se_an) > invnormal(0.975) if !missing(se_an)
gen byte rej_bs = abs(b_bs/se_bs) > invnormal(0.975) if !missing(se_bs)
qui sum b_an
dis "P05 TRUE_SD_OF_b = " %12.8f r(sd) "  (n=" r(N) ")"
qui sum se_an
dis "P05 MEAN_ANALYTIC_SE = " %12.8f r(mean)
qui sum se_bs
dis "P05 MEAN_BOOTSTRAP_SE = " %12.8f r(mean)
qui sum rej_an
dis "P05 ANALYTIC REJECTION RATE = " %8.4f r(mean) "  n=" r(N)
qui sum rej_bs
dis "P05 BOOTSTRAP REJECTION RATE = " %8.4f r(mean) "  n=" r(N)
save "p05_results.dta", replace
dis "ALLDONE"
