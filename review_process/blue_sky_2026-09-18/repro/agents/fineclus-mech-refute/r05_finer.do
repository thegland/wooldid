clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
*----------------------------------------------------------------------
* R05.  NEW PROBE (not seeded).  Is the keep-i regime CALIBRATED where it
* is the design's INTENDED rule -- a cluster genuinely FINER than i?
* i = cohort group (12), cluster = firm (96, 8 per group).  Treatment
* effects are iid ACROSS FIRMS with no group component, so the firm-level
* cluster bootstrap can in principle capture all of the sampling variance
* of att_i.  Compare mean bootstrap se(att_i) with the Monte Carlo SD.
*----------------------------------------------------------------------
prog def mkdat2
  syntax , [seed(int 1) ng(int 12) nf(int 8) tt(int 8)]
  clear
  set seed `seed'
  set sortseed `seed'
  set obs `=`ng'*`nf''
  gen long firm = _n
  gen long cgrp = ceil(_n/`nf')
  gen int onset = .
  replace onset = 2004 if inrange(cgrp,1,3)
  replace onset = 2006 if inrange(cgrp,4,6)
  replace onset = 2008 if inrange(cgrp,7,9)
  gen double tau = 1 + rnormal()*0.6
  gen double ai  = rnormal()
  expand `tt'
  bysort firm (cgrp): gen int year = 2000 + _n
  gen double e = rnormal()*0.3
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = ai + 0.05*(year-2001) + tau*D + e
  sort firm year
end
prog def getatti
  tempname b
  mat `b' = e(b)
  local nm : colnames `b'
  local K = colsof(`b')
  local v = .
  forvalues q=1/`K' {
    local this : word `q' of `nm'
    if "`this'"=="main_att_i" local v = `b'[1,`q']
  }
  c_local ATTI = `v'
end

* confirm the REGIME the code would pick (cluster finer than i -> keep-i)
mkdat2 , seed(1)
tempvar mn mx
qui gegen double `mn' = min(firm), by(cgrp)
qui gegen double `mx' = max(firm), by(cgrp)
cap assert `mn' == `mx'
dis "FINER_PREDNEST=" (_rc==0)

*--- Monte Carlo sampling SD of att_i ---
tempname M
postfile `M' int m double atti double taubar using "mc_finer.dta", replace
forvalues m = 1/300 {
  mkdat2 , seed(`=7000+`m'')
  preserve
    qui bysort firm: keep if _n==1
    qui sum tau if !missing(onset), meanonly
    local tb = r(mean)
  restore
  cap qui wooldid y cgrp year onset , cluster(firm) att att_i fastnose
  if _rc == 0 {
    getatti
    post `M' (`m') (`ATTI') (`tb')
  }
  if mod(`m',100)==0 dis "MCF_AT=`m'"
}
postclose `M'
use "mc_finer.dta", clear
qui count
dis "MCF_N=" r(N)
gen double dsup = atti - 1
gen double dfin = atti - taubar
qui sum dsup
dis "FINER_SD_SUPER=" r(sd)
qui sum dfin
dis "FINER_SD_FINITE=" r(sd)
qui sum taubar
dis "FINER_SD_TAUBAR=" r(sd)

*--- bootstrap SE on 12 fresh datasets ---
tempname B
postfile `B' int d double se_boot double se_ana double atti using "bs_finer.dta", replace
forvalues d = 1/12 {
  mkdat2 , seed(`=7000+`d'')
  cap qui wooldid y cgrp year onset , cluster(firm) att att_i
  local sa = .
  if _rc == 0 {
     tempname b V
     mat `b' = e(b)
     mat `V' = e(V)
     local nm : colnames `b'
     local K = colsof(`b')
     forvalues q=1/`K' {
        local this : word `q' of `nm'
        if "`this'"=="main_att_i" local sa = sqrt(`V'[`q',`q'])
     }
  }
  cap qui wooldid y cgrp year onset , cluster(firm) att att_i bootstrap(standard) reps(150) rseed(`=500+`d'') dots(0)
  local sb = .
  local ai = .
  if _rc == 0 {
     tempname b2 V2
     mat `b2' = e(b)
     mat `V2' = e(V)
     local nm : colnames `b2'
     local K = colsof(`b2')
     forvalues q=1/`K' {
        local this : word `q' of `nm'
        if "`this'"=="main_att_i" {
           local sb = sqrt(`V2'[`q',`q'])
           local ai = `b2'[1,`q']
        }
     }
  }
  post `B' (`d') (`sb') (`sa') (`ai')
  dis "BSF_d`d'_seboot=`sb'_seana=`sa'"
}
postclose `B'
use "bs_finer.dta", clear
qui sum se_boot
dis "FINER_MEAN_SEBOOT=" r(mean) " n=" r(N)
qui sum se_ana
dis "FINER_MEAN_SEANA=" r(mean)
dis "R05_DONE"
