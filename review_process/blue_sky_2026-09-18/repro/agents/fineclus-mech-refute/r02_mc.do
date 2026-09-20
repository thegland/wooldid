clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
*----------------------------------------------------------------------
* R02.  CALIBRATION.  What IS the sampling variability of att_i in this
* DGP?  Two truths:
*   TRUTH_SUPER  = E[ptau] = 1        (counties are iid draws; this is what
*                                      a pairs-CLUSTER bootstrap targets)
*   TRUTH_FINITE = mean of ptau over this dataset's treated plants
* SD over MC draws of (atti - 1) is the superpopulation sampling SD.
* SD over MC draws of (atti - taubar) is pure estimation noise.
* Compare both with the two bootstrap regimes measured in r01.
*----------------------------------------------------------------------
prog def mkdat
  syntax , [seed(int 101) nc(int 24) np(int 4) tt(int 8)]
  clear
  set seed `seed'
  set sortseed `seed'
  set obs `=`nc'*`np''
  gen long plant  = _n
  gen long county = ceil(_n/`np')
  gen int onset = .
  replace onset = 2004 if mod(plant,4)==1
  replace onset = 2006 if mod(plant,4)==2
  replace onset = 2008 if mod(plant,4)==3
  by county, sort: gen double cshk = rnormal()          if _n==1
  by county:       gen double ctau = 1 + rnormal()*0.8  if _n==1
  by county:       replace cshk = cshk[1]
  by county:       replace ctau = ctau[1]
  gen double ptau = ctau + rnormal()*0.4
  gen double ai = rnormal()
  expand `tt'
  bysort plant (county): gen int year = 2000 + _n
  gen double e = rnormal()*0.3
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = ai + cshk + 0.05*(year-2001) + ptau*D + e
  sort plant year
end

* --- does fastnose post main_att_i ? ---
mkdat , seed(1)
cap noi qui wooldid y plant year onset , cluster(county) att att_i fastnose
dis "FN_RC=" _rc
mat b = e(b)
local nm : colnames b
dis "FN_COLS=`nm'"

tempname M
postfile `M' int draw double atti double attv double taubar double attihat_est using "mcres.dta", replace
local M0 400
forvalues m = 1/`M0' {
  mkdat , seed(`=2000+`m'')
  * finite-sample truth: simple mean of ptau over treated plants
  preserve
    qui bysort plant: keep if _n==1
    qui sum ptau if !missing(onset), meanonly
    local tb = r(mean)
  restore
  cap qui wooldid y plant year onset , cluster(county) att att_i fastnose
  if _rc == 0 {
    tempname b
    mat `b' = e(b)
    local nm : colnames `b'
    local K = colsof(`b')
    local ai = .
    local aa = .
    forvalues q=1/`K' {
      local this : word `q' of `nm'
      if "`this'"=="main_att_i" local ai = `b'[1,`q']
      if "`this'"=="main_att"   local aa = `b'[1,`q']
    }
    post `M' (`m') (`ai') (`aa') (`tb') (.)
  }
  if mod(`m',50)==0 dis "MC_AT=`m'"
}
postclose `M'

use "mcres.dta", clear
gen double dsuper = atti - 1
gen double dfinite = atti - taubar
gen double dsuperA = attv - 1
qui count
dis "MC_N=" r(N)
qui sum dsuper
dis "SD_SUPER_ATTI=" %21x r(sd) "  mean=" %21x r(mean) "  (dec " r(sd) ")"
qui sum dfinite
dis "SD_FINITE_ATTI=" %21x r(sd) "  mean=" %21x r(mean) "  (dec " r(sd) ")"
qui sum taubar
dis "SD_TAUBAR=" r(sd)
qui sum dsuperA
dis "SD_SUPER_ATT=" r(sd)
qui sum atti
dis "MEAN_ATTI=" r(mean) " SD_ATTI=" r(sd)
dis "R02_DONE"
