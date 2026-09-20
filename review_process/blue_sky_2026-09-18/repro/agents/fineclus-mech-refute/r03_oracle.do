clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
*----------------------------------------------------------------------
* R03.  INDEPENDENT hand-rolled pairs-CLUSTER bootstrap of att_i, done
* entirely OUTSIDE wooldid's bootstrap engine.  Two arms:
*   ORACLE_NEST : re-identify plants by (draw copy, plant)  [nesting rule]
*   ORACLE_KEEP : keep the original plant codes             [keep-i rule]
* Same data as r01 ARM_A / ARM_B (plant 1's y missing from 2005).
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

mkdat , seed(101)
qui replace y = . if plant==1 & year>=2005
tempfile base
qui save "`base'"

* empirical between-county dispersion of the plant-level truth, for reference
preserve
  qui bysort plant: keep if _n==1
  qui collapse (mean) ptau if !missing(onset), by(county)
  qui sum ptau
  dis "EMP_SD_COUNTYMEAN_PTAU=" r(sd) "  ->  implied SE = " r(sd)/sqrt(24)
restore

qui wooldid y plant year onset , cluster(county) att att_i fastnose
getatti
dis "OBS_ATTI=" %21x `ATTI'

set seed 313131
set sortseed 313131
local R 250
tempname P
postfile `P' int rep double a_nest double a_keep using "oracle.dta", replace
forvalues r = 1/`R' {
  * ---- draw 24 counties with replacement ----
  clear
  qui set obs 24
  gen long drawid = _n
  gen long county = ceil(runiform()*24)
  tempfile dr
  qui save "`dr'"
  * ---- NEST arm ----
  use "`base'", clear
  qui joinby county using "`dr'"
  qui gegen double newplant = group(drawid plant)
  local an = .
  cap qui wooldid y newplant year onset , att_i fastnose
  if _rc == 0 {
    getatti
    local an = `ATTI'
  }
  * ---- KEEP arm : same drawn rows, original plant codes ----
  local ak = .
  cap qui wooldid y plant year onset , att_i fastnose
  if _rc == 0 {
    getatti
    local ak = `ATTI'
  }
  post `P' (`r') (`an') (`ak')
  if mod(`r',50)==0 dis "ORA_AT=`r'"
}
postclose `P'
use "oracle.dta", clear
qui sum a_nest
dis "ORACLE_NEST_SD=" r(sd) " n=" r(N) " mean=" r(mean)
qui sum a_keep
dis "ORACLE_KEEP_SD=" r(sd) " n=" r(N) " mean=" r(mean)
dis "R03_DONE"
