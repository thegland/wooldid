clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
*----------------------------------------------------------------------
* R01. (a) REGIME ISOLATION: d0m vs d1 have bit-identical estimation
*          samples, point estimates and analytic SEs (proved in r00).
*          Only the regime flag differs.  Does the bootstrap SE move?
*      (b) DOSE: d0 (0 straddlers) vs d2 (1 real straddler) vs
*          k=8 and k=24 straddlers.
*      (c) DEAD CLUSTER: e(df_r) vs e(N_clust) under bootdf(clusters).
*----------------------------------------------------------------------
prog def mkdat
  syntax , [seed(int 101) nc(int 24) np(int 4) tt(int 8) k(int 0)]
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
  if `k' > 0 {
    forvalues q = 1/`k' {
      qui replace county = cond(`q'==`nc', 1, `q'+1) if plant==`=(`q'-1)*`np'+1' & year>=2005
    }
  }
  sort plant year
end
prog def show
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  forvalues q=1/`K' {
     local this : word `q' of `nm'
     foreach w in main_att main_att_i main_att_it main_att_itime {
       if "`this'"=="`w'" dis "`tag'|`w'|b=" %21x `b'[1,`q'] "|se=" %21x sqrt(`V'[`q',`q'])
     }
  }
  dis "`tag'|N=" e(N) "|NCLUST=" e(N_clust) "|DFR=" e(df_r) "|DONE=" e(boot_completed)
end

local RP 250
local RS 606060

* ---- (a) regime isolation ----
mkdat , seed(101)
qui replace y = . if plant==1 & year>=2005
tempfile A
qui save "`A'"
cap noi wooldid y plant year onset , cluster(county) att att_i att_it att_itime bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
dis "ARMA_RC=" _rc
show , tag(ARM_A_NEST)

use "`A'", clear
qui replace county = 2 if plant==1 & year>=2005
cap noi wooldid y plant year onset , cluster(county) att att_i att_it att_itime bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
dis "ARMB_RC=" _rc
show , tag(ARM_B_KEEPI)

* ---- (b) dose ----
foreach k in 0 1 8 24 {
  mkdat , seed(101) k(`k')
  cap noi wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(`RP') rseed(`RS') dots(0)
  dis "DOSE`k'_RC=" _rc
  show , tag(DOSE_k`k')
}

* ---- (c) dead cluster + bootdf(clusters) ----
mkdat , seed(101)
qui replace y = . if county==1
cap noi wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(60) rseed(`RS') bootdf(clusters) dots(0)
dis "DEAD_RC=" _rc
show , tag(DEAD)
dis "DEAD_BOOTDF=" e(boot_df)

* control: physically drop the dead cluster's rows
mkdat , seed(101)
qui drop if county==1
cap noi wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(60) rseed(`RS') bootdf(clusters) dots(0)
dis "DEADDROP_RC=" _rc
show , tag(DEADDROP)
dis "DEADDROP_BOOTDF=" e(boot_df)

dis "R01_DONE"
