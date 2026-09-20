*======================================================================
* SELF-CONTAINED REPRO -- FCM-2 (skeptic's independent version)
*   "the regime predicate runs on `smp', which is NOT the estimation
*    sample, so cluster values on rows the regression never sees change
*    the reported att_i bootstrap SE"
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FCM2.do
* ~5 minutes.
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
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
prog def grab
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  forvalues q=1/`K' {
     local this : word `q' of `nm'
     foreach w in main_att main_att_it main_att_i main_att_itime {
       if "`this'"=="`w'" dis "`tag'|`w'|b=" %21x `b'[1,`q'] "|se=" %12.10f sqrt(`V'[`q',`q'])
     }
  }
  dis "`tag'|N=" e(N) "|NCLUST=" e(N_clust)
end

*----------------------------------------------------------------------
* DOOR 1 -- a MISSING OUTCOME.  Plant 1 shuts down after 2004: y is
* missing from 2005 on.  In arm B the registry ALSO records it in county 2
* on exactly those (already-excluded) rows.
*----------------------------------------------------------------------
foreach arm in A B {
  mkdat , seed(101)
  qui replace y = . if plant==1 & year>=2005
  if "`arm'"=="B" qui replace county = 2 if plant==1 & year>=2005
  * what the estimator sees, and what ado:8986-8995 sees
  tempvar mn mx mn2 mx2
  * what ado:8986-8995 computes: min/max of the cluster over `smp' = all rows
  qui gegen double `mn' = min(county), by(plant)
  qui gegen double `mx' = max(county), by(plant)
  cap assert `mn'==`mx'
  dis "ARM`arm'  NEST on smp (all rows)         = " (_rc==0)
  * what the ESTIMATION SAMPLE says: min/max computed over rows with y only
  qui gegen double `mn2' = min(county) if !missing(y), by(plant)
  qui gegen double `mx2' = max(county) if !missing(y), by(plant)
  cap assert `mn2'==`mx2' if !missing(y)
  dis "ARM`arm'  NEST on e(sample) rows only    = " (_rc==0)
  drop `mn' `mx' `mn2' `mx2'
  qui wooldid y plant year onset , cluster(county) att att_i att_it att_itime
  grab , tag(ANALYTIC_ARM`arm')
  qui wooldid y plant year onset , cluster(county) att att_i att_it att_itime bootstrap(standard) reps(250) rseed(606060) dots(0)
  grab , tag(BOOT_ARM`arm')
}

*----------------------------------------------------------------------
* DOOR 2 -- a MISSING WEIGHT.  Same story with [aw=wt]; wt is missing on
* the rows that carry the phantom move.
*----------------------------------------------------------------------
foreach arm in A B {
  mkdat , seed(101)
  gen double wt = 1
  qui replace wt = . if plant==1 & year>=2005
  if "`arm'"=="B" qui replace county = 2 if plant==1 & year>=2005
  qui wooldid y plant year onset [aw=wt] , cluster(county) att att_i bootstrap(standard) reps(150) rseed(4321) dots(0)
  grab , tag(WBOOT_ARM`arm')
}
dis "REPRO_FCM2_DONE"
* EXPECTED READING
*  ARMA and ARMB: same N (764), same b to every hex digit, same analytic se.
*  NEST on smp (all rows): A=1, B=0.   NEST on e(sample) rows only: A=1, B=1.
*  BOOT att   se: identical to 3 ulps (the draws are the same).
*  BOOT att_i se: A ~ .1257   B ~ .0962   ratio .765
*  WBOOT att_i se: A ~ .1264  B ~ .1053   (the same phantom via a missing WEIGHT)
