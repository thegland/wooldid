*======================================================================
* SELF-CONTAINED REPRO -- FCM-1 (skeptic's independent version)
*   "one straddling unit out of 96 flips the GLOBAL bootstrap unit rule
*    and cuts se(att_i) by ~23%, silently"
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FCM1.do
* Reads repro_FCM1.log.  ~6 minutes.
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

* 24 counties x 4 plants, T=8.  cluster(county) is COARSER than i=plant,
* so the design's NESTING rule is the one that applies.  k of the 96 plants
* relocate to the next county halfway through the panel -- an ordinary
* administrative reclassification.  county enters ONLY through cluster(),
* so the point estimate cannot move: the estimation data is identical.
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
prog def grab
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  forvalues q=1/`K' {
     local this : word `q' of `nm'
     if "`this'"=="main_att"   dis "`tag'|att  |b=" %12.10f `b'[1,`q'] "|se=" %12.10f sqrt(`V'[`q',`q'])
     if "`this'"=="main_att_i" dis "`tag'|att_i|b=" %12.10f `b'[1,`q'] "|se=" %12.10f sqrt(`V'[`q',`q'])
  }
end

* ---- how many plants actually straddle, and what regime does ado:8986-8995 pick? ----
prog def regpred
  syntax varlist(min=2 max=2)
  tokenize `varlist'
  tempvar smp mn mx bad t
  qui gen byte `smp' = !missing(`2')
  qui gegen double `mn' = min(`2') if `smp', by(`1')
  qui gegen double `mx' = max(`2') if `smp', by(`1')
  cap assert `mn' == `mx' if `smp'
  local p = (_rc==0)
  qui gen byte `bad' = (`mn' != `mx')
  qui gegen byte `t' = tag(`1')
  qui count if `t'==1 & `bad'==1
  dis "  PREDNEST=`p'   straddling units = " r(N)
end

* ---- ANALYTIC (no bootstrap): the same four datasets ----
foreach k in 0 1 8 24 {
  mkdat , seed(101) k(`k')
  dis "=== k=`k' ==="
  regpred plant county
  qui wooldid y plant year onset , cluster(county) att att_i
  grab , tag(ANALYTIC_k`k')
}

* ---- BOOTSTRAP: identical seed, identical reps ----
foreach k in 0 1 8 24 {
  mkdat , seed(101) k(`k')
  qui wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(250) rseed(606060) dots(0)
  grab , tag(BOOT_k`k')
  dis "BOOT_k`k'|N=" e(N) "|NCLUST=" e(N_clust) "|completed=" e(boot_completed)
}
dis "REPRO_FCM1_DONE"
* EXPECTED READING OF THE LOG
*   ANALYTIC att_i b is IDENTICAL for k=0,1,8,24 (county only enters cluster()).
*   BOOT att_i se:  k=0 ~ .1259   k=1 ~ .0965   k=8 ~ .0976   k=24 ~ .0942
*   One straddler out of 96 buys the whole 23% drop; the next 23 buy nothing.
*   BOOT att se moves the OTHER way (up), as reassigning rows should.
