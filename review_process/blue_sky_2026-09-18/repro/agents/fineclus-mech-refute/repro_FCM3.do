*======================================================================
* SELF-CONTAINED REPRO -- FCM-3 : nothing records which bootstrap unit
* regime ran.  Two runs with the SAME estimation sample, SAME point
* estimates and DIFFERENT standard errors; scan every e() name and the
* whole verbose display for any trace of the regime.
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FCM3.do
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
prog def mkdat
  syntax , [seed(int 101)]
  clear
  set seed `seed'
  set sortseed `seed'
  set obs 96
  gen long plant  = _n
  gen long county = ceil(_n/4)
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
  expand 8
  bysort plant (county): gen int year = 2000 + _n
  gen double e = rnormal()*0.3
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = ai + cshk + 0.05*(year-2001) + ptau*D + e
  sort plant year
end
foreach arm in A B {
  mkdat
  qui replace y = . if plant==1 & year>=2005
  if "`arm'"=="B" qui replace county = 2 if plant==1 & year>=2005
  qui wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(40) rseed(9) dots(0)
  dis "---- ARM `arm' : every e() macro name ----"
  dis "`:e(macros)'"
  dis "---- ARM `arm' : every e() scalar name ----"
  dis "`:e(scalars)'"
  local hits 0
  foreach nm in `:e(macros)' `:e(scalars)' {
    if strpos("`nm'","nest") | strpos("`nm'","reid") | strpos("`nm'","regime") local ++hits
  }
  dis "ARM`arm'_REGIME_RECORD_HITS=`hits'"
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  forvalues q=1/`=colsof(`b')' {
     local this : word `q' of `nm'
     if "`this'"=="main_att_i" dis "ARM`arm'|att_i|b=" %21x `b'[1,`q'] "|se=" %12.10f sqrt(`V'[`q',`q'])
  }
}
* and with verbose, in full
mkdat
qui replace y = . if plant==1 & year>=2005
qui replace county = 2 if plant==1 & year>=2005
dis "==== VERBOSE RUN BEGINS ===="
wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(8) rseed(9) dots(0) verbose
dis "==== VERBOSE RUN ENDS ===="
dis "REPRO_FCM3_DONE"
