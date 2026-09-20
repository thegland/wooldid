clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
*----------------------------------------------------------------------
* R04.  (1) FCM-3: is the regime recorded ANYWHERE in e() or on screen?
*       (2) FCM-4 severity: MANY dead clusters -> how wrong does
*           bootdf(clusters) get, and does the few-cluster caution fire?
*       (3) FCM-2 extra doors: phantom straddle carried on zero-weight
*           and missing-weight rows.
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
prog def escan
  syntax , tag(string)
  local hits 0
  foreach nm in `:e(macros)' {
    if strpos("`nm'","nest") | strpos("`nm'","reid") | strpos("`nm'","regime") {
      dis "`tag'_EMAC_HIT=`nm'"
      local ++hits
    }
  }
  foreach nm in `:e(scalars)' {
    if strpos("`nm'","nest") | strpos("`nm'","reid") | strpos("`nm'","regime") {
      dis "`tag'_ESCA_HIT=`nm'"
      local ++hits
    }
  }
  dis "`tag'_RECORD_HITS=`hits'"
  dis "`tag'_NCLUST=" e(N_clust) " DFR=" e(df_r) " BOOTDF=" e(boot_df) " N=" e(N)
end

*=== (1) FCM-3 : record of the regime ===
mkdat , seed(101)
qui replace y = . if plant==1 & year>=2005
cap noi qui wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(20) rseed(9) dots(0)
dis "NESTRUN_RC=" _rc
escan , tag(NESTRUN)
mat bA = e(b)
mat VA = e(V)

mkdat , seed(101)
qui replace y = . if plant==1 & year>=2005
qui replace county = 2 if plant==1 & year>=2005
cap noi qui wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(20) rseed(9) dots(0)
dis "KEEPRUN_RC=" _rc
escan , tag(KEEPRUN)

* verbose: does anything print about the regime?
mkdat , seed(101)
qui replace y = . if plant==1 & year>=2005
qui replace county = 2 if plant==1 & year>=2005
dis "VERBOSE_BEGIN"
cap noi wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(6) rseed(9) dots(0) verbose
dis "VERBOSE_END_RC=" _rc

*=== (2) FCM-4 : many dead clusters ===
* 20 counties x 4 plants; the outcome is only collected in 5 counties.
mkdat , seed(202) nc(20) np(4)
qui replace y = . if county > 5
cap noi wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(40) rseed(11) bootdf(clusters) dots(0)
dis "DEAD15_RC=" _rc
escan , tag(DEAD15)

* control: physically drop the 15 dead counties
mkdat , seed(202) nc(20) np(4)
qui drop if county > 5
cap noi wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(40) rseed(11) bootdf(clusters) dots(0)
dis "DEAD15DROP_RC=" _rc
escan , tag(DEAD15DROP)

*=== (3) FCM-2 extra doors: weights ===
* (a) zero aweight on the rows that carry the straddle
mkdat , seed(101)
gen double wt = 1
qui replace wt = 0    if plant==1 & year>=2005
qui replace county = 2 if plant==1 & year>=2005
cap noi qui wooldid y plant year onset [aw=wt] , cluster(county) att att_i
dis "W0_RC=" _rc " N=" e(N)
* (b) MISSING aweight on those rows
mkdat , seed(101)
gen double wt = 1
qui replace wt = .    if plant==1 & year>=2005
qui replace county = 2 if plant==1 & year>=2005
cap noi qui wooldid y plant year onset [aw=wt] , cluster(county) att att_i
dis "WM_RC=" _rc " N=" e(N)
* control (b0): missing aweight, no move
mkdat , seed(101)
gen double wt = 1
qui replace wt = .    if plant==1 & year>=2005
cap noi qui wooldid y plant year onset [aw=wt] , cluster(county) att att_i
dis "WM0_RC=" _rc " N=" e(N)

* now the bootstrap pair for the missing-weight door
foreach v in WMB0 WMB1 {
  mkdat , seed(101)
  gen double wt = 1
  qui replace wt = . if plant==1 & year>=2005
  if "`v'"=="WMB1" qui replace county = 2 if plant==1 & year>=2005
  cap noi qui wooldid y plant year onset [aw=wt] , cluster(county) att att_i bootstrap(standard) reps(150) rseed(4321) dots(0)
  dis "`v'_RC=" _rc
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  forvalues q=1/`K' {
    local this : word `q' of `nm'
    if "`this'"=="main_att_i" dis "`v'|att_i|b=" %21x `b'[1,`q'] "|se=" %21x sqrt(`V'[`q',`q'])
    if "`this'"=="main_att"   dis "`v'|att|b="   %21x `b'[1,`q'] "|se=" %21x sqrt(`V'[`q',`q'])
  }
  dis "`v'|N=" e(N)
}
dis "R04_DONE"
