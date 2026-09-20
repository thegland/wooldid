clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150

which wooldid

*---------------------------------------------------------------
* MY OWN generator.  24 counties x 4 plants = 96 plants, T=8.
* cluster(county) is COARSER than i=plant, so the nesting regime
* is the one the design says applies.
*---------------------------------------------------------------
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
  * county-level effect heterogeneity + plant-level heterogeneity
  by county, sort: gen double cshk = rnormal()      if _n==1
  by county:       gen double ctau = 1 + rnormal()*0.8 if _n==1
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

prog def regpred
  * replicate ado:8986-8995 EXACTLY on the smp concept the code uses
  syntax varlist(min=2 max=2) , [ifr(string asis)]
  tokenize `varlist'
  local iv `1'
  local cv `2'
  tempvar smp mn mx
  if `"`ifr'"' != "" {
    qui gen byte `smp' = ((`ifr')==1)
    qui replace `smp' = 0 if `smp' >= .
  }
  else qui gen byte `smp' = 1
  qui replace `smp' = 0 if missing(`cv')
  qui gegen double `mn' = min(`cv') if `smp', by(`iv')
  qui gegen double `mx' = max(`cv') if `smp', by(`iv')
  cap assert `mn' == `mx' if `smp'
  c_local PRED = (_rc==0)
end

*-------------------------------------------------
* D0 : clean nesting
*-------------------------------------------------
mkdat , seed(101)
regpred plant county
dis "D0_PRED=`PRED'"
save "d0.dta", replace

*-------------------------------------------------
* D1 : PHANTOM straddle.  plant 1 stops reporting y after 2004 (it shuts
*      down) and the registry moves it to county 2 on exactly those rows.
*      Rows with missing y are NOT in e(sample) but ARE in `smp'.
*-------------------------------------------------
use "d0.dta", clear
replace y = . if plant==1 & year>=2005
replace county = 2 if plant==1 & year>=2005
regpred plant county
dis "D1_PRED=`PRED'"
save "d1.dta", replace

*-------------------------------------------------
* D0m : same missing y, NO county move.  Identical estimation sample to D1.
*-------------------------------------------------
use "d0.dta", clear
replace y = . if plant==1 & year>=2005
regpred plant county
dis "D0m_PRED=`PRED'"
save "d0m.dta", replace

*-------------------------------------------------
* D2 : a REAL straddler -- plant 1 genuinely moves to county 2 in 2005,
*      y observed throughout.  One plant of 96.
*-------------------------------------------------
use "d0.dta", clear
replace county = 2 if plant==1 & year>=2005
regpred plant county
dis "D2_PRED=`PRED'"
save "d2.dta", replace

*-------------------------------------------------
* D3 : dead cluster -- county 1's four plants never report y.
*-------------------------------------------------
use "d0.dta", clear
replace y = . if county==1
regpred plant county
dis "D3_PRED=`PRED'"
save "d3.dta", replace

*-------------------------------------------------
* D4 : STRING cluster variable (state FIPS as a string -- very common)
*-------------------------------------------------
use "d0.dta", clear
gen str6 cstr = string(county,"%03.0f")
save "d4.dta", replace

*=================================================
* STRUCTURAL CHECKS (no bootstrap yet)
*=================================================
prog def show
  syntax , tag(string)
  tempname b V
  cap mat `b' = e(b)
  if _rc {
     dis "`tag'_NOB=1"
     exit
  }
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  forvalues q=1/`K' {
     local this : word `q' of `nm'
     if "`this'"=="main_att"    dis "`tag'_ATT=" %21x `b'[1,`q'] "  SE=" %21x sqrt(`V'[`q',`q'])
     if "`this'"=="main_att_i"  dis "`tag'_ATTI=" %21x `b'[1,`q'] "  SE=" %21x sqrt(`V'[`q',`q'])
  }
  dis "`tag'_N=" e(N) " NCLUST=" e(N_clust) " DFR=" e(df_r)
end

foreach d in d0 d1 d0m d2 d3 {
  use "`d'.dta", clear
  cap noi qui wooldid y plant year onset , cluster(county) att att_i
  dis "`d'_RC=" _rc
  show , tag(`d')
}

* ---- string cluster variable, plain (no bootstrap) ----
use "d4.dta", clear
cap noi qui wooldid y plant year onset , cluster(cstr) att att_i
dis "STRPLAIN_RC=" _rc
show , tag(STRPLAIN)

* ---- string cluster variable, bootstrap(standard) + att_i  -> reid leg ----
use "d4.dta", clear
cap noi wooldid y plant year onset , cluster(cstr) att att_i bootstrap(standard) reps(4) rseed(77) dots(0)
dis "STRBOOT_RC=" _rc

* ---- string cluster variable, bootstrap(standard) att only (count path) ----
use "d4.dta", clear
cap noi qui wooldid y plant year onset , cluster(cstr) att bootstrap(standard) reps(4) rseed(77) dots(0)
dis "STRBOOTATT_RC=" _rc

* ---- missing i values (rows the estimator drops) with a county move on them ----
use "d0.dta", clear
gen long plant2 = plant
replace plant2 = . if plant==1 & year>=2005
replace county = 2 if plant==1 & year>=2005
regpred plant2 county
dis "D5_PRED=`PRED'"
cap noi qui wooldid y plant2 year onset , cluster(county) att att_i
dis "D5_RC=" _rc
show , tag(D5)

dis "R00_DONE"
