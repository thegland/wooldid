clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
* R06: does a dead cluster as a resampling unit move the bootstrap SE?
* 500 reps, three seeds each arm, same data.
prog def mkdat
  syntax , [seed(int 202) nc(int 20) np(int 4) tt(int 8)]
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
prog def rep2
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  forvalues q=1/`=colsof(`b')' {
     local this : word `q' of `nm'
     if "`this'"=="main_att_i" dis "`tag'|att_i|se=" %12.10f sqrt(`V'[`q',`q'])
     if "`this'"=="main_att"   dis "`tag'|att  |se=" %12.10f sqrt(`V'[`q',`q'])
  }
  dis "`tag'|done=" e(boot_completed) "|failed=" e(boot_failed) "|df_r=" e(df_r)
end
foreach s in 11 22 33 {
  mkdat
  qui replace y = . if county > 5
  cap qui wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(500) rseed(`s') dots(0)
  dis "PAD_s`s'_RC=" _rc
  rep2 , tag(PAD_s`s')
  mkdat
  qui drop if county > 5
  cap qui wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(500) rseed(`s') dots(0)
  dis "DRP_s`s'_RC=" _rc
  rep2 , tag(DRP_s`s')
}
dis "R06_DONE"
