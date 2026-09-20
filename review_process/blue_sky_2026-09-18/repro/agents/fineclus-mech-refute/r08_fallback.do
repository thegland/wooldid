clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
* R08: r07 showed wooldid's ANALYTIC cluster(county) SE is the SAME number
* (to a few ulps) for 1, 2, 8 and 24 straddling plants and for different
* straddling plants, while reghdfe's moves 12%.  What is it actually clustering on?
prog def mk
  syntax , [seed(int 101) k(int 0) who(int 0)]
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
  if `k' > 0 {
    forvalues q = 1/`k' {
      qui replace county = cond(`q'==24, 1, `q'+1) if plant==`=(`q'-1)*4+1' & year>=2005
    }
  }
  if `who' > 0 qui replace county = county + 1 if plant==`who' & year>=2005 & county < 24
  sort plant year
end
prog def gsee
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  forvalues q=1/`=colsof(`b')' {
     local this : word `q' of `nm'
     if "`this'"=="main_att"   dis "`tag'|att  |se=" %21x sqrt(`V'[`q',`q'])
     if "`this'"=="main_att_i" dis "`tag'|att_i|se=" %21x sqrt(`V'[`q',`q'])
  }
  dis "`tag'|NCLUST=" e(N_clust) "|DFR=" e(df_r) "|vcetype=" e(vce) "|clustvar=" e(clustvar)
end
* --- k=0 baseline, three VCEs ---
mk , k(0)
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(K0_CLcounty)
mk , k(0)
qui wooldid y plant year onset , cluster(plant) att att_i
gsee , tag(K0_CLplant)
mk , k(0)
qui wooldid y plant year onset , robust att att_i
gsee , tag(K0_ROBUST)
* --- k=1 ---
mk , k(1)
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(K1_CLcounty)
mk , k(1)
qui wooldid y plant year onset , cluster(plant) att att_i
gsee , tag(K1_CLplant)
* --- k=24 ---
mk , k(24)
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(K24_CLcounty)
mk , k(24)
qui wooldid y plant year onset , cluster(plant) att att_i
gsee , tag(K24_CLplant)
* --- a cluster variable that is CONSTANT within plant but re-labels the counties
*     (a pure relabelling: the partition is identical to k=0's) ---
mk , k(0)
qui gen long county2 = 25 - county
qui wooldid y plant year onset , cluster(county2) att att_i
gsee , tag(K0_RELABEL)
* --- a genuinely different coarser partition: 12 regions of 2 counties ---
mk , k(0)
qui gen long region = ceil(county/2)
qui wooldid y plant year onset , cluster(region) att att_i
gsee , tag(K0_REGION)
mk , k(1)
qui gen long region = ceil(county/2)
qui wooldid y plant year onset , cluster(region) att att_i
gsee , tag(K1_REGION)
dis "R08_DONE"
