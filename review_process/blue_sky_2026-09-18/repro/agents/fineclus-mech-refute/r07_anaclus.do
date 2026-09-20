clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
* R07: repro_FCM1 showed the ANALYTIC clustered se identical to 10 decimals for
* k=1, k=8 and k=24 straddlers but different from k=0.  Does wooldid's analytic
* clustered VCE see WHICH rows belong to which cluster, or only a binary flag?
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
  dis "`tag'|N=" e(N) "|NCLUST=" e(N_clust) "|DFR=" e(df_r)
end
* how many distinct (plant,county) pairs, i.e. how much the cluster map changed
prog def pairs
  tempvar t
  qui gegen double `t' = group(plant county)
  qui sum `t'
  dis "  distinct (plant,county) pairs = " r(max)
end
foreach k in 0 1 2 8 24 {
  mk , k(`k')
  dis "=== k=`k' ==="
  pairs
  qui wooldid y plant year onset , cluster(county) att att_i
  gsee , tag(K`k')
  * independent benchmark: static TWFE, reghdfe, same clustering
  qui reghdfe y D , absorb(plant year) cluster(county)
  dis "K`k'|reghdfe_D|b=" %21x _b[D] "|se=" %21x _se[D]
}
* same k=1 but a DIFFERENT plant straddles
foreach w in 1 45 93 {
  mk , who(`w')
  dis "=== single straddler = plant `w' ==="
  pairs
  qui wooldid y plant year onset , cluster(county) att att_i
  gsee , tag(W`w')
  qui reghdfe y D , absorb(plant year) cluster(county)
  dis "W`w'|reghdfe_D|b=" %21x _b[D] "|se=" %21x _se[D]
}
dis "R07_DONE"
