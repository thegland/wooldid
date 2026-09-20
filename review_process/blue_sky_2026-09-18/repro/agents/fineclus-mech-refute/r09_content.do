clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
* R09: does the analytic clustered SE respond to the CONTENT of the straddling,
* or only to its existence?  y is fixed in every variant (county never enters y).
prog def mk
  clear
  set seed 101
  set sortseed 101
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
prog def gsee
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  forvalues q=1/`=colsof(`b')' {
     local this : word `q' of `nm'
     if "`this'"=="main_att"   dis "`tag'|att  |b=" %21x `b'[1,`q'] "|se=" %21x sqrt(`V'[`q',`q'])
     if "`this'"=="main_att_i" dis "`tag'|att_i|se=" %21x sqrt(`V'[`q',`q'])
  }
  dis "`tag'|NCLUST=" e(N_clust) "|DFR=" e(df_r)
end
* (a) plant 1 -> county 2 from 2005
mk
qui replace county = 2 if plant==1 & year>=2005
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(A_p1c2_2005)
* (b) plant 1 -> county 2 from 2003 (6 rows move, not 4)
mk
qui replace county = 2 if plant==1 & year>=2003
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(B_p1c2_2003)
* (c) plant 1 -> county 13 from 2005 (a far county)
mk
qui replace county = 13 if plant==1 & year>=2005
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(C_p1c13_2005)
* (d) plant 2 -> county 3 from 2005 (a different plant, different cohort)
mk
qui replace county = 3 if plant==2 & year>=2005
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(D_p2c3_2005)
* (e) plant 1 -> a BRAND NEW county 25 from 2005
mk
qui replace county = 25 if plant==1 & year>=2005
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(E_p1c25_2005)
* (f) HALF the plants straddle
mk
qui replace county = county + 1 if mod(plant,2)==1 & year>=2005 & county < 24
qui wooldid y plant year onset , cluster(county) att att_i
gsee , tag(F_half)
* (g) cluster BELOW i for contrast: cluster on plant-halves
mk
gen long half = plant*2 - (year<2005)
qui wooldid y plant year onset , cluster(half) att att_i
gsee , tag(G_finer)
dis "R09_DONE"
