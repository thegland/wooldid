*======================================================================
* SELF-CONTAINED REPRO -- FCM-4 and the finding beside it.
* A cluster with NO rows in e(sample) is still a resampling unit:
*   (i)  bootdf(clusters) uses it, so e(df_r) is the wrong t;
*   (ii) the few-cluster and few-treated-cluster CAUTIONS are silenced,
*        so the run that most needs them is the one that does not get them.
* Control: the same estimation sample with the dead clusters physically
* dropped.  Point estimates are identical; everything else moves.
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FCM4.do
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
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
prog def rep
  syntax , tag(string)
  dis "`tag'|N=" e(N) "|N_clust=" e(N_clust) "|df_r=" e(df_r) "|boot_df=" e(boot_df)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  forvalues q=1/`=colsof(`b')' {
     local this : word `q' of `nm'
     if "`this'"=="main_att_i" {
        local se = sqrt(`V'[`q',`q'])
        dis "`tag'|att_i|b=" %12.10f `b'[1,`q'] "|se=" %12.10f `se'
        if e(boot_df) < . {
          dis "`tag'|t crit at e(boot_df)=" e(boot_df) " is " %8.6f invttail(e(boot_df),.025) ///
              "; at N_clust-1=" e(N_clust)-1 " it is " %8.6f invttail(e(N_clust)-1,.025) ///
              "; CI half-width ratio = " %6.4f invttail(e(boot_df),.025)/invttail(e(N_clust)-1,.025)
        }
     }
  }
end

* ===== A. 20 counties in the file; the outcome is collected in only 5 =====
dis "===== A: 15 of 20 clusters have no outcome at all ====="
mkdat
qui replace y = . if county > 5
wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(60) rseed(11) bootdf(clusters) bootdisplay(wald) dots(0)
rep , tag(DEAD15)

* ===== B. control: same estimation sample, dead clusters removed =====
dis "===== B: control, the 15 dead clusters physically dropped ====="
mkdat
qui drop if county > 5
wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(60) rseed(11) bootdf(clusters) bootdisplay(wald) dots(0)
rep , tag(DEAD15DROP)

* ===== C. the milder, very common version: ONE dead cluster of 24 =====
dis "===== C: one dead cluster of 24 ====="
mkdat , nc(24)
qui replace y = . if county == 1
wooldid y plant year onset , cluster(county) att att_i bootstrap(standard) reps(60) rseed(11) bootdf(clusters) dots(0)
rep , tag(DEAD1)
dis "REPRO_FCM4_DONE"
* EXPECTED READING
*  A: header prints "N clusters = 5", e(df_r)=e(boot_df)=19, and NEITHER
*     few-cluster caution is printed.
*  B: identical point estimates, e(df_r)=e(boot_df)=4, and BOTH cautions print.
*     t(19)=2.093 vs t(4)=2.776 -- A's wald interval is 25% too narrow.
*  C: e(N_clust)=23 but e(df_r)=23 (should be 22).
