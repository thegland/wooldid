* r02: SCOPE of bb-01.  The prober measured only cluster() == i with SINGLETON treated
* cells.  Two designs they did not test:
*   B: i = onset cohort (5 groups), cluster(state) - cluster FINER than i, cells hold 8 rows
*   C: i = state, cluster(region)               - cluster COARSER than i, cells singleton
* If bayes tracks standard in B, bb-01 scopes to "cluster nests i AND cells are singletons".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

*=================== DESIGN B: i = cohort group, cluster(state) finer than i
local k = 0
foreach het in 0 1.5 {
  local ++k
  set seed 909090
  set sortseed 909090
  clear
  set obs 40
  gen int st = _n
  gen int grp = ceil(_n/8)
  gen double a_s = rnormal()*1.2
  gen double shk = rnormal()
  gen int onset = .
  replace onset = 2006 if grp==1
  replace onset = 2008 if grp==2
  replace onset = 2010 if grp==3
  replace onset = 2012 if grp==4
  expand 14
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_s + 0.04*(year-2000) + (0.8 + `het'*shk)*post + rnormal()*0.6
  tempfile dd
  qui save `dd'
  qui gunique grp year if post
  dis "B[het=`het'] treated (i,t) cells=" r(unique) "   rows=" _N

  cap noi qui wooldid y grp year onset, att att_it cluster(st)
  dis "B_A_RC_`het'=" _rc
  matrix A = e(wdidmainresults)
  dis "B[het=`het'] ANALYTIC att b=" %21x A[1,1] " se=" %21x A[1,2] "  att_it b=" %21x A[2,1] " se=" %21x A[2,2]

  use `dd', clear
  cap noi qui wooldid y grp year onset, att att_it cluster(st) bootstrap(standard) reps(150) rseed(77777) dots(0)
  dis "B_S_RC_`het'=" _rc
  matrix VS = e(V)
  use `dd', clear
  cap noi qui wooldid y grp year onset, att att_it cluster(st) bootstrap(bayes) reps(150) rseed(77777) dots(0) safety(off)
  dis "B_B_RC_`het'=" _rc
  matrix VB = e(V)
  dis "B[het=`het'] STANDARD se_att=" %12.8f sqrt(VS[1,1]) " se_attit=" %12.8f sqrt(VS[2,2])
  dis "B[het=`het'] BAYES    se_att=" %12.8f sqrt(VB[1,1]) " se_attit=" %12.8f sqrt(VB[2,2])
  dis "B[het=`het'] RATIO bayes/std: att=" %8.4f sqrt(VB[1,1]/VS[1,1]) "  att_it=" %8.4f sqrt(VB[2,2]/VS[2,2])
  global BSBIT`k' = sqrt(VB[2,2])
  global BSSIT`k' = sqrt(VS[2,2])
}
dis "B INVARIANCE: bayes se(att_it) het0=" %21x $BSBIT1 "  het1.5=" %21x $BSBIT2
dis "B relmove bayes se(att_it)=" %12.4e ($BSBIT2-$BSBIT1)/$BSBIT1 "   std=" %12.4e ($BSSIT2-$BSSIT1)/$BSSIT1

*=================== DESIGN C: i = state, cluster(region) coarser than i
local k = 0
foreach het in 0 1.5 {
  local ++k
  set seed 313131
  set sortseed 313131
  clear
  set obs 28
  gen int st = _n
  gen int reg = ceil(_n/4)
  gen double a_i = rnormal()*1.2
  gen double shk = rnormal()
  gen int onset = .
  replace onset = 2006 if mod(st,4)==1
  replace onset = 2009 if mod(st,4)==2
  replace onset = 2012 if mod(st,4)==3
  expand 14
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.04*(year-2000) + (0.8 + `het'*shk)*post + rnormal()*0.6
  tempfile dd
  qui save `dd'

  cap noi qui wooldid y st year onset, att att_it cluster(reg)
  dis "C_A_RC_`het'=" _rc
  matrix A = e(wdidmainresults)
  dis "C[het=`het'] ANALYTIC att b=" %21x A[1,1] " se=" %21x A[1,2] "  att_it b=" %21x A[2,1] " se=" %21x A[2,2]

  use `dd', clear
  cap noi qui wooldid y st year onset, att att_it cluster(reg) bootstrap(standard) reps(150) rseed(77777) dots(0)
  dis "C_S_RC_`het'=" _rc
  matrix VS = e(V)
  use `dd', clear
  cap noi qui wooldid y st year onset, att att_it cluster(reg) bootstrap(bayes) reps(150) rseed(77777) dots(0) safety(off)
  dis "C_B_RC_`het'=" _rc
  matrix VB = e(V)
  dis "C[het=`het'] STANDARD se_att=" %12.8f sqrt(VS[1,1]) " se_attit=" %12.8f sqrt(VS[2,2])
  dis "C[het=`het'] BAYES    se_att=" %12.8f sqrt(VB[1,1]) " se_attit=" %12.8f sqrt(VB[2,2])
  dis "C[het=`het'] RATIO bayes/std: att=" %8.4f sqrt(VB[1,1]/VS[1,1]) "  att_it=" %8.4f sqrt(VB[2,2]/VS[2,2])
  global CSBIT`k' = sqrt(VB[2,2])
  global CSSIT`k' = sqrt(VS[2,2])
}
dis "C INVARIANCE: bayes se(att_it) het0=" %21x $CSBIT1 "  het1.5=" %21x $CSBIT2
dis "C relmove bayes se(att_it)=" %12.4e ($CSBIT2-$CSBIT1)/$CSBIT1 "   std=" %12.4e ($CSSIT2-$CSSIT1)/$CSSIT1
dis "R02_DONE"
