* r05: two controls on bb-01's mechanism.
*  D: i = state, cluster(state), but THREE rows per (state,year) cell (individual-level
*     administrative data with a state-level policy) - cells are no longer singletons,
*     but the Dirichlet multiplier is still constant within a cell.
*  E: the same design as r01 with only the IDIOSYNCRATIC noise changed (0.6 -> 1.5) and
*     no cluster-level heterogeneity at all.  If bayes se(att_it) moves there, the
*     invariance in r01 is specific to the cluster-level component and is not a frozen
*     constant.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

local k = 0
foreach het in 0 1.5 {
  local ++k
  set seed 717171
  set sortseed 717171
  clear
  set obs 24
  gen int st = _n
  gen double a_i = rnormal()*1.2
  gen double shk = rnormal()
  gen int onset = .
  replace onset = 2006 if st <=  6
  replace onset = 2009 if st >  6 & st <= 12
  replace onset = 2012 if st > 12 & st <= 17
  expand 12
  bysort st: gen int year = 2000 + _n
  expand 3
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.04*(year-2000) + (0.8 + `het'*shk)*post + rnormal()*0.8
  tempfile dd
  qui save `dd'
  dis "D[het=`het'] rows=" _N
  cap noi qui wooldid y st year onset, att att_it cluster(st)
  dis "D_A_RC_`het'=" _rc
  matrix A = e(wdidmainresults)
  dis "D[het=`het'] ANALYTIC se(att)=" %21x A[1,2] "  se(att_it)=" %21x A[2,2]
  use `dd', clear
  cap noi qui wooldid y st year onset, att att_it cluster(st) bootstrap(standard) reps(150) rseed(5555) dots(0)
  dis "D_S_RC_`het'=" _rc
  matrix VS = e(V)
  use `dd', clear
  cap noi qui wooldid y st year onset, att att_it cluster(st) bootstrap(bayes) reps(150) rseed(5555) dots(0) safety(off)
  dis "D_B_RC_`het'=" _rc
  matrix VB = e(V)
  dis "D[het=`het'] STANDARD se_att=" %12.8f sqrt(VS[1,1]) " se_attit=" %12.8f sqrt(VS[2,2])
  dis "D[het=`het'] BAYES    se_att=" %12.8f sqrt(VB[1,1]) " se_attit=" %12.8f sqrt(VB[2,2])
  dis "D[het=`het'] RATIO bayes/std: att=" %8.4f sqrt(VB[1,1]/VS[1,1]) " att_it=" %8.4f sqrt(VB[2,2]/VS[2,2])
  global DSBIT`k' = sqrt(VB[2,2])
  global DAN`k'   = A[2,2]
}
dis "D INVARIANCE bayes se(att_it): het0=" %21x $DSBIT1 " het1.5=" %21x $DSBIT2 "  relmove=" %12.4e ($DSBIT2-$DSBIT1)/$DSBIT1
dis "D INVARIANCE analytic se(att_it): het0=" %21x $DAN1 " het1.5=" %21x $DAN2 "  relmove=" %12.4e ($DAN2-$DAN1)/$DAN1

*========== E: idiosyncratic-noise sensitivity, no cluster heterogeneity
local k = 0
foreach sig in 0.6 1.5 {
  local ++k
  set seed 424242
  set sortseed 424242
  clear
  set obs 28
  gen int st = _n
  gen double a_i = rnormal()*1.2
  gen double shk = rnormal()
  gen int onset = .
  replace onset = 2006 if st <=  7
  replace onset = 2009 if st >  7 & st <= 14
  replace onset = 2012 if st > 14 & st <= 20
  expand 14
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.04*(year-2000) + 0.8*post + rnormal()*`sig'
  cap noi qui wooldid y st year onset, att att_it cluster(st) bootstrap(bayes) reps(200) rseed(77777) dots(0) safety(off)
  dis "E_RC_`sig'=" _rc
  matrix VB = e(V)
  dis "E[sig=`sig'] BAYES se_attit=" %12.8f sqrt(VB[2,2]) "  hex=" %21x sqrt(VB[2,2])
  global ESB`k' = sqrt(VB[2,2])
}
dis "E relmove bayes se(att_it) with idiosyncratic sd 0.6 -> 1.5 = " %12.4e ($ESB2-$ESB1)/$ESB1
dis "R05_DONE"
