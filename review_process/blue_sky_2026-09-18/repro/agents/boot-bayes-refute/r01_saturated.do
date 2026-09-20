* r01: REFUTE/CONFIRM bb-01 on MY OWN data, their structural shape.
* i = unit = cluster, one row per (i,t) -> every treated (i,relyr) cell is a SINGLETON
* fitted by its own dummy.  Two arms: homogeneous TE, and a cluster-constant shock on
* the post rows.  Does bayes se(att_it) really not move at all with the heterogeneity?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

local k = 0
foreach het in 0 1.5 {
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
  * st 21..28 never treated
  expand 14
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.04*(year-2000) + (0.8 + `het'*shk)*post + rnormal()*0.6
  tempfile dd
  qui save `dd'

  cap noi qui wooldid y st year onset, att att_it cluster(st)
  dis "A_RC_`het'=" _rc
  matrix A = e(wdidmainresults)
  dis "rownames A: " "`: rowfullnames A'"
  local bat  = A[1,1]
  local sat  = A[1,2]
  local bit  = A[2,1]
  local sit  = A[2,2]
  dis "[het=`het'] ANALYTIC att    b=" %21x `bat' "  se=" %21x `sat'
  dis "[het=`het'] ANALYTIC att_it b=" %21x `bit' "  se=" %21x `sit'
  dis "[het=`het'] identical_b=" (`bat'==`bit') "  identical_se=" (`sat'==`sit')

  use `dd', clear
  cap noi qui wooldid y st year onset, att att_it cluster(st) bootstrap(standard) reps(200) rseed(77777) dots(0)
  dis "S_RC_`het'=" _rc
  matrix VS = e(V)
  local ss_att = sqrt(VS[1,1])
  local ss_it  = sqrt(VS[2,2])

  use `dd', clear
  cap noi qui wooldid y st year onset, att att_it cluster(st) bootstrap(bayes) reps(200) rseed(77777) dots(0) safety(off)
  dis "B_RC_`het'=" _rc
  matrix VB = e(V)
  local sb_att = sqrt(VB[1,1])
  local sb_it  = sqrt(VB[2,2])

  dis "[het=`het'] STANDARD se_att=" %12.8f `ss_att' "  se_attit=" %12.8f `ss_it' "  hexit=" %21x `ss_it'
  dis "[het=`het'] BAYES    se_att=" %12.8f `sb_att' "  se_attit=" %12.8f `sb_it' "  hexit=" %21x `sb_it'
  dis "[het=`het'] RATIO bayes/std: att=" %8.4f `sb_att'/`ss_att' "  att_it=" %8.4f `sb_it'/`ss_it'
  global SBIT`k' = `sb_it'
  global SSIT`k' = `ss_it'
  global SBAT`k' = `sb_att'
}
dis "INVARIANCE CHECK  bayes se(att_it): het0=" %21x $SBIT1 "   het1.5=" %21x $SBIT2
dis "relative move bayes se(att_it) = " %12.4e (${SBIT2}-${SBIT1})/${SBIT1}
dis "relative move std   se(att_it) = " %12.4e (${SSIT2}-${SSIT1})/${SSIT1}
dis "relative move bayes se(att)    = " %12.4e (${SBAT2}-${SBAT1})/${SBAT1}
dis "R01_DONE"
