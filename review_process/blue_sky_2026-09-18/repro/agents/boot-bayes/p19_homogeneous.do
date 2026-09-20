* p19: the refusal message conditions the understatement on treatment effects that
* "vary across clusters" (ado:6439, sthlp:241). Does it also fire with HOMOGENEOUS
* effects?  Same design as p15 with the cluster-level effect shock switched off.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
foreach het in 0 1.2 {
  set seed 20260917
  set sortseed 20260917
  clear
  set obs 24
  gen int st = _n
  gen double a_i = rnormal()
  gen int onset = .
  replace onset = 2005 if st <= 6
  replace onset = 2008 if st > 6 & st <= 12
  replace onset = 2011 if st > 12 & st <= 16
  gen double shk = rnormal()
  expand 12
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.03*(year-2000) + (1.0 + `het'*shk)*post + rnormal()*0.5
  tempfile dd
  qui save `dd'
  qui wooldid y st year onset, att att_it cluster(st)
  matrix A = e(wdidmainresults)
  use `dd', clear
  qui wooldid y st year onset, att att_it cluster(st) bootstrap(standard) reps(400) rseed(31337) dots(0)
  matrix VS = e(V)
  use `dd', clear
  qui wooldid y st year onset, att att_it cluster(st) bootstrap(bayes) reps(400) rseed(31337) dots(0) safety(off)
  matrix VB = e(V)
  dis "[het=`het'] analytic se(att)=" %10.6f A[1,2] "  identical b: " (A[1,1]==A[2,1]) " identical se: " (A[1,2]==A[2,2])
  dis "[het=`het'] standard se(att)=" %10.6f sqrt(VS[1,1]) "  se(att_it)=" %10.6f sqrt(VS[2,2]) "  ratio=" %7.4f sqrt(VS[2,2]/VS[1,1])
  dis "[het=`het'] bayes    se(att)=" %10.6f sqrt(VB[1,1]) "  se(att_it)=" %10.6f sqrt(VB[2,2]) "  ratio=" %7.4f sqrt(VB[2,2]/VB[1,1])
  dis "[het=`het'] bayes/standard: att=" %7.4f sqrt(VB[1,1]/VS[1,1]) "  att_it=" %7.4f sqrt(VB[2,2]/VS[2,2])
}
dis "PROBE_DONE"
