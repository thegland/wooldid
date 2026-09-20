* r09 (own hunt): how general is the ANALYTIC-se blindness?
* Four arms of the same panel, identical except for the treatment-effect heterogeneity:
*   0  homogeneous
*   1  cluster-constant effect shock, sd 1.5
*   2  cluster-constant effect shock, sd 4.0   (the truth is 2.7x more variable)
*   3  cluster-AND-relative-time varying effect shock, sd 1.5  (dynamic heterogeneity)
* i = state = cluster.  The reported ANALYTIC clustered se, e(V) and the printed CI are
* compared bitwise across the arms, and bootstrap(standard) is run at arms 0 and 2 to
* show what a cell-resampling variance does over the same range.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

forvalues arm = 0/3 {
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
  gen double dyn = rnormal()
  local eff 0.8
  if `arm' == 1 local eff (0.8 + 1.5*shk)
  if `arm' == 2 local eff (0.8 + 4.0*shk)
  if `arm' == 3 local eff (0.8 + 1.5*dyn)
  gen double y = a_i + 0.04*(year-2000) + `eff'*post + rnormal()*0.6
  tempfile dd
  qui save `dd'
  cap noi qui wooldid y st year onset, att att_it cluster(st)
  dis "ARM`arm'_RC=" _rc
  matrix A = e(wdidmainresults)
  dis "ARM`arm' b(att)=" %21x A[1,1] "  se(att)=" %21x A[1,2] "  CI=[" %12.8f A[1,5] "," %12.8f A[1,6] "]"
  dis "ARM`arm' se(att_it)=" %21x A[2,2]
  global SE`arm' = A[1,2]
  global B`arm'  = A[1,1]
  if `arm' == 0 | `arm' == 2 {
    use `dd', clear
    cap noi qui wooldid y st year onset, att att_it cluster(st) bootstrap(standard) reps(200) rseed(77777) dots(0)
    dis "ARM`arm'_STD_RC=" _rc
    matrix VS = e(V)
    dis "ARM`arm' STANDARD se(att)=" %12.8f sqrt(VS[1,1]) "  se(att_it)=" %12.8f sqrt(VS[2,2])
    global STD`arm' = sqrt(VS[1,1])
  }
}
dis "analytic se(att): arm0=" %21x $SE0 " arm1=" %21x $SE1 " arm2=" %21x $SE2 " arm3=" %21x $SE3
dis "bitwise equal 0-1: " ($SE0==$SE1) "   0-2: " ($SE0==$SE2) "   0-3: " ($SE0==$SE3)
dis "point estimates:  arm0=" %10.6f $B0 " arm1=" %10.6f $B1 " arm2=" %10.6f $B2 " arm3=" %10.6f $B3
dis "bootstrap(standard) se(att): arm0=" %10.6f $STD0 "  arm2=" %10.6f $STD2 "   ratio=" %8.4f $STD2/$STD0
dis "analytic/standard: arm0=" %8.4f $SE0/$STD0 "  arm2=" %8.4f $SE2/$STD2
dis "R09_DONE"
