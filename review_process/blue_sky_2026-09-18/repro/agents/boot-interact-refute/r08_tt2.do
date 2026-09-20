* r08_tt2.do -- SELF-CONTAINED.  BOOTX-4 retry: hunt the sthlp:1059 geometry
* (timetrends + cluster() ON the cohort identifier i giving an ALL-ZERO analytic
* covariance), then check what bootstrap() does on it.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

cap prog drop mk
prog def mk
  syntax , ns(int) np(int)
  clear
  set obs `ns'
  gen state = _n
  gen onset = .
  replace onset = ceil(`np'/3) + mod(_n,4) if mod(_n,2)==0
  expand `np'
  bysort state: gen t = _n
  bysort state (t): gen double ueff = rnormal(0,.4) if _n==1
  bysort state (t): replace ueff = ueff[_n-1] if _n>1
  bysort state (t): gen double utr = rnormal(0,.05) if _n==1
  bysort state (t): replace utr = utr[_n-1] if _n>1
  gen double e = rnormal(0,1)
  gen double rel = t - onset
  gen double te = cond(onset<. & rel>=0, .8 + .1*rel, 0)
  gen double y = 2 + .05*t + ueff + utr*t + te + e
  drop rel te
end

cap prog drop st8
prog def st8
  args tag
  cap confirm matrix e(V)
  local rcv = _rc
  cap confirm matrix e(V_es_att)
  local rcves = _rc
  cap confirm matrix e(V_boot)
  local rcvb = _rc
  dis "TT|`tag'|rc=" _rc "|vcvstate=`e(wdid_vcvstate)'|confirmV=`rcv'|confirmVes=`rcves'|confirmVboot=`rcvb'|vposted=" e(boot_v_posted) "|completed=" e(boot_completed)
end

foreach NS in 40 24 {
  foreach NP in 8 10 12 {
    mk , ns(`NS') np(`NP')
    local L = min(3, floor(`NP'/4))
    qui cap noi wooldid y state t onset, cluster(state) esfixedbaseperiod esprelength(`L') espostlength(`L') timetrends
    dis "A_RC_`NS'_`NP'=" _rc
    st8 ANA_ns`NS'_np`NP'
    qui cap noi wooldid y state t onset, cluster(state) timetrends
    dis "A2_RC_`NS'_`NP'=" _rc
    st8 ANAnoES_ns`NS'_np`NP'
  }
}
dis "======== bootstrap arm on whichever shape came back 'absent' (all shapes tried)"
foreach NS in 40 24 {
  foreach NP in 8 10 12 {
    mk , ns(`NS') np(`NP')
    local L = min(3, floor(`NP'/4))
    qui cap noi wooldid y state t onset, cluster(state) esfixedbaseperiod esprelength(`L') espostlength(`L') timetrends bootstrap(standard) reps(30) rseed(4242)
    dis "B_RC_`NS'_`NP'=" _rc
    st8 BOOT_ns`NS'_np`NP'
  }
}
dis "DONE"
