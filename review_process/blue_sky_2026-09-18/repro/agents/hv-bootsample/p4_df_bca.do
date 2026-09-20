* p4_df_bca.do
*  A) TWO live clusters padded to twelve by dead ones: bootdf(clusters) and both
*     few-cluster cautions, full vs oracle.
*  B) bootalpha(equal) with WHOLLY dead clusters at reps(600): is the SE itself right?
*  C) bootdisplay(bca): the jackknife loops 1..wdb_C, so dead clusters are
*     zero-influence points in the acceleration.
*  D) event-study window as a complete-data trigger.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootsample"

cap prog drop mkpanel
prog def mkpanel
  syntax , nstates(integer) ntreat(integer) nunits(integer) nyears(integer)
  clear
  set obs `nstates'
  gen int state = _n
  gen double sfx = rnormal()*0.4
  gen int onset = .
  replace onset = 2003 + mod(_n-1,4) if _n <= `ntreat'
  expand `nunits'
  bysort state: gen int unit = (state-1)*`nunits' + _n
  gen double ufx = rnormal()*0.3
  expand `nyears'
  bysort unit: gen int year = 1999 + _n
  gen byte post = (onset < . & year >= onset)
  gen double y = 1 + .10*(year-1999) + sfx + ufx + .5*post + rnormal()*0.6
end

*=================================================== A: two live clusters of twelve
mkpanel , nstates(12) ntreat(6) nunits(3) nyears(10)
tempfile a0
save `a0', replace
* kill everything except state 6 (treated) and state 12 (never treated)
replace y = . if !inlist(state,6,12)
qui levelsof state if !missing(y), local(lv)
dis "A_LIVECLUS=" `: word count `lv''
dis "=== A1: FULL DATA (12 clusters in the frame, 2 in the estimation sample) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(150) bootdf(clusters) rseed(11) dots(0)
dis "A1_RC=" _rc
dis "A1_N=" e(N)
dis "A1_BOOTDF=" e(boot_df)
cap mat va = e(V)
cap mat ba = e(b)
cap dis "A1_B=" ba[1,1]
cap dis "A1_SE=" sqrt(va[1,1])
cap mat li e(boot_ci_percentile)
cap mat ca = e(boot_ci_wald)
cap mat li ca

dis "=== A2: DEAD ROWS DROPPED (oracle) ==="
use `a0', clear
keep if inlist(state,6,12)
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(150) bootdf(clusters) rseed(11) dots(0)
dis "A2_RC=" _rc
dis "A2_N=" e(N)
dis "A2_BOOTDF=" e(boot_df)
cap mat vb = e(V)
cap mat bb = e(b)
cap dis "A2_B=" bb[1,1]
cap dis "A2_SE=" sqrt(vb[1,1])
cap mat cb = e(boot_ci_wald)
cap mat li cb

*=================================================== B: equal + wholly dead clusters
set seed 2026
mkpanel , nstates(18) ntreat(8) nunits(2) nyears(10)
tempfile b0
save `b0', replace
* five never-treated states report no outcome at all
replace y = . if inlist(state,9,10,11,12,13)
dis "=== B1: FULL, bootalpha(equal), reps(600) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(600) bootdf(normal) rseed(606) dots(0)
dis "B1_RC=" _rc
dis "B1_N=" e(N)
mat vb1 = e(V)
mat bb1 = e(b)
dis "B1_B=" bb1[1,1]
dis "B1_SE=" sqrt(vb1[1,1])
use `b0', clear
drop if inlist(state,9,10,11,12,13)
dis "=== B2: ORACLE, bootalpha(equal), reps(600) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(600) bootdf(normal) rseed(606) dots(0)
dis "B2_RC=" _rc
dis "B2_N=" e(N)
mat vb2 = e(V)
mat bb2 = e(b)
dis "B2_B=" bb2[1,1]
dis "B2_SE=" sqrt(vb2[1,1])
dis "B_SE_RATIO=" sqrt(vb1[1,1])/sqrt(vb2[1,1])
dis "B_MCNOISE_PCT=" 100/sqrt(2*600)

*=================================================== C: BCa acceleration
use `b0', clear
replace y = . if inlist(state,9,10,11,12,13)
dis "=== C1: FULL, bootdisplay(bca) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(standard) reps(300) bootdisplay(bca) rseed(88) dots(0)
dis "C1_RC=" _rc
cap mat li e(boot_ci_bca)
cap mat li e(boot_accel)
cap mat li e(boot_a)
use `b0', clear
drop if inlist(state,9,10,11,12,13)
dis "=== C2: ORACLE, bootdisplay(bca) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(standard) reps(300) bootdisplay(bca) rseed(88) dots(0)
dis "C2_RC=" _rc
cap mat li e(boot_ci_bca)
cap mat li e(boot_accel)
cap mat li e(boot_a)

*=================================================== D: event-study window trigger
use `b0', clear
dis "=== D0: no window ==="
cap noi wooldid y unit year onset, cluster(state) fastnose
dis "D0_N=" e(N)
dis "=== D1: espostlength(2) esprelength(2) ==="
cap noi wooldid y unit year onset, cluster(state) esprelength(2) espostlength(2) fastnose
dis "D1_N=" e(N)
dis "P4_ALLDONE"
