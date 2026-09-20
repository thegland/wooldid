* p1_guards.do -- does the bootstrap engine count clusters that contribute ZERO rows
* to e(sample)?  Oracle: physically drop the dead rows and re-run identically.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid
dis "BUILDCHECK_DONE"

*---------------------------------------------------------------- data generator
cap prog drop mkpanel
prog def mkpanel
  syntax , nstates(integer) ntreat(integer) [deadto(integer 0)]
  clear
  set obs `nstates'
  gen int state = _n
  gen double sfx = rnormal()*0.4
  gen int onset = .
  replace onset = 2003 + mod(_n-1,5) if _n <= `ntreat'
  expand 2
  bysort state: gen int unit = (state-1)*2 + _n
  gen double ufx = rnormal()*0.3
  expand 10
  bysort unit: gen int year = 1999 + _n
  gen byte post = (onset < . & year >= onset)
  gen double y = 1 + .10*(year-1999) + sfx + ufx + .5*post + rnormal()*0.6
  * "dead" states: the outcome is simply not reported there (very common in admin data)
  if `deadto' > 0 replace y = . if state <= `deadto'
  label var y "outcome"
end

*================================================================ DESIGN 1
* 20 states, 8 treated. Kill states 1-4 (all treated) -> 16 live clusters,
* 4 LIVE treated clusters. The "<6 treated clusters" caution should fire.
mkpanel , nstates(20) ntreat(8) deadto(4)
tempfile d1
save `d1', replace

qui count
dis "D1_ROWS_FULL=" r(N)
qui count if !missing(y)
dis "D1_ROWS_LIVE=" r(N)
qui levelsof state if !missing(y), local(ls)
dis "D1_LIVECLUS=" `: word count `ls''

dis "=== D1A: FULL DATA (dead states present) ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(100) bootdf(clusters) rseed(777) dots(0)
dis "D1A_RC=" _rc
dis "D1A_N=" e(N)
dis "D1A_BOOTDF=" e(boot_df)
dis "D1A_COMPLETED=" e(boot_completed)
dis "D1A_FAILED=" e(boot_failed)
mat b1a = e(b)
mat v1a = e(V)
dis "D1A_B=" b1a[1,1]
dis "D1A_SE=" sqrt(v1a[1,1])

dis "=== D1B: DEAD ROWS PHYSICALLY DROPPED (the oracle) ==="
use `d1', clear
drop if missing(y)
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(100) bootdf(clusters) rseed(777) dots(0)
dis "D1B_RC=" _rc
dis "D1B_N=" e(N)
dis "D1B_BOOTDF=" e(boot_df)
dis "D1B_COMPLETED=" e(boot_completed)
dis "D1B_FAILED=" e(boot_failed)
mat b1b = e(b)
mat v1b = e(V)
dis "D1B_B=" b1b[1,1]
dis "D1B_SE=" sqrt(v1b[1,1])

*================================================================ DESIGN 2
* 12 states, 6 treated. Kill 4 -> 8 live clusters. The "<10 clusters" caution
* should fire on the live sample and must not be suppressed by the dead ones.
set seed 4242
mkpanel , nstates(12) ntreat(6) deadto(4)
tempfile d2
save `d2', replace
qui levelsof state if !missing(y), local(ls2)
dis "D2_LIVECLUS=" `: word count `ls2''

dis "=== D2A: FULL DATA ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(standard) reps(100) bootdf(clusters) rseed(777) dots(0)
dis "D2A_RC=" _rc
dis "D2A_N=" e(N)
dis "D2A_BOOTDF=" e(boot_df)
dis "D2A_COMPLETED=" e(boot_completed)
dis "D2A_FAILED=" e(boot_failed)
mat v2a = e(V)
dis "D2A_SE=" sqrt(v2a[1,1])

dis "=== D2B: DEAD ROWS DROPPED ==="
use `d2', clear
drop if missing(y)
cap noi wooldid y unit year onset, cluster(state) bootstrap(standard) reps(100) bootdf(clusters) rseed(777) dots(0)
dis "D2B_RC=" _rc
dis "D2B_N=" e(N)
dis "D2B_BOOTDF=" e(boot_df)
dis "D2B_COMPLETED=" e(boot_completed)
dis "D2B_FAILED=" e(boot_failed)
mat v2b = e(V)
dis "D2B_SE=" sqrt(v2b[1,1])

*================================================================ DESIGN 3 (degenerate)
* 6 states, 5 dead -> ONE live cluster. The r(459) message claims to count
* clusters "in the estimation sample".
set seed 99
mkpanel , nstates(6) ntreat(3) deadto(5)
tempfile d3
save `d3', replace
* make the single live state treated so the run is estimable at all
dis "=== D3A: FULL DATA, 1 live cluster ==="
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(60) bootdf(clusters) rseed(777) dots(0)
dis "D3A_RC=" _rc
dis "D3A_N=" e(N)
dis "D3A_BOOTDF=" e(boot_df)
cap mat v3a = e(V)
cap dis "D3A_SE=" sqrt(v3a[1,1])

dis "=== D3B: DEAD ROWS DROPPED ==="
use `d3', clear
drop if missing(y)
cap noi wooldid y unit year onset, cluster(state) bootstrap(bayes) reps(60) bootdf(clusters) rseed(777) dots(0)
dis "D3B_RC=" _rc

dis "P1_ALLDONE"
