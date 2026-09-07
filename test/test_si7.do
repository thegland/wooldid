* test_si7.do — certification of the SI-7 fix (author-ruled 2026-08-27): when margins
* marks a subgroup-by-horizon margin as not estimable, that horizon is dropped from
* THAT subgroup's event study with a caution naming it, the subgroup's remaining
* horizons and every other subgroup are unaffected, and the run completes (rc 0).
* Before the fix the omission marker survived the sanitizer and hit an arithmetic
* evaluator, aborting the whole run with a bare "5o. invalid name" (r198).
* Run from test/:  "<stata>" -e do test_si7.do   -> expect SI7_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 7717
set sortseed 7717

local REPO ".."
do "`REPO'/wooldidinstall/wooldid.ado"

global pf 0

cap program drop scanlog
program define scanlog, rclass
  args lf needle
  tempname fh
  local found 0
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') != 0 local found 1
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
end

* rows of a subgroup's result matrix, and whether a given row label is present
cap program drop sgrows
program define sgrows, rclass
  args sg lbl
  tempname S
  cap matrix `S' = e(wdidsgresults_`sg')
  if _rc != 0 {
    return scalar n = -1
    return scalar has = -1
    exit
  }
  local rn : rownames `S'
  local has 0
  foreach r of local rn {
    if "`r'" == "`lbl'" local has 1
  }
  return scalar n = rowsof(`S')
  return scalar has = `has'
end

* ---- T1: the pinned config, where a subgroup pre-period margin is not estimable ----
* On this fixture the non-estimability is intrinsic to the data (the documented
* knife-edge design), not caused by fe() - the recorded bisection shows dropping
* fe() does NOT avoid it. Expect: rc 0, the caution names pre2 and sg1, sg1 loses
* ONLY that horizon, sg2 keeps all of its. fullprecisionslow because safety(off)
* with the default cell collapse keeps the collapsed rows (documented), on which the
* horizon happens to be estimable (bug sweep 5c-1).
use suite/data/survey_rcs.dta, clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
qui log using si7_t1.txt, text replace name(si7)
cap noi wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) ///
    fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) ///
    safety(off) adversarial fullprecisionslow
local rc1 = _rc
qui log close si7
scanlog si7_t1.txt "are not estimable for subgroup 1"
local caut = r(found)
scanlog si7_t1.txt "invalid name"
local crash = r(found)
* the adversarial line test must be REFUSED for a subgroup with a dropped horizon:
* the null is that the coefficients lie on a line in relative time, which an
* incomplete (possibly interior-holed) grid cannot express
scanlog si7_t1.txt "adversarial pre-trend test was skipped for subgroup 1"
local advskip = r(found)
sgrows 1 pre2
local n1 = r(n)
local has1 = r(has)
sgrows 2 pre2
local n2 = r(n)
local has2 = r(has)
cap erase si7_t1.txt
if `rc1' == 0 & `caut' == 1 & `crash' == 0 & `has1' == 0 & `has2' == 1 & `n1' < `n2' & `n1' > 0 & `advskip' == 1 {
  dis as result "T1 PASS: rc 0, caution named sg1/pre2, sg1 dropped it (`n1' rows) and sg2 kept it (`n2' rows), adversarial refused for sg1"
}
else {
  dis as error "T1 FAIL: rc=`rc1' caution=`caut' crashtext=`crash' sg1rows=`n1' sg1haspre2=`has1' sg2rows=`n2' sg2haspre2=`has2' advskip=`advskip'"
  global pf = $pf + 1
}

* ---- T2: a clean synthetic panel where every subgroup horizon IS estimable ----
* the guard must stay completely dormant: no caution, every subgroup keeps pre2
clear
set seed 7717
qui set obs 16
gen i = _n
gen sg = mod(_n, 2) + 1
gen ttre = cond(mod(_n, 4) == 0 | mod(_n, 4) == 1, 6, .)
qui expand 6
bysort i: gen unit = _n
gen unitid = (i-1)*6 + unit
qui expand 10
bysort i unitid: gen t = _n
gen y = 0.05*i + 0.04*t + 0.5*(ttre < . & t >= ttre) + rnormal()
sort i unitid t
ereturn clear
qui log using si7_t2.txt, text replace name(si7)
cap noi wooldid y i t ttre, att subgroups(sg) esfixedbaseperiod esprelength(1) espostlength(1) robust
local rc2 = _rc
qui log close si7
scanlog si7_t2.txt "are not estimable for subgroup"
local caut2 = r(found)
sgrows 1 pre2
local h1b = r(has)
sgrows 2 pre2
local h2b = r(has)
cap erase si7_t2.txt
if `rc2' == 0 & `caut2' == 0 & `h1b' == 1 & `h2b' == 1 {
  dis as result "T2 PASS: healthy panel -> guard dormant, both subgroups keep pre2"
}
else {
  dis as error "T2 FAIL: rc=`rc2' caution=`caut2' sg1haspre2=`h1b' sg2haspre2=`h2b'"
  global pf = $pf + 1
}


* ---- T3 (bug sweep BS-26): fast mode keeps the SI-7 drop ----
* fastnose runs every other margins call with nose, which suppresses the omission
* marker the SI-7 scan reads; the subgroup-by-horizon call is exempt, so the drop and
* its caution survive fast mode and sg1 keeps the same row set as the analytic run
use suite/data/survey_rcs.dta, clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
qui log using si7_t3.txt, text replace name(si7)
cap noi wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) ///
    fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) ///
    fastnose fullprecisionslow
local rc3 = _rc
qui log close si7
scanlog si7_t3.txt "are not estimable for subgroup 1"
local caut3 = r(found)
sgrows 1 pre2
local n3 = r(n)
local has3 = r(has)
cap erase si7_t3.txt
if `rc3' == 0 & `caut3' == 1 & `has3' == 0 & `n3' == `n1' {
  dis as result "T3 PASS: fastnose keeps the SI-7 drop (caution printed, sg1 has `n3' rows, no pre2)"
}
else {
  dis as error "T3 FAIL: rc=`rc3' caution=`caut3' sg1rows=`n3' (analytic `n1') sg1haspre2=`has3'"
  global pf = $pf + 1
}

* ---- T4-T7 (BS-26): the replication engines on the same fixture ----
* the observed run and the replications must agree on the column set: before the fix
* the reps (run under fastnose) kept the non-estimable horizon while the observed
* analytic run dropped it, and bootstrap()/permute() permanalyticse aborted r(459) at
* rep 0, while plain permute() reported and ranked the structural-zero row
local k 3
foreach leg in "permute(shuffle) permreps(20) rseed(11)" "bootstrap(standard) reps(10) rseed(11)" "bootstrap(bayes) reps(10) rseed(11)" "permute(shuffle) permreps(20) rseed(11) permanalyticse" {
  local ++k
  use suite/data/survey_rcs.dta, clear
  qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
  cap noi wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) ///
      fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state) `leg'
  local rck = _rc
  sgrows 1 pre2
  local nk = r(n)
  local hask = r(has)
  if `rck' == 0 & `hask' == 0 & `nk' == `n1' {
    dis as result "T`k' PASS: `leg' -> rc 0, sg1 has `nk' rows and no pre2"
  }
  else {
    dis as error "T`k' FAIL: `leg' -> rc=`rck' sg1rows=`nk' (analytic `n1') sg1haspre2=`hask'"
    global pf = $pf + 1
  }
}

if $pf == 0 dis as result "SI7_ALL_PASS"
else dis as error "SI7_FAILURES: $pf"
exit, clear
