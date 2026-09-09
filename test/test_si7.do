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
  * Stata wraps a long display line in the log and continues it after "> ". A needle that
  * spans the wrap is then on no single line, which is how T10 failed while its shorter
  * twin passed. Rejoin continuations before searching.
  tempname fh
  local found 0
  local acc
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if substr(`"`macval(line)'"', 1, 2) == "> " {
      local acc `"`macval(acc)'`=substr(`"`macval(line)'"', 3, .)'"'
    }
    else {
      if strpos(`"`macval(acc)'"', `"`needle'"') != 0 local found 1
      local acc `"`macval(line)'"'
    }
    file read `fh' line
  }
  if strpos(`"`macval(acc)'"', `"`needle'"') != 0 local found 1
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


* ---- T8-T10 (bug sweep BS-8): the contreat subgroup event study under an SI-7 drop ----
* the continuous-treatment (AME) horizon call was labelled from the binary call's
* already-filtered relyr column, so any SI-7 drop was a bare conformability error
* r(503) right after the caution. The AME call is now scanned and labelled on its own.
* T8: the sweep's fixture (dose constant within cohort-period cells, so no AME horizon
*     is identified): rc 0, the binary caution, both "skipped entirely" cautions and
*     empty ES_ame blocks. T9/T10: a dose with within-cell variation everywhere except
*     subgroup 2's pre2 cells, so exactly that AME horizon is dropped: sg2 ES_ame (and
*     ES_ame_it) = post0 with the caution naming it, sg1 keeps pre2 post0, the binary
*     rows are untouched, and no label repeats (Stata pads a short namelist by repeating
*     its last name, so a duplicate is the signature of the old mislabelling)
cap program drop ameblock
program define ameblock, rclass
  args sg pfx
  tempname S
  cap matrix `S' = e(wdidsgresults_`sg')
  if _rc != 0 {
    return scalar n = -1
    return scalar dup = -1
    exit
  }
  local rn : rowfullnames `S'
  local n 0
  local dup 0
  local seen
  local rows
  foreach r of local rn {
    if strpos("`r'", "`pfx':") == 1 {
      local ++n
      if strpos(" `seen' ", " `r' ") > 0 local dup 1
      local seen `seen' `r'
      local rows `rows' `=subinstr("`r'", "`pfx':", "", .)'
    }
  }
  return scalar n = `n'
  return scalar dup = `dup'
  return local rows `rows'
end
use suite/data/survey_rcs.dta, clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
qui gen double dose = 1 + 0.2*mod(state,7) + 0.05*t
qui log using si7_t8.txt, text replace name(si7)
cap noi wooldid y_rcs state_group t ttre [pw=survey_weight], att subgroups(group) contreat(dose) ///
    contreatcontrols(basic) fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state)
local rc8 = _rc
qui log close si7
scanlog si7_t8.txt "the relative period effect(s) pre2 are not estimable for subgroup 2"
local caut8 = r(found)
scanlog si7_t8.txt "no continuous-treatment relative period effect is estimable for subgroup 1"
local skip8a = r(found)
scanlog si7_t8.txt "no continuous-treatment relative period effect is estimable for subgroup 2"
local skip8b = r(found)
scanlog si7_t8.txt "conformability"
local conf8 = r(found)
cap erase si7_t8.txt
ameblock 1 ES_ame
local n81 = r(n)
ameblock 2 ES_ame
local n82 = r(n)
if `rc8' == 0 & `caut8' == 1 & `skip8a' == 1 & `skip8b' == 1 & `conf8' == 0 & `n81' == 0 & `n82' == 0 {
  dis as result "T8 PASS: the sweep's fixture runs rc 0 (no AME horizon identified: both contreat event studies skipped with a caution, binary caution intact)"
}
else {
  dis as error "T8 FAIL: rc=`rc8' binarycaution=`caut8' skip1=`skip8a' skip2=`skip8b' conformability=`conf8' ES_ame n1=`n81' n2=`n82'"
  global pf = $pf + 1
}

use suite/data/survey_rcs.dta, clear
qui replace y_rcs = y_rcs - rcs_cell_shock if y_rcs < .
set seed 8008
qui gen double dose = 1 + 0.2*mod(state,7) + 0.05*t + rnormal()
qui replace dose = 1 + 0.2*mod(state,7) + 0.05*t if group == 2 & ttre < . & t - ttre == -2
tempfile fx9
save `fx9'
local k 8
foreach leg in "att" "att att_it" {
  local ++k
  use `fx9', clear
  qui log using si7_t9.txt, text replace name(si7)
  cap noi wooldid y_rcs state_group t ttre [pw=survey_weight], `leg' subgroups(group) contreat(dose) ///
      contreatcontrols(basic) fe(i.group##i.t) esfixedbaseperiod esprelength(1) espostlength(1) cluster(state)
  local rck = _rc
  qui log close si7
  scanlog si7_t9.txt "the continuous-treatment relative period effect(s) pre2 are not estimable for subgroup 2"
  local cautk = r(found)
  scanlog si7_t9.txt "not estimable for subgroup 1"
  local caut1k = r(found)
  scanlog si7_t9.txt "conformability"
  local confk = r(found)
  cap erase si7_t9.txt
  local okk 1
  local pfxs ES_ame
  if "`leg'" == "att att_it" local pfxs ES_ame ES_ame_it
  local rep
  foreach pfx of local pfxs {
    ameblock 1 `pfx'
    local r1 `r(rows)'
    local d1 = r(dup)
    ameblock 2 `pfx'
    local r2 `r(rows)'
    local d2 = r(dup)
    local rep `rep' `pfx' sg1=[`r1'] sg2=[`r2']
    if "`r1'" != "pre2 post0" | "`r2'" != "post0" | `d1' != 0 | `d2' != 0 local okk 0
  }
  ameblock 1 ES_att
  local b1 `r(rows)'
  ameblock 2 ES_att
  local b2 `r(rows)'
  if "`b1'" != "pre2 post0" | "`b2'" != "pre2 post0" local okk 0
  if `rck' == 0 & `cautk' == 1 & `caut1k' == 0 & `confk' == 0 & `okk' {
    dis as result "T`k' PASS: contreat + `leg' -> rc 0, sg2's pre2 AME dropped with its caution;`rep'; ES_att sg1=[`b1'] sg2=[`b2']"
  }
  else {
    dis as error "T`k' FAIL: contreat + `leg' -> rc=`rck' caution2=`cautk' caution1=`caut1k' conformability=`confk';`rep'; ES_att sg1=[`b1'] sg2=[`b2']"
    global pf = $pf + 1
  }
}

if $pf == 0 dis as result "SI7_ALL_PASS"
else dis as error "SI7_FAILURES: $pf"
exit, clear
