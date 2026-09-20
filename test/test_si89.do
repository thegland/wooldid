* test_si89.do — certification of SI-8 (granular-i caution + maxvar preflight,
* author-ruled 2026-08-25) and, since 2026-09-20, the retirement of the options
* SI-9 used to police. poisexpresults, semielasticity and the eyex/eydx forms of
* contreatelasticitytype were removed (author ruling 2026-09-20); T1 and T2 pin
* that they are now rejected by name. This file is the cert HOST for that pin:
* test/suite/configs_ext.csv cannot carry it, because coverage_check.py treats an
* option token it cannot resolve against the syntax line as FATAL (author ruling
* 2026-09-11), so a retired NAME cannot appear in a config row at all.
* Run from test/:  "<stata>" -e do test_si89.do   -> expect SI89_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 813
set sortseed 813

local REPO ".."
do "`REPO'/wooldidinstall/wooldid.ado"

global pf 0

cap program drop mkgran
program define mkgran
  args n
  clear
  set seed 813
  qui set obs `n'
  gen i = _n
  gen ttre = cond(mod(_n, 2) == 0, 6, .)
  qui expand 10
  bysort i: gen t = _n
  gen dose = exp(0.3*runiform())
  gen y = 0.02*i + 0.05*t + 0.4*(ttre < . & t >= ttre) + rnormal()
  sort i t
  ereturn clear
end

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

* ---- T1: the retired poisson-scale options are rejected by name (rc 198) ----
* Stata's own syntax parser raises 198 once the option leaves the syntax line; this is
* the clean-break signal the fastcollapse retirement set as precedent
* (review_process/option_consolidation_2026-09-16.md). Both spellings are checked with
* AND without poisson, since SI-9's old contract screen distinguished those two cases
* and the retirement must collapse them to the same answer.
mkgran 30
local t1bad 0
foreach o in "poisexpresults" "poisson poisexpresults" "semielasticity" "poisson semielasticity" {
  cap wooldid y i t ttre, att robust `o'
  if _rc != 198 {
    dis as error "T1 FAIL: `o' -> rc " _rc " (wanted 198)"
    local t1bad = `t1bad' + 1
  }
}
if `t1bad' == 0 dis as result "T1 PASS: retired poisexpresults/semielasticity -> r(198) in all 4 forms"
else global pf = $pf + 1

* ---- T2: the retired contreatelasticitytype VALUES are refused, survivors still run ----
* Here the option NAME survives, so the refusal is wooldid's own domain screen rather
* than Stata's parser, and it must name the removed forms rather than collapsing them
* silently to dydx (the pre-2026-09-06 behaviour).
use suite/data/continuous_control_compact.dta, clear
local t2bad 0
foreach v in eyex eydx {
  cap wooldid y_positive_common cohort t ttre, att contreat(dose_common) contreatcontrols(basic) ///
      contreatelasticitytype(`v') robust safety(off)
  if _rc != 198 {
    dis as error "T2 FAIL: contreatelasticitytype(`v') -> rc " _rc " (wanted 198)"
    local t2bad = `t2bad' + 1
  }
}
foreach v in dydx dyex {
  ereturn clear
  cap wooldid y_positive_common cohort t ttre, att contreat(dose_common) contreatcontrols(basic) ///
      contreatelasticitytype(`v') robust safety(off)
  if _rc != 0 {
    dis as error "T2 FAIL: contreatelasticitytype(`v') -> rc " _rc " (wanted 0)"
    local t2bad = `t2bad' + 1
  }
}
if `t2bad' == 0 dis as result "T2 PASS: eyex/eydx refused r(198); dydx/dyex still estimate"
else global pf = $pf + 1

* ---- T3: >100 distinct i values -> granular-i caution, estimation proceeds ----
mkgran 150
cap log close si89
qui log using si89_t3.txt, text replace name(si89)
cap noi wooldid y i t ttre, att robust
qui log close si89
local rc3 = _rc
scanlog si89_t3.txt "distinct values. wooldid builds"
if `rc3' == 0 & r(found) == 1 dis as result "T3 PASS: 150-level i -> caution shown, run completes"
else {
  dis as error "T3 FAIL: rc = `rc3', caution found = " r(found)
  global pf = $pf + 1
}

* ---- T4: safety(off) suppresses the granular-i caution ----
mkgran 150
qui log using si89_t4.txt, text replace name(si89)
cap noi wooldid y i t ttre, att robust safety(off)
qui log close si89
local rc4 = _rc
scanlog si89_t4.txt "distinct values. wooldid builds"
if `rc4' == 0 & r(found) == 0 dis as result "T4 PASS: safety(off) -> no caution, run completes"
else {
  dis as error "T4 FAIL: rc = `rc4', caution found = " r(found)
  global pf = $pf + 1
}

* ---- T5: maxvar preflight -> early r(907) with the fe()/coarser-i steer ----
* set maxvar refuses to run while data/matrices/stored estimates linger
clear
ereturn clear
estimates clear
matrix drop _all
clear mata
set maxvar 2500
* 200 treated cohorts x (5 post cells + 1 pooled-reference cell) = 1200 cells;
* demand = 1200*2 + 400 + 10 + 200 = 3010 > 2500 with real margin
mkgran 400
qui log using si89_t5.txt, text replace name(si89)
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(basic)
qui log close si89
local rc5 = _rc
scanlog si89_t5.txt "so estimation would fail with r(907)"
local f5a = r(found)
scanlog si89_t5.txt "Consider specifying unit-level fixed effects via fe()"
if `rc5' == 907 & `f5a' == 1 & r(found) == 1 dis as result "T5 PASS: maxvar preflight -> early r(907) + steer"
else {
  dis as error "T5 FAIL: rc = `rc5' (wanted 907), demand line = `f5a', steer line = " r(found)
  global pf = $pf + 1
}
clear
ereturn clear
estimates clear
matrix drop _all
clear mata
set maxvar 20000

* ---- verdict ----
cap erase si89_t3.txt
cap erase si89_t4.txt
cap erase si89_t5.txt
if $pf == 0 dis as result "SI89_ALL_PASS"
else dis as error "SI89_FAILURES: $pf"
exit, clear
