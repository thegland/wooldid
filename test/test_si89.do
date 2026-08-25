* test_si89.do — certification of SI-8 (granular-i caution + maxvar preflight,
* author-ruled 2026-08-25) and SI-9 (poisexpresults requires poisson, rc 198).
* The poisexpresults+poisson POSITIVE path is pinned by the oldnew suite
* (ext_bin_ext_poisson, ext_cont_ext_poisson) and re-checked here as T2.
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

* ---- T1: poisexpresults without poisson is an option-contract error (rc 198) ----
mkgran 30
cap wooldid y i t ttre, att robust poisexpresults
if _rc == 198 dis as result "T1 PASS: poisexpresults without poisson -> r(198)"
else {
  dis as error "T1 FAIL: rc = " _rc " (wanted 198)"
  global pf = $pf + 1
}

* ---- T2: poisexpresults WITH poisson still estimates (suite-pinned positive) ----
use suite/data/binary_balanced.dta, clear
ereturn clear
cap wooldid y_count_zero cohort t ttre, att poisson poisexpresults customvce(robust) esprelength(2) espostlength(2) safety(off)
if _rc == 0 & _b[main_att] < . dis as result "T2 PASS: poisexpresults+poisson estimates (att = " %9.6f _b[main_att] ")"
else {
  dis as error "T2 FAIL: rc = " _rc
  global pf = $pf + 1
}

* ---- T3: >100 distinct i values -> granular-i caution, estimation proceeds ----
mkgran 150
cap log close si89
qui log using si89_t3.txt, text replace name(si89)
cap noi wooldid y i t ttre, att robust
qui log close si89
local rc3 = _rc
scanlog si89_t3.txt "distinct values. wooldid builds cohort-by-period"
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
scanlog si89_t4.txt "distinct values. wooldid builds cohort-by-period"
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
