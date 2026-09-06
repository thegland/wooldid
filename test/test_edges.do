* test_edges.do - edge-of-panel shapes (bug sweep BS-33, 2026-09-06).
* Run from test/: stata-mp -e do test_edges.do -> expect EDGES_ALL_PASS
*
* E1 panel A (the onset period dropped for the treated cohort, a common convention):
*    bare att -> rc 0 with att = 1.5 (the noise-free DGP effect); att espostlength(2) -> rc 0,
*    att = 1.5, and the only ES row is relative period 1 (0 is not realized); no pre-att row
* E2 panel B (biennial t = 1,3,..,11 with onset 6, so relative period 0 never occurs):
*    bare att -> rc 0, att = 1.5; espostlength(2) -> rc 0 with the single ES row at 1
* E3 panel B espostlength(1): the requested window is relative period 0 only, which nothing
*    realizes -> the informative error 459, not the bare r(198) from matrix ()
* E4 panel A esprelength(1) espostlength(2): rc 0 (relative period -1 survives; 4 rows)
* E5 regression pin: dta_binary att cluster(i) esprelength(2) espostlength(2) is unmoved
*    (att .8235285, se .0629643, six rows: att, pre att, pre2 pre1 post0 post1; e(wdid_espre) 2)
clear all
set more off
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop rowk
program rowk, rclass
  args mname rname
  local rn : rownames `mname'
  local k : list posof "`rname'" in rn
  return scalar k = `k'
end

* panel A: 12 units x 6 periods, onset 4 for units 1-4, onset row dropped
clear
set obs 12
gen i = _n
gen ttre = cond(i <= 4, 4, .)
expand 6
bysort i: gen t = _n
gen y = 1.5*(t >= ttre & ttre < .) + 0.3*i + 0.2*t
drop if i <= 4 & t == 4
tempfile pa
save `pa'

* panel B: biennial periods, onset 6 (never a realized period)
clear
set obs 12
gen i = _n
gen ttre = cond(i <= 4, 6, .)
expand 6
bysort i: gen t = 2*_n - 1
gen y = 1.5*(t >= ttre & ttre < .) + 0.3*i + 0.2*t
tempfile pb
save `pb'

* ---- E1
use `pa', clear
cap noi wooldid y i t ttre, att
local rc1a = _rc
local att1a = .
if `rc1a' == 0 {
  matrix E1A = e(wdidmainresults)
  rowk E1A att
  local att1a = E1A[r(k), 1]
}
use `pa', clear
cap noi wooldid y i t ttre, att espostlength(2)
local rc1b = _rc
local att1b = .
local rn1b
if `rc1b' == 0 {
  matrix E1B = e(wdidmainresults)
  rowk E1B att
  local att1b = E1B[r(k), 1]
  local rn1b : rownames E1B
}
di "E1 bare rc=`rc1a' att=`att1a' | espostlength(2) rc=`rc1b' att=`att1b' rows=[`rn1b']"
if `rc1a' != 0 | abs(`att1a' - 1.5) > 1e-6 | `rc1b' != 0 | abs(`att1b' - 1.5) > 1e-6 | rowsof(E1B) != 2 {
  di as error "E1 FAIL: the onset-period gap still breaks a bare att call"
  local ++fails
}

* ---- E2
use `pb', clear
cap noi wooldid y i t ttre, att
local rc2a = _rc
local att2a = .
if `rc2a' == 0 {
  matrix E2A = e(wdidmainresults)
  rowk E2A att
  local att2a = E2A[r(k), 1]
}
use `pb', clear
cap noi wooldid y i t ttre, att espostlength(2)
local rc2b = _rc
local att2b = .
local rn2b
if `rc2b' == 0 {
  matrix E2B = e(wdidmainresults)
  rowk E2B att
  local att2b = E2B[r(k), 1]
  local rn2b : rownames E2B
}
di "E2 bare rc=`rc2a' att=`att2a' | espostlength(2) rc=`rc2b' att=`att2b' rows=[`rn2b']"
if `rc2a' != 0 | abs(`att2a' - 1.5) > 1e-6 | `rc2b' != 0 | abs(`att2b' - 1.5) > 1e-6 | rowsof(E2B) != 2 {
  di as error "E2 FAIL: the off-grid onset still breaks att"
  local ++fails
}

* ---- E3
use `pb', clear
cap noi wooldid y i t ttre, att espostlength(1)
local rc3 = _rc
di "E3 espostlength(1) on the biennial panel rc=`rc3' (expect the informative 459)"
if `rc3' != 459 {
  di as error "E3 FAIL: an unrealizable window should give the informative error 459, got `rc3'"
  local ++fails
}

* ---- E4
use `pa', clear
cap noi wooldid y i t ttre, att esprelength(1) espostlength(2)
local rc4 = _rc
local att4 = .
local rn4
if `rc4' == 0 {
  matrix E4 = e(wdidmainresults)
  rowk E4 att
  local att4 = E4[r(k), 1]
  local rn4 : rownames E4
}
di "E4 rc=`rc4' att=`att4' rows=[`rn4']"
if `rc4' != 0 | abs(`att4' - 1.5) > 1e-6 | rowsof(E4) != 4 {
  di as error "E4 FAIL: the surviving pre period should still run"
  local ++fails
}

* ---- E5
use data/dta_binary.dta, clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
local rc5 = _rc
local b5 = .
local se5 = .
local rn5
if `rc5' == 0 {
  matrix E5 = e(wdidmainresults)
  rowk E5 att
  local b5 = E5[r(k), 1]
  local se5 = E5[r(k), 2]
  local rn5 : rownames E5
  local espre5 = e(wdid_espre)
}
di "E5 rc=`rc5' att=" %10.7f `b5' " se=" %10.7f `se5' " rows=[`rn5'] e(wdid_espre)=`espre5'"
if `rc5' != 0 | abs(`b5' - .8235285) > 5e-8 | abs(`se5' - .0629643) > 5e-8 | rowsof(E5) != 6 | `espre5' != 2 {
  di as error "E5 FAIL: dta_binary moved"
  local ++fails
}

if `fails' == 0 di as result "EDGES_ALL_PASS"
else di as error "EDGES_FAILS=`fails'"
