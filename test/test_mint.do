* Targeted suite for the mint-staleness fix: the no-pre-period drop is unit-level by
* construction and mint/mintp1 are recomputed post-drop, so sentinel-coded ttre under
* safety(off) can no longer strand surviving units without their reference period.
* Run from test/:  "<stata>" -e do test_mint.do   -> expect MINT_ALL_PASS
* On the pre-fix build, M1/M4 return silently biased estimates instead of the bench values.
clear all
set more off
set maxvar 20000
set graphics off
local fails 0
do "../wooldidinstall/wooldid.ado"

* panel: 10 treated (ttre=5) + 20 never; unit 1 optionally carries sentinel-coded ttre
* (0 on not-yet-treated rows), the common raw-data convention F13 guards against
cap program drop mkmint
program define mkmint
  args sentinel
  clear
  set seed 4242
  set sortseed 4242
  qui set obs 30
  qui gen i = _n
  qui gen ttre = cond(_n <= 10, 5, .)
  qui expand 10
  qui bysort i: gen t = _n
  qui gen te = 2*(t >= ttre & ttre < .)
  qui gen y = te + 0.05*i + rnormal()
  qui drop te
  if "`sentinel'" == "one" qui replace ttre = 0 if i == 1 & t < 5
  if "`sentinel'" == "all" qui replace ttre = 0 if ttre == 5 & t < 5
end

* ---- M1: one sentinel-coded unit + safety(off): equals physically deleting that unit ----
mkmint one
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "M1 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname b1
  matrix `b1' = e(b)
  mkmint one
  qui drop if i == 1
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(off)
  if _rc != 0 | mreldif(`b1', e(b)) != 0 {
    dis as error "M1 FAIL: bench rc=" _rc " mreldif=" mreldif(`b1', e(b))
    local ++fails
  }
}

* ---- M2: every treated unit sentinel-coded + safety(off): informative refusal, not a sign flip ----
mkmint all
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(off)
if _rc != 459 {
  dis as error "M2 FAIL: rc=" _rc " (expected 459: all treated units lose their pre-period)"
  local ++fails
}

* ---- M3: clean data + safety(off): estimates match the safety-on run exactly ----
mkmint none
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "M3 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname b3
  matrix `b3' = e(b)
  mkmint none
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
  if _rc != 0 | mreldif(`b3', e(b)) != 0 {
    dis as error "M3 FAIL: safety-on rc=" _rc " mreldif=" mreldif(`b3', e(b))
    local ++fails
  }
}

* ---- M4: the timetrends (mintp1) twin of M1 ----
mkmint one
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends safety(off)
if _rc != 0 {
  dis as error "M4 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname b4
  matrix `b4' = e(b)
  mkmint one
  qui drop if i == 1
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends safety(off)
  if _rc != 0 | mreldif(`b4', e(b)) != 0 {
    dis as error "M4 FAIL: bench rc=" _rc " mreldif=" mreldif(`b4', e(b))
    local ++fails
  }
}

* ---- M5: a 0/1 dummy passed where an onset period is expected + safety(off): refuses ----
mkmint none
qui replace ttre = (ttre < .)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) safety(off)
if _rc == 0 {
  dis as error "M5 FAIL: rc 0 on a 0/1 dummy passed as the onset variable"
  local ++fails
}

if `fails' == 0 dis "MINT_ALL_PASS"
else dis as error "MINT_FAILS=`fails'"
exit, clear
