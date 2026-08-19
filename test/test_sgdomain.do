* Targeted suite for the subgroup sign-screen domain harmonization (VAL-G15A): the
* negative-sg check moved from wooldid_validate (raw if-domain, no markout knowledge)
* into the worker beside the integer screen, on `ifstatement & internal-ttre < .`.
* Run from test/: "<stata>" -e do test_sgdomain.do -> SGDOM_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 24680
set sortseed 24680
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop mksg
program define mksg
  use "data/dta_binary", clear
end

* ---- S1: false-positive cure - sg=-1 confined to one missing-y treated-period row:
* was r(411); now rc 0 and bit-identical to the run with that row deleted ----
tempname B1 B2
mksg
qui su i if ttre < . & t >= ttre, meanonly
local i0 = r(min)
qui su t if ttre < . & t >= ttre & i == `i0', meanonly
local t0 = r(min)
qui replace y = . if i == `i0' & t == `t0'
qui replace sg = -1 if i == `i0' & t == `t0'
qui count if sg == -1
local npoison = r(N)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 | `npoison' == 0 {
  dis as error "S1 FAIL: rc=" _rc " npoison=" `npoison' " (expected rc 0 on poisoned markout row)"
  local ++fails
}
else {
  matrix `B1' = e(b)
  mksg
  qui drop if i == `i0' & t == `t0'
  qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
  matrix `B2' = e(b)
  if mreldif(`B1', `B2') != 0 {
    dis as error "S1 FAIL: poisoned-markout-row twin not bit-identical, mreldif=" mreldif(`B1', `B2')
    local ++fails
  }
}

* ---- S2: sg=-1 on a live treated-post row still errors r(411) ----
mksg
qui replace sg = -1 if ttre < . & t >= ttre & _n < 50
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 411 {
  dis as error "S2 FAIL: live treated-post negative sg rc=" _rc " (expected 411)"
  local ++fails
}

* ---- S3: sg=-1 on a live PRE-onset row of a treated unit still errors r(411)
* (the SG pre-trend margins consume pre-onset rows; the t>=ttre domain would miss this) ----
mksg
qui replace sg = -1 if ttre < . & t < ttre & _n < 300
qui count if sg == -1
local npre = r(N)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 411 | `npre' == 0 {
  dis as error "S3 FAIL: pre-onset negative sg rc=" _rc " npre=" `npre' " (expected 411)"
  local ++fails
}

* ---- S4: newly-tolerated class certification - sg=-1 confined to never-treated rows
* runs rc 0 bit-identical to the same run with sg=. on those rows ----
tempname B4 B5
mksg
qui replace sg = -1 if ttre >= . & mod(i, 7) == 0
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
matrix `B4' = e(b)
mksg
qui replace sg = . if ttre >= . & mod(i, 7) == 0
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
matrix `B5' = e(b)
if mreldif(`B4', `B5') != 0 {
  dis as error "S4 FAIL: never-treated sg=-1 vs sg=. twins differ, mreldif=" mreldif(`B4', `B5')
  local ++fails
}

* ---- S5: doubly-bad value (sg = -1.5) on a live treated-post row: negative message
* fires first (sign screen precedes the integer screen) ----
mksg
qui log using sgdom_s5.txt, text replace name(sgd)
qui replace sg = -1.5 if ttre < . & t >= ttre & _n < 50
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
local rc5 = _rc
qui log close sgd
if `rc5' != 411 {
  dis as error "S5 FAIL: doubly-bad sg rc=" `rc5' " (expected 411)"
  local ++fails
}
else if strpos(fileread("sgdom_s5.txt"), "contains negative values") == 0 {
  dis as error "S5 FAIL: negative-message-first precedence lost for doubly-bad values"
  local ++fails
}

* ---- S6: clean sg run unchanged (screen is domain-relaxation only) ----
mksg
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 {
  dis as error "S6 FAIL: clean sg run rc=" _rc
  local ++fails
}

if `fails' == 0 dis as result "SGDOM_ALL_PASS"
else dis as error "SGDOM_FAILS=`fails'"
