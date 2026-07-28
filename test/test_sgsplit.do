* Targeted suite for the subgroup-split gate (author directive 2026-07-28): a unit
* split across subgroups over time is refused under default safety, and proceeds with
* share weighting plus a caution under safety(off). Also pins the numeric-only
* cluster() syntax. Run from test/: "<stata>" -e do test_sgsplit.do -> SGSPLIT_ALL_PASS
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

* panel with one unit whose subgroup changes over time
cap program drop mksplit
program define mksplit
  clear
  set seed 24680
  set sortseed 24680
  qui set obs 30
  qui gen i = _n
  qui gen ttre = cond(_n <= 10, 5, .)
  qui expand 10
  qui bysort i: gen t = _n
  qui gen y = 2*(t >= ttre & ttre < .) + 0.1*i + rnormal()
  qui gen sg = mod(i, 2)
  qui replace sg = 1 - sg if i == 3 & t >= 6
end

* ---- G1: split unit + default safety -> informative r(459) ----
mksplit
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 459 {
  dis as error "G1 FAIL: rc=" _rc " (expected 459)"
  local ++fails
}

* ---- G2: split unit + safety(off) -> runs with the share-weighting caution ----
mksplit
qui log using sgsp_g2.txt, text replace name(sgsp)
cap noi wooldid y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) subgroups(sg) safety(off)
local rc2 = _rc
qui log close sgsp
if `rc2' != 0 {
  dis as error "G2 FAIL: rc=" `rc2'
  local ++fails
}
else if strpos(fileread("sgsp_g2.txt"), "splits some units") == 0 | strpos(fileread("sgsp_g2.txt"), "share weighting") == 0 {
  dis as error "G2 FAIL: share-weighting caution missing under safety(off)"
  local ++fails
}

* ---- G3: unit-constant subgroups unaffected, both safety settings ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 {
  dis as error "G3 FAIL: clean sg run rc=" _rc
  local ++fails
}
tempname b3
matrix `b3' = e(b)
use "`DATA'/dta_binary", clear
qui log using sgsp_g3.txt, text replace name(sgsp)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) safety(off)
qui log close sgsp
if _rc != 0 | mreldif(`b3', e(b)) != 0 {
  dis as error "G3 FAIL: safety(off) clean sg run rc=" _rc " mreldif=" mreldif(`b3', e(b))
  local ++fails
}
if strpos(fileread("sgsp_g3.txt"), "splits some units") != 0 {
  dis as error "G3 FAIL: split caution fired on unit-constant subgroups"
  local ++fails
}

* ---- G4: string cluster variables are refused at parse (numeric-only syntax) ----
use "`DATA'/dta_binary", clear
qui gen str8 sclus = "c" + string(i)
cap wooldid y i t ttre, att cluster(sclus)
if _rc != 109 {
  dis as error "G4 FAIL: rc=" _rc " (expected 109 from varlist numeric)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i)
if _rc != 0 {
  dis as error "G4 FAIL: numeric cluster run rc=" _rc
  local ++fails
}

if `fails' == 0 dis "SGSPLIT_ALL_PASS"
else dis as error "SGSPLIT_FAILS=`fails'"
exit, clear
