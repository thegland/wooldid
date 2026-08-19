* Targeted suite for the cfxlattice() screen domain narrowing (VAL-G12A/VAL-G12B):
* both content screens now tolerate rows with a missing recorded onset (never-treated
* rows feed no executed CFX subpop), curing false rejections / false cautions, while
* keeping full protection on nonmissing-onset rows.
* Run from test/: "<stata>" -e do test_cfxlat.do -> CFXLAT_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 13579
set sortseed 13579
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop mkcfx
program define mkcfx
  use "data/dta_cont", clear
  qui gen byte clat = mod(i, 3)
end

* ---- L1: bad value confined to one never-treated row -> rc 0 with the CFX block
* actually executing (was r411 at parse) ----
mkcfx
qui replace clat = 0.5 if ttre >= . & _n == _N
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
if _rc != 0 {
  dis as error "L1 FAIL: bad-value-on-never-treated rc=" _rc " (expected 0)"
  local ++fails
}
else {
  cap confirm matrix e(cfxatt_cfxpost)
  if _rc != 0 {
    dis as error "L1 FAIL: rc 0 but no e(cfxatt_cfxpost) - CFX block did not execute (vacuous pass)"
    local ++fails
  }
}

* ---- L2: twin parity - the poisoned never-treated row vs a benign integer there;
* returned cfx matrices must be bit-identical ----
tempname T1 T2
mkcfx
qui replace clat = 0.5 if ttre >= . & _n == _N
qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
matrix `T1' = e(cfxatt_cfxpost)
mkcfx
qui replace clat = 1 if ttre >= . & _n == _N
qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
matrix `T2' = e(cfxatt_cfxpost)
if mreldif(`T1', `T2') != 0 {
  dis as error "L2 FAIL: never-treated-row value moved cfx results, mreldif=" mreldif(`T1', `T2')
  local ++fails
}

* ---- L3: bad value on a treated-unit row still errors r(411) ----
mkcfx
qui replace clat = 0.5 if ttre < . & _n == 1
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
if _rc != 411 {
  dis as error "L3 FAIL: bad-value-on-treated rc=" _rc " (expected 411)"
  local ++fails
}

* ---- L4: missing on never-treated only -> NO caution; missing on a treated row -> caution ----
mkcfx
qui replace clat = . if ttre >= . & _n == _N
qui log using cfxlat_l4a.txt, text replace name(cfxl)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
qui log close cfxl
if _rc != 0 {
  dis as error "L4a FAIL: rc=" _rc
  local ++fails
}
else if strpos(fileread("cfxlat_l4a.txt"), "cfxlattice() variable is missing") != 0 {
  dis as error "L4a FAIL: caution fired for a never-treated-only missing"
  local ++fails
}
mkcfx
qui replace clat = . if ttre < . & t >= ttre & _n <= 400
qui count if clat >= .
local nmiss = r(N)
qui log using cfxlat_l4b.txt, text replace name(cfxl)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
qui log close cfxl
if `nmiss' == 0 {
  dis as error "L4b FAIL: fixture built no missing treated rows"
  local ++fails
}
else if strpos(fileread("cfxlat_l4b.txt"), "observation(s) with a nonmissing treatment onset") == 0 {
  dis as error "L4b FAIL: onset-scoped caution missing for treated-row missings"
  local ++fails
}

* ---- L5: safety(off) still bypasses the value screen ----
mkcfx
qui replace clat = 0.5 if ttre >= . & _n == _N
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat) safety(off)
if _rc != 0 {
  dis as error "L5 FAIL: safety(off) leg rc=" _rc
  local ++fails
}

if `fails' == 0 dis as result "CFXLAT_ALL_PASS"
else dis as error "CFXLAT_FAILS=`fails'"
