* Targeted suite for the ttre/t double-storage fix: large calendar codes (YYYYMMDD
* class), half-integers, and sparse grids must give results identical to the same
* panel with small recoded time values. On the pre-fix build, F1 silently MERGES
* cohorts whose onsets float-round together and F4 aborts r(461).
* Run from test/:  "<stata>" -e do test_ttrefloat.do   -> expect TTREFLOAT_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
local fails 0
do "../wooldidinstall/wooldid.ado"

* panel maker: 30 units, 10 periods on a parameterized time grid; two close odd onsets
cap program drop mkbig
program define mkbig
  args base step o1 o2
  clear
  set seed 2468
  set sortseed 2468
  qui set obs 30
  qui gen i = _n
  * the user's own data must hold large codes in double - that is a data-authoring
  * requirement, not wooldid's; these tests target wooldid's INTERNAL copies
  qui gen double ttre = cond(_n <= 8, `o1', cond(_n <= 16, `o2', .))
  qui expand 10
  qui bysort i: gen double t = `base' + (_n - 1)*`step'
  qui gen te = 2*(t >= ttre & ttre < .)
  qui gen y = te + 0.05*i + rnormal()
  qui drop te
  qui gen double trecode = cond(ttre < ., (ttre - `base')/`step' + 1, .)
  qui gen double t2 = (t - `base')/`step' + 1
end

* helper: run config on raw grid and on the recode, assert identical e(b)
cap program drop cmprun
program define cmprun
  args label opts
  tempname braw
  cap noi wooldid y i t ttre, `opts'
  if _rc != 0 {
    dis as error "`label' FAIL: raw-grid rc=" _rc
    global s6f = $s6f + 1
    exit
  }
  matrix `braw' = e(b)
  cap noi wooldid y i t2 trecode, `opts'
  if _rc != 0 | mreldif(`braw', e(b)) != 0 {
    dis as error "`label' FAIL: recode rc=" _rc " mreldif=" mreldif(`braw', e(b))
    global s6f = $s6f + 1
  }
end
global s6f 0

* ---- F1: dense daily-style grid, two odd onsets one code apart (float-merges pre-fix) ----
mkbig 20200101 1 20200107 20200109
cmprun F1 "att cluster(i) esprelength(2) espostlength(2)"

* ---- F2: odd integers below 2^24 (float-exact class; provable no-harm pin) ----
mkbig 16000001 1 16000004 16000006
cmprun F2 "att cluster(i) esprelength(2) espostlength(2)"

* ---- F3: half-integer codes in [2^23, 2^24) (silently shifted pre-fix) ----
mkbig 8500000.5 1 8500003.5 8500005.5
cmprun F3 "att cluster(i) esprelength(2) espostlength(2)"

* ---- F4: sparse step-100 large grid (r(461) pre-fix); no ES window, since window
* lengths are in t-units and would differ between the raw grid and the recode ----
mkbig 20200100 100 20200400 20200600
cmprun F4 "att cluster(i)"

* ---- F5: esfixedbaseperiod on the dense large grid ----
mkbig 20200101 1 20200107 20200109
cmprun F5 "att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(2) espostlength(2)"

* ---- F6: anticipation(1) on the dense large grid (shifted onsets stay exact) ----
mkbig 20200101 1 20200107 20200109
cmprun F6 "att cluster(i) esprelength(2) espostlength(2) anticipation(1)"

* ---- F7: gteffects carries the exact large onset/caltime codes ----
mkbig 20200101 1 20200107 20200109
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "F7 FAIL: rc=" _rc
  global s6f = $s6f + 1
}
else {
  tempname G
  matrix `G' = e(gteffects)
  local o7 0
  local o9 0
  local calbad 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',2] == 20200107 local ++o7
    if `G'[`r',2] == 20200109 local ++o9
    if `G'[`r',3] != `G'[`r',2] + `G'[`r',4] local ++calbad
  }
  if `o7' == 0 | `o9' == 0 | `calbad' != 0 {
    dis as error "F7 FAIL: onset codes not exact (o7=`o7' o9=`o9' calbad=`calbad')"
    global s6f = $s6f + 1
  }
}

* ---- F8: timetrends on the dense large grid == its recode (trends no longer inert) ----
mkbig 20200101 1 20200107 20200109
cmprun F8 "att cluster(i) esprelength(2) espostlength(2) timetrends"

* ---- F9: trends do real work at large t (tt vs no-tt must differ materially) ----
mkbig 20200101 1 20200107 20200109
qui replace y = y + 0.3*(t - 20200101)*(ttre < .)
tempname bnt
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "F9 FAIL: no-tt rc=" _rc
  global s6f = $s6f + 1
}
else {
  matrix `bnt' = e(b)
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
  if _rc != 0 {
    dis as error "F9 FAIL: tt rc=" _rc
    global s6f = $s6f + 1
  }
  else if reldif(`bnt'[1,1], e(b)[1,1]) < 1e-4 {
    dis as error "F9 FAIL: timetrends numerically inert at large t (reldif=" reldif(`bnt'[1,1], e(b)[1,1]) ")"
    global s6f = $s6f + 1
  }
}

* ---- F5 (bug sweep BS-34): float-stored large codes draw a caution ----
* a 20-period panel whose t is an 8-digit calendar code. Stored as float, consecutive
* codes above 2^24 round together and adjacent periods merge with no other symptom;
* wooldid cannot recover the values but must say so. Expect: the caution on the float
* run (rc 0), none on the double twin (att = 1.5, the noise-free DGP effect), none on a
* small float t, and none under safety(off)
cap program drop scanlog34
program scanlog34, rclass
  args fname needle
  tempname fh
  local found 0
  file open `fh' using "`fname'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "`needle'") > 0 local found 1
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
end
cap program drop mk34
program define mk34
  args ttype base
  clear
  qui set obs 12
  qui gen i = _n
  qui gen double onset = cond(_n <= 6, 8, .)
  qui expand 20
  qui bysort i: gen double k = _n
  qui gen `ttype' t = `base' + k
  qui gen `ttype' ttre = cond(onset < ., `base' + onset, .)
  qui gen double y = 1.5*(k >= onset & onset < .) + 0.3*i + 0.2*k
end
set linesize 200
local f5needle stored as float and holds values above
local f5 0
mk34 float 20150000
qui log using ttf_f5a.txt, text replace name(ttf)
cap noi wooldid y i t ttre, att
local rc5a = _rc
qui log close ttf
scanlog34 ttf_f5a.txt "`f5needle'"
local c5a = r(found)
local att5a = cond(`rc5a' == 0, e(wdidmainresults)[1,1], .)
mk34 double 20150000
qui log using ttf_f5b.txt, text replace name(ttf)
cap noi wooldid y i t ttre, att
local rc5b = _rc
qui log close ttf
scanlog34 ttf_f5b.txt "`f5needle'"
local c5b = r(found)
local att5b = cond(`rc5b' == 0, e(wdidmainresults)[1,1], .)
mk34 float 0
qui log using ttf_f5c.txt, text replace name(ttf)
cap noi wooldid y i t ttre, att
local rc5c = _rc
qui log close ttf
scanlog34 ttf_f5c.txt "`f5needle'"
local c5c = r(found)
mk34 float 20150000
qui log using ttf_f5d.txt, text replace name(ttf)
cap noi wooldid y i t ttre, att safety(off)
local rc5d = _rc
qui log close ttf
scanlog34 ttf_f5d.txt "`f5needle'"
local c5d = r(found)
cap erase ttf_f5a.txt
cap erase ttf_f5b.txt
cap erase ttf_f5c.txt
cap erase ttf_f5d.txt
dis "F5 float big: rc=`rc5a' caution=`c5a' att=`att5a' | double: rc=`rc5b' caution=`c5b' att=`att5b' | float small: rc=`rc5c' caution=`c5c' | float big safety(off): rc=`rc5d' caution=`c5d'"
if `rc5a' == 0 & `c5a' == 1 & `rc5b' == 0 & `c5b' == 0 & abs(`att5b' - 1.5) < 1e-6 & `rc5c' == 0 & `c5c' == 0 & `rc5d' == 0 & `c5d' == 0 {
  dis as result "F5 PASS: the float-magnitude caution fires only on the large float codes at default safety"
}
else {
  dis as error "F5 FAIL: float-magnitude caution"
  global s6f = $s6f + 1
}

* ---- F6 (ruling D3, 2026-09-23): a treated unit whose i code exceeds 2,147,483,647, the
* factor-variable cap, is refused in validate with wooldid's own error (r(411)) at any safety(),
* where it used to die r(452) after the checks naming a tempvar; a never-treated unit with such
* a code still runs, and its att equals the run on the original codes
use "data/dta_binary", clear
qui sum i if ttre < .
local f6tr = r(max)
qui sum i if ttre >= .
local f6nt = r(max)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
local f6att0 = cond(_rc == 0, e(wdidmainresults)[1,1], .)
qui gen double ib = i
qui replace ib = 3000000000 + i if i == `f6tr'
qui log using ttf_f6a.txt, text replace name(ttf)
cap noi wooldid y ib t ttre, att cluster(ib) esprelength(1) espostlength(1)
local rc6a = _rc
cap noi wooldid y ib t ttre, att cluster(ib) esprelength(1) espostlength(1) safety(off)
local rc6b = _rc
qui log close ttf
local c6a = strpos(fileread("ttf_f6a.txt"), "Error: ib must not exceed 2,147,483,647 on treated units") > 0
cap erase ttf_f6a.txt
qui replace ib = i
qui replace ib = 3000000000 + i if i == `f6nt'
cap noi wooldid y ib t ttre, att cluster(ib) esprelength(1) espostlength(1)
local rc6c = _rc
local f6att1 = cond(`rc6c' == 0, e(wdidmainresults)[1,1], .)
dis "F6 i cap: treated big rc=`rc6a' (safety(off) rc=`rc6b') message=`c6a' | never-treated big rc=`rc6c' att=`f6att1' vs original `f6att0'"
if `rc6a' == 411 & `rc6b' == 411 & `c6a' == 1 & `rc6c' == 0 & reldif(`f6att1', `f6att0') < 1e-10 {
  dis as result "F6 PASS: a treated i past the factor-variable cap is refused by name; a never-treated one runs unchanged"
}
else {
  dis as error "F6 FAIL: treated-i factor-variable cap"
  global s6f = $s6f + 1
}

if $s6f == 0 dis "TTREFLOAT_ALL_PASS"
else dis as error "TTREFLOAT_FAILS=$s6f"
global s6f
exit, clear
