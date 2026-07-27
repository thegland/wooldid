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

if $s6f == 0 dis "TTREFLOAT_ALL_PASS"
else dis as error "TTREFLOAT_FAILS=$s6f"
global s6f
exit, clear
