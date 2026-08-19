* Targeted suite for the esf support-first reorder (VAL-G10B): under esfixedbaseperiod
* the calendar-support drop now runs BEFORE the per-group reference checks, so a group
* whose only reference row sits in an unsupported period is dropped gracefully instead
* of passing the reference check and dying opaquely downstream (r322 default safety;
* silently wrong on some shapes under safety(off)).
* Run from test/: "<stata>" -e do test_esfsupport.do -> ESFSUP_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 31415
set sortseed 31415
do "../wooldidinstall/wooldid.ado"
local fails 0

* fixture A: 6 control units observed t=1..6 EXCEPT t==3; treated cohorts onset 4
* (reference t=3 unsupported) and onset 5 (reference t=4 supported)
cap program drop mkfixa
program define mkfixa
  clear
  set seed 31415
  set sortseed 31415
  qui set obs 14
  qui gen i = _n
  qui gen ttre = cond(_n <= 4, 4, cond(_n <= 8, 5, .))
  qui expand 6
  qui bysort i: gen t = _n
  qui gen y = 1.5*(t >= ttre & ttre < .) + 0.07*i + rnormal()
  qui drop if ttre >= . & t == 3
end

* ---- E1: bug cure - fixture A completes and matches the hand-deleted-cohort benchmark ----
tempname BA BB
mkfixa
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "E1 FAIL: fixture A rc=" _rc " (expected 0 post-reorder)"
  local ++fails
}
else {
  matrix `BA' = e(b)
  mkfixa
  qui drop if ttre == 4
  qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2)
  matrix `BB' = e(b)
  if mreldif(`BA', `BB') != 0 {
    dis as error "E1 FAIL: fixture A != hand-deleted-cohort benchmark, mreldif=" mreldif(`BA', `BB')
    local ++fails
  }
}

* ---- E2: the support Caution announces the dropped period on fixture A ----
mkfixa
qui log using esfsup_e2.txt, text replace name(esfs)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2)
qui log close esfs
if strpos(fileread("esfsup_e2.txt"), "time periods exist in sample without a control unit") == 0 {
  dis as error "E2 FAIL: support Caution missing on fixture A"
  local ++fails
}

* ---- E3: safety(off) leg is silent-correct (matches the same benchmark; the support
* Caution and geometry refusal are safety-gated by design) ----
tempname BC
mkfixa
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) safety(off)
if _rc != 0 {
  dis as error "E3 FAIL: fixture A safety(off) rc=" _rc
  local ++fails
}
else {
  matrix `BC' = e(b)
  if mreldif(`BC', `BB') != 0 {
    dis as error "E3 FAIL: safety(off) leg != benchmark, mreldif=" mreldif(`BC', `BB')
    local ++fails
  }
}

* ---- E4: fixture B - support drop fires benignly (no reference row affected): runs,
* caution fires, results match hand-deleting the unsupported period ----
tempname BD BE
clear
qui set obs 14
qui gen i = _n
qui gen ttre = cond(_n <= 4, 4, cond(_n <= 8, 5, .))
qui expand 6
qui bysort i: gen t = _n
qui gen y = 1.5*(t >= ttre & ttre < .) + 0.07*i + 0.5*sin(i*t)
qui drop if ttre >= . & t == 6
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(1)
if _rc != 0 {
  dis as error "E4 FAIL: fixture B rc=" _rc
  local ++fails
}
else {
  matrix `BD' = e(b)
  clear
  qui set obs 14
  qui gen i = _n
  qui gen ttre = cond(_n <= 4, 4, cond(_n <= 8, 5, .))
  qui expand 6
  qui bysort i: gen t = _n
  qui gen y = 1.5*(t >= ttre & ttre < .) + 0.07*i + 0.5*sin(i*t)
  qui drop if t == 6
  qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(1)
  matrix `BE' = e(b)
  if mreldif(`BD', `BE') != 0 {
    dis as error "E4 FAIL: fixture B != hand-deleted-period twin, mreldif=" mreldif(`BD', `BE')
    local ++fails
  }
}

* ---- E5: anticipation attribution cell - shifted reference lands on a support-dropped
* period; the anticesf Caution must fire with the amended (truthful) wording ----
clear
qui set obs 14
qui gen i = _n
qui gen ttre = cond(_n <= 4, 5, cond(_n <= 8, 6, .))
qui expand 7
qui bysort i: gen t = _n
qui gen y = 1.5*(t >= ttre & ttre < .) + 0.07*i + rnormal()
qui drop if ttre >= . & t == 3
qui log using esfsup_e5.txt, text replace name(esfs)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(1) espostlength(1) anticipation(1)
local rc5 = _rc
qui log close esfs
if strpos(fileread("esfsup_e5.txt"), "a period unusable for them") == 0 {
  dis as error "E5 FAIL: amended anticesf caution not observed (rc=" `rc5' ")"
  local ++fails
}

if `fails' == 0 dis as result "ESFSUP_ALL_PASS"
else dis as error "ESFSUP_FAILS=`fails'"
