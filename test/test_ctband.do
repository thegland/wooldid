* Targeted suite for the commensurability tier's timetrends band fix: under timetrends
* the worker drops treated units with onsets at or below the SECOND observed period, so
* the validate tier's lower bound lifts to the second-smallest distinct t there.
* Run from test/:  "<stata>" -e do test_ctband.do   -> expect CTBAND_ALL_PASS
* On the pre-fix build, T1/T4 abort with r(459) instead of running.
clear all
set more off
set maxvar 20000
set graphics off
set seed 777
set sortseed 777
local fails 0
do "../wooldidinstall/wooldid.ado"

* annual panel, one fractional onset in the band (min(t), second t]
cap program drop mkpanel
program define mkpanel
  args spacing bandonset
  clear
  set seed 777
  set sortseed 777
  qui set obs 40
  qui gen i = _n
  qui gen ttre = .
  qui replace ttre = 2004 if _n <= 8
  qui replace ttre = 2005 if _n > 8 & _n <= 16
  qui replace ttre = `bandonset' if _n == 40
  qui expand 8
  qui bysort i: gen t = 2000 + (_n - 1)*`spacing'
  qui gen te = 2*(t >= ttre & ttre < .)
  qui gen y = te + 0.05*i + 0.1*t - 195 + rnormal()
  qui drop te
end

* ---- T1: annual band case + timetrends -> runs, and e(b) equals the offender-deleted run ----
mkpanel 1 2000.5
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
if _rc != 0 {
  dis as error "T1 FAIL: rc=" _rc " (expected the band onset to be dropped by the worker, not an error)"
  local ++fails
}
else {
  tempname b1
  matrix `b1' = e(b)
  mkpanel 1 2000.5
  qui drop if ttre == 2000.5
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
  if _rc != 0 | mreldif(`b1', e(b)) != 0 {
    dis as error "T1 FAIL: offender-deleted run rc=" _rc " mreldif=" mreldif(`b1', e(b))
    local ++fails
  }
}

* ---- T2: same panel WITHOUT timetrends -> the tier still errors (bound unchanged there) ----
mkpanel 1 2000.5
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
if _rc != 459 {
  dis as error "T2 FAIL: rc=" _rc " (expected 459: non-timetrends bound must not move)"
  local ++fails
}

* ---- T3: onset ABOVE the second period + timetrends -> still errors (bound is exact) ----
mkpanel 1 2001.5
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
if _rc != 459 {
  dis as error "T3 FAIL: rc=" _rc " (expected 459: 2001.5 sits above the lifted bound)"
  local ++fails
}

* ---- T4: biennial panel: band spans two calendar years; both offsets pass to the worker ----
* (the once-proposed +1 bound would still error on 2001.5 here; second-smallest t = 2002)
foreach bo in 2000.5 2001.5 {
  mkpanel 2 `bo'
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
  if _rc != 0 {
    dis as error "T4 FAIL: rc=" _rc " at band onset `bo' on the biennial grid"
    local ++fails
  }
}

* ---- T5: clean integer data + timetrends unaffected ----
use "data/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
if _rc != 0 {
  dis as error "T5 FAIL: rc=" _rc
  local ++fails
}

if `fails' == 0 dis "CTBAND_ALL_PASS"
else dis as error "CTBAND_FAILS=`fails'"
exit, clear
