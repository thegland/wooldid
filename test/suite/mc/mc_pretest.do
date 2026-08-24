* mc_pretest.do — null-fidelity pretest (PREREG_NEWFEAT_2026-08-24 §M2): before any
* planted-line MC fires, verify on a deterministic large-N tiny-noise panel that the
* REPORTED ES coefficients lie on the planted line (guards against endcap pooling
* silently breaking the composite null). Also reports (info) the pooling behavior of
* the 12-period A1 window shape.
* args: 1=adopath.  Expect NULLFID_ALL_PASS at column 1.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
global pf 0

* ---- the planted-line shape used by A2/A3b/ptp arms: G=40, t=1..9, onset 6 ----
clear
set seed 771001
qui set obs 40
gen i = _n
gen ttre = cond(mod(_n, 2) == 0, 6, .)
qui expand 25
bysort i: gen unit = _n
qui expand 9
bysort i unit: gen t = _n
gen y = .05*t + .15*(t - 6 + 1)*(ttre < .) + .0001*rnormal()
sort i unit t
cap noi wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "NULLFID FAIL: 9-period run rc=" _rc
  global pf = $pf + 1
}
else {
  tempname B R
  matrix `B' = e(b_es_att)
  local lbl : colnames `B'
  local want "pre5 pre4 pre3 pre2 post0 post1 post2 post3"
  if "`lbl'" != "`want'" {
    dis as error "NULLFID FAIL: 9-period labels [`lbl'] != [`want']"
    global pf = $pf + 1
  }
  local rys "-5 -4 -3 -2 0 1 2 3"
  forval j = 1/8 {
    local ry : word `j' of `rys'
    local truth = .15*(`ry' + 1)
    local got = `B'[1, `j']
    if abs(`got' - `truth') > .002 {
      dis as error "NULLFID FAIL: horizon `: word `j' of `lbl'' coef " %12.6f `got' " vs planted line " %12.6f `truth'
      global pf = $pf + 1
    }
  }
  dis "NULLFID info: 9-period shape max |coef - line| check complete"
}

* ---- info only: the 12-period shape (A1/A3a/A5/A6 use it with a FLAT null, where
* pooling is harmless; this documents what post3 does with relyr 3..6 present) ----
clear
set seed 771002
qui set obs 40
gen i = _n
gen ttre = cond(mod(_n, 2) == 0, 6, .)
qui expand 25
bysort i: gen unit = _n
qui expand 12
bysort i unit: gen t = _n
gen y = .05*t + .15*(t - 6 + 1)*(ttre < .) + .0001*rnormal()
sort i unit t
cap noi wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
if _rc != 0 dis "NULLFID info: 12-period run rc=" _rc
else {
  tempname B2
  matrix `B2' = e(b_es_att)
  dis "NULLFID info: 12-period labels: `: colnames `B2''"
  dis "NULLFID info: post3 coef = " %12.6f `B2'[1, colsof(`B2')] " (pure relyr-3 line value would be .60; a pooled 3..6 average would exceed it)"
}

if $pf == 0 dis "NULLFID_ALL_PASS"
else dis as error "NULLFID_FAILS=$pf"
exit, clear
