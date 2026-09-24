* Targeted suite for F3 gteffects: the cohort-by-relative-period building-block table.
* Run from test/:  "<stata>" -e do test_gteffects.do   -> expect GTEFFECTS_ALL_PASS
* Core invariant: the weighted mean of the estimate column under each estimand's weight
* column reproduces that aggregate estimate (the cells ARE the building blocks).
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* ---- G1: basic binary run - shape, names, aggregation identity, sane SEs ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G1 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  if colsof(`G') != 15 {
    dis as error "G1 FAIL: colsof=" colsof(`G')
    local ++fails
  }
  local cn : colnames `G'
  if "`cn'" != "i onset caltime relyr estimate se p lb_95ci ub_95ci inwindow nobs w_att w_att_it w_att_i w_att_itime" {
    dis as error "G1 FAIL: colnames=[`cn']"
    local ++fails
  }
  local num 0
  local den 0
  local badse 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',4] >= 0 {
      local num = `num' + `G'[`r',12]*`G'[`r',5]
      local den = `den' + `G'[`r',12]
    }
    if `G'[`r',6] != . & `G'[`r',6] <= 0 local ++badse
  }
  if reldif(`num'/`den', _b[main_att]) > 1e-10 {
    dis as error "G1 FAIL: att aggregation identity reldif=" reldif(`num'/`den', _b[main_att])
    local ++fails
  }
  if `badse' != 0 {
    dis as error "G1 FAIL: `badse' nonpositive SEs"
    local ++fails
  }
}

* ---- G2: aw weights - att vs att_it differ, both identities hold, horizon cross-check ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G2 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G BE
  matrix `G' = e(gteffects)
  matrix `BE' = e(b_es_att)
  local num 0
  local den 0
  local numit 0
  local denit 0
  local h1 0
  local hw 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',4] >= 0 {
      local num = `num' + `G'[`r',12]*`G'[`r',5]
      local den = `den' + `G'[`r',12]
      local numit = `numit' + `G'[`r',13]*`G'[`r',5]
      local denit = `denit' + `G'[`r',13]
    }
    if `G'[`r',4] == 1 & `G'[`r',10] == 1 {
      local h1 = `h1' + `G'[`r',12]*`G'[`r',5]
      local hw = `hw' + `G'[`r',12]
    }
  }
  if reldif(`num'/`den', _b[main_att]) > 1e-10 | reldif(`numit'/`denit', _b[main_att_it]) > 1e-10 {
    dis as error "G2 FAIL: identities att=" reldif(`num'/`den', _b[main_att]) " att_it=" reldif(`numit'/`denit', _b[main_att_it])
    local ++fails
  }
  if reldif(_b[main_att], _b[main_att_it]) <= 1e-12 {
    dis as error "G2 FAIL: att and att_it coincide under aw (test has no power)"
    local ++fails
  }
  local c1 = colnumb(`BE', "post1")
  if reldif(`h1'/`hw', `BE'[1,`c1']) > 1e-10 {
    dis as error "G2 FAIL: post1 horizon cross-check reldif=" reldif(`h1'/`hw', `BE'[1,`c1'])
    local ++fails
  }
}

* ---- G3: esfixedbaseperiod - relyr mapping, caltime identity, window flag ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G3 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G BE
  matrix `G' = e(gteffects)
  matrix `BE' = e(b_es_att)
  local h0 0
  local hw0 0
  local calbad 0
  local inw 0
  local outw 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',4] == 0 & `G'[`r',10] == 1 {
      local h0 = `h0' + `G'[`r',12]*`G'[`r',5]
      local hw0 = `hw0' + `G'[`r',12]
    }
    if `G'[`r',3] != `G'[`r',2] + `G'[`r',4] local ++calbad
    if `G'[`r',10] == 1 local ++inw
    else local ++outw
  }
  local c0 = colnumb(`BE', "post0")
  if reldif(`h0'/`hw0', `BE'[1,`c0']) > 1e-10 {
    dis as error "G3 FAIL: post0 horizon cross-check reldif=" reldif(`h0'/`hw0', `BE'[1,`c0'])
    local ++fails
  }
  if `calbad' != 0 {
    dis as error "G3 FAIL: caltime != onset + relyr on `calbad' rows"
    local ++fails
  }
  if `inw' == 0 | `outw' == 0 {
    dis as error "G3 FAIL: inwindow split degenerate (inw=`inw' outw=`outw')"
    local ++fails
  }
}

* ---- G4: contreat() is refused at parse ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) gteffects
if _rc != 198 {
  dis as error "G4 FAIL: rc=" _rc " (expected 198)"
  local ++fails
}

* ---- G5: fastnose synergy - full SEs in the table while e(V) is absent ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fastnose gteffects
if _rc != 0 {
  dis as error "G5 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  local misse 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',6] == . local ++misse
  }
  cap confirm matrix e(V)
  if `misse' != 0 | _rc == 0 {
    dis as error "G5 FAIL: missing SEs=`misse', e(V) present rc=" _rc
    local ++fails
  }
}

* ---- G6: timetrends - degenerate variances report as missing, never zero ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends gteffects
if _rc != 0 {
  dis as error "G6 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  local zse 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',6] == 0 local ++zse
  }
  if `zse' != 0 {
    dis as error "G6 FAIL: `zse' exact-zero SEs shipped"
    local ++fails
  }
}

* ---- G7: poisson - runs; the xb-scale note under the table was removed by the author's 2026-09-21 message review ----
use "`DATA'/dta_binary", clear
qui log using gt_t7.txt, text replace name(gt)
cap noi wooldid ycount i t ttre, att cluster(i) poisson gteffects
local rc7 = _rc
qui log close gt
if `rc7' != 0 | strpos(fileread("gt_t7.txt"), "linear-index (xb) scale") != 0 {
  dis as error "G7 FAIL: rc=" `rc7' " xb note printed but was removed"
  local ++fails
}

* ---- G8: subgroups - one full-sample table ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) gteffects
cap confirm matrix e(gteffects)
if _rc != 0 {
  dis as error "G8 FAIL: no table under subgroups"
  local ++fails
}

* ---- G9: 9-digit administrative ids round-trip exactly ----
use "`DATA'/dta_binary", clear
qui gen double bigi = 100000000 + i*7919
cap noi wooldid y bigi t ttre, att cluster(bigi) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G9 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  qui glevelsof bigi if ttre != ., clean local(bl)
  local bad 0
  forval r = 1/`=rowsof(`G')' {
    local v = `G'[`r',1]
    if strpos(" `bl' ", " `=string(`v', "%12.0f")' ") == 0 local ++bad
  }
  if `bad' != 0 {
    dis as error "G9 FAIL: `bad' id values do not round-trip"
    local ++fails
  }
}

* ---- G10: replay reprints the table; G11: off-state leaves e() clean ----
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) gteffects
qui log using gt_t10.txt, text replace name(gt)
cap noi wooldid
qui log close gt
if _rc != 0 | strpos(fileread("gt_t10.txt"), "Cohort-by-Relative-Period Effects") == 0 {
  dis as error "G10 FAIL: replay rc=" _rc " or table not reprinted"
  local ++fails
}
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
cap confirm matrix e(gteffects)
if _rc == 0 {
  dis as error "G11 FAIL: e(gteffects) posted without the option"
  local ++fails
}

* ---- G12 (change audit 2026-09-22, gt-containment-customvce-arm): the containment count and
* its caution under a heteroskedasticity-robust customvce(). every treated cell of
* dta_path_onepercell is one observation, so robust counts all of them; customvce(robust) is the
* same meat and must count the same cells (it used to post no e(gt_ncontained) and print nothing) ----
set linesize 255
use "`DATA'/dta_path_onepercell", clear
cap noi wooldid y i t ttre, att robust esprelength(1) espostlength(1) gteffects
local rc12a = _rc
local nc12a = e(gt_ncontained)
local nr12a = cond(`rc12a' == 0, rowsof(e(gteffects)), .)
use "`DATA'/dta_path_onepercell", clear
qui log using gt_g12.txt, text replace name(gt12)
cap noi wooldid y i t ttre, att customvce(robust) esprelength(1) espostlength(1) gteffects
local rc12b = _rc
qui log close gt12
local nc12b = e(gt_ncontained)
local cau12b = strpos(fileread("gt_g12.txt"), "e(gteffects) cell(s) lie wholly inside a single observation") > 0
cap erase gt_g12.txt
dis "G12: robust rc=`rc12a' ncontained=`nc12a' rows=`nr12a' | customvce(robust) rc=`rc12b' ncontained=`nc12b' caution=`cau12b'"
if `rc12a' != 0 | `rc12b' != 0 | `nc12a' >= . | `nc12a' < 1 | `nc12b' != `nc12a' | `cau12b' != 1 {
  dis as error "G12 FAIL: customvce(robust) must count the one-observation cells as robust does"
  local ++fails
}

if `fails' == 0 dis "GTEFFECTS_ALL_PASS"
else dis as error "GTEFFECTS_FAILS=`fails'"
exit, clear
