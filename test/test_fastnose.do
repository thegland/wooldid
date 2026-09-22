* Targeted suite for fastnose (F9): point estimates only, "." cells, no tests.
* Run from test/:  "<stata>" -e do test_fastnose.do   -> expect FASTNOSE_ALL_PASS
* Fast-vs-full point estimates agree to ~1e-12 (margins moves in the last ulp when
* it skips the VCE; author-accepted); shapes/absences per the documented contract.
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* ---- T1: fast-vs-full e(b) equality across representative configs ----
local cfgs `""att cluster(i) esprelength(2) espostlength(2)""'
local cfgs `"`cfgs' "att att_it att_i att_itime cluster(i)""'
local cfgs `"`cfgs' "att cluster(i) esprelength(2) espostlength(2) subgroups(sg)""'
local cfgs `"`cfgs' "att robust esprelength(2) espostlength(2)""'
local cfgs `"`cfgs' "att customvce(robust) esprelength(2) espostlength(2)""'
local cfgs `"`cfgs' "att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(2)""'
local n 0
foreach c of local cfgs {
  local ++n
  use "`DATA'/dta_binary", clear
  cap noi wooldid y i t ttre [aw=w], `c'
  if _rc != 0 {
    dis as error "T1.`n' FAIL: full run rc=" _rc
    local ++fails
    continue
  }
  tempname bfull
  matrix `bfull' = e(b)
  use "`DATA'/dta_binary", clear
  cap noi wooldid y i t ttre [aw=w], `c' fastnose
  if _rc != 0 {
    dis as error "T1.`n' FAIL: fast run rc=" _rc
    local ++fails
    continue
  }
  if mreldif(`bfull', e(b)) > 1e-10 {
    dis as error "T1.`n' FAIL: mreldif(e(b)) = " mreldif(`bfull', e(b))
    local ++fails
  }
  if "`e(properties)'" != "b" {
    dis as error "T1.`n' FAIL: e(properties) = `e(properties)'"
    local ++fails
  }
}
* continuous twin
use "`DATA'/dta_cont", clear
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2)
tempname bfullc
matrix `bfullc' = e(b)
use "`DATA'/dta_cont", clear
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) fastnose
if _rc != 0 | mreldif(`bfullc', e(b)) > 1e-10 {
  dis as error "T1.c FAIL: rc=" _rc " mreldif=" mreldif(`bfullc', e(b))
  local ++fails
}

* ---- T2: shape/absence contract on a subgroups fast run ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg) summarystats fastnose
if _rc != 0 {
  dis as error "T2 FAIL: rc=" _rc
  local ++fails
}
else {
  cap confirm matrix e(V)
  if _rc == 0 {
    dis as error "T2 FAIL: e(V) posted under fastnose"
    local ++fails
  }
  * 2026-09-13 author ruling: fast mode posts the event-study COEFFICIENT vector and
  * e(numpre_*) - sthlp:1295 always said it did - and still no covariance. Both halves
  * are pinned here, on the main family and on a subgroup twin
  cap confirm matrix e(b_es_att)
  if _rc != 0 {
    dis as error "T2 FAIL: e(b_es_att) not posted under fastnose"
    local ++fails
  }
  cap confirm matrix e(b_es_att_s1)
  if _rc != 0 {
    dis as error "T2 FAIL: e(b_es_att_s1) not posted under fastnose"
    local ++fails
  }
  if e(numpre_att) != 2 {
    dis as error "T2 FAIL: e(numpre_att) = " e(numpre_att) " under fastnose, expected 2"
    local ++fails
  }
  foreach t2v in V_es_att V_es_att_it V_es_att_s1 V_es_att_it_s1 {
    cap confirm matrix e(`t2v')
    if _rc == 0 {
      dis as error "T2 FAIL: e(`t2v') posted under fastnose"
      local ++fails
    }
  }
  if "`e(joint_att_allsgssame)'" != "" | "`e(joint_horiz)'" != "" {
    dis as error "T2 FAIL: joint-test results posted under fastnose"
    local ++fails
  }
  if "`e(wdid_fastnose)'" != "fastnose" {
    dis as error "T2 FAIL: e(wdid_fastnose) missing"
    local ++fails
  }
  tempname W M
  matrix `W' = e(pwsgs)
  matrix `M' = e(wdidmainresults)
  if colsof(`W') != 6 | !missing(`W'[1,2]) | missing(`W'[1,1]) {
    dis as error "T2 FAIL: e(pwsgs) shape/cells wrong under fastnose"
    local ++fails
  }
  if !missing(`M'[1,2]) | missing(`M'[1,1]) {
    dis as error "T2 FAIL: e(wdidmainresults) SE cell not missing / estimate missing"
    local ++fails
  }
  * escomp and summarystats are margins-free and must match the full run bit-for-bit
  tempname E1m
  matrix `E1m' = e(escomp)
  use "`DATA'/dta_binary", clear
  qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg) summarystats
  if mreldif(`E1m', e(escomp)) != 0 {
    dis as error "T2 FAIL: e(escomp) differs between fast and full"
    local ++fails
  }
}

* ---- T3: cfx matrices keep the cfxplotnose 1-column shape ----
use "`DATA'/dta_cont", clear
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) makecfxplots fastnose
if _rc != 0 {
  dis as error "T3 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname C
  cap matrix `C' = e(cfxatt_cfxpost)
  if _rc != 0 | colsof(`C') != 1 {
    dis as error "T3 FAIL: e(cfxatt_cfxpost) not 1-column under fastnose"
    local ++fails
  }
}

* ---- T4: the four remaining incompatibility gates ----
* impvalues left this list on 2026-09-21: the author ruled that fastnose alongside impvalues
* cautions and is dropped rather than refused, because impvalues already runs every margins
* call with nose. T4b below pins the new behaviour.
foreach g in "unconditionalse robust" "adversarial" "jointtests" "pretrendpower(80) esfixedbaseperiod esprelength(2)" {
  use "`DATA'/dta_binary", clear
  cap wooldid y i t ttre, att cluster(i) fastnose `g'
  if _rc != 198 {
    dis as error "T4 FAIL: fastnose + `g' gave rc=" _rc " (expected 198)"
    local ++fails
  }
}

* ---- T4b: fastnose + impvalues runs, cautions, and drops fast mode ----
use "`DATA'/dta_binary", clear
qui log using fn_t4b.txt, text replace name(fn4b)
cap noi wooldid y i t ttre, att cluster(i) fastnose impvalues
local rc4b = _rc
qui log close fn4b
if `rc4b' != 0 {
  dis as error "T4b FAIL: fastnose + impvalues gave rc=" `rc4b' " (expected 0)"
  local ++fails
}
else if strpos(fileread("fn_t4b.txt"), "fastnose is unnecessary with impvalues") == 0 {
  dis as error "T4b FAIL: the redundant-fastnose caution did not print"
  local ++fails
}
else if "`e(wdid_fastnose)'" != "" {
  dis as error "T4b FAIL: fastnose was carried into the run rather than dropped"
  local ++fails
}

* ---- T5: replay of a fast run is byte-identical from the banner onward ----
use "`DATA'/dta_binary", clear
qui log using fn_est.txt, text replace name(fn)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) fastnose
qui log close fn
local rc1 = _rc
qui log using fn_rep.txt, text replace name(fn)
cap noi wooldid
qui log close fn
if `rc1' != 0 | _rc != 0 {
  dis as error "T5 FAIL: rc est=" `rc1' " rep=" _rc
  local ++fails
}
* the actual replay comparison: the two logs must be LINE-IDENTICAL from the
* "Wooldid Estimation" banner to EOF. The only benign difference is the pre-banner
* region, which echoes the estimation command line in fn_est.txt but the bare
* `wooldid` replay call in fn_rep.txt (verified byte-identical from the banner onward
* on the certified build; both logs are written by the same session, so the wrap
* width is identical too).
local t5same 0
local t5n 0
cap confirm file fn_est.txt
local frc1 = _rc
cap confirm file fn_rep.txt
if `frc1' != 0 | _rc != 0 {
  dis as error "T5 FAIL: replay log files missing (rc `frc1'/" _rc ")"
  local ++fails
}
else {
  mata: a_ = cat("fn_est.txt") \ ""; b_ = cat("fn_rep.txt") \ ""
  mata: ia_ = selectindex(strpos(a_, "Wooldid Estimation") :== 1); ib_ = selectindex(strpos(b_, "Wooldid Estimation") :== 1)
  mata: hb_ = (rows(ia_) > 0) & (rows(ib_) > 0)
  mata: sa_ = a_[|(ia_ \ 1)[1] \ rows(a_)|]; sb_ = b_[|(ib_ \ 1)[1] \ rows(b_)|]
  mata: n_ = rows(sa_)
  mata: st_local("t5same", strofreal(hb_ & (n_ > 20) & (sa_ == sb_))); st_local("t5n", strofreal(n_))
  if `t5same' != 1 {
    dis as error "T5 FAIL: replay log differs from the estimation log from the banner onward (banner-to-EOF lines in fn_est: `t5n')"
    local ++fails
  }
  else dis "T5 info: replay byte-identical from the banner onward (`t5n' lines)"
}

* ---- T6: makeplots under fastnose - plots draw without bands, no rarea failure ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) makeplots saveplots(fnplot) replace fastnose
if _rc != 0 {
  dis as error "T6 FAIL: rc=" _rc
  local ++fails
}
else {
  qui use fnplot_ES.dta, clear
  qui count if !missing(lb_95ci) & rownamesvar != "reference"
  if r(N) != 0 {
    dis as error "T6 FAIL: non-reference CI bounds present under fastnose"
    local ++fails
  }
  qui count if missing(estimate)
  if r(N) != 0 {
    dis as error "T6 FAIL: missing point estimates in the saved ES data"
    local ++fails
  }
}

if `fails' == 0 dis "FASTNOSE_ALL_PASS"
else dis as error "FASTNOSE_FAILS=`fails'"
exit, clear
