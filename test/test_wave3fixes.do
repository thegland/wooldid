* Targeted suite for the final issue-file batch (NFCOMP-02 cwm sign screen; NFES-02
* degenerate-V_es caution + e(es_ndegen_*); NFRPL-03 replay existence guard; NFFAST-02
* stable CFX row keys; NFCOMP-01/04 + NFES-01 + NFRPL-01/02 are docs-only).
* Run from test/: "<stata>" -e do test_wave3fixes.do -> WAVE3FIX_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 4321
set sortseed 4321
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

* ---- W1: negative cwm on a live in-sample row -> r(411), was silent rc 0 ----
use "`DATA'/dta_binary", clear
qui gen cwm = 1
qui replace cwm = -1 in 5
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) customweightmultiplier(cwm)
if _rc != 411 {
  dis as error "W1 FAIL: negative cwm rc=" _rc " (expected 411)"
  local ++fails
}

* ---- W2: negative cwm confined to a deterministically-discarded row (missing y) ->
* rc 0, bit-identical to the clean twin (estimation-aware domain; G15A precedent) ----
tempname B1 B2
use "`DATA'/dta_binary", clear
qui gen cwm = 1
qui replace y = . in 5
qui replace cwm = -1 in 5
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) customweightmultiplier(cwm)
if _rc != 0 {
  dis as error "W2 FAIL: markout-row negative cwm rc=" _rc " (expected 0)"
  local ++fails
}
else {
  matrix `B1' = e(b)
  use "`DATA'/dta_binary", clear
  qui gen cwm = 1
  qui replace y = . in 5
  qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) customweightmultiplier(cwm)
  matrix `B2' = e(b)
  if mreldif(`B1', `B2') != 0 {
    dis as error "W2 FAIL: twin mreldif=" mreldif(`B1', `B2')
    local ++fails
  }
}

* ---- W3: safety(off) bypasses the cwm screen (documented escape hatch) ----
use "`DATA'/dta_binary", clear
qui gen cwm = 1
qui replace cwm = -1 in 5
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) customweightmultiplier(cwm) safety(off)
if _rc != 0 {
  dis as error "W3 FAIL: safety(off) bypass rc=" _rc
  local ++fails
}

* ---- W4: partially degenerate V_es (untreatedZero: pre rows structurally zero) ->
* caution fires naming the export, e(es_ndegen_ame) counts the degenerate rows ----
use "suite/data/continuous_untreated_zero", clear
qui log using w3f_w4.txt, text replace name(w3f)
cap noi wooldid y_cont_zero cohort t ttre, att contreat(dose_common_zero) contreatcontrols(untreatedZero) esprelength(2) espostlength(2) robust safety(off)
local rc4 = _rc
qui log close w3f
if `rc4' != 0 {
  dis as error "W4 FAIL: rc=" `rc4'
  local ++fails
}
else {
  cap confirm matrix e(V_es_ame)
  if _rc != 0 {
    dis as text "W4 NOTE: no V_es_ame on this fixture; degenerate-row leg vacuous - check fixture varnames"
    dis as error "W4 FAIL: expected V_es_ame to post"
    local ++fails
  }
  else {
    if strpos(fileread("w3f_w4.txt"), "structurally degenerate rows") == 0 {
      dis as error "W4 FAIL: degenerate-V_es caution missing"
      local ++fails
    }
    local ndeg = e(es_ndegen_ame)
    if missing(`ndeg') | `ndeg' < 1 {
      dis as error "W4 FAIL: e(es_ndegen_ame)=" `ndeg' " (expected >= 1)"
      local ++fails
    }
  }
}

* ---- W5: clean config posts NO es_ndegen_* and no caution ----
use "`DATA'/dta_binary", clear
qui log using w3f_w5.txt, text replace name(w3f)
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2)
qui log close w3f
if "`e(es_ndegen_att)'" != "" {
  dis as error "W5 FAIL: es_ndegen_att posted on a clean run"
  local ++fails
}
if strpos(fileread("w3f_w5.txt"), "structurally degenerate rows") != 0 {
  dis as error "W5 FAIL: degenerate caution false-fired on a clean run"
  local ++fails
}

* ---- W6: NFRPL-03 replay guard - forged store without e(wdidmainresults) -> loud 301 ----
macro drop WOOLDID_LASTRC
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
local wrepv `e(wdid_replayv)'
tempname zb zv
matrix `zb' = e(b)
matrix `zv' = e(V)
tempvar zsmp
qui gen byte `zsmp' = e(sample)
cap program drop forge_nomain
program define forge_nomain, eclass
  args bb vv repv smp
  ereturn post `bb' `vv', esample(`smp')
  ereturn local cmd wooldid
  ereturn hidden local wdid_replayv `repv'
end
forge_nomain `zb' `zv' "`wrepv'" `zsmp'
qui log using w3f_w6.txt, text replace name(w3f)
cap noi wooldid
local rc6 = _rc
qui log close w3f
if `rc6' != 301 | strpos(fileread("w3f_w6.txt"), "compatible stored results not found") == 0 {
  dis as error "W6 FAIL: forged store rc=" `rc6' " (expected loud 301 with the targeted message)"
  local ++fails
}

* ---- W7: genuine store still replays after the guard ----
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
qui log using w3f_w7.txt, text replace name(w3f)
cap noi wooldid
local rc7 = _rc
qui log close w3f
if `rc7' != 0 | strpos(fileread("w3f_w7.txt"), "Results (Full Estimation Sample)") == 0 {
  dis as error "W7 FAIL: genuine replay rc=" `rc7'
  local ++fails
}

* ---- W8: NFFAST-02 stable CFX row keys - exact lattice_# sequence (incl. a multi-digit
* and a 0 level), null row equations, and run-to-run stripe stability ----
use "`DATA'/dta_cont", clear
qui gen clat = mod(i, 3)
qui su i if ttre < ., meanonly
qui replace clat = 12 if i == r(min)
qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
local rn1 : rowfullnames e(cfxatt_cfxpost)
foreach tok of local rn1 {
  if !regexm("`tok'", "^lattice_[0-9]+$") {
    dis as error "W8 FAIL: unexpected row key '`tok''"
    local ++fails
  }
}
if strpos("`rn1'", "lattice_0") == 0 | strpos("`rn1'", "lattice_12") == 0 {
  dis as error "W8 FAIL: expected lattice_0 and lattice_12 keys, got [`rn1']"
  local ++fails
}
local re1 : roweq e(cfxatt_cfxpost)
foreach tok of local re1 {
  if "`tok'" != "_" {
    dis as error "W8 FAIL: non-null row equation '`tok''"
    local ++fails
  }
}
use "`DATA'/dta_cont", clear
qui gen clat = mod(i, 3)
qui su i if ttre < ., meanonly
qui replace clat = 12 if i == r(min)
qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
local rn2 : rowfullnames e(cfxatt_cfxpost)
if "`rn1'" != "`rn2'" {
  dis as error "W8 FAIL: row keys differ across identical runs: [`rn1'] vs [`rn2']"
  local ++fails
}

if `fails' == 0 dis as result "WAVE3FIX_ALL_PASS"
else dis as error "WAVE3FIX_FAILS=`fails'"
