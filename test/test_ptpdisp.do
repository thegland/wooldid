* Targeted suite for the PTP/display medium batch (NFPTP-03/04/06 + DOC-03 + NFADV-07):
* - NFPTP-03: PTP bias block names its estimand; surrogate note for att_i/att_itime
* - NFPTP-04: below-size NIS power levels report mdt_nis/bias_nis = 0 (boundary answer)
* - DOC-03: escomp coverage scalars displayed under summarystats
* - NFADV-07: altmethod_adversarial without adversarial -> r(198), was a silent no-op
* - T7: the pre-test note names the normal reference, live and on replay
* Run from test/: "<stata>" -e do test_ptpdisp.do -> PTPDISP_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 8642
set sortseed 8642
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

* ---- T1: att_it-only request - no surrogate note; block names att_it ----
use "`DATA'/dta_binary", clear
qui log using ptpd_t1.txt, text replace name(ptpd)
cap noi wooldid y i t ttre, att_it cluster(i) esfixedbaseperiod esprelength(3) espostlength(3) pretrendpower(80)
local rc1 = _rc
qui log close ptpd
if `rc1' != 0 {
  dis as error "T1 FAIL: rc=" `rc1'
  local ++fails
}
else {
  if strpos(fileread("ptpd_t1.txt"), "att_it pooled estimand if it continued") == 0 {
    dis as error "T1 FAIL: estimand-named bias line missing"
    local ++fails
  }
  if strpos(fileread("ptpd_t1.txt"), "stands in for the requested") != 0 {
    dis as error "T1 FAIL: surrogate note fired on an att_it-only request"
    local ++fails
  }
}

* ---- T2: att_i-only request - surrogate note present ----
use "`DATA'/dta_binary", clear
qui log using ptpd_t2.txt, text replace name(ptpd)
cap noi wooldid y i t ttre, att_i cluster(i) esfixedbaseperiod esprelength(3) espostlength(3) pretrendpower(80)
local rc2 = _rc
qui log close ptpd
if `rc2' != 0 {
  dis as error "T2 FAIL: rc=" `rc2'
  local ++fails
}
else if strpos(fileread("ptpd_t2.txt"), "stands in for the requested att_i estimand(s)") == 0 {
  dis as error "T2 FAIL: surrogate note missing for att_i-only request"
  local ++fails
}

* ---- T3: live vs replay - roth block byte-parity (incl. the new note lines) ----
use "`DATA'/dta_binary", clear
qui log using ptpd_t3a.txt, text replace name(ptpd)
cap noi wooldid y i t ttre, att_i cluster(i) esfixedbaseperiod esprelength(3) espostlength(3) pretrendpower(80)
qui log close ptpd
qui log using ptpd_t3b.txt, text replace name(ptpd)
cap noi wooldid
qui log close ptpd
preserve
clear
local bad 0
foreach ff in ptpd_t3a ptpd_t3b {
  qui infix str300 line 1-300 using "`ff'.txt", clear
  qui keep if strpos(line, "Pre-Trend Power") > 0 | strpos(line, "stands in for") > 0 | strpos(line, "pooled estimand if it continued") > 0
  qui count
  local n_`ff' = r(N)
  qui gen n = _n
  qui tempfile tf_`ff'
  qui save `tf_`ff''
}
if `n_ptpd_t3a' != `n_ptpd_t3b' | `n_ptpd_t3a' == 0 {
  local bad 1
}
else {
  qui use `tf_ptpd_t3a', clear
  qui rename line linea
  qui merge 1:1 n using `tf_ptpd_t3b', nogen
  qui count if linea != line
  if r(N) > 0 local bad 1
}
restore
if `bad' {
  dis as error "T3 FAIL: live and replay roth note lines differ"
  local ++fails
}

* ---- T4: NFPTP-04 - below-size level reports mdt_nis = bias_nis = 0 ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(3) espostlength(3) pretrendpower(6 80)
if _rc != 0 {
  dis as error "T4 FAIL: rc=" _rc
  local ++fails
}
else {
  local nsize = e(roth_nissize_att)
  tempname RP
  matrix `RP' = e(roth_power_att)
  if `nsize' <= 0.06 {
    dis as text "T4 NOTE: NIS size `nsize' <= 0.06 on this fixture; boundary cell not exercised (skipping)"
  }
  else {
    if `RP'[1,3] != 0 | `RP'[1,5] != 0 {
      dis as error "T4 FAIL: p6 mdt_nis=" `RP'[1,3] " bias_nis=" `RP'[1,5] " (expected exact 0; NIS size=" `nsize' ")"
      local ++fails
    }
    if `RP'[2,3] == 0 | missing(`RP'[2,3]) {
      dis as error "T4 FAIL: p80 mdt_nis=" `RP'[2,3] " (expected positive)"
      local ++fails
    }
  }
}

* ---- T5: DOC-03 - escomp coverage scalars displayed under summarystats ----
use "`DATA'/dta_binary", clear
qui log using ptpd_t5.txt, text replace name(ptpd)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) summarystats
qui log close ptpd
if strpos(fileread("ptpd_t5.txt"), "Reference-normalization pool: N =") == 0 | strpos(fileread("ptpd_t5.txt"), "Entire treated estimation sample: N =") == 0 {
  dis as error "T5 FAIL: coverage-scalar display lines missing under summarystats"
  local ++fails
}

* ---- T6: NFADV-07 - altmethod alone r(198); the explicit pair still runs ----
use "`DATA'/dta_binary", clear
cap wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) altmethod_adversarial
if _rc != 198 {
  dis as error "T6 FAIL: altmethod-alone rc=" _rc " (expected 198, was a silent no-op)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) adversarial altmethod_adversarial
if _rc != 0 {
  dis as error "T6b FAIL: adversarial altmethod pair rc=" _rc
  local ++fails
}
else {
  * NFADV-01/02: the adversarial p-value is esf-scoped; this pooled-reference pair
  * still runs (the NFADV-07 contract), and since 2026-09-23 (author ruling: no line-only
  * results) it returns neither a p-value nor the line
  cap confirm scalar e(adversarial_pval_att)
  if _rc == 0 {
    dis as error "T6c FAIL: adversarial p returned on the pooled-reference pair (esf-scoped since NFADV-01/02)"
    local ++fails
  }
  cap confirm matrix e(adversarial_line_att)
  if _rc == 0 {
    dis as error "T6c FAIL: adversarial line returned without a p-value on the pooled-reference pair"
    local ++fails
  }
}

* ---- T7 (author ruling 2026-09-23): the pre-test note says both pre-tests are under the
* normal reference, in the live block (wooldid_ptp) and on replay (wooldid_display). Reads
* T3's two logs, written at linesize 255 so the line is not wrapped ----
local t7line "Both pre-tests are at the 5% level under the normal reference, whatever level() is in force; level() moves the slope CI only."
foreach ff in ptpd_t3a ptpd_t3b {
  if strpos(fileread("`ff'.txt"), "`t7line'") == 0 {
    dis as error "T7 FAIL: `ff'.txt (" cond("`ff'" == "ptpd_t3a", "live", "replay") ") lacks the pre-test note naming the normal reference"
    local ++fails
  }
}

if `fails' == 0 dis as result "PTPDISP_ALL_PASS"
else dis as error "PTPDISP_FAILS=`fails'"
