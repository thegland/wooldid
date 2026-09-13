* ============================================================================
* REPRO_ADV: adversarial line-test family, one build per Stata process.
* BUILD = mid  (adopath dir <workdir>/mid)
* Configs (exact option strings from test/suite/configs*.csv):
*   esf_adv, adv_gapb, adv_altb, cov_advpool, ext_b1_h11_f1  (value cells)
*   adv_gapa (old r198), adv_r2001 (old r2001), trig_sgadv_deeppre (old r2000) - crash configs LAST
* Dumps: <workdir>/mid/out/<cfg>.csv in %21x (compare with test/compare.py or hexdec.py)
*
* ARGUMENTS (positional, both required):
*   `1' repo    - the wooldid_dev checkout root (test/data and test/suite/data)
*   `2' workdir - scratch directory holding old/, mid/, new/ (each an extracted wooldid.ado)
* EXACT INVOCATION:
*   stata-mp -e do <repro>/adversarial/mid/repro_adv_mid.do <repo> <workdir>
* ============================================================================
args repo W_
local W "`W_'"
local BUILD "mid"
if "`repo'" == "" | "`W'" == "" {
  di as error "usage: do repro_adv_mid.do <repo-root> <workdir>"
  exit 198
}
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
adopath ++ "`W'/`BUILD'"
cap mkdir "`W'/`BUILD'/out"
local D "`repo'/test/data"
local S "`repo'/test/suite/data"

cap program drop dumpres
program define dumpres
  args cfg rc
  tempname fh
  file open `fh' using "`W'/`BUILD'/out/`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    if _rc != 0 continue
    local R = rowsof(`M')
    local C = colsof(`M')
    local cn : colnames `M'
    file write `fh' "matnames|`m'|`R'x`C'|`cn'" _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  file close `fh'
end

* print the decisive adversarial objects for one config
cap program drop showadv
program define showadv
  args cfg rc
  dis as result "ADV[`cfg'][mid] rc=`rc'"
  foreach s in adversarial_pval_att adversarial_pval_att_s0 adversarial_pval_att_s1 adversarial_pval_att_s2 {
    cap confirm scalar e(`s')
    if _rc == 0 dis as result "ADV[`cfg'][mid] e(`s') = " %21x e(`s') " = " %20.12e e(`s')
    else dis as result "ADV[`cfg'][mid] e(`s') NOT POSTED"
  }
  foreach m in adversarial_line_att adversarial_line_att_s0 adversarial_line_att_s1 {
    cap confirm matrix e(`m')
    if _rc == 0 {
      tempname L
      matrix `L' = e(`m')
      local cn : colnames `L'
      dis as result "ADV[`cfg'][mid] e(`m') cols=[`cn'] [1,1]=" %21x `L'[1,1] " = " %20.12e `L'[1,1] "  [1,2]=" %21x `L'[1,2] " = " %20.12e `L'[1,2]
    }
    else dis as result "ADV[`cfg'][mid] e(`m') NOT POSTED"
  }
end

* full state reset between configs: a crashed old-build config can leave the session
* inside its adversarial frame with stored estimates pending
cap program drop freshstart
program define freshstart
  args dtafile
  cap frame change default
  cap est drop _all
  ereturn clear
  use "`dtafile'", clear
end

* ---- 1. esf_adv (core; dta_binary) ----
freshstart "`D'/dta_binary.dta"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2) adversarial
local rc = _rc
showadv esf_adv `rc'
dumpres esf_adv `rc'

* ---- 2. adv_gapb (core; dta_path_relyrgap) ----
freshstart "`D'/dta_path_relyrgap.dta"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-4) esprelength(3) espostlength(3) adversarial
local rc = _rc
showadv adv_gapb `rc'
dumpres adv_gapb `rc'

* ---- 3. adv_altb (core; dta_binary) ----
freshstart "`D'/dta_binary.dta"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-3) esprelength(2) espostlength(2) adversarial altmethod_adversarial
local rc = _rc
showadv adv_altb `rc'
dumpres adv_altb `rc'

* ---- 4. cov_advpool (covars; dta_s_covars) - pooled reference ----
freshstart "`S'/dta_s_covars.dta"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) adversarial
local rc = _rc
showadv cov_advpool `rc'
dumpres cov_advpool `rc'

* ---- 5. ext_b1_h11_f1 (extgrid; binary_balanced) - two real horizons ----
freshstart "`S'/binary_balanced.dta"
cap noi wooldid y_single cohort t ttre_single, att esprelength(1) espostlength(1) esfixedbaseperiod esrelativeto(-1) jointtests adversarial makeplots safety(off)
local rc = _rc
showadv ext_b1_h11_f1 `rc'
dumpres ext_b1_h11_f1 `rc'

* ---- 6. adv_gapa (core; dta_path_relyrgap) - old r198 ----
freshstart "`D'/dta_path_relyrgap.dta"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(4) espostlength(4) adversarial
local rc = _rc
showadv adv_gapa `rc'
dumpres adv_gapa `rc'

* ---- 7. adv_r2001 (core; dta_binary) - old r2001 ----
freshstart "`D'/dta_binary.dta"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(1) espostlength(2) adversarial altmethod_adversarial
local rc = _rc
showadv adv_r2001 `rc'
dumpres adv_r2001 `rc'

* ---- 8. trig_sgadv_deeppre (trig; dta_s_deeppre) - old r2000 ----
freshstart "`S'/dta_s_deeppre.dta"
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esfixedbaseperiod esrelativeto(-3) esprelength(2) adversarial altmethod_adversarial
local rc = _rc
showadv trig_sgadv_deeppre `rc'
dumpres trig_sgadv_deeppre `rc'

dis "REPRO_B16_DONE"
dis "REPRO_ADVPOOL_DONE"
dis "REPRO_NFADV_DONE"
dis "REPRO_BS12_DONE"
dis "REPRO_ADV_DONE"
exit, clear
