clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"

*==================== A. retired option names ====================
use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att fastcollapse
dis "A1_fastcollapse_RC=" _rc
cap noisily wooldid y cohort t ttre, att fullprecisionslow
dis "A2_fullprecisionslow_RC=" _rc
cap noisily wooldid y cohort t ttre, att nofastcollapse
dis "A3_nofastcollapse_RC=" _rc

*==================== B. hidden results contract ====================
use "`DD'/dC.dta", clear
qui wooldid y cohort t ttre, att esprelength(2) espostlength(2)
dis "B_fcol=" e(wdid_fcol) " fcncell=" e(wdid_fcncell) " fcskip=<" e(wdid_fcskip) "> fcguard=<" e(wdid_fcguard) ">"
qui wooldid y cohort t ttre, att esprelength(2) espostlength(2) nofastcollapse
dis "B_NOFC_fcol=" e(wdid_fcol) " fcncell=" e(wdid_fcncell) " fcskip=<" e(wdid_fcskip) "> fcguard=<" e(wdid_fcguard) ">"
* ineligible skip should CAUTION and set fcskip
qui wooldid y cohort t ttre, att esprelength(2) espostlength(2) unconditionalse
dis "B_UNC_fcol=" e(wdid_fcol) " fcskip=<" e(wdid_fcskip) ">"

*==================== C. data state: signature / sortorder / N / k ====================
foreach arm in "" "nofastcollapse" {
  use "`DD'/dC.dta", clear
  qui datasignature set, reset
  local s0 = r(datasignature)
  local n0 = _N
  local k0 = c(k)
  local sb0 : sortedby
  qui wooldid y cohort t ttre, att att_i att_it att_itime esprelength(3) espostlength(3) jointtests `arm'
  qui datasignature
  local s1 = r(datasignature)
  local sb1 : sortedby
  dis "C|arm=<`arm'>|sigsame=" ("`s0'"=="`s1'") "|N=`n0'->" _N "|k=`k0'->" c(k) "|sortedby=<`sb0'>-><`sb1'>"
}

*==================== D. user's own preserve around wooldid ====================
use "`DD'/dC.dta", clear
local n0 = _N
preserve
qui drop if t > 2005
qui count
dis "D_inner_N=" r(N)
cap noisily wooldid y cohort t ttre, att esprelength(2) espostlength(2)
dis "D1_RC=" _rc
cap noisily restore
dis "D1_restore_RC=" _rc
dis "D1_after_N=" _N " expected=`n0'"

*==================== E. error raised INSIDE the collapse window ====================
* graph save to a nonexistent directory, from wooldid_hist, which runs inside the window
use "`DD'/dC.dta", clear
qui datasignature set, reset
local s0 = r(datasignature)
local n0 = _N
cap noisily wooldid y cohort t ttre, att att_i esprelength(2) espostlength(2) histogramcohorteffects saveplots(/no/such/dir/zz)
dis "E1_RC=" _rc
qui datasignature
dis "E1_sigsame=" ("`s0'"=="`=r(datasignature)'") " N=" _N " expected=`n0'" " k=" c(k)
qui describe, short
cap confirm variable y
dis "E1_y_present_RC=" _rc

*==================== F. same, but the user had preserved first ====================
use "`DD'/dC.dta", clear
local n0 = _N
preserve
qui keep if cohort != 2
cap noisily wooldid y cohort t ttre, att esprelength(2) espostlength(2) histogramcohorteffects saveplots(/no/such/dir/zz)
dis "F1_RC=" _rc
dis "F1_N_now=" _N
cap noisily restore
dis "F1_restore_RC=" _rc
dis "F1_after_N=" _N " expected=`n0'"

*==================== G. wooldid twice in a row after an error ====================
use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att histogramcohorteffects saveplots(/no/such/dir/zz)
cap noisily wooldid y cohort t ttre, att esprelength(2) espostlength(2)
dis "G1_RC=" _rc
cap dis "G1_b_att=" %21x el(e(b),1,1)
dis "G1_disp_rc=" _rc
dis "P04_DONE"
