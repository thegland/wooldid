* ============================================================================
* repro_inference.do — inference-convention families, one build per Stata process.
* Run from INSIDE the build's own subdirectory so logs/dumps/plots do not collide:
*   cd <workdir>/old && stata-mp -e do <repro>/inference/repro_inference.do old <repo> <workdir>
*   cd <workdir>/new && stata-mp -e do <repro>/inference/repro_inference.do new <repo> <workdir>
* ARGUMENTS (positional, all three required):
*   `1' build   - "old" or "new"; sources <workdir>/<build>/wooldid.ado
*   `2' repo    - the wooldid_dev checkout root (test/data and test/suite/data)
*   `3' workdir - scratch directory holding old/ and new/
* Blocks: BS10+B4 (impv), A18 (trig_uncond_jt), A16 (trig_imdegen),
*         A17 (trig_imdeg_part), TM2+CS60 (trig_tm2_sgjoint), DDG
*         (ext_cov10_rcs_pweight_ddd_c1), BS15 (adm_svpsg, slow, last).
* Each block prints REPRO_<id>_DONE. Dumps in %21x go to <cfg>.csv in cwd.
* ============================================================================
args build repo W_
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
local W "`W_'"
if "`build'" == "" | "`repo'" == "" | "`W'" == "" {
  di as error "usage: do repro_inference.do <old|new> <repo-root> <workdir>"
  exit 198
}
local D "`repo'/test/data"
local S "`repo'/test/suite/data"
do "`W'/`build'/wooldid.ado"
dis "BUILD=`build'"

* ---- dumpres: every e() scalar/macro/matrix cell in %21x (harness template) ----
cap program drop dumpres
program define dumpres
  args cfg rc
  tempname fh
  file open `fh' using "`cfg'.csv", write replace
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

* ============================ BS-10 + B4: impv ==============================
ereturn clear
use "`D'/dta_binary.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) impvalues
dis "BS10_RC=" _rc
dumpres impv `=_rc'
dis "BS10_df_r=" e(df_r)
dis "BS10_df_r_analytic=" e(df_r_analytic)
dis "BS10_imdf_min=" e(imdf_min) " imdf_max=" e(imdf_max)
dis "B4_ll=" %21.12g e(ll)
matrix WM = e(wdidmainresults)
local wmr : rownames WM
dis "BS10_wdidmainresults_rows=[`wmr']"
dis "BS10_table_pretreatment_p=" %14.10g WM[2,4] " lb=" %14.10g WM[2,5] " ub=" %14.10g WM[2,6]
cap noi lincom pretreatment_att
dis "BS10_lincom_rc=" _rc
dis "BS10_lincom_p=" %14.10g r(p) " lb=" %14.10g r(lb) " ub=" %14.10g r(ub) " df=" r(df)
dis "BS10_2ttail_at_e_df_r=" %14.10g 2*ttail(e(df_r), abs(WM[2,1]/WM[2,2])) " at_t49=" %14.10g 2*ttail(49, abs(WM[2,1]/WM[2,2])) " at_t79=" %14.10g 2*ttail(79, abs(WM[2,1]/WM[2,2]))
dis "REPRO_BS10_DONE"
dis "REPRO_B4_DONE"

* ============================ A18: trig_uncond_jt ===========================
ereturn clear
use "`D'/dta_binary.dta", clear
cap noi wooldid y i t ttre [aw=w], att robust unconditionalse esprelength(2) espostlength(2) jointtests
dis "A18_RC=" _rc
dumpres trig_uncond_jt `=_rc'
dis "A18_joint_horiz=" %14.10g e(joint_horiz)
dis "REPRO_A18_DONE"

* ============================ A16: trig_imdegen =============================
ereturn clear
use "`S'/dta_s_imdegen.dta", clear
cap noi wooldid y i t ttre, att cluster(clus) esprelength(1) espostlength(1) impvalues
dis "A16_RC=" _rc
dumpres trig_imdegen `=_rc'
dis "A16_properties=[`e(properties)']"
matrix WM = e(wdidmainresults)
dis "A16_main_att_est=" %14.10g WM[1,1] " se=" WM[1,2] " t=" WM[1,3] " p=" WM[1,4] " lb=" WM[1,5] " ub=" WM[1,6]
dis "A16_ll=" %14.8g e(ll)
cap confirm matrix e(V)
dis "A16_eV_confirm_rc=" _rc
if _rc == 0 dis "A16_eV_11=" e(V)[1,1]
dis "REPRO_A16_DONE"

* ============================ A17: trig_imdeg_part ==========================
ereturn clear
use "`D'/dta_path_imbalanced.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(9) impvalues
dis "A17_RC=" _rc
dumpres trig_imdeg_part `=_rc'
dis "A17_properties=[`e(properties)']"
dis "A17_df_r=" e(df_r) " df_r_analytic=" e(df_r_analytic) " imdf_min=" e(imdf_min) " imdf_max=" e(imdf_max)
cap confirm matrix e(V)
dis "A17_eV_confirm_rc=" _rc
if _rc == 0 {
  matrix VD = vecdiag(e(V))
  local vn : colnames VD
  local k 0
  foreach n of local vn {
    local ++k
    dis "A17_Vdiag " %-22s "`n'" " = " VD[1,`k']
  }
}
matrix WM = e(wdidmainresults)
local wmr : rownames WM
local k 0
foreach n of local wmr {
  local ++k
  dis "A17_table " %-22s "`n'" " est=" %10.6g WM[`k',1] " se=" %10.6g WM[`k',2] " p=" %10.4g WM[`k',4]
}
dis "REPRO_A17_DONE"

* ===================== TM2 + CS-60: trig_tm2_sgjoint ========================
ereturn clear
drop _all
set obs 15
gen i = _n
gen ttre = cond(_n <= 4, 8, cond(_n <= 8, 10, cond(_n <= 12, 12, .)))
expand 6
bysort i: gen unit = _n
gen sg = mod(i, 3) + 1
gen dose = exp(0.11*mod(i*7 + unit*3, 10) + 0.05*mod(unit, 4))
expand 12
bysort i unit: gen t = _n
gen rel = t - ttre
gen y = 0.05*i + 0.04*t + 0.35*dose*(rel >= 0 & rel < .) + 0.06*(sg == 2)*dose*cond(rel < 0 & rel < ., rel + 8, 0) + 0.10*(sg == 3)*dose*cond(rel < 0 & rel < ., rel + 8, 0) + 0.3*sin(2.7*i + 1.3*t + 0.9*unit)
sort i unit t
ereturn clear
cap noi wooldid y i t ttre, att robust contreat(dose) contreatcontrols(basic) esprelength(3) espostlength(3) subgroups(sg) jointtests
dis "TM2_RC=" _rc
dumpres trig_tm2_sgjoint `=_rc'
dis "TM2_joint_Choriz_pooled=" e(joint_Choriz)
dis "TM2_joint_horiz_pooled=" e(joint_horiz)
forval s = 1/3 {
  dis "TM2_joint_Choriz_s`s'_macroform=[`e(joint_Choriz_s`s')']"
  dis "TM2_joint_Choriz_s`s'_exprform=" e(joint_Choriz_s`s')
  dis "CS60_joint_horiz_s`s'_macroform=[`e(joint_horiz_s`s')']"
  dis "CS60_joint_horiz_s`s'_exprform=" e(joint_horiz_s`s')
}
local esc : e(scalars)
local emc : e(macros)
dis "CS60_scalars_named_joint_horiz_s: " cond(strpos("`esc'", "joint_horiz_s1") > 0, "YES", "no")
dis "CS60_macros_named_joint_horiz_s:  " cond(strpos("`emc'", "joint_horiz_s1") > 0, "YES", "no")
cap noi dis "CS60_numeric_compare_result=" (e(joint_horiz_s1) < .05)
dis "CS60_numeric_compare_rc=" _rc
dis "REPRO_TM2_DONE"
dis "REPRO_CS60_DONE"

* ===================== DDG: ext_cov10_rcs_pweight_ddd_c1 ====================
ereturn clear
use "`S'/survey_rcs.dta", clear
cap noi wooldid y_rcs state_group t ttre [pw=survey_weight], att att_i att_it att_itime subgroups(group) fe(group_time) esfixedbaseperiod esrelativeto(-1) esprelength(1) espostlength(1) cluster(state) safety(off)
dis "DDG_RC=" _rc
dumpres ext_cov10_rcs_pweight_ddd_c1 `=_rc'
dis "DDG_pwsgs:"
matlist e(pwsgs), format(%12.4g)
foreach it in att att_it att_i att_itime {
  dis "DDG_joint_`it'_allsgssame=" e(joint_`it'_allsgssame) "   joint_pre_`it'_allsgssame=" e(joint_pre_`it'_allsgssame)
}
dis "DDG_sgdegen=" e(wdid_sgdegen) " act=[`e(wdid_sgdegenact)'] r=[`e(wdid_sgdegenr)'] fcguard=[`e(wdid_fcguard)']"
dis "REPRO_DDG_DONE"

* ===================== BS-15: adm_svpsg (slow, last) ========================
ereturn clear
use "`S'/dta_s_admin.dta", clear
qui keep if mod(fips - 1000*floor(fips/1000), 7) == 1
cap noi wooldid y fips year ttre, att cluster(state) subgroups(region) esprelength(1) espostlength(1) makeplots saveplots(svpadm)
dis "BS15_RC=" _rc
dumpres adm_svpsg `=_rc'
dis "BS15_pwsgs:"
matlist e(pwsgs), format(%12.4g)
dis "BS15_joint_att_allsgssame=" e(joint_att_allsgssame)
dis "BS15_joint_pre_att_allsgssame=" e(joint_pre_att_allsgssame)
dis "BS15_sgdegen=" e(wdid_sgdegen) " act=[`e(wdid_sgdegenact)'] r=[`e(wdid_sgdegenr)']"
* implied correlation per main-att pairwise row, from the level SE columns already in e():
*   rho = (V_a + V_b - V_ab) / (2*sqrt(V_a*V_b)); the level SE is wdidsgresults_<sg>[main:att, se]
matrix PW = e(pwsgs)
local prn : rowfullnames PW
local r 0
foreach rn of local prn {
  local ++r
  if strpos("`rn'", "att:") != 1 continue
  local nm = subinstr("`rn'", "att:", "", 1)
  if !regexm("`nm'", "^([0-9]+)(bn)?vs([0-9]+)(bn)?[.]") continue
  local a = regexs(1)
  local b = regexs(3)
  matrix WSA = e(wdidsgresults_`a')
  matrix WSB = e(wdidsgresults_`b')
  local sea = WSA[1,2]
  local seb = WSB[1,2]
  local seab = PW[`r',2]
  local rho = (`sea'^2 + `seb'^2 - `seab'^2) / (2*sqrt(`sea'^2*`seb'^2))
  dis "BS15_row " %-18s "`rn'" " est=" %9.5f PW[`r',1] " se_contrast=" %10.6f `seab' " se_a=" %9.6f `sea' " se_b=" %9.6f `seb' " rho=" %8.5f `rho' " p=" %10.3g PW[`r',4]
}
dis "REPRO_BS15_DONE"
dis "REPRO_ALL_DONE"
exit, clear
