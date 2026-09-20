* ================================================================================
* THE THIRD POST-FAILURE e() STATE: a PARTIAL wooldid surface with e(cmd)=wooldid
* recon F1 enumerated two states (prior e() untouched; reghdfe's internal model).
* There is a third: a failure AFTER wooldid's own `ereturn post' and AFTER
* `ereturn local cmd "wooldid"' (ado:4051 / 4198) leaves a wooldid-labelled e()
* that is missing whatever the run had not yet added.
* ================================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid

cap prog drop DUMPNAMES
prog def DUMPNAMES
  args tag
  tempname fh
  file open `fh' using "partial_names.txt", write text append
  file write `fh' "T|`tag'|CMD|`e(cmd)'" _n
  foreach s in `: e(scalars)' {
    file write `fh' "T|`tag'|SCA|`s'" _n
  }
  foreach s in `: e(macros)' {
    file write `fh' "T|`tag'|MAC|`s'" _n
  }
  foreach s in `: e(matrices)' {
    file write `fh' "T|`tag'|MAT|`s'" _n
  }
  file close `fh'
end
cap erase partial_names.txt

use D1.dta, clear
gen long sg7 = 1000000 + mod(id,2)
gen long sg8 = 10000000 + mod(id,2)

* ---------- 1. the reference: the same spec with 7-digit values, which fits ----------
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg7) jointtests adversarial cluster(clu)
dis "GOOD_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b))
DUMPNAMES GOOD
local ngood : word count `: e(scalars)'
dis "GOOD_NSCALARS=`ngood'"

* ---------- 2. the failing twin: 8-digit values overflow adversarial_*_s<v> ----------
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg8) jointtests adversarial cluster(clu)
local frc = _rc
dis "FAIL_RC=`frc' cmd=`e(cmd)' K=" colsof(e(b)) " sglist=`e(wdid_sglist)'"
DUMPNAMES FAIL
dis "FAIL_NSCALARS=" `: word count `: e(scalars)''
dis "LASTRC=${WOOLDID_LASTRC}"
dis "cmdline=`e(cmdline)'"
* does the surface LOOK valid to a consumer?
cap confirm matrix e(wdidmainresults)
dis "HAS_MAINRESULTS=" (_rc==0)
if _rc == 0 {
  tempname MM
  matrix `MM' = e(wdidmainresults)
  dis "MAINRESULTS DIM=" rowsof(`MM') "x" colsof(`MM')
  matlist `MM', format(%9.4f)
}
cap noi ereturn display
dis "EREDISP_RC=" _rc
cap noi lincom main_att
dis "LINCOM_RC=" _rc " est=" r(estimate) " se=" r(se)
cap noi estimates store PARTIAL
dis "ESTSTORE_RC=" _rc
cap noi wooldid
dis "BAREREPLAY_RC=" _rc
* the documented defeat of the replay guard
macro drop WOOLDID_LASTRC
cap noi wooldid
dis "BAREREPLAY_AFTER_MACRODROP_RC=" _rc

* ---------- 3. which joint/adversarial scalars survived? ----------
foreach s in joint_horiz joint_horiz_i joint_att_allsgssame joint_pre_att_allsgssame adversarial_pval_att adversarial_line_att {
  cap confirm scalar e(`s')
  local a = (_rc==0)
  cap confirm matrix e(`s')
  local b = (_rc==0)
  dis "SURV|`s'|scalar=`a'|matrix=`b'"
}
foreach v in 10000000 10000001 {
  foreach s in joint_horiz_s`v' adversarial_pval_att_s`v' {
    cap confirm scalar e(`s')
    dis "SURVSG|`s'|exists=" (_rc==0) "|len=" strlen("`s'")
  }
}

* ---------- 4. the reghdfe state, for contrast ----------
use D1.dta, clear
gen long sg8b = 10000000 + mod(id,2)
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod subgroups(sg8b)
dis "REGHDFE_RC=" _rc " cmd=`e(cmd)' K=" colsof(e(b)) " props=`e(properties)'"
DUMPNAMES REGHDFESTATE
cap noi predict double ph
dis "PREDICT_AFTER_REGHDFE_RC=" _rc
cap sum ph
dis "PREDICT_N=" r(N) " mean=" r(mean)
cap drop ph
cap noi estimates store REGH
dis "ESTSTORE_AFTER_REGHDFE_RC=" _rc
cap noi wooldid
dis "BAREREPLAY_AFTER_REGHDFE_RC=" _rc

dis "PARTIAL_DONE"
