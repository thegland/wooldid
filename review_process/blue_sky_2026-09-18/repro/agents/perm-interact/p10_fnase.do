*==============================================================================
* p10: fastnose x permanalyticse under permute().
*
* wooldid_permctl REFUSES fastnose together with permtstat/permtci ("permtstat/permtci
* randomize on t-statistics, which fastnose does not compute").  It does NOT refuse
* fastnose together with permanalyticse, which asks for exactly the same object.
* ado:5448 only ADDS fastnose in the fast-mode default; it never removes a user-typed
* one.  So permanalyticse is expected to be a SILENT NO-OP here, while the header and
* the footer (ado:8214 / ado:11679, keyed on e(perm_analyticse)) still assert that the
* analytic standard errors and e(V)/e(V_es_*) are there.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear

cap prog drop fchk
prog def fchk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0) `opts'
  dis "FN|`tag'|RC=" _rc
  if _rc != 0 exit
  dis "FN|`tag'|perm_analyticse=" e(perm_analyticse) "|vcvstate=" "`e(wdid_vcvstate)'"
  cap confirm matrix e(V)
  dis "FN|`tag'|HASV=" (_rc==0)
  cap confirm matrix e(V_es_att)
  dis "FN|`tag'|HASVES=" (_rc==0)
  cap confirm matrix e(b_es_att)
  dis "FN|`tag'|HASBES=" (_rc==0)
  * the se_analytic column of the main results matrix (col 2)
  tempname M
  matrix `M' = e(wdidmainresults)
  dis "FN|`tag'|SE_ANALYTIC_col2_row1=" `M'[1,2]
  dis "FN|`tag'|header=" "`e(wdid_vcespec)'"
end

dis "===== BASELINE: fast-mode default ====="
fchk fastdefault
fchk ts_only , opts(permtstat)
dis "===== permanalyticse alone (should restore SEs) ====="
fchk ase , opts(permanalyticse)
dis "===== fastnose alone ====="
fchk fn , opts(fastnose)
dis "===== fastnose + permanalyticse  <-- the contradiction ====="
fchk fn_ase , opts(fastnose permanalyticse)
dis "===== for contrast: fastnose + permtstat (refused?) ====="
fchk fn_ts , opts(fastnose permtstat)
dis "===== fastnose + permtci (refused?) ====="
fchk fn_tci , opts(fastnose permtci)

dis "===== what does the FOOTER say under fastnose+permanalyticse? ====="
cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0) fastnose permanalyticse
dis "FOOTERCHK_RC=" _rc
dis "P10_DONE"
