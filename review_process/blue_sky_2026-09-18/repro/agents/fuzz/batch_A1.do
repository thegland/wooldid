* generic fuzz runner.  Expects globals: FZSPEC (spec file), FZOUT (result csv)
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/assert.do"

file open FS using "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/spec_A1.txt", read text
file open FO using "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/res_A1.csv", write replace
file write FO "id|ds|rc|secs|eN|nviol|viol|cmd" _n
file read FS line
while r(eof) == 0 {
  local L `"`line'"'
  * fields: id | ds | yvar | ivar | ifclause | wclause | levtxt | eng | opts
  local f1 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f2 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f3 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f4 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f5 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f6 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f7 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local L  = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)
  local f8 = substr(`"`L'"',1,strpos(`"`L'"',"|")-1)
  local f9 = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)

  qui use "`D'/ds_`f2'.dta", clear
  wdfz_snap
  local n0 `r(nobs)'
  local v0 `r(nvar)'
  local s0 `r(sort)'
  local g0 `r(sig)'
  dis "@@BEGIN|`f1'|`f2'|wooldid `f3' `f4' t ttre `f5' `f6', `f9'"
  timer clear 1
  timer on 1
  ereturn clear
  set seed 4242
  set sortseed 4242
  capture noisily wooldid `f3' `f4' t ttre `f5' `f6', `f9'
  local rc = _rc
  timer off 1
  qui timer list 1
  local secs = r(t1)
  local viol ""
  if `rc' == 0 {
    capture wdfz_assert, id(`f1') levtxt(`f7') eng(`f8') nobs0(`n0') nvar0(`v0') sort0(`s0') sig0(`g0')
    if _rc local viol "ASSERTFAIL(`=_rc')"
    else local viol `"`r(viol)'"'
  }
  else {
    capture wdfz_snap
    if "`r(sig)'" != "`g0'" local viol "DATA_SIG_AFTER_ERROR"
  }
  local nv : word count `viol'
  dis "@@END|`f1'|rc=`rc'|nviol=`nv'|`viol'"
  file write FO "`f1'|`f2'|`rc'|`secs'|`=e(N)'|`nv'|" `"`viol'"' "|wooldid `f3' `f4' t ttre `f5' `f6', `f9'" _n
  file read FS line
}
file close FS
file close FO
dis "RUNNER_DONE"
