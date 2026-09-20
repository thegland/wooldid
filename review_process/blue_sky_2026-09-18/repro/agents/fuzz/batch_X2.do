* exactness runner.  spec fields: id|ds|y|i|ifclause|wclause|wvar|opts
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/truth.do"

file open FS using "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/spec_X2.txt", read text
file open FO using "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/res_X2.csv", write replace
file write FO "id|ds|rc|worst|detail|cmd" _n
file read FS line
while r(eof) == 0 {
  local L `"`line'"'
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
  local f8 = substr(`"`L'"',strpos(`"`L'"',"|")+1,.)

  qui use "`D'/ds_`f2'.dta", clear
  dis "@@BEGIN|`f1'|`f2'|wooldid `f3' `f4' t ttre `f5' `f6', att att_it att_i att_itime `f8'"
  ereturn clear
  set seed 4242
  set sortseed 4242
  capture noisily wooldid `f3' `f4' t ttre `f5' `f6', att att_it att_i att_itime `f8'
  local rc = _rc
  local worst 0
  local det ""
  if `rc' == 0 {
    capture {
      qui gen byte __insmp = e(sample)
      matrix WM = e(wdidmainresults)
      local reqs : roweq WM
      local rnms : rownames WM
      local R = rowsof(WM)
      wdfz_truth , smp(__insmp) ivar(`f4') wvar(`f7')
      foreach e in att att_it att_i att_itime {
        local TR_`e' = r(`e')
      }
      forvalues r = 1/`R' {
        local eq : word `r' of `reqs'
        local nm : word `r' of `rnms'
        if "`eq'" != "main" continue
        if !inlist("`nm'","att","att_it","att_i","att_itime") continue
        local got = WM[`r',1]
        local tru = `TR_`nm''
        local d = abs(`got' - `tru') / max(1, abs(`tru'))
        if `d' > `worst' local worst = `d'
        if `d' > 1e-8 local det "`det' `nm'(got=`got' true=`tru' rel=`d')"
      }
      drop __insmp
    }
    if _rc local det "TRUTHFAIL(`=_rc')"
  }
  dis "@@END|`f1'|rc=`rc'|worst=`worst'|`det'"
  file write FO "`f1'|`f2'|`rc'|`worst'|" `"`det'"' "|wooldid `f3' `f4' t ttre `f5' `f6', att att_it att_i att_itime `f8'" _n
  file read FS line
}
file close FS
file close FO
dis "EXACT_DONE"
