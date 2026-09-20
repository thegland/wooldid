* metamorphic pair runner.  spec: id|ds|y|i|ifclause|wclause|transform|opts
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_cmp
program define wdfz_cmp, rclass
  * compares matrices WM1 and WM2 under a transform rule
  syntax , RULE(string) [ SCALE(real 1) ]
  local bad ""
  local worst = 0
  local R1 = rowsof(WM1)
  local R2 = rowsof(WM2)
  if `R1' != `R2' {
    return local bad "DIM(`R1'vs`R2')"
    return scalar worst = 99
    exit
  }
  local rn1 : rownames WM1
  local rn2 : rownames WM2
  if "`rn1'" != "`rn2'" local bad "`bad' ROWNAMES"
  forvalues r = 1/`R1' {
    local nm : word `r' of `rn1'
    forvalues c = 1/6 {
      local a = WM1[`r',`c']
      local b = WM2[`r',`c']
      local tgt = .
      if "`rule'" == "same"   local tgt = `a'
      if "`rule'" == "scale" {
        if inlist(`c',1,2,5,6) local tgt = `scale'*`a'
        else                   local tgt = `a'
      }
      if "`rule'" == "negate" {
        if inlist(`c',1,3)  local tgt = -`a'
        if inlist(`c',2,4)  local tgt =  `a'
        if `c' == 5         local tgt = -WM1[`r',6]
        if `c' == 6         local tgt = -WM1[`r',5]
      }
      if "`rule'" == "pointonly" {
        if `c' != 1 continue
        local tgt = `a'
      }
      if missing(`tgt') & missing(`b') continue
      if missing(`tgt') | missing(`b') {
        local bad "`bad' MISSMATCH(`nm'c`c':`a'/`b')"
        continue
      }
      local d = abs(`b' - `tgt') / max(1, abs(`tgt'))
      if `d' > `worst' local worst = `d'
      if `d' > 1e-9 local bad "`bad' DIFF(`nm'c`c':`a'->`b' exp`tgt' rel`d')"
    }
  }
  return local bad "`bad'"
  return scalar worst = `worst'
end

file open FS using "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/spec_P1.txt", read text
file open FO using "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/res_P1.csv", write replace
file write FO "id|ds|xform|rc1|rc2|worst|bad|cmd" _n
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

  local ivar `f4'
  local ifc  `f5'
  local ifc2 `f5'
  local opt1 `f8'
  local opt2 `f8'
  local rule same
  local scl  1

  dis "@@BEGIN|`f1'|`f2'|`f7'|wooldid `f3' `f4' t ttre `f5' `f6', `f8'"
  qui use "`D'/ds_`f2'.dta", clear
  ereturn clear
  set seed 4242
  set sortseed 4242
  capture noisily wooldid `f3' `f4' t ttre `f5' `f6', `f8'
  local rc1 = _rc
  if `rc1' == 0 matrix WM1 = e(wdidmainresults)

  * ---------------- transform ----------------
  if "`f7'" == "shuffle" {
    set seed 771
    qui gen double __rr = runiform()
    sort __rr
    qui drop __rr
    local rule same
  }
  else if "`f7'" == "tshift" {
    qui replace t = t + 137
    qui replace ttre = ttre + 137
    local rule same
  }
  else if "`f7'" == "idperm" {
    qui egen long __g = group(`f4')
    qui summarize __g, meanonly
    qui gen long __ni = 500000 - __g
    local ivar __ni
    local rule same
  }
  else if "`f7'" == "yscale" {
    qui replace `f3' = 2.5 * `f3'
    local rule scale
    local scl 2.5
  }
  else if "`f7'" == "yshift" {
    qui replace `f3' = `f3' + 137
    local rule same
  }
  else if "`f7'" == "ynegate" {
    qui replace `f3' = -`f3'
    local rule negate
  }
  else if "`f7'" == "wscale" {
    qui replace w = w * 1000
    local rule same
  }
  else if "`f7'" == "ifkeep" {
    qui keep if (`=subinstr("`f5'","if ","",1)') == 1
    local ifc2 ""
    local rule same
  }
  else if "`f7'" == "collin" {
    local opt2 = subinstr("`f8'","controls(x1)","controls(x1 xcollin)",1)
    local rule same
  }
  else if "`f7'" == "rowdup" {
    local rule pointonly
    qui expand 2
    local opt2 "`f8'"
  }
  else if "`f7'" == "sortrev" {
    gsort -`f4' -t
    local rule same
  }

  ereturn clear
  set seed 4242
  set sortseed 4242
  capture noisily wooldid `f3' `ivar' t ttre `ifc2' `f6', `opt2'
  local rc2 = _rc
  local bad ""
  local worst 0
  if `rc1' == 0 & `rc2' == 0 {
    matrix WM2 = e(wdidmainresults)
    capture wdfz_cmp , rule(`rule') scale(`scl')
    if _rc local bad "CMPFAIL(`=_rc')"
    else {
      local bad `"`r(bad)'"'
      local worst = r(worst)
    }
  }
  else if `rc1' != `rc2' local bad "RCDIFF(`rc1'/`rc2')"
  dis "@@END|`f1'|rc=`rc1'|rc2=`rc2'|worst=`worst'|`bad'"
  file write FO "`f1'|`f2'|`f7'|`rc1'|`rc2'|`worst'|" `"`bad'"' "|wooldid `f3' `f4' t ttre `f5' `f6', `f8'" _n
  file read FS line
}
file close FS
file close FO
dis "PAIR_DONE"
