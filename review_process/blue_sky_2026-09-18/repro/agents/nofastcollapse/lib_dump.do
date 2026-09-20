* ---- differential-testing library for the nofastcollapse probe ----
cap prog drop wddump
prog def wddump
  syntax , FILE(string)
  tempname fh
  file open `fh' using "`file'", write text replace
  local ss : e(scalars)
  foreach s of local ss {
    local v = string(e(`s'), "%21x")
    file write `fh' `"S|`s'|`v'"' _n
  }
  local ms : e(macros)
  foreach m of local ms {
    local v `"`e(`m')'"'
    file write `fh' `"M|`m'|`v'"' _n
  }
  local mm : e(matrices)
  foreach m of local mm {
    tempname M
    matrix `M' = e(`m')
    local r = rowsof(`M')
    local c = colsof(`M')
    local rn : rownames `M'
    local cn : colnames `M'
    local re : roweq `M'
    local ce : coleq `M'
    file write `fh' `"X|`m'|`r'|`c'|`rn'|`cn'|`re'|`ce'"' _n
    forvalues ii = 1/`r' {
      local line
      forvalues jj = 1/`c' {
        local q = string(`M'[`ii',`jj'], "%21x")
        local line `line' `q'
      }
      file write `fh' `"V|`line'"' _n
    }
  }
  file close `fh'
end

* run one cell: arm A (collapse default) vs arm B (nofastcollapse)
cap prog drop wdcell
prog def wdcell
  syntax , DTA(string) TAG(string) BASE(string) OUT(string) [OPTS(string)]
  * ---- arm A : collapse ON by default ----
  qui use "`dta'", clear
  qui datasignature set, reset
  local sigA = r(datasignature)
  local nA = _N
  local kA = c(k)
  cap noisily `base', `opts'
  local rcA = _rc
  local fcolA "`e(wdid_fcol)'"
  local fcskipA "`e(wdid_fcskip)'"
  local fcguardA "`e(wdid_fcguard)'"
  if `rcA' == 0 wddump , file("`out'/`tag'.A.txt")
  qui datasignature
  local sigA2 = r(datasignature)
  local nA2 = _N
  local kA2 = c(k)
  local leak 0
  if "`sigA'" != "`sigA2'" local leak 1
  * ---- arm B : nofastcollapse ----
  qui use "`dta'", clear
  qui datasignature set, reset
  local sigB = r(datasignature)
  cap noisily `base', `opts' nofastcollapse
  local rcB = _rc
  local fcolB "`e(wdid_fcol)'"
  if `rcB' == 0 wddump , file("`out'/`tag'.B.txt")
  qui datasignature
  local sigB2 = r(datasignature)
  local nB2 = _N
  local kB2 = c(k)
  local leakB 0
  if "`sigB'" != "`sigB2'" local leakB 1
  dis "WDCELL|`tag'|rcA=`rcA'|rcB=`rcB'|fcolA=`fcolA'|fcolB=`fcolB'|guardA=`fcguardA'|leakA=`leak'|leakB=`leakB'|skipA=`fcskipA'"
  dis "WDCELLSIG|`tag'|preA=`sigA'|postA=`sigA2'|postB=`sigB2'|NA=`nA'->`nA2'|kA=`kA'->`kA2'|NB2=`nB2'|kB2=`kB2'"
end
