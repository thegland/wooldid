* Is the P0016 yscale mismatch the reghdfe alternating-projections floor?
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
foreach rt in 9 12 13 {
  use "`D'/ds_canon.dta", clear
  qui wooldid y unit t ttre if !missing(x1), att cluster(unit region) esprelength(1) espostlength(3) anticipation(2) level(99) regtol(`rt')
  matrix A = e(wdidmainresults)
  use "`D'/ds_canon.dta", clear
  qui replace y = 2.5*y
  qui wooldid y unit t ttre if !missing(x1), att cluster(unit region) esprelength(1) espostlength(3) anticipation(2) level(99) regtol(`rt')
  matrix B = e(wdidmainresults)
  local R = rowsof(A)
  local worst1 = 0
  local worst3 = 0
  forvalues r = 1/`R' {
    local a1 = 2.5*A[`r',1]
    local b1 = B[`r',1]
    if !missing(`a1') & !missing(`b1') {
      local d = abs(`a1'-`b1')/max(1,abs(`a1'))
      if `d' > `worst1' local worst1 = `d'
    }
    local a3 = A[`r',3]
    local b3 = B[`r',3]
    if !missing(`a3') & !missing(`b3') {
      local d = abs(`a3'-`b3')/max(1,abs(`a3'))
      if `d' > `worst3' local worst3 = `d'
    }
  }
  dis "REGTOL|`rt'|worst_estimate=`worst1'|worst_t=`worst3'"
}
dis "REGTOL_DONE"
