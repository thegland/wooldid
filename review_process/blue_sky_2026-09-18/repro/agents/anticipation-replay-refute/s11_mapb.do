clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916
* MAP B, measured with TYPED literals: for every legal 2-dp level L, does
*   set level L   ->   global g = c(level)   ->   set level $g
* round-trip?  (this is exactly what wooldid.ado:575 + ado:69 do)
tempname fh
file open `fh' using "s11_mapb.csv", write replace
file write `fh' "L,rc_set,rendered,rc_restore" _n
local nbad = 0
local nset = 0
forvalues k = 1000/9999 {
  local ip = floor(`k'/100)
  local fp : display %02.0f mod(`k',100)
  local L "`ip'.`fp'"
  cap set level `L'
  local rcs = _rc
  if `rcs' != 0 {
    local ++nset
    file write `fh' "`L',`rcs',,"  _n
    continue
  }
  local g = c(level)
  cap set level `g'
  local rcr = _rc
  if `rcr' != 0 local ++nbad
  file write `fh' "`L',0,`g',`rcr'" _n
}
file close `fh'
set level 95
dis "@@MAPB set_refused=`nset' restore_failed=`nbad' of 9000"
dis "S11_DONE"
