clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916

* exact message from set level on a 16-digit rendering
cap noi set level 99.90000000000001
dis "MSG_RC=" _rc
set level 95

* for every 2-dp level in a grid, render it the way `syntax , Level(real 0)` does
* (i.e. macro-assign the double) and ask set level whether it accepts that rendering.
cap prog drop s2p
prog def s2p
  syntax [, Level(real 0)]
  c_local rendered "`level'"
end
tempname fh
file open `fh' using "s2_map.csv", write replace
file write `fh' "L,rendered,rc_direct,rc_rendered" _n
forvalues k = 1000/9999 {
  local L = `k'/100
  * only probe every 7th to keep it fast, plus all of 99.xx and the .x0 family
  local doit = 0
  if mod(`k',7)==0 local doit = 1
  if `k' >= 9900 local doit = 1
  if mod(`k',10)==0 local doit = 1
  if `doit' == 0 continue
  cap set level `L'
  local rcd = _rc
  s2p , level(`L')
  cap set level `rendered'
  local rcr = _rc
  file write `fh' "`L',`rendered',`rcd',`rcr'" _n
}
file close `fh'
set level 95
dis "S2_DONE"
