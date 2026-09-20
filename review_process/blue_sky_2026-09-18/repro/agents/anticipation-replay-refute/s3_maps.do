clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916

cap prog drop s3p
prog def s3p
  syntax [, Level(real 0)]
  c_local rendered "`level'"
end

* exact hex of the two representations for a handful
foreach L in 99.9 99.91 99.99 95 92.5 {
  s3p , level(`L')
  local plit = `rendered'
  cap set level `L'
  local rcset = _rc
  local cl = c(level)
  dis "HEX L=`L' parsed=" %21x `plit' " clevel=" %21x `cl' " render_parsed=[`rendered'] render_clevel=[`cl'] rcset=`rcset'"
}
set level 95

* FULL enumeration of both maps over all 9000 legal 2-dp levels
tempname fh
file open `fh' using "s3_full.csv", write replace
file write `fh' "L,rA_optrender,rB_clevelrender" _n
forvalues k = 1000/9999 {
  local L = `k'/100
  * map A: the option-path rendering (syntax local) fed to set level
  s3p , level(`L')
  cap set level `rendered'
  local rA = _rc
  * map B: session level L, then global = c(level), then restore
  cap set level `L'
  local rcs = _rc
  if `rcs' != 0 {
    local rB = -1
  }
  else {
    local cl = c(level)
    cap set level `cl'
    local rB = _rc
  }
  file write `fh' "`L',`rA',`rB'" _n
}
file close `fh'
set level 95
dis "S3_DONE"
