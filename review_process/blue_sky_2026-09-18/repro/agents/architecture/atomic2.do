clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local BD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/architecture"

dis "=== P1: poisson, plain ==="
use "`BD'/base.dta", clear
gen double ycount = ceil(exp(y))
sort tt id
cap noi wooldid ycount id tt ttre [pw=wt] , att poisson
dis "ATOM2 P1 rc=" _rc " N=" e(N)
cap confirm variable _ppmlhdfe_d
dis "ATOM2 P1 ppmld_present=" (_rc==0) " nvars=" c(k) " sortedby=`: sortedby'"
qui est dir
dis "ATOM2 P1 estdir=[`r(names)']"

dis "=== P2: poisson, verbose ==="
use "`BD'/base.dta", clear
gen double ycount = ceil(exp(y))
sort tt id
cap noi wooldid ycount id tt ttre [pw=wt] , att poisson verbose
dis "ATOM2 P2 rc=" _rc
cap confirm variable _ppmlhdfe_d
dis "ATOM2 P2 ppmld_present=" (_rc==0) " nvars=" c(k) " sortedby=`: sortedby'"
qui est dir
dis "ATOM2 P2 estdir=[`r(names)']"

dis "=== P3: poisson, y<0 on a ZERO-WEIGHT row (pw) ==="
use "`BD'/base.dta", clear
gen double ycount = ceil(exp(y))
replace ycount = -1 in 1
replace wt = 0 in 1
cap noi wooldid ycount id tt ttre [pw=wt] , att poisson
dis "ATOM2 P3 rc=" _rc
cap noi wooldid ycount id tt ttre [pw=wt] , att poisson safety(off)
dis "ATOM2 P3safetyoff rc=" _rc

dis "=== P4: binary, id==0 on zero-weight row, safety(off) ==="
use "`BD'/base.dta", clear
set obs `=_N+1'
replace id = 0 in L
replace tt = 3 in L
replace y = 1 in L
replace wt = 0 in L
cap noi wooldid y id tt ttre [aw=wt] , att safety(off)
dis "ATOM2 P4 rc=" _rc

dis "=== W: size ceiling with 100-char lines ==="
foreach k in 18 19 {
  cap prog drop _all
  clear mata
  cap noi run "`BD'/sz/wide_`k'.ado"
  dis "ATOM2 WIDE kb=`k' rc=" _rc
}
dis "ATOM2DONE"
