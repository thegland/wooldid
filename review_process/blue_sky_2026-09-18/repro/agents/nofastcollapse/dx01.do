* DIAGNOSTIC dx01: is margins' delta-method Jacobian entry for a CONTROL coefficient
* nonzero after reghdfe, and does collapsing the rows change it?
* Analytically the entry is EXACTLY 0: the control contributes identically to
* xb(tf=1) and xb(tf=0).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4711
qui set obs 400
gen long u = _n
gen byte g = mod(_n-1,4)+1
qui expand 10
qui bysort u: gen int t = 2000 + _n - 1
gen byte tf = (g<4 & t>=2004)
gen double x = 52000 + 9000*rnormal()
gen double y = g + 0.2*(t-2000) + 0.00005*x + 1.3*tf + rnormal()*0.3
gen double wgt = 1
qui reghdfe y i.g#i.t#i.tf x , absorb(i.g i.t)
tempvar es
qui gen byte `es' = e(sample)
qui margins [iw=wgt], dydx(tf) noestimcheck
matrix J1 = r(Jacobian)
matrix V1 = r(V)
local cn : colnames J1
local k : word count `cn'
local jx 0
forvalues j = 1/`k' {
  local nm : word `j' of `cn'
  if "`nm'" == "x" local jx `j'
}
dis "FULL  DATA : N=" _N "  J[x]=" %21x J1[1,`jx'] "   V=" %21x V1[1,1]
* now collapse to one row per (g,t,tf) cell carrying the summed weight, as wooldid does
preserve
qui keep if `es' == 1
qui gen long ord = _n
qui gcollapse (sum) wgt (min) pick = ord, by(g t tf) merge replace
qui keep if ord == pick
qui drop ord pick
qui sort g t tf
qui margins [iw=wgt], dydx(tf) noestimcheck noesample
matrix J2 = r(Jacobian)
matrix V2 = r(V)
dis "COLLAPSED  : N=" _N "  J[x]=" %21x J2[1,`jx'] "   V=" %21x V2[1,1]
restore
dis "relative move in the margins VCE: " %9.2e (abs(V1[1,1]-V2[1,1])/max(abs(V1[1,1]),abs(V2[1,1])))
dis "DX01_DONE"
