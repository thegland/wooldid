clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260917
set sortseed 20260917

local G 15
local U 4
local T 7
clear
set obs `G'
gen long cty = _n
gen double a_g = rnormal(0,.4)
gen onset = .
replace onset = 4 if cty <= 5
replace onset = 6 if cty > 5 & cty <= 10
expand `U'
sort cty
gen long unit = _n
gen double b_i = rnormal(0,.4)
gen double uf  = rnormal(0,1)
expand `T'
bysort unit (cty): gen int t = _n
gen double tau = 1 + a_g + b_i
gen byte D = (onset < . & t >= onset)
gen double y = 2 + uf + .1*t + tau*D + rnormal(0,1)
keep y unit cty t onset
compress
save "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-cellwt-refute/smoke.dta", replace
tab onset cty if t==1

timer clear
timer on 1
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(49) rseed(777) dots(0)
dis "SMOKE_RC=" _rc
timer off 1
timer list 1
dis "N=" e(N)
mat b = e(b)
mat V = e(V)
mat list b
mat list V
