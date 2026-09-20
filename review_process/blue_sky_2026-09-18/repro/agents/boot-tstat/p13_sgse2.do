clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid
* 30 states x 4 units, i = UNIT (120 cohorts is too slow) -> use i = a
* district id with 2 districts per state, so cluster(state) is STRICTLY ABOVE i.
clear
set obs 30
gen int state = _n
gen int onset = .
forval s = 1/22 {
  local w = mod(`s'-1, 3) + 1
  local o : word `w' of 2002 2003 2004
  qui replace onset = `o' in `s'
}
expand 2
bys state: gen int dist = _n + 100*state
gen double a_i = rnormal()
gen byte reg2 = mod(dist, 2)
expand 3
bys dist: gen int unit = _n
gen double u_i = rnormal()*0.5
expand 6
bys dist unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (dist unit): replace tshock = rnormal() if _n==1
qui bys t (dist unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t - onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
dis "NOBS=" _N "  NDIST=" 60

dis "=== W1: i = dist (60), cluster(state) STRICTLY ABOVE i ==="
cap noi wooldid y dist t ttre , cluster(state) att esprelength(1) espostlength(1) subgroups(reg2) safety(off)
dis "W1_RC=" _rc
mat W1M = e(wdidmainresults)
mat W1A = e(wdidsgresults_0)
mat W1B = e(wdidsgresults_1)
dis "W1 main att b=" %12.9f W1M[1,1] " se=" %21x W1M[1,2]
dis "W1 sg0  att b=" %12.9f W1A[1,1] " se=" %21x W1A[1,2]
dis "W1 sg1  att b=" %12.9f W1B[1,1] " se=" %21x W1B[1,2]
dis "W1_RATIO0=" %14.11f W1A[1,2]/W1M[1,2] "  W1_RATIO1=" %14.11f W1B[1,2]/W1M[1,2]

dis "=== W2: same but cluster(dist) == i ==="
cap noi wooldid y dist t ttre , cluster(dist) att esprelength(1) espostlength(1) subgroups(reg2) safety(off)
dis "W2_RC=" _rc
mat W2M = e(wdidmainresults)
mat W2A = e(wdidsgresults_0)
mat W2B = e(wdidsgresults_1)
dis "W2 main att se=" %21x W2M[1,2]
dis "W2 sg0  att se=" %21x W2A[1,2]
dis "W2 sg1  att se=" %21x W2B[1,2]
dis "W2_RATIO0=" %14.11f W2A[1,2]/W2M[1,2] "  W2_RATIO1=" %14.11f W2B[1,2]/W2M[1,2]
dis "P13_DONE"
