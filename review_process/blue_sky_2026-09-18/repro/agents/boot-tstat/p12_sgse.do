clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid
local gG 24
local gUnits 3
local gT 7
local gnever 6
include "`WD'/gen_small.do"
gen byte reg2 = mod(state,2)

dis "=== V1: analytic (no bootstrap), cluster(state) ==="
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(reg2)
dis "V1_RC=" _rc
mat V1M = e(wdidmainresults)
mat V1A = e(wdidsgresults_0)
mat V1B = e(wdidsgresults_1)
dis "V1 main att  b=" %12.9f V1M[1,1] "  se=" %21x V1M[1,2]
dis "V1 sg0  att  b=" %12.9f V1A[1,1] "  se=" %21x V1A[1,2]
dis "V1 sg1  att  b=" %12.9f V1B[1,1] "  se=" %21x V1B[1,2]
dis "V1_SE_RATIO_sg0_main=" %12.9f V1A[1,2]/V1M[1,2]
dis "V1_SE_RATIO_sg1_main=" %12.9f V1B[1,2]/V1M[1,2]

dis "=== V2: same but robust (no cluster) ==="
cap noi wooldid y state t ttre , robust att esprelength(2) espostlength(2) subgroups(reg2)
dis "V2_RC=" _rc
mat V2M = e(wdidmainresults)
mat V2A = e(wdidsgresults_0)
mat V2B = e(wdidsgresults_1)
dis "V2 main att  se=" %21x V2M[1,2]
dis "V2 sg0  att  se=" %21x V2A[1,2]
dis "V2 sg1  att  se=" %21x V2B[1,2]
dis "V2_SE_RATIO_sg0_main=" %12.9f V2A[1,2]/V2M[1,2]
dis "V2_SE_RATIO_sg1_main=" %12.9f V2B[1,2]/V2M[1,2]

dis "=== V3: separate runs on each subgroup (the honest comparison) ==="
cap noi wooldid y state t ttre if reg2==0 , cluster(state) att esprelength(2) espostlength(2)
dis "V3a_RC=" _rc
mat V3A = e(wdidmainresults)
dis "V3 reg2==0 alone: b=" %12.9f V3A[1,1] "  se=" %12.9f V3A[1,2]
cap noi wooldid y state t ttre if reg2==1 , cluster(state) att esprelength(2) espostlength(2)
dis "V3b_RC=" _rc
mat V3B = e(wdidmainresults)
dis "V3 reg2==1 alone: b=" %12.9f V3B[1,1] "  se=" %12.9f V3B[1,2]
dis "P12_DONE"
