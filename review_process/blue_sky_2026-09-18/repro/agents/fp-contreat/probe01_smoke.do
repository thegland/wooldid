clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/gen_ct.do" basic

dis "NOBS=" _N
tab i t if ttre<. & t>=ttre

cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) robust regtol(14)
dis "SMOKE_RC=" _rc
ereturn list
matrix list e(wdidmainresults)

* independent truth: simple mean of betatrue over treated post rows
qui sum betatrue if ttre<. & t>=ttre
dis "TRUTH_AME_unweighted = " %21x r(mean)
