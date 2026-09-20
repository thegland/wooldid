clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"
which wooldid

use "dataA.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(permA) replace
dis "SMOKE_RC=" _rc
dis "PERM_IREBUILT=" e(perm_irebuilt)
dis "PERM_M=" e(perm_M)
mat list e(b), format(%21x)
cap mat list e(perm_comparator_b), format(%21x)
ereturn list

use "permA.dta", clear
desc
list in 1/5
