clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"

which wooldid
dis "BUILDCHECK done"

do mkdata_R.do
use panelR.dta, clear
dis "N=" _N
timer clear 1
timer on 1
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(S_BASE, replace)
dis "SMOKE_RC=" _rc
timer off 1
timer list 1
ereturn list
mat B = e(b)
dis "BOBS %21x = " %21x B[1,1]
dis "COLNAME1 = " abbrev("`: word 1 of `: colnames B''",40)
