* p02b_ctw_map.do -- contreatwithin + contreat() rewrites the WORKING safety value to
* off (ado:4829).  Which guards does that silently disable, and does a number move?
* SELF-CONTAINED.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

local gN 90
local gT 10
local gt0 2000
local gonsets 2004 2005 2006
local gdose unitconst
local gysd 0
local gaisd 0
local gtrendsd 0
include "`GEN'/gen_core.do"
gen long cohort0 = cohort - 1          // a PURE RELABEL: 0-indexed cohort codes
save "p02b.dta", replace

prog def dumpres
  syntax , tag(string)
  cap mat drop M
  cap mat M = e(wdidmainresults)
  cap dis "`tag'_N=" e(N)
  cap mat list M, format(%21x)
  dis "`tag'_END"
end

dis "=== D. contreat+contreatwithin, NO safety() typed, 1-indexed i (reference) ==="
capture noisily wooldid ync cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) contreatwithin
dis "D_RC=" _rc
dumpres, tag(D)

dis "=== E. same run, 0-indexed i.  ado:4959 would have refused this under safety(on) ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0) contreat(dose) contreatcontrols(basic) contreatwithin
dis "E_RC=" _rc
dumpres, tag(E)

dis "=== E2. same 0-indexed run with an EXPLICIT safety(off): documented bypass ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0) contreat(dose) contreatcontrols(basic) contreatwithin safety(off)
dis "E2_RC=" _rc
dumpres, tag(E2)

dis "=== A. binary, default safety, 0-indexed i : the guard DOES fire here ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0)
dis "A_RC=" _rc

dis "=== B. binary, safety(off), 0-indexed i vs C. binary default, 1-indexed i ==="
capture noisily wooldid ync cohort0 t ttre, att cluster(cohort0) safety(off)
dis "B_RC=" _rc
dumpres, tag(B)
capture noisily wooldid ync cohort t ttre, att cluster(cohort)
dis "C_RC=" _rc
dumpres, tag(C)

dis "=== G. non-integer i under contreatwithin ==="
gen double cohorthalf = cohort/2
capture noisily wooldid ync cohorthalf t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) contreatwithin
dis "G_RC=" _rc

dis "=== H. the float-2^24 caution (ado:4973), with and without contreatwithin ==="
use "p02b.dta", clear
gen float tdate    = 20000101 + (t - 2000)*10000
gen float ttredate = 20000101 + (ttre - 2000)*10000
capture noisily wooldid ync cohort tdate ttredate, att cluster(cohort)
dis "H1_RC=" _rc
capture noisily wooldid ync cohort tdate ttredate, att cluster(cohort) contreat(dose) contreatcontrols(basic) contreatwithin
dis "H2_RC=" _rc
dis "P02B_DONE"
