* =================================================================
* probe_sort.do -- a probe NOBODY suggested: is bootstrap() reproducible
* at a fixed rseed() when the SAME data arrive in a different ROW ORDER?
* The count leg draws on a map frame it sorts itself (ado:9240 `qui sort
* clusv`), so it should be order-invariant. The physical / re-identification
* legs (ado:9246, ado:9344) `frame copy` the user's frame and bsample it
* unsorted, so they may not be.
* A user who sorts by (year unit) instead of (unit year), or who merged in a
* different order, has exactly this.
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_sort.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 313131
set sortseed 313131

set obs 30
gen long unit = _n
gen int ttre = .
replace ttre = 2004 if unit <= 10
replace ttre = 2006 if unit > 10 & unit <= 20
gen double ufe = rnormal()
expand 8
bysort unit: gen int year = 2000 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.5, 0))
gen double y  = 1 + ufe + 0.1*(year-2000) + te + rnormal()*0.7
sort unit year
save "sk_sortP.dta", replace          /* panel order */
gen double u = runiform()
sort u
drop u
save "sk_sortS.dta", replace          /* same rows, shuffled order */

local SEEDV 777001
local NREP 80

* --- COUNT leg (att only), clustered ---
foreach o in P S {
  use "sk_sort`o'.dta", clear
  cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(`NREP') rseed(`SEEDV') dots(0)
  matrix K`o' = e(bootdist)
  dis "S_COUNTCL_`o'_B="  %21x el(e(b),1,1)
  dis "S_COUNTCL_`o'_SE=" %21x sqrt(el(e(V),1,1))
}
dis "S_COUNTCL_MRELDIF=" mreldif(KP, KS)

* --- RE-ID leg (att att_i), clustered ---
foreach o in P S {
  use "sk_sort`o'.dta", clear
  cap noi wooldid y unit year ttre, att att_i cluster(unit) bootstrap(standard) reps(`NREP') rseed(`SEEDV') dots(0)
  matrix R`o' = e(bootdist)
  dis "S_REIDCL_`o'_B="  %21x el(e(b),1,1)
  dis "S_REIDCL_`o'_SE=" %21x sqrt(el(e(V),1,1))
}
matrix RPa = RP[1...,1]
matrix RSa = RS[1...,1]
dis "S_REIDCL_MRELDIF=" mreldif(RPa, RSa)

* --- COUNT leg, UNCLUSTERED ---
foreach o in P S {
  use "sk_sort`o'.dta", clear
  cap noi wooldid y unit year ttre, att bootstrap(standard) reps(`NREP') rseed(`SEEDV') dots(0)
  matrix U`o' = e(bootdist)
  dis "S_COUNTUN_`o'_SE=" %21x sqrt(el(e(V),1,1))
}
dis "S_COUNTUN_MRELDIF=" mreldif(UP, US)

* --- RE-ID leg, UNCLUSTERED ---
foreach o in P S {
  use "sk_sort`o'.dta", clear
  cap noi wooldid y unit year ttre, att att_i bootstrap(standard) reps(`NREP') rseed(`SEEDV') dots(0)
  matrix V`o' = e(bootdist)
  dis "S_REIDUN_`o'_SE=" %21x sqrt(el(e(V),1,1))
}
matrix VPa = VP[1...,1]
matrix VSa = VS[1...,1]
dis "S_REIDUN_MRELDIF=" mreldif(VPa, VSa)

dis "PROBE_SORT_DONE"
