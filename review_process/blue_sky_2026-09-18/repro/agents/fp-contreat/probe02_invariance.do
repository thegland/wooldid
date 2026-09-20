clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/gen_ct.do"

* --- 0. how does Stata expand i.a#c.x vs i.a##c.x ? -------------------------
sysuse auto, clear
qui reg price i.foreign##c.mpg
mat lb2 = e(b)
local n2 : colnames lb2
dis "COLNAMES_DD=`n2'"
qui reg price i.foreign#c.mpg
mat lb1 = e(b)
local n1 : colnames lb1
dis "COLNAMES_SG=`n1'"

* --- 1. invariance battery -------------------------------------------------
foreach nz in none basic byTime byCohort byTreated byTtre {
 foreach cc in "basic" "byTreated" "byTtre" "byCohort" "byTime" "basic byTime" "basic byCohort" "basic byTreated" "basic byTtre" {
   qui do "`GEN'" `nz'
   qui sum betatrue if ttre<. & t>=ttre
   local tru = r(mean)
   cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(`cc') robust regtol(14)
   local rc = _rc
   if `rc' == 0 {
     mat M = e(wdidmainresults)
     local a = M[2,1]
     local rd = abs(`a'-`tru')/max(abs(`tru'),1e-12)
     local r2 = e(r2)
     dis "RES|`nz'|`cc'|ame=" %12.8f `a' "|truth=" %12.8f `tru' "|reld=" %10.3e `rd' "|r2=" %12.10f `r2'
   }
   else dis "RES|`nz'|`cc'|RC=`rc'"
 }
}
dis "DONE"
