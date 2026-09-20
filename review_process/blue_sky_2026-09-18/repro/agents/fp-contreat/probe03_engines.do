clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/gen_ct.do"

* the nuisance planted in y is ALWAYS the one the contreatcontrols() asked for,
* so a correct fit is exact and the AME must equal mean(betatrue) on treated post
foreach cc in basic byTime byCohort byTreated byTtre {
  foreach path in "DEFAULT" "ccnoabsorb" "unconditionalse" "customvce(robust)" {
    qui do "`GEN'" `cc'
    qui sum betatrue if ttre<. & t>=ttre
    local tru = r(mean)
    local opt
    if "`path'" != "DEFAULT" local opt `path'
    local vce robust
    if "`path'" == "unconditionalse" | "`path'" == "customvce(robust)" local vce
    cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(`cc') `vce' `opt' regtol(14)
    local rc = _rc
    if `rc' == 0 {
      mat M = e(wdidmainresults)
      local a = M[2,1]
      local rd = abs(`a'-`tru')/abs(`tru')
      dis "ENG|`cc'|`path'|ame=" %14.10f `a' "|truth=" %14.10f `tru' "|reld=" %10.3e `rd' "|r2=" %14.12f e(r2)
    }
    else dis "ENG|`cc'|`path'|RC=`rc'"
  }
}
dis "DONE"
