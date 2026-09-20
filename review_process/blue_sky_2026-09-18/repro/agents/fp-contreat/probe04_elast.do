clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/gen_ct.do"

* why did unconditionalse refuse?
qui do "`GEN'" basic
cap noisily wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) unconditionalse regtol(14)
dis "UNCSE_RC=" _rc

* ---- elasticity family: exact fit, so yhat should be y -------------------
foreach shift in 0 1000 {
 foreach path in "DEFAULT" "customvce(robust)" {
  foreach et in dydx dyex eydx eyex {
    qui do "`GEN'" basic
    qui replace y = y + `shift'
    * truths, computed on the treated post rows with yhat == y (exact fit)
    qui sum betatrue if ttre<. & t>=ttre
    local T_dydx = r(mean)
    tempvar q
    qui gen double `q' = betatrue*dose
    qui sum `q' if ttre<. & t>=ttre
    local T_dyex = r(mean)
    qui replace `q' = betatrue/y
    qui sum `q' if ttre<. & t>=ttre
    local T_eydx = r(mean)
    qui replace `q' = betatrue*dose/y
    qui sum `q' if ttre<. & t>=ttre
    local T_eyex = r(mean)
    local vce robust
    local opt
    if "`path'" != "DEFAULT" {
      local vce
      local opt `path'
    }
    cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') `vce' `opt' regtol(14)
    local rc = _rc
    if `rc' == 0 {
      mat M = e(wdidmainresults)
      local a = M[2,1]
      local rd = abs(`a'-`T_`et'')/abs(`T_`et'')
      dis "EL|shift=`shift'|`path'|`et'|ame=" %16.10f `a' "|truth=" %16.10f `T_`et'' "|reld=" %10.3e `rd'
    }
    else dis "EL|shift=`shift'|`path'|`et'|RC=`rc'"
  }
 }
}
dis "DONE"
