clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/gen_ct.do"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat"

qui do "`GEN'" basic
qui save "`D'/fix1.dta", replace

* truths on the exact-fit data (yhat == y)
use "`D'/fix1.dta", clear
tempvar q
qui gen double `q' = betatrue
qui sum `q' if ttre<. & t>=ttre
local T_dydx = r(mean)
qui replace `q' = betatrue*dose
qui sum `q' if ttre<. & t>=ttre
local T_dyex = r(mean)
qui replace `q' = betatrue/y
qui sum `q' if ttre<. & t>=ttre
local T_eydx = r(mean)
qui replace `q' = betatrue*dose/y
qui sum `q' if ttre<. & t>=ttre
local T_eyex = r(mean)
dis "TRUTH dydx=" %18.12f `T_dydx' " dyex=" %18.12f `T_dyex' " eydx=" %18.12f `T_eydx' " eyex=" %18.12f `T_eyex'

foreach et in dydx dyex eydx eyex {
  foreach path in "robust" "customvce(robust)" "unconditionalse robust" {
    use "`D'/fix1.dta", clear
    cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) contreatelasticitytype(`et') `path' regtol(14)
    local rc = _rc
    if `rc'==0 {
      mat M = e(wdidmainresults)
      dis "SAME|`et'|`path'|ame=" %18.12f M[2,1] "|truth=" %18.12f `T_`et'' "|reld=" %10.3e abs(M[2,1]-`T_`et'')/abs(`T_`et'')
    }
    else dis "SAME|`et'|`path'|RC=`rc'"
  }
}

* ---- what prediction is margins using after reghdfe? ---------------------
use "`D'/fix1.dta", clear
qui gen byte post = (ttre<. & t>=ttre)
qui egen long cell = group(i t) if post
qui replace cell = 0 if cell==.
qui gen double ringer = dose
qui reghdfe y i.cell##c.dose , absorb(i.i i.t i.cell#c.ringer) tol(1e-14)
qui predict double xbhat , xb
qui gen double diffy = y - xbhat
qui sum diffy
dis "PRED xb-after-reghdfe: mean(y-xb)=" %12.6f r(mean) "  sd=" %12.6f r(sd)
qui sum xbhat if post
dis "PRED mean xb on post=" %12.6f r(mean)
qui sum y if post
dis "PRED mean y  on post=" %12.6f r(mean)
dis "DONE"
