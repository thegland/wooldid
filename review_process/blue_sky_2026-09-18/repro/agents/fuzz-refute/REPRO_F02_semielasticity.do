clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ============ FUZZ-02: semi-elasticity on a zero-straddling outcome ============
* realistic: the outcome is an annual growth rate
clear
set seed 606
set sortseed 606
qui set obs 60
gen long unit = _n
gen double a_i = rnormal(0,0.5)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.4 & u <= 0.6
qui replace ttre = 2005 if u > 0.6 & u <= 0.8
qui replace ttre = 2007 if u > 0.8
qui expand 10
bysort unit: gen int t = 1999 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double g = a_i + 0.3*rnormal() + 0.6*D + 0.2*x1
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i
tempfile gd
qui save "`gd'", replace
qui sum g
dis "outcome range " r(min) " .. " r(max)
cap noi qui wooldid g i t ttre , att cluster(i)
dis "CONTROL (no semielasticity) rc=" _rc
* every configuration fails the same way on this data
cap noi wooldid g i t ttre , att cluster(i) semielasticity
dis "SEMI i=unit cluster(i) rc=" _rc
cap noi qui wooldid g icoh t ttre , att cluster(icoh) semielasticity
dis "SEMI i=cohort rc=" _rc
cap noi qui wooldid g icoh t ttre , att cluster(icoh) semielasticity nofastcollapse
dis "SEMI i=cohort nofastcollapse rc=" _rc
cap noi qui wooldid g i t ttre , att robust semielasticity safety(off)
dis "SEMI robust safety(off) rc=" _rc
* under permute() the message names the wrong cause
cap noi wooldid g i t ttre , att cluster(i) semielasticity permute(shuffle) permreps(24) rseed(3) dots(0)
dis "SEMI + permute rc=" _rc
* ATTRIBUTION: plain margins gives the clean r(459); wooldid gives r(3598)+tempvar
use "`gd'", clear
qui reghdfe g x1 i.icoh#i.t#i.D, absorb(i t) vce(cluster i)
cap noi margins, dydx(x1)
dis "plain margins dydx rc=" _rc
cap noi margins, eydx(x1)
dis "plain margins eydx rc=" _rc "   <- 459 is genuine margins pass-through"
* TRIGGER: the minimum SATURATED CELL MEAN, not the raw outcome minimum
foreach sh in 0 1 2 {
  use "`gd'", clear
  qui replace g = g + `sh'
  qui egen double cm = mean(g), by(icoh t D)
  qui sum cm
  local cmin = r(min)
  qui sum g
  local ymin = r(min)
  cap noi qui wooldid g i t ttre , att cluster(i) semielasticity
  dis "shift `sh': ymin=" %9.4f `ymin' "  min cell mean=" %9.4f `cmin' "  rc=" _rc
}
