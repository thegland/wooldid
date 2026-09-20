clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- self-contained staggered-DiD panel (60 units, T=10, 3 onsets) ----
clear
set seed 4242
set sortseed 4242
qui set obs 60
gen long unit = _n
gen double a_i = rnormal(0,1)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.35 & u <= 0.5667
qui replace ttre = 2005 if u > 0.5667 & u <= 0.7833
qui replace ttre = 2007 if u > 0.7833
qui expand 10
bysort unit: gen int t = 1999 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double y = a_i + 0.4*(t-2000) + 2*D + 0.2*x1 + 0.3*rnormal()
gen double dose = D*(0.5 + runiform())
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i

* ============ FUZZ-03 + the wider level() defect ============
* (a) plain Stata accepts these levels
foreach L in 99.9 99.99 95.4 99.1 90.1 {
  cap noi set level `L'
  dis "PLAIN set level `L' rc=" _rc " c(level)=" c(level)
}
set level 95
* (b) wooldid refuses them
foreach L in 95 99.5 97.5 99.95 99.9 99.99 95.4 99.1 90.1 {
  cap noi qui wooldid y i t ttre , att robust safety(off) level(`L')
  dis "WOOLDID level(`L') rc=" _rc
}
* (c) MECHANISM: Stata renders macro `=' results at 16 significant digits
cap prog drop lvlecho
prog def lvlecho
  syntax , [Level(real 0)]
  dis "   syntax real echo of `level' -> [`level']"
end
foreach L in 99.5 99.9 99.95 99.99 95.4 {
  lvlecho , level(`L')
}
set level 99.9
local w = c(level)
dis "   local w = c(level) -> [`w']"
set level 95
* (d) NEW: a session level of 99.9 with NO level() option bricks every call
set level 99.9
foreach o in "att robust" "att robust safety(off)" "att unconditionalse" "att robust fastnose" {
  cap noi qui wooldid y i t ttre , `o'
  dis "SESSION99.9  wooldid , `o'   rc=" _rc
}
* ... while plain Stata commands are fine at the same session level
cap noi qui regress y x1
dis "SESSION99.9 regress rc=" _rc
cap noi qui reghdfe y x1, absorb(i t)
dis "SESSION99.9 reghdfe rc=" _rc
set level 95
