*=====================================================================
* f3_dots.do -- SELF-CONTAINED repro.
* dots() is range-checked in wooldid_bootctl (ado:6272) which returns
* early at its stray-companion screen whenever bootstrap() is absent
* (ado:6225).  So dots(-5) is r(198) under bootstrap() and accepted
* SILENTLY under permute(), where it then behaves as dots(0).
* The ado's own comment at ado:6268-6271 admits this.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 60
gen long unit = _n
gen byte arm  = mod(_n,3)
gen int ttre  = .
replace ttre  = 2004 if arm==1
replace ttre  = 2006 if arm==2
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double a_i = rnormal()
gen double y   = a_i + 0.10*(t-2000) + 1.50*(t>=ttre & ttre<.) + rnormal()
keep unit t ttre y
tempfile base
save `base'

use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(20) dots(-5)
dis "BOOT_DOTSNEG_RC=" _rc

use `base', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) dots(-5)
dis "PERM_DOTSNEG_RC=" _rc
dis "F3_DONE"
