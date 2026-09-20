*=====================================================================
* pH2: FUZZ-02 controls -- strictly positive outcome, contreat forms,
*      the permute() rep-0 message, and whether nofastcollapse helps.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

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
gen double dose = D*(0.5 + runiform())
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i
tempfile gd
qui save "`gd'", replace
qui sum g
dis "H2_GMIN=" r(min) " H2_GMAX=" r(max)

* (1) strictly positive outcome, same design, same options
use "`gd'", clear
qui replace g = g + 50
qui sum g
dis "H2_POSMIN=" r(min)
cap noi qui wooldid g i t ttre , att cluster(i) semielasticity
dis "H2_POS_SEMI_i_RC=" _rc
cap noi qui wooldid g icoh t ttre , att cluster(icoh) semielasticity
dis "H2_POS_SEMI_coh_RC=" _rc
cap noi qui wooldid g i t ttre , att robust semielasticity safety(off)
dis "H2_POS_SEMI_rob_RC=" _rc

* (2) straddling zero: contreat elasticity forms
use "`gd'", clear
foreach f in dydx dyex eydx eyex {
  cap noi qui wooldid g i t ttre , att cluster(i) contreat(dose) contreatcontrols(basic) contreatelasticitytype(`f')
  dis "H2_CET_`f'_RC=" _rc
}

* (3) permute() leg
use "`gd'", clear
cap noi wooldid g i t ttre , att cluster(i) semielasticity permute(shuffle) permreps(12) rseed(3) dots(0)
dis "H2_PERM_SEMI_RC=" _rc

* (4) how far does a sign-straddling outcome have to be shifted before it works?
foreach sh in 0 1 2 3 5 10 {
  use "`gd'", clear
  qui replace g = g + `sh'
  qui sum g
  local mn = r(min)
  cap noi qui wooldid g i t ttre , att cluster(i) semielasticity
  dis "H2_SHIFT_`sh'_RC=" _rc "  ymin=" %9.4f `mn'
}
dis "=================== DONE pH2 ==================="
