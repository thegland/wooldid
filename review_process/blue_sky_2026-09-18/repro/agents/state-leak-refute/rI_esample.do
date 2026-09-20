* rI_esample.do -- SKEPTIC new-probe: is the POSTED e(sample) usable after a wooldid run?
* sthlp:1426 promises: "e(sample): A function characterizing the estimation sample used;
* useful for identifying which observations are dropped by wooldid during its sample
* formation and sanitation."
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*--------------------------------------------------- UNIT-LEVEL panel (collapse is a no-op)
clear
set obs 30
gen int st = _n
gen int onset = .
replace onset = 2008 if st <= 8
replace onset = 2011 if st > 8 & st <= 16
replace onset = 2013 if st > 16 & st <= 22
gen double sfe = rnormal()
expand 10
bys st: gen int yr = 2004 + _n
gen double y = sfe + 0.02*(yr-2004) + 0.35*((yr>=onset) & onset<.) + rnormal()*0.5
gen byte keepme = (st != 5)
save "rI_unit.dta", replace

*--------------------------------------------------- MICRODATA (many rows per (i,t) cell)
clear
set obs 30
gen int st = _n
gen int onset = .
replace onset = 2008 if st <= 8
replace onset = 2011 if st > 8 & st <= 16
replace onset = 2013 if st > 16 & st <= 22
gen double sfe = rnormal()
expand 10
bys st: gen int yr = 2004 + _n
expand 40
gen double y = sfe + 0.02*(yr-2004) + 0.35*((yr>=onset) & onset<.) + rnormal()
gen byte keepme = (st != 5)
save "rI_micro.dta", replace

cap prog drop esamptest
prog def esamptest
  syntax , tag(string)
  dis "ES|`tag'|eN|"      e(N)
  dis "ES|`tag'|dataN|"   c(N)
  cap count if e(sample)
  dis "ES|`tag'|countrc|" _rc
  if _rc == 0 dis "ES|`tag'|countN|" r(N)
  else        dis "ES|`tag'|countN|" -1
  cap su y if e(sample), meanonly
  dis "ES|`tag'|surc|"    _rc
  if _rc == 0 dis "ES|`tag'|suN|" r(N)
  else        dis "ES|`tag'|suN|" -1
end

*=========== I1: unit-level, default
use "rI_unit.dta", clear
qui wooldid y st yr onset, cluster(st) att
esamptest , tag(I1_unit_default)

*=========== I2: microdata, default (fast collapse active)
use "rI_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att
esamptest , tag(I2_micro_default)

*=========== I3: microdata, nofastcollapse
use "rI_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att nofastcollapse
esamptest , tag(I3_micro_nofast)

*=========== I4: microdata with real drops -- can e(sample) identify them?
use "rI_micro.dta", clear
replace y = . if st == 7 & yr == 2009
qui wooldid y st yr onset, cluster(st) att
esamptest , tag(I4_micro_drops)
dis "I4_NDROP=" e(N_dropped)

*=========== I5: unit-level with an if() restriction
use "rI_unit.dta", clear
qui wooldid y st yr onset if keepme == 1, cluster(st) att
esamptest , tag(I5_unit_if)
dis "I5_NDROP=" e(N_dropped)

*=========== I6: microdata with an if() restriction
use "rI_micro.dta", clear
qui wooldid y st yr onset if keepme == 1, cluster(st) att
esamptest , tag(I6_micro_if)
dis "I6_NDROP=" e(N_dropped)

*=========== I7: reference behaviour -- what reghdfe does on the same microdata
use "rI_micro.dta", clear
gen byte post = (yr>=onset) & onset<.
qui reghdfe y post, absorb(st yr) cluster(st)
esamptest , tag(I7_reghdfe_ref)

*=========== I8: permute(shuffle) state residue + row order (properly pinned)
use "rI_unit.dta", clear
gen long myrow = _n
sort myrow
qui datasignature
local d0 "`r(datasignature)'"
cap noi wooldid y st yr onset, cluster(st) att permute(shuffle) permreps(20)
dis "I8_RC=" _rc
qui datasignature
dis "I8_DSIG_UNCHANGED=" ("`r(datasignature)'"=="`d0'")
cap assert myrow == _n
dis "I8_ROWORDER=" (_rc==0)
dis "I8_N=" c(N) " K=" c(k) " SORTED=" "`: sortedby'"
qui frames dir
dis "I8_FRAMES=" "`r(frames)'"
local gl : all globals
dis "I8_GLOBALS=" "`gl'"

*=========== I9: bootstrap state residue + row order (properly pinned)
use "rI_unit.dta", clear
gen long myrow = _n
sort myrow
qui datasignature
local d0 "`r(datasignature)'"
cap noi wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
dis "I9_RC=" _rc
qui datasignature
dis "I9_DSIG_UNCHANGED=" ("`r(datasignature)'"=="`d0'")
cap assert myrow == _n
dis "I9_ROWORDER=" (_rc==0)
qui frames dir
dis "I9_FRAMES=" "`r(frames)'"

dis "ALLDONE"
