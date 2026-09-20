* rC_state_catchall.do -- SKEPTIC: session-state catch-all. Settings, globals, RNG,
* sort order, row order, frames, estimates, on success AND on the error path.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

cap prog drop snap
prog def snap
  syntax , tag(string)
  dis "SNAP|`tag'|level|"        c(level)
  dis "SNAP|`tag'|emptycells|"   "`c(emptycells)'"
  dis "SNAP|`tag'|coeftabres|"   "`c(coeftabresults)'"
  dis "SNAP|`tag'|type|"         "`c(type)'"
  dis "SNAP|`tag'|varabbrev|"    "`c(varabbrev)'"
  dis "SNAP|`tag'|matsize|"      c(matsize)
  dis "SNAP|`tag'|dp|"           "`c(dp)'"
  dis "SNAP|`tag'|cformat|"      "`c(cformat)'"
  dis "SNAP|`tag'|pformat|"      "`c(pformat)'"
  dis "SNAP|`tag'|sformat|"      "`c(sformat)'"
  dis "SNAP|`tag'|showbase|"     "`c(showbaselevels)'"
  dis "SNAP|`tag'|showempty|"    "`c(showemptycells)'"
  dis "SNAP|`tag'|showomit|"     "`c(showomitted)'"
  dis "SNAP|`tag'|fvlabel|"      "`c(fvlabel)'"
  dis "SNAP|`tag'|linesize|"     c(linesize)
  dis "SNAP|`tag'|k|"            c(k)
  dis "SNAP|`tag'|N|"            c(N)
  dis "SNAP|`tag'|sortlist|"     "`: sortedby'"
  dis "SNAP|`tag'|rngstate|"     substr("`c(rngstate)'",1,40)
  dis "SNAP|`tag'|rng|"          "`c(rng)'"
  dis "SNAP|`tag'|frame|"        "`c(frame)'"
  qui datasignature
  dis "SNAP|`tag'|dsig|"         "`r(datasignature)'"
  qui estimates dir
  dis "SNAP|`tag'|estdir|"       "`r(names)'"
  local gl : all globals
  dis "SNAP|`tag'|globals|"      "`gl'"
  local mt : all matrices
  dis "SNAP|`tag'|matrices|"     "`mt'"
  local sc : all scalars
  dis "SNAP|`tag'|scalars|"      "`sc'"
  qui frames dir
  dis "SNAP|`tag'|frames|"       "`r(frames)'"
end

*---------------------------------------------------------------- own DGP
clear
set obs 40
gen int st = _n
gen int onset = .
replace onset = 2009 if st <= 10
replace onset = 2011 if st > 10 & st <= 20
replace onset = 2013 if st > 20 & st <= 30
gen double sfe = rnormal()
expand 12
bys st: gen int yr = 2004 + _n
gen byte post = (yr >= onset) & onset < .
gen double y = sfe + 0.03*(yr-2004) + 0.40*post + rnormal()*0.5
gen long myrow = _n
save "rC_base.dta", replace

*================== ARM 1: plain run, non-default session settings set first
use "rC_base.dta", clear
set emptycells keep
set level 80
set varabbrev off
set cformat %9.4f
sort myrow
snap , tag(A1PRE)
qui wooldid y st yr onset, cluster(st) att
snap , tag(A1POST)
cap assert myrow == _n
dis "A1_ROWORDER_PRESERVED=" (_rc==0)

*================== ARM 2: level(90) explicit, caller at 80
use "rC_base.dta", clear
set level 80
snap , tag(A2PRE)
qui wooldid y st yr onset, cluster(st) att level(90)
snap , tag(A2POST)
dis "A2_LEVEL_AFTER=" c(level)

*================== ARM 3: error path (r459) with level(90) and non-default emptycells
use "rC_base.dta", clear
set level 80
set emptycells keep
drop if onset >= .
snap , tag(A3PRE)
cap noi wooldid y st yr onset, cluster(st) att level(90)
dis "A3_RC=" _rc
snap , tag(A3POST)
dis "A3_LEVEL_AFTER=" c(level) "  A3_EMPTYCELLS_AFTER=" "`c(emptycells)'"

*================== ARM 4: parse-stage error with level(90)
use "rC_base.dta", clear
set level 80
snap , tag(A4PRE)
cap noi wooldid y st yr onset, cluster(st) att level(90) bootstrap(nonsense)
dis "A4_RC=" _rc
snap , tag(A4POST)
dis "A4_LEVEL_AFTER=" c(level)

*================== ARM 5: RNG -- does a PLAIN run move the stream?
use "rC_base.dta", clear
set seed 12345
local rs0 = c(rngstate)
qui wooldid y st yr onset, cluster(st) att
local rs1 = c(rngstate)
dis "A5_RNG_UNCHANGED_PLAIN=" ("`rs0'"=="`rs1'")
set seed 12345
qui wooldid y st yr onset, cluster(st) att gteffects makeplots
local rs2 = c(rngstate)
dis "A5_RNG_UNCHANGED_GTE=" ("`rs0'"=="`rs2'")

*================== ARM 6: bootstrap with rseed() -- RNG afterwards, and repeatability
use "rC_base.dta", clear
set seed 777
local rsb0 = c(rngstate)
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
local rsb1 = c(rngstate)
dis "A6_RNG_UNCHANGED_BOOT=" ("`rsb0'"=="`rsb1'")
dis "A6_BOOTSE=" %21x e(wdid_att_bootse)
snap , tag(A6POST)

*================== ARM 7: back-to-back determinism under a moved RNG stream
use "rC_base.dta", clear
set seed 111
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
matrix R1 = e(wdidmainresults)
use "rC_base.dta", clear
set seed 999
qui runiform()
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
matrix R2 = e(wdidmainresults)
mata: st_numscalar("D7", max(abs(st_matrix("R1"):-st_matrix("R2"))))
dis "A7_RSEED_REPRODUCIBLE_MAXDIFF=" %21x D7

*================== ARM 8: bootstrap WITHOUT rseed -- does set seed alone pin it?
use "rC_base.dta", clear
set seed 2468
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20)
matrix S1 = e(wdidmainresults)
use "rC_base.dta", clear
set seed 2468
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20)
matrix S2 = e(wdidmainresults)
mata: st_numscalar("D8", max(abs(st_matrix("S1"):-st_matrix("S2"))))
dis "A8_SETSEED_REPRODUCIBLE_MAXDIFF=" %21x D8

dis "ALLDONE"
