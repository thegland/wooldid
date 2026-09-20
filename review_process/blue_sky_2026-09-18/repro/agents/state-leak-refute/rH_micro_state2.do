* rH_micro_state2.do -- SKEPTIC new-probe part 2: xtset survival, caller preserve,
* the verbose marker under the fast collapse, and bootstrap/permute state residue.
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
  dis "SNAP|`tag'|level|"      c(level)
  dis "SNAP|`tag'|emptycells|" "`c(emptycells)'"
  dis "SNAP|`tag'|k|"          c(k)
  dis "SNAP|`tag'|N|"          c(N)
  dis "SNAP|`tag'|sortlist|"   "`: sortedby'"
  dis "SNAP|`tag'|frame|"      "`c(frame)'"
  qui datasignature
  dis "SNAP|`tag'|dsig|"       "`r(datasignature)'"
  qui estimates dir
  dis "SNAP|`tag'|estdir|"     "`r(names)'"
  local gl : all globals
  dis "SNAP|`tag'|globals|"    "`gl'"
  qui frames dir
  dis "SNAP|`tag'|frames|"     "`r(frames)'"
  qui ds
  dis "SNAP|`tag'|vars|"       "`r(varlist)'"
end

*---------------------------------------------------------------- MICRODATA DGP (same as rG)
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
gen double tfe = 0.02*(yr-2004)
expand 40
gen byte post = (yr >= onset) & onset < .
gen double y = sfe + tfe + 0.35*post + rnormal()
gen double indwt = 1 + runiform()
gen long myrow = _n
bys st yr: gen int pid = _n
egen long uid = group(st pid)
order myrow
save "rH_micro.dta", replace

*================== H1: xtset survival across a wooldid call
use "rH_micro.dta", clear
qui xtset uid yr
dis "H1_BEFORE_PANEL=" "`r(panelvar)'" " TIME=" "`r(timevar)'"
dis "H1_BEFORE_CHAR_IIS=" "`: char _dta[iis]'" " TIS=" "`: char _dta[tis]'"
qui wooldid y st yr onset, cluster(st) att
cap noi xtset
dis "H1_AFTER_RC=" _rc
dis "H1_AFTER_PANEL=" "`r(panelvar)'" " TIME=" "`r(timevar)'"
dis "H1_AFTER_CHAR_IIS=" "`: char _dta[iis]'" " TIS=" "`: char _dta[tis]'"

*================== H2: the caller's own preserve straddling a wooldid call
use "rH_micro.dta", clear
qui datasignature
local d0 "`r(datasignature)'"
preserve
qui drop if yr > 2011
dis "=====H2 BEGIN====="
cap noi wooldid y st yr onset, cluster(st) att
dis "H2_INNER_RC=" _rc " INNER_N=" c(N)
restore
dis "=====H2 END====="
qui datasignature
dis "H2_DSIG_RESTORED=" ("`r(datasignature)'"=="`d0'") "  N=" c(N)
cap assert myrow == _n
dis "H2_ROWORDER=" (_rc==0)

*================== H3: verbose on MICRODATA -- does the retained marker survive the collapse,
*                       and what scale does it mark?
use "rH_micro.dta", clear
dis "H3_DATA_N=" c(N)
qui wooldid y st yr onset, cluster(st) att verbose
dis "H3_POST_DATA_N=" c(N) " K=" c(k)
cap confirm variable _est____wooldidfullmodel
local mex = (_rc==0)
dis "H3_MARKER_EXISTS=" `mex'
if `mex' == 1 {
  qui count if _est____wooldidfullmodel == 1
  dis "H3_MARKER_ONES=" r(N) "  DATA_N=" c(N) "  eN=" e(N)
}
qui est restore ___wooldidfullmodel
cap noi count if e(sample)
dis "H3_RESTORED_ESAMPLE_RC=" _rc " N=" r(N)

*================== H4: bootstrap residue, PRE and POST
use "rH_micro.dta", clear
set seed 4321
snap , tag(H4PRE)
cap noi wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
dis "H4_RC=" _rc
snap , tag(H4POST)
cap assert myrow == _n
dis "H4_ROWORDER=" (_rc==0)

*================== H5: permute residue, PRE and POST
use "rH_micro.dta", clear
set seed 4321
snap , tag(H5PRE)
cap noi wooldid y st yr onset, cluster(st) att permute(onset) permreps(20)
dis "H5_RC=" _rc
snap , tag(H5POST)
cap assert myrow == _n
dis "H5_ROWORDER=" (_rc==0)

*================== H6: rseed() reproducibility across a moved RNG stream
use "rH_micro.dta", clear
set seed 111
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
matrix R1 = e(wdidmainresults)
use "rH_micro.dta", clear
set seed 999
qui gen double zz = runiform()
drop zz
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
matrix R2 = e(wdidmainresults)
mata: st_numscalar("D6", max(abs(st_matrix("R1"):-st_matrix("R2"))))
dis "H6_RSEED_REPRO_MAXDIFF=" %21x D6

*================== H7: no rseed -- does `set seed` alone pin the bootstrap?
use "rH_micro.dta", clear
set seed 2468
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20)
matrix S1 = e(wdidmainresults)
use "rH_micro.dta", clear
set seed 2468
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20)
matrix S2 = e(wdidmainresults)
mata: st_numscalar("D7", max(abs(st_matrix("S1"):-st_matrix("S2"))))
dis "H7_SETSEED_REPRO_MAXDIFF=" %21x D7

*================== H8: two back-to-back runs in one session vs each alone
use "rH_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att_i att_it att_itime att gteffects
matrix T1 = e(wdidmainresults)
use "rH_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att
matrix T2 = e(wdidmainresults)
dis "H8_T2_ROWS=" rowsof(T2) " COLS=" colsof(T2)
matrix list T2, format(%21x)
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use "rH_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att
matrix T3 = e(wdidmainresults)
mata: st_numscalar("D8", max(abs(st_matrix("T2"):-st_matrix("T3"))))
dis "H8_BACKTOBACK_MAXDIFF=" %21x D8

dis "ALLDONE"
