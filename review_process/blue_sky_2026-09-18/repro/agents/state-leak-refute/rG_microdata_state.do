* rG_microdata_state.do -- SKEPTIC new-probe: does wooldid leak state on the FAST-COLLAPSE
* path (microdata: many rows per (i,t) cell), and does it survive xtset / a caller preserve?
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

*---------------------------------------------------------------- MICRODATA DGP
* 30 states x 10 years x 40 individuals = 12,000 rows; i = state (many rows per cell)
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
order myrow
save "rG_micro.dta", replace

*================== G1: default (fast collapse) -- full state snapshot
use "rG_micro.dta", clear
sort myrow
snap , tag(G1PRE)
qui wooldid y st yr onset, cluster(st) att
snap , tag(G1POST)
cap assert myrow == _n
dis "G1_ROWORDER_PRESERVED=" (_rc==0)
dis "G1_eN=" e(N) "  G1_NDROP=" e(N_dropped)
matrix G1 = e(wdidmainresults)
matrix list G1, format(%21x)

*================== G2: nofastcollapse -- same answer? same state?
use "rG_micro.dta", clear
sort myrow
snap , tag(G2PRE)
qui wooldid y st yr onset, cluster(st) att nofastcollapse
snap , tag(G2POST)
cap assert myrow == _n
dis "G2_ROWORDER_PRESERVED=" (_rc==0)
dis "G2_eN=" e(N) "  G2_NDROP=" e(N_dropped)
matrix G2 = e(wdidmainresults)
matrix list G2, format(%21x)
mata: st_numscalar("DG", max(abs(st_matrix("G1"):-st_matrix("G2"))))
mata: st_numscalar("RG", max(abs((st_matrix("G1"):-st_matrix("G2")):/editmissing(abs(st_matrix("G1")),1))))
dis "G_FASTvsNOFAST_MAXABS=" %21x DG
dis "G_FASTvsNOFAST_MAXREL=" %21x RG

*================== G3: xtset survival
use "rG_micro.dta", clear
bys st yr: gen int pid = _n
qui xtset st yr
dis "G3_XTSET_BEFORE_PANEL=" "`r(panelvar)'" " TIME=" "`r(timevar)'"
qui wooldid y st yr onset, cluster(st) att
cap noi xtset
dis "G3_XTSET_AFTER_RC=" _rc
dis "G3_XTSET_AFTER_PANEL=" "`r(panelvar)'" " TIME=" "`r(timevar)'"
char li _dta[]
dis "G3_DTA_IIS=" "`: char _dta[iis]'" " TIS=" "`: char _dta[tis]'"

*================== G4: caller's own preserve straddling a wooldid call
use "rG_micro.dta", clear
qui datasignature
local d0 "`r(datasignature)'"
preserve
qui drop if yr > 2010
dis "=====G4 BEGIN====="
cap noi wooldid y st yr onset, cluster(st) att
dis "G4_INNER_RC=" _rc
restore
dis "=====G4 END====="
qui datasignature
dis "G4_DSIG_RESTORED=" ("`r(datasignature)'"=="`d0'") "  N=" c(N)

*================== G5: verbose on microdata -- what does the retained marker mark?
use "rG_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att verbose
cap confirm variable _est____wooldidfullmodel
dis "G5_MARKER_EXISTS=" (_rc==0)
if _rc == 0 {
  qui count if _est____wooldidfullmodel == 1
  dis "G5_MARKER_ONES=" r(N) "  DATA_N=" c(N) "  eN=" e(N)
  qui su _est____wooldidfullmodel, meanonly
  dis "G5_MARKER_MIN=" r(min) " MAX=" r(max)
}
qui est restore ___wooldidfullmodel
qui count if e(sample)
dis "G5_RESTORED_ESAMPLE_N=" r(N)

*================== G6: bootstrap / permute state, with PRE snapshots this time
use "rG_micro.dta", clear
set seed 4321
snap , tag(G6PRE)
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
snap , tag(G6POST)
cap assert myrow == _n
dis "G6_ROWORDER_PRESERVED=" (_rc==0)

use "rG_micro.dta", clear
set seed 4321
snap , tag(G7PRE)
cap noi wooldid y st yr onset, cluster(st) att permute(onset) permreps(20)
dis "G7_RC=" _rc
snap , tag(G7POST)
cap assert myrow == _n
dis "G7_ROWORDER_PRESERVED=" (_rc==0)

*================== G8: rseed() reproducibility across a moved RNG stream
use "rG_micro.dta", clear
set seed 111
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
matrix R1 = e(wdidmainresults)
use "rG_micro.dta", clear
set seed 999
qui gen double zz = runiform()
drop zz
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20) rseed(4242)
matrix R2 = e(wdidmainresults)
mata: st_numscalar("D8", max(abs(st_matrix("R1"):-st_matrix("R2"))))
dis "G8_RSEED_REPRO_MAXDIFF=" %21x D8

*================== G9: no rseed -- does set seed alone pin the bootstrap?
use "rG_micro.dta", clear
set seed 2468
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20)
matrix S1 = e(wdidmainresults)
use "rG_micro.dta", clear
set seed 2468
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(20)
matrix S2 = e(wdidmainresults)
mata: st_numscalar("D9", max(abs(st_matrix("S1"):-st_matrix("S2"))))
dis "G9_SETSEED_REPRO_MAXDIFF=" %21x D9

dis "ALLDONE"
