* rB_verbose_and_remedy.do -- SKEPTIC: (1) finding 02 verbose _est_ variable,
*                                       (2) finding 01 follow-ups: exact-hex A3-vs-A5 and the
*                                           "modify the sample window" remedy the message names.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*---------------------------------------------------------------- own DGP (same shape as rA)
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
gen double tfe = 0.03*(yr-2004)
gen byte post = (yr >= onset) & onset < .
gen double y = sfe + tfe + 0.40*post + rnormal()*0.5
gen double wt = 1 + runiform()
save "rB_base.dta", replace

***********************************************************************
*** PART 1 -- finding 01 follow-ups
***********************************************************************
* B1: y missing on all never-treated, safety(off) -- exact hex
use "rB_base.dta", clear
replace y = . if onset >= .
qui wooldid y st yr onset, cluster(st) att safety(off)
matrix B1 = e(wdidmainresults)
matrix list B1, format(%21x)
dis "B1_N=" e(N) " B1_NDROP=" e(N_dropped)

* B2: never-treated rows physically deleted, safety(off) -- exact hex
use "rB_base.dta", clear
drop if onset >= .
qui wooldid y st yr onset, cluster(st) att safety(off)
matrix B2 = e(wdidmainresults)
matrix list B2, format(%21x)
dis "B2_N=" e(N) " B2_NDROP=" e(N_dropped)

mata: st_numscalar("SAMEDIM", (rows(st_matrix("B1"))==rows(st_matrix("B2")))*(cols(st_matrix("B1"))==cols(st_matrix("B2"))))
dis "SAMEDIM=" SAMEDIM
if SAMEDIM == 1 {
  mata: st_numscalar("MAXABSDIFF", max(abs(st_matrix("B1"):-st_matrix("B2"))))
  dis "MAXABSDIFF_B1_B2=" %21x MAXABSDIFF
  dis "BITIDENTICAL=" (MAXABSDIFF==0)
}

* B3: the SECOND remedy the r(459) message names -- "modify the sample window manually".
*     Truncate so the last cohort becomes never-treated. Data still has y missing on the
*     original never-treated units.
use "rB_base.dta", clear
replace y = . if onset >= .
dis "=====B3 BEGIN====="
cap noi wooldid y st yr onset if yr <= 2012, cluster(st) att
dis "B3_RC=" _rc
dis "=====B3 END====="
dis "B3_N=" e(N) " B3_NDROP=" e(N_dropped)

* B4: same truncation via a manually restricted dataset
use "rB_base.dta", clear
replace y = . if onset >= .
keep if yr <= 2012
dis "=====B4 BEGIN====="
cap noi wooldid y st yr onset, cluster(st) att
dis "B4_RC=" _rc
dis "=====B4 END====="

***********************************************************************
*** PART 2 -- finding 02: the verbose _est_ marker variable
***********************************************************************
use "rB_base.dta", clear
local k0 = c(k)
qui wooldid y st yr onset, cluster(st) att
local k1 = c(k)
dis "PLAIN_K_BEFORE=" `k0' " PLAIN_K_AFTER=" `k1' " PLAIN_ADDED=" (`k1'-`k0')
cap confirm variable _est____wooldidfullmodel
dis "PLAIN_ESTVAR_EXISTS=" (_rc==0)
qui estimates dir
dis "PLAIN_ESTDIR=" "`r(names)'"

* verbose arm
use "rB_base.dta", clear
local k0 = c(k)
qui wooldid y st yr onset, cluster(st) att verbose
local k1 = c(k)
dis "VERB_K_BEFORE=" `k0' " VERB_K_AFTER=" `k1' " VERB_ADDED=" (`k1'-`k0')
cap confirm variable _est____wooldidfullmodel
dis "VERB_ESTVAR_EXISTS=" (_rc==0)
qui estimates dir
dis "VERB_ESTDIR=" "`r(names)'"

* --- 2a: does the DOCUMENTED cleanup action remove the variable too?
cap noi est drop ___wooldidfullmodel
dis "ESTDROP_RC=" _rc
cap confirm variable _est____wooldidfullmodel
dis "AFTER_ESTDROP_ESTVAR_EXISTS=" (_rc==0)
dis "AFTER_ESTDROP_K=" c(k)

* --- 2b: does the prober's proposed remedy (drop the marker var, keep the estimate) work?
use "rB_base.dta", clear
qui wooldid y st yr onset, cluster(st) att verbose
qui count if _est____wooldidfullmodel == 1
dis "ESAMPLE_ROWS_IN_MARKER=" r(N)
* simulate the proposed remedy
drop _est____wooldidfullmodel
dis "=====REMEDY BEGIN====="
cap noi est restore ___wooldidfullmodel
dis "RESTORE_RC=" _rc
cap noi count if e(sample)
dis "ESAMPLE_COUNT_RC=" _rc " N=" r(N)
cap noi estimates describe ___wooldidfullmodel
dis "ESTDESC_RC=" _rc
dis "=====REMEDY END====="

* --- 2c: does a subsequent PLAIN run clean up a prior verbose run's leftovers?
use "rB_base.dta", clear
qui wooldid y st yr onset, cluster(st) att verbose
cap confirm variable _est____wooldidfullmodel
dis "PRE2C_ESTVAR=" (_rc==0) " K=" c(k)
qui wooldid y st yr onset, cluster(st) att
cap confirm variable _est____wooldidfullmodel
dis "POST2C_ESTVAR=" (_rc==0) " K=" c(k)

* --- 2d: back-to-back verbose runs -- does the count grow without bound?
use "rB_base.dta", clear
qui wooldid y st yr onset, cluster(st) att verbose
local ka = c(k)
qui wooldid y st yr onset, cluster(st) att verbose
local kb = c(k)
qui wooldid y st yr onset, cluster(st) att verbose
local kc = c(k)
dis "VERB_REPEAT_K=" `ka' " " `kb' " " `kc'

dis "ALLDONE"
