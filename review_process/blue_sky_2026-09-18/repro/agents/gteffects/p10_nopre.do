* p10_nopre.do -- designs in which a treated cohort-by-period dummy is collinear with
* the fixed effects, so the fitted model carries an OMITTED (o.) cell.  Two realistic
* shapes: (K) a cohort treated in the FIRST observed period (no pre-treatment rows at
* all), and (L) a control sample that stops being collected before the treated one does.
* Question: is the cell silently missing from e(gteffects), and does the documented
* weighted-mean identity then fail?
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime

*=========== K: one cohort treated in the first observed period ==============
clear
set obs 90
gen long unit = _n
gen byte cohp = 1
qui replace cohp = 2 if _n > 40 & _n <= 65
qui replace cohp = 3 if _n > 65
gen double ttre = .
qui replace ttre = 2000 if cohp == 2
qui replace ttre = 2004 if cohp == 3
gen double ai = rnormal()
qui expand 9
bys unit: gen int year = 1999 + _n
tempvar tsh
bys year: gen double `tsh' = rnormal() if _n==1
bys year: replace `tsh' = `tsh'[1]
gen double te = 0
qui replace te = 1 + 0.2*(year-ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`tsh' + te + rnormal()*0.3
tempfile kd
qui save `kd'
qui tab cohp if year == 2000

dis "--- K1: default safety ---"
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2)
dis "K1_RC=" _rc

dis "--- K2: safety(off) ---"
use `kd', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) safety(off)
dis "K2_RC=" _rc
gtcheck , tag(K2_nopre_off) pre es
cap noi matrix list e(gteffects), format(%12.6g)
* how many treated (cohort, calendar year) cells does the estimation sample realise?
qui gen byte esx = e(sample)
qui gegen byte tgc = tag(cohp year) if esx & cohp > 1
qui count if tgc == 1
dis "K2_REALISED_TREATED_CELLS=" r(N)

dis "--- K3: verbose so the omitted coefficient is visible ---"
use `kd', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) safety(off) verbose
dis "K3_RC=" _rc
qui est restore ___wooldidfullmodel
local bn : colnames e(b)
local nomit 0
foreach v of local bn {
  if strpos("`v'", "o.") != 0 & strpos("`v'", "0b.") == 0 local ++nomit
}
dis "K3_OMITTED_TERMS=" `nomit'
foreach v of local bn {
  if strpos("`v'", "o.") != 0 & strpos("`v'", "0b.") == 0 dis "K3_OMIT: `v'"
}

*=========== L: control sample ends before the treated sample does ============
clear
set obs 70
gen long unit = _n
gen byte cohp = cond(_n<=45, 1, 2)
gen double ttre = cond(cohp==2, 2004, .)
gen double ai = rnormal()
qui expand 13
bys unit: gen int year = 1999 + _n
qui drop if cohp == 1 & year > 2008
tempvar tsh2
bys year: gen double `tsh2' = rnormal() if _n==1
bys year: replace `tsh2' = `tsh2'[1]
gen double te = 0
qui replace te = 1 + 0.2*(year-ttre) if year >= ttre & ttre < .
gen double y = 2 + ai + 0.5*`tsh2' + te + rnormal()*0.3
tempfile ld
qui save `ld'

dis "--- L1: default safety ---"
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2)
dis "L1_RC=" _rc

dis "--- L2: safety(off) ---"
use `ld', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) safety(off)
dis "L2_RC=" _rc
gtcheck , tag(L2_ctrlend_off) pre es
cap noi matrix list e(gteffects), format(%12.6g)
qui gen byte esy = e(sample)
qui gegen byte tgc2 = tag(cohp year) if esy & cohp > 1
qui count if tgc2 == 1
dis "L2_REALISED_TREATED_CELLS=" r(N)

dis "DONE_P10"
