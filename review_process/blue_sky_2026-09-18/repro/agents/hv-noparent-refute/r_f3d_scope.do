* ============================================================================
* SKEPTIC FOLLOW-UP 3 for hv-noparent-F3: how far does the undefended domain go?
*  (a) event study: does a straddling outcome give a PARTIALLY contaminated
*      horizon table -- some rows fine, some not -- with full inference and no
*      message?
*  (b) the cell-weighted estimands (att_it) on the same data.
*  (c) contreatelasticitytype(eydx) and (eyex), which the finding lumps into the
*      same family by assertion and never measured.
* The panel is the pole-test panel shifted by -8, where the treated-post fitted
* values straddle zero. The unshifted panel is the control.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f3d_scope.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31337
set sortseed 31337

which wooldid

clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = 8 + rnormal()*0.5
gen double dose0 = 0.5 + runiform()*2
expand 10
bysort unit: gen int t = _n
gen double ybase = ufe + 0.05*t + rnormal()*0.3
replace ybase = ybase + 0.8*(t>=ttre) if ttre<.
gen double dose = dose0*(t>=ttre & ttre<.)
tempfile P
save `P', replace

* ---- (a) event study, control (no shift) -----------------------------------
use `P', clear
gen double y = ybase
capture noisily wooldid y unit t ttre, att robust semielasticity esprelength(2) espostlength(2) jointtests
dis "A_CTRL_RC=" _rc

* ---- (a) event study, shifted so fitted values straddle zero ---------------
use `P', clear
gen double y = ybase - 8
capture noisily wooldid y unit t ttre, att robust semielasticity esprelength(2) espostlength(2) jointtests
dis "A_SHIFT_RC=" _rc

* ---- (b) cell-weighted estimand on the same straddling data ----------------
use `P', clear
gen double y = ybase - 8
capture noisily wooldid y unit t ttre, att_it robust semielasticity
dis "B_ATTIT_RC=" _rc

* ---- (c) contreatelasticitytype on a straddling outcome --------------------
use `P', clear
gen double y = ybase
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(eydx)
dis "C_CTRL_EYDX_RC=" _rc

use `P', clear
gen double y = ybase - 8
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(eydx)
dis "C_SHIFT_EYDX_RC=" _rc

use `P', clear
gen double y = ybase - 8
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(untreatedZero) contreatelasticitytype(eyex)
dis "C_SHIFT_EYEX_RC=" _rc

dis "DONE_F3D"
