* ============================================================================
* SKEPTIC REPRO for hv-noparent-F3
* CLAIM UNDER TEST: semielasticity (margins eydx) has no domain screen, and on
*   an outcome whose FITTED values cross zero it returns a complete inference
*   row with no caution (a "silently wrong number"), or dies r(459), or dies
*   inside Mata r(3598), depending on data and weights.
* This probe separates two very different things:
*   (1) y merely takes NEGATIVE VALUES but the fitted values in the evaluated
*       subpopulation are strictly positive  -> eydx is DEFINED; a large number
*       there is arithmetically correct, not wrong.
*   (2) the FITTED values themselves straddle zero in the evaluated
*       subpopulation -> eydx is genuinely undefined.
* It also runs the decisive functional test: eydx must be INVARIANT to a
*   multiplicative rescaling of y and must MOVE under an additive shift.
* Self-contained.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f3_semielas.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 990211
set sortseed 990211

which wooldid

* ============================ PANEL P: well-behaved ==========================
* unit effects all well above zero, tiny noise -> fitted values >> 0
clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = 8 + rnormal()*0.5
expand 10
bysort unit: gen int t = _n
gen double e = rnormal()*0.3
gen double ybase = ufe + 0.05*t + e
replace ybase = ybase + 0.8*(t>=ttre) if ttre<.
gen double w = 1 + runiform()
tempfile P
save `P', replace

* ---- P1: level ATT ----------------------------------------------------------
use `P', clear
gen double y = ybase
capture noisily wooldid y unit t ttre, att robust
dis "P1_RC=" _rc
dis "P1_ATT=" %21x el(e(wdidmainresults),1,1)

* ---- P2: semielasticity on y (min y well above 0) ---------------------------
use `P', clear
gen double y = ybase
qui su y
dis "P2_MINY=" r(min)
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "P2_RC=" _rc
capture dis "P2_SEMI=" el(e(wdidmainresults),1,1)

* ---- P3: SCALE INVARIANCE. y*10. eydx must be IDENTICAL to P2. --------------
use `P', clear
gen double y = ybase*10
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "P3_RC=" _rc
capture dis "P3_SEMI_x10=" el(e(wdidmainresults),1,1)

* ---- P4: ADDITIVE SHIFT. y-7. level ATT must be unchanged, eydx must MOVE. --
* min(y) goes negative but fitted values in the treated-post subpop stay >0.
use `P', clear
gen double y = ybase - 7
qui su y
dis "P4_MINY=" r(min)
qui su y if t>=ttre & ttre<.
dis "P4_MINY_TREATEDPOST=" r(min)
capture noisily wooldid y unit t ttre, att robust
dis "P4_LEVEL_RC=" _rc
capture dis "P4_LEVELATT=" %21x el(e(wdidmainresults),1,1)
use `P', clear
gen double y = ybase - 7
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "P4_RC=" _rc
capture dis "P4_SEMI=" el(e(wdidmainresults),1,1)

* ---- P5: the SAME shifted outcome, scaled by 10. eydx must equal P4. -------
use `P', clear
gen double y = (ybase - 7)*10
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "P5_RC=" _rc
capture dis "P5_SEMI=" el(e(wdidmainresults),1,1)

* ============================ PANEL Q: fitted values STRADDLE zero ===========
* unit effects spread across zero with tiny noise, so yhat itself crosses 0
* inside the treated-post region that margins evaluates.
clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = -4 + 8*runiform()
expand 10
bysort unit: gen int t = _n
gen double y = ufe + 0.05*t + rnormal()*0.2
replace y = y + 0.8*(t>=ttre) if ttre<.
gen double w = 1 + runiform()
tempfile Q
save `Q', replace

use `Q', clear
qui su y if t>=ttre & ttre<.
dis "Q_MINY_TREATEDPOST=" r(min)
dis "Q_MAXY_TREATEDPOST=" r(max)
capture noisily wooldid y unit t ttre, att robust
dis "Q0_LEVEL_RC=" _rc
capture dis "Q0_LEVELATT=" el(e(wdidmainresults),1,1)

use `Q', clear
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "Q1_RC=" _rc
capture dis "Q1_SEMI=" el(e(wdidmainresults),1,1)
capture dis "Q1_SE=" el(e(wdidmainresults),1,2)

use `Q', clear
capture noisily wooldid y unit t ttre, att robust semielasticity [aw=w]
dis "Q2_AW_RC=" _rc
capture dis "Q2_SEMI=" el(e(wdidmainresults),1,1)

use `Q', clear
capture noisily wooldid y unit t ttre, att robust semielasticity [pw=w]
dis "Q3_PW_RC=" _rc

* control: same data, same weights, NO semielasticity
use `Q', clear
capture noisily wooldid y unit t ttre, att robust [aw=w]
dis "Q4_CONTROL_AW_RC=" _rc

* ============================ PANEL R: fitted value exactly near zero =======
* a tighter straddle: unit effects centred on 0 with small spread
clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = rnormal()*0.4
expand 10
bysort unit: gen int t = _n
gen double y = ufe + 0.02*t + rnormal()*0.1
replace y = y + 0.3*(t>=ttre) if ttre<.
gen double w = 1 + runiform()
tempfile R
save `R', replace
use `R', clear
qui su y if t>=ttre & ttre<.
dis "R_MINY_TREATEDPOST=" r(min)
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "R1_RC=" _rc
capture dis "R1_SEMI=" el(e(wdidmainresults),1,1)
use `R', clear
capture noisily wooldid y unit t ttre, att robust semielasticity [aw=w]
dis "R2_AW_RC=" _rc
capture dis "R2_SEMI=" el(e(wdidmainresults),1,1)

dis "DONE_F3"
