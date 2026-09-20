* ============================================================================
* SKEPTIC REPRO for hv-noparent-F2
* CLAIM UNDER TEST: contreatelasticitytype(eyex|eydx|dyex) given WITHOUT
*   contreat() leaves every returned number bit-identical but changes the
*   printed results heading to an elasticity heading, and sets e(wdid_cetype).
* Also tests the prober's SURFACE claims: event-study table header, ES graph
*   titles, histogramcohorteffects labels, subgroup headers.
* Self-contained. Generates its own binary staggered panel.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f2_cet.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 71135
set sortseed 71135

which wooldid

* ---------- my own data: 90 units, 12 periods, 3 treated cohorts + never ----
clear
set obs 90
gen long unit = _n
gen byte grp = mod(_n,4)            /* 0 = never treated */
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen sgv = mod(_n,2)                 /* a subgroup variable */
gen double ufe = rnormal()
expand 12
bysort unit: gen int t = _n
gen double y = 2 + ufe + 0.15*t + rnormal()
replace y = y + 0.8*(t>=ttre) if ttre<.
gen double w = 1 + runiform()
tempfile dat
save `dat', replace

* ---------- leg A: baseline, no contreatelasticitytype -----------------------
use `dat', clear
capture noisily wooldid y unit t ttre, att robust
dis "A_RC=" _rc
scalar bA = el(e(wdidmainresults),1,1)
dis "A_ATT_HEX=" %21x bA
dis "A_CETYPE=[" e(wdid_cetype) "]"
matrix MA = e(wdidmainresults)

* ---------- legs B,C,D: each elasticity type, no contreat() ------------------
foreach ct in eyex eydx dyex {
  use `dat', clear
  capture noisily wooldid y unit t ttre, att robust contreatelasticitytype(`ct')
  dis "`ct'_RC=" _rc
  scalar b`ct' = el(e(wdidmainresults),1,1)
  dis "`ct'_ATT_HEX=" %21x b`ct'
  dis "`ct'_CETYPE=[" e(wdid_cetype) "]"
  matrix M`ct' = e(wdidmainresults)
  dis "`ct'_MRELDIF=" mreldif(M`ct', MA)
}

* ---------- leg E: does it reach the EVENT-STUDY table + graph titles? -------
use `dat', clear
capture noisily wooldid y unit t ttre, att robust esprelength(2) espostlength(2) ///
   makeplots saveplots("f2pfx") replace
dis "E_RC=" _rc
use `dat', clear
capture noisily wooldid y unit t ttre, att robust esprelength(2) espostlength(2) ///
   makeplots saveplots("f2pfxE") replace contreatelasticitytype(eyex)
dis "E2_RC=" _rc

* ---------- leg F: subgroup header --------------------------------------------
use `dat', clear
capture noisily wooldid y unit t ttre, att robust subgroups(sgv) contreatelasticitytype(eyex)
dis "F_RC=" _rc

* ---------- leg G: histogramcohorteffects labels ------------------------------
use `dat', clear
capture noisily wooldid y unit t ttre, att_i robust histogramcohorteffects contreatelasticitytype(eyex)
dis "G_RC=" _rc

* ---------- leg H: POISSON without poisexpresults, no contreat() -------------
* the prober could not measure this (separation). Build a clean count DGP.
clear
set obs 90
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = rnormal()*0.2
expand 12
bysort unit: gen int t = _n
gen double lam = exp(2.2 + ufe + 0.02*t + 0.3*(t>=ttre & ttre<.))
gen ycount = rpoisson(lam)
tempfile pd
save `pd', replace

use `pd', clear
capture noisily wooldid ycount unit t ttre, att robust poisson
dis "H0_RC=" _rc
use `pd', clear
capture noisily wooldid ycount unit t ttre, att robust poisson contreatelasticitytype(eyex)
dis "H_RC=" _rc
dis "H_CETYPE=[" e(wdid_cetype) "] H_CETYPERAW=[" e(wdid_cetyperaw) "]"

dis "DONE_F2"
