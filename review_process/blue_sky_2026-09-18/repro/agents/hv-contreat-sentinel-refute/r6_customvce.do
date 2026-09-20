* =============================================================================
* R6 -- SKEPTIC, NEW FINDING under test:
*   customvce() / unconditionalse (the `regress'-based estimation path) return
*   DIFFERENT POINT ESTIMATES from the default reghdfe path on a continuous
*   treatment.  The help documents customvce() as changing the variance only
*   (sthlp:264, 270-273: "passes the user specified text without modification
*   into the vce() option"); nothing says the estimates move.
*
*   r5_ameoracle.do established that the DEFAULT path's main_ame reproduces an
*   independent hand-built saturated-regression oracle to 10 decimal places
*   (0.9127941639 both), so the default path is the correct side.
*
* This file is self-contained: it rebuilds its own panel.
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r6_customvce.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
which wooldid

clear
set obs 12
gen byte state = _n
gen int expyr = .
replace expyr = 2012 if inlist(state,1,2,3,4)
replace expyr = 2014 if inlist(state,5,6,7)
replace expyr = 2016 if inlist(state,8,9)
expand 60
bys state: gen int hosp = _n
gen long hid = state*1000 + hosp
gen double expo = 0.2 + 1.6*runiform()
gen double hfe  = rnormal(0,1)
expand 12
bys hid: gen int year = 2007 + _n
gen double dose = expo*(expyr<. & year>=expyr)
gen double y    = 4 + hfe + 0.03*(year-2007) + 0.8*dose + rnormal(0,0.05)
gen double doseinv = expo + 0.0*dose          /* time-invariant dose */
gen double yinv = 4 + hfe + 0.03*(year-2007) + 0.8*doseinv*(expyr<. & year>=expyr) + rnormal(0,0.05)
compress
save "r6_panel.dta", replace

cap prog drop brow
prog def brow
  syntax , tag(string)
  matrix BB = e(b)
  local bn : colfullnames BB
  local j 0
  foreach nm of local bn {
    local ++j
    dis "B_`tag'|" %-22s "`nm'" "|" %16.10f BB[1,`j']
  }
  dis "N_`tag'=" e(N) "  r2=" e(r2) "  dfr=" e(df_r)
end

* ============ 1. continuous treatment, untreatedZero, five VCE paths ========
foreach CFG in DEFAULT ROBUSTOPT CVROBUST CVOLS CVCLUS UNCOND {
  use "r6_panel.dta", clear
  local opt ""
  if "`CFG'" == "ROBUSTOPT" local opt robust
  if "`CFG'" == "CVROBUST"  local opt customvce(robust)
  if "`CFG'" == "CVOLS"     local opt customvce(ols)
  if "`CFG'" == "CVCLUS"    local opt customvce(cluster hid)
  if "`CFG'" == "UNCOND"    local opt unconditionalse cluster(hid)
  dis "############### 1/`CFG' ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(untreatedZero) `opt'
  dis "RC1_`CFG'=" _rc
  if _rc == 0 brow , tag(1`CFG')
}

* ============ 2. BINARY treatment: does the same split appear? ==============
foreach CFG in DEFAULT CVROBUST {
  use "r6_panel.dta", clear
  local opt ""
  if "`CFG'" == "CVROBUST" local opt customvce(robust)
  dis "############### 2/`CFG' (binary) ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) `opt'
  dis "RC2_`CFG'=" _rc
  if _rc == 0 brow , tag(2`CFG')
}

* ============ 3. contreatcontrols(basic) ===================================
foreach CFG in DEFAULT CVROBUST {
  use "r6_panel.dta", clear
  local opt ""
  if "`CFG'" == "CVROBUST" local opt customvce(robust)
  dis "############### 3/`CFG' (basic) ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(basic) `opt'
  dis "RC3_`CFG'=" _rc
  if _rc == 0 brow , tag(3`CFG')
}

* ============ 4. dose defined in EVERY period (no structural zeros) ========
foreach CFG in DEFAULT CVROBUST {
  use "r6_panel.dta", clear
  local opt ""
  if "`CFG'" == "CVROBUST" local opt customvce(robust)
  dis "############### 4/`CFG' (time-invariant dose, basic) ###############"
  cap noisily wooldid yinv state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(doseinv) contreatcontrols(basic) `opt'
  dis "RC4_`CFG'=" _rc
  if _rc == 0 brow , tag(4`CFG')
}

* ============ 5. smaller, fully-in-window design ===========================
* every treated cohort's whole post path fits inside espostlength(), so there
* are no out-of-window treated cells at all.
use "r6_panel.dta", clear
keep if year >= 2010 & year <= 2014
foreach CFG in DEFAULT CVROBUST {
  preserve
  local opt ""
  if "`CFG'" == "CVROBUST" local opt customvce(robust)
  dis "############### 5/`CFG' (2010-2014 window) ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(2) espostlength(2) ///
      contreat(dose) contreatcontrols(untreatedZero) `opt'
  dis "RC5_`CFG'=" _rc
  if _rc == 0 brow , tag(5`CFG')
  restore
}

dis "R6_DONE=1"
