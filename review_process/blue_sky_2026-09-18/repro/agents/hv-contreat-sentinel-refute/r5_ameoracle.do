* =============================================================================
* R5 -- SKEPTIC, NEW LEAD: is the ame POINT ESTIMATE right on an exposure x post
*       design?  In r4 section (E) a low-noise panel whose true average marginal
*       effect of the dose is EXACTLY 0.8 came back as main_ame = -0.3095.
*       Here I build an independent oracle: the saturated Wooldridge regression
*       written out by hand, and the cell-weighted mean of its own dose slopes.
*
* Design: 12 states (i), 60 hospitals each, T = 2008..2019.  Onsets 2012/2014/2016,
*         3 never-expanding states.  dose = expo_h * 1{year >= expyr}.
*         y = 4 + hfe_h + 0.03*(year-2007) + 0.8*dose + N(0, sd).
*         TRUE ame = 0.8 at every horizon.  TRUE dose-zero intercept = 0.
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r5_ameoracle.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
which wooldid

* ------------------------------------------------------------------- DGP ---
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
gen byte trt    = (expyr<. & year>=expyr)
compress
save "r5_panel.dta", replace
count
sum dose if trt==1

* ============================ ORACLE 1 ======================================
* pooled slope with hospital fixed effects: must recover 0.8 tightly
dis "=== ORACLE1: reghdfe y dose, absorb(hid year) ==="
cap noisily reghdfe y dose, absorb(hid year)
dis "ORACLE1_b=" %14.10f _b[dose] "  se=" %12.8f _se[dose]

dis "=== ORACLE1b: reghdfe y dose, absorb(state year) (wooldid's own FE set) ==="
cap noisily reghdfe y dose, absorb(state year)
dis "ORACLE1b_b=" %14.10f _b[dose] "  se=" %12.8f _se[dose]

* ============================ ORACLE 2 ======================================
* the saturated Wooldridge/untreatedZero regression, written out by hand:
*   y = a_state + d_year + sum_cells [ g_c * 1{cell} + b_c * 1{cell} * dose ]
* then ame = obs-weighted mean of b_c over treated cells
*      att = obs-weighted mean of g_c over treated cells (dose-zero intercept)
use "r5_panel.dta", clear
egen long cellg = group(state year) if trt == 1
replace cellg = 0 if missing(cellg)
qui levelsof cellg if cellg > 0, local(cells)
dis "ORACLE2_ncells=" wordcount("`cells'")
dis "=== ORACLE2: saturated cell intercepts + cell dose slopes ==="
qui reghdfe y ib0.cellg ib0.cellg#c.dose, absorb(state year)
dis "ORACLE2_rc=" _rc
local sb 0
local sg 0
local sw 0
local ndrop 0
foreach c of local cells {
  qui count if cellg == `c'
  local w = r(N)
  local bc = _b[`c'.cellg#c.dose]
  local gc = _b[`c'.cellg]
  if missing(`bc') | `bc' == 0 local ++ndrop
  local sb = `sb' + `w'*`bc'
  local sg = `sg' + `w'*`gc'
  local sw = `sw' + `w'
}
dis "ORACLE2_ndroppedslopes=" `ndrop'
dis "ORACLE2_ame=" %14.10f (`sb'/`sw')
dis "ORACLE2_att=" %14.10f (`sg'/`sw')
dis "ORACLE2_totobs=" `sw'

* ============================== WOOLDID =====================================
foreach CFG in DEFAULT ROBUST CLHOSP FEHID {
  use "r5_panel.dta", clear
  local opt ""
  if "`CFG'" == "ROBUST" local opt customvce(robust)
  if "`CFG'" == "CLHOSP" local opt cluster(hid)
  if "`CFG'" == "FEHID"  local opt fe(hid)
  dis "############### WOOLDID `CFG' ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(untreatedZero) `opt'
  dis "WRC_`CFG'=" _rc
  if _rc continue
  matrix BW = e(b)
  local bn : colfullnames BW
  local j 0
  foreach nm of local bn {
    local ++j
    dis "W_`CFG'|" %-22s "`nm'" "|b=" %14.10f BW[1,`j']
  }
  dis "W_`CFG'_N=" e(N)
}

* =========== does a SHORTER panel (no out-of-window horizons) fix it? =======
use "r5_panel.dta", clear
keep if year <= 2014
dis "############### WOOLDID SHORTPANEL (2008-2014, onset 2012 only treated) ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(2) ///
    contreat(dose) contreatcontrols(untreatedZero)
dis "WRC_SHORT=" _rc
if _rc == 0 {
  matrix BS = e(b)
  local bn : colfullnames BS
  local j 0
  foreach nm of local bn {
    local ++j
    dis "W_SHORT|" %-22s "`nm'" "|b=" %14.10f BS[1,`j']
  }
}

dis "R5_DONE=1"
