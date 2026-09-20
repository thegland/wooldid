* ================================================================
* p08_simE.do -- DESIGN E: small-G state-level policy panel, the
*   canonical applied DiD shape, on a clustering geometry wooldid
*   does NOT caution about (cluster strictly finer than i).
*     16 states, 10 years, state-level outcome
*     i = onset-cohort code (4 codes), fe(state), cluster(state)
*     G = 16  ->  df_r = 15,  invttail(15,.025) = 2.131 vs 1.96
*   This is where the t-critical-value-inside-a-normal-probability
*   hybrid in mdt_nis / roth_nissize should bite, and where a small-G
*   F test can be oversized without any of wooldid's cautions firing.
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/pretrendpower/ptplib.do"

capture program drop drawy
program define drawy
  args delta
  qui replace e = rnormal()
  forvalues k = 2/10 {
    qui replace e = 0.6*e[_n-1] + sqrt(0.64)*rnormal() if t == `k'
  }
  qui replace y = f_t + 1.0*D + `delta'*trendterm + e
end

set obs 16
gen long st = _n
gen byte coh = 4   // never-treated cohort code; wooldid requires positive i
replace coh = 1 if st <= 4
replace coh = 2 if st > 4 & st <= 8
replace coh = 3 if st > 8 & st <= 12
gen int ttre = .
replace ttre = 5 if coh == 1
replace ttre = 7 if coh == 2
replace ttre = 9 if coh == 3
expand 10
bysort st: gen int t = _n
sort st t
gen double f_t = 0.1*t + 0.3*sin(t)
gen byte D = (ttre<. & t>=ttre)
gen double relyr = t - ttre
gen double trendterm = cond(ttre<., relyr + 1, 0)
gen double e = 0
gen double y = 0

global WDCALL "y coh t ttre"
local EST "fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)"

drawy 0
cap noi wooldid $WDCALL , `EST' pretrendpower(50 80)
dis "E_PILOT_RC=" _rc
dis "E df_r = " e(df_r)

* pilot: the NIS quadrature is expensive, so take the median of
* e(roth_nissize_att) and mdt_nis(80) over a dozen null draws
file open ph using "E_pilot.csv", write text replace
file write ph "rep,nissize,mdtnis80,mdtwald80" _n
file close ph
forvalues r = 1/12 {
  set seed `= 20260000 + `r''
  drawy 0
  cap qui wooldid $WDCALL , `EST' pretrendpower(80)
  local rc = _rc
  local ns = .
  local mn = .
  local mw = .
  if `rc' == 0 {
    local ns = e(roth_nissize_att)
    mat RR = e(roth_power_att)
    local mn = RR[1,3]
    local mw = RR[1,2]
  }
  file open ph using "E_pilot.csv", write text append
  file write ph "`r',`ns',`mn',`mw'" _n
  file close ph
}
preserve
  qui import delimited using "E_pilot.csv", clear varnames(1)
  qui sum nissize, detail
  global ENIS = r(p50)
  qui sum mdtnis80, detail
  global EDNIS = r(p50)
restore
dis "E MEDIAN reported nissize = $ENIS"
dis "E MEDIAN reported mdt_nis(80) = $EDNIS"

timer clear
timer on 1
runarm , delta(0) reps(400) out("E0_null.csv") est(`"`EST'"')
timer off 1
dis "E_ARM0_DONE"
timer list 1

preserve
  qui import delimited using "E0_null.csv", clear varnames(1)
  qui sum mdt80, detail
  global D80 = r(p50)
restore
dis "E MEDIAN_MDT80 = $D80"

timer on 2
runarm , delta($D80) reps(300) out("E1_pow80.csv") est(`"`EST'"')
timer off 2
dis "E_ARM1_DONE"
timer list 2

timer on 3
runarm , delta($EDNIS) reps(300) out("E3_nis80.csv") est(`"`EST'"')
timer off 3
dis "E_ARM3_DONE"
timer list 3
dis "E_ALL_DONE"
