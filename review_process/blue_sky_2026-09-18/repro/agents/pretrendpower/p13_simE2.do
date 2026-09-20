* ================================================================
* p13_simE2.do -- DESIGN E, extra power arms.
*   E1 measured 75.3% nominal / 58.0% size-adjusted power at the
*   REPORTED mdt_wald(80) = 0.5463.  Extrapolating the noncentral
*   chi-squared, an honest (size-corrected) 5% test needs about
*   0.6456 for 80% power.  Verify that directly.
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
gen byte coh = 4
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

timer clear
timer on 1
runarm , delta(0.6456) reps(300) out("E4_d6456.csv") est(`"`EST'"')
timer off 1
dis "E4_DONE"
timer list 1
