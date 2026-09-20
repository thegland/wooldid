* ================================================================
* p09_simB.do -- DESIGN B: i = unit AND cluster(unit), i.e. clustering
*   AT the level of i.  wooldid prints a caution here and the project
*   has measured the joint pre-trend test as oversized on this shape.
*   Nobody has measured what that does to pretrendpower()'s MDT.
*   Report the (size, power) PAIR and the ratio mdt_reported/mdt_true.
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
  qui replace y = a_i + f_t + 1.0*D + `delta'*trendterm + e
end

set obs 80
gen long unit = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 5 if mod(unit,5)==1
replace onset = 7 if mod(unit,5)==2
replace onset = 9 if mod(unit,5)==3
expand 10
bysort unit: gen int t = _n
sort unit t
gen double f_t = 0.1*t + 0.3*sin(t)
gen byte D = (onset<. & t>=onset)
gen double relyr = t - onset
gen double trendterm = cond(onset<., relyr + 1, 0)
gen double e = 0
gen double y = 0

global WDCALL "y unit t onset"
local EST "cluster(unit) esfixedbaseperiod esprelength(3) espostlength(2)"

drawy 0
cap noi wooldid $WDCALL , `EST' pretrendpower(50 80)
dis "B_PILOT_RC=" _rc
dis "B df_r = " e(df_r)

timer clear
timer on 1
runarm , delta(0) reps(250) out("B0_null.csv") est(`"`EST'"')
timer off 1
dis "B_ARM0_DONE"
timer list 1

preserve
  qui import delimited using "B0_null.csv", clear varnames(1)
  qui sum mdt80, detail
  global D80 = r(p50)
restore
dis "B MEDIAN_MDT80 = $D80"

timer on 2
runarm , delta($D80) reps(220) out("B1_pow80.csv") est(`"`EST'"')
timer off 2
dis "B_ARM1_DONE"
timer list 2
dis "B_ALL_DONE"
