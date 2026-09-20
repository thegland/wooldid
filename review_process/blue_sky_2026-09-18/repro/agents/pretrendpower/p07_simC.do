* ================================================================
* p07_simC.do -- DESIGN C: the geometry the help itself recommends
*   for pretrendpower() (sthlp:1334: "a coarser i (cohorts, or
*   treatment-onset groups) with unit fixed effects supplied through
*   fe()"), with cluster() STRICTLY FINER than i -- the clustering
*   ado:1443 names as giving a valid joint test.
*     i = onset-cohort code (4 codes), fe(unit), cluster(unit), G=80
*     AR(1) errors within unit, so the clustering is doing real work.
*   arm0 delta=0; arm1 delta=median(mdt80); arm2 delta=median(mdt50)
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

set obs 80
gen long unit = _n
gen byte coh = 4   // never-treated cohort code; wooldid requires positive i
replace coh = 1 if mod(unit,5)==1
replace coh = 2 if mod(unit,5)==2
replace coh = 3 if mod(unit,5)==3
gen int ttre = .
replace ttre = 5 if coh==1
replace ttre = 7 if coh==2
replace ttre = 9 if coh==3
expand 10
bysort unit: gen int t = _n
sort unit t
gen double f_t = 0.1*t + 0.3*sin(t)
gen byte D = (ttre<. & t>=ttre)
gen double relyr = t - ttre
gen double trendterm = cond(ttre<., relyr + 1, 0)
gen double e = 0
gen double y = 0

global WDCALL "y coh t ttre"
local EST "fe(unit) cluster(unit) esfixedbaseperiod esprelength(3) espostlength(2)"

* one loud run so the log shows the header, the df and whether any
* at-or-above-i caution fires on this geometry
drawy 0
cap noi wooldid $WDCALL , `EST' pretrendpower(50 80)
dis "C_PILOT_RC=" _rc
dis "C df_r = " e(df_r)

timer clear
timer on 1
runarm , delta(0) reps(220) out("C0_null.csv") est(`"`EST'"')
timer off 1
dis "C_ARM0_DONE"
timer list 1

preserve
  qui import delimited using "C0_null.csv", clear varnames(1)
  qui sum mdt80, detail
  global D80 = r(p50)
  qui sum mdt50, detail
  global D50 = r(p50)
restore
dis "C MEDIAN_MDT80 = $D80"
dis "C MEDIAN_MDT50 = $D50"

timer on 2
runarm , delta($D80) reps(220) out("C1_pow80.csv") est(`"`EST'"')
timer off 2
dis "C_ARM1_DONE"
timer list 2

dis "C_ALL_DONE"
