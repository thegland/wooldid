*=====================================================================
* r5_mechanism.do -- the margins-level mechanism the prober left inconclusive.
*
* wooldid's SI-7 estimability screen reads margins' not-estimable marker.  ado:1973
* sets `marginsdof' to `nose' under fastnose; ado:1981 gives the BINARY
* subgroup-by-horizon call back its df(`df_r') (BS-26) but the CONTINUOUS twin at
* ado:3290 still threads `marginsdof'.  Question: does `nose' really suppress the
* marker that df() leaves visible?
*
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r5_mechanism.do
*=====================================================================
clear all
set more off
set processors 2
set seed 4711
set sortseed 4711

clear
set obs 600
gen int t = 1 + mod(_n-1, 6)
gen double x = 0
replace x = runiform() if t >= 4          // no dose variation at t = 1,2,3
gen double y = 0.3*t + 1.2*x + rnormal()
qui reg y i.t##c.x

dis _n "===== margins with df() (variance computed) ====="
qui margins , dydx(x) over(t) df(500) noestimcheck
matrix A = r(table)
dis "DF_cols: `: colnames A'"
dis "DF_rows: `: rownames A'"
matrix list A, format(%12.0g)

dis _n "===== margins with nose (no variance) ====="
qui margins , dydx(x) over(t) nose noestimcheck
matrix B = r(table)
dis "NOSE_cols: `: colnames B'"
dis "NOSE_rows: `: rownames B'"
matrix list B, format(%12.0g)

dis _n "===== r(b) colnames both ways ====="
qui margins , dydx(x) over(t) df(500) noestimcheck
dis "DF_b_cols: `: colnames r(b)'"
dis "DF_b_n=" colsof(r(b))
qui margins , dydx(x) over(t) nose noestimcheck
dis "NOSE_b_cols: `: colnames r(b)'"
dis "NOSE_b_n=" colsof(r(b))

dis "DONE_R5"
