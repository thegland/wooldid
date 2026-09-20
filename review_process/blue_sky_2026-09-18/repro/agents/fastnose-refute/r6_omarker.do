*=====================================================================
* r6_omarker.do -- does `nose' suppress margins' "o." omission marker?
*
* wooldid_sglbl (ado:12575) drops an event-study horizon iff its margins colname
* carries the "o." prefix.  ado:1975 asserts nose suppresses that marker.  Here is
* the direct demonstration, with no wooldid involved: a dose-slope design in which
* three period-by-dose terms are omitted for collinearity.
*
* RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r6_omarker.do
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
replace x = runiform() if t >= 4          // NO dose variation at t = 1,2,3
gen double y = 0.3*t + 1.2*x + rnormal()
* dose slopes are period-specific with NO main effect -- wooldid's actual shape
qui reg y i.t i.t#c.x

dis _n "=== df() arm ==="
qui margins , dydx(x) over(t) df(500) noestimcheck
matrix A = r(table)
local an : colnames A
dis "DF_table_cols: `an'"
dis "DF_b_cols: `: colnames r(b)'"
local nd = 0
foreach w of local an {
  if strpos("`w'", "o.") != 0 local ++nd
}
dis "DF_n_omitted_marked=" `nd'

dis _n "=== nose arm ==="
qui margins , dydx(x) over(t) nose noestimcheck
matrix B = r(table)
local bn : colnames B
dis "NOSE_table_cols: `bn'"
dis "NOSE_b_cols: `: colnames r(b)'"
local nd2 = 0
foreach w of local bn {
  if strpos("`w'", "o.") != 0 local ++nd2
}
dis "NOSE_n_omitted_marked=" `nd2'

dis "DONE_R6"
