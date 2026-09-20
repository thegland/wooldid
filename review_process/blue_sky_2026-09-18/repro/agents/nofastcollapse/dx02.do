* DIAGNOSTIC dx02: wooldid_fcol keeps the cell's first row with
*     qui gen long `ord' = _n
*     qui gcollapse (sum) `wv' (min) `pick' = `ord', by(`key') merge replace
*     qui keep if `ord' == `pick'
* If gcollapse types `pick' as FLOAT, then for an estimation sample larger than
* 2^24 = 16,777,216 rows the comparison `ord' == `pick' is a long-vs-float
* comparison and can keep the wrong number of rows per cell.
* This file only measures the TYPE and the comparison behaviour - it does not
* need a 17-million-row panel.
clear all
set more off
set processors 2
dis "default c(type) = " c(type)

* ---- 1. what type does gcollapse give (min) of a long? ----
set obs 1000
gen long ord = _n
gen byte g = mod(_n,10)
qui gcollapse (min) pick = ord, by(g) merge
dis "T1 type of pick = " "`: type pick'"
qui count if ord == pick
dis "T1 rows kept (10 cells, must be 10): " r(N)

* ---- 2. the same, with ord values above 2^24 ----
clear
set obs 200
gen long ord = 16777200 + _n
gen byte g = 1
qui gcollapse (min) pick = ord, by(g) merge
dis "T2 type of pick = " "`: type pick'"
qui count if ord == pick
dis "T2 rows kept in ONE cell (must be exactly 1): " r(N)

* ---- 3. and with the double weight variable alongside, as wooldid does ----
clear
set obs 200
gen long ord = 16777200 + _n
gen byte g = 1
gen double w = 1
qui gcollapse (sum) w (min) pick = ord, by(g) merge replace
dis "T3 type of pick = " "`: type pick'"  "   type of w = " "`: type w'"
qui count if ord == pick
dis "T3 rows kept in ONE cell (must be exactly 1): " r(N)
dis "T3 summed weight (must be 200): " w[1]

* ---- 4. under set type double, in case the user changed it ----
clear
set type double
set obs 200
gen long ord = 16777200 + _n
gen byte g = 1
qui gcollapse (min) pick = ord, by(g) merge
dis "T4 type of pick = " "`: type pick'"
qui count if ord == pick
dis "T4 rows kept in ONE cell (must be exactly 1): " r(N)
set type float
dis "DX02_DONE"
